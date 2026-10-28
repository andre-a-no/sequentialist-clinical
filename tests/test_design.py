# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Novikov Laboratories LLC
# Commercial licenses for use outside the AGPL: see COMMERCIAL.md
"""Designs: boundaries against an independent computation, maximum information against closed formulas."""
import math
import unittest
from statistics import NormalDist

from sequentialist.clinical import TrialSpec, check_trial, plan_report

from .helpers import references, spec

Z = NormalDist().inv_cdf


class BoundaryTest(unittest.TestCase):
    def test_boundaries_match_the_multivariate_normal_computation(self):
        # tests/data/references.json: scipy's multivariate normal (Genz), see tools/make_references.py
        for ref in references()["boundaries"]:
            with self.subTest(spending=ref["spending"], analyses=ref["analyses"]):
                design = {"alpha": ref["alpha"], "analyses": ref["analyses"], "spending": ref["spending"],
                          "rho": ref["rho"], "two_sided": ref["two_sided"], "maximum": 300, "power": None}
                designed = spec(design=design).data["designed"]
                for ours, theirs in zip(designed["bounds"], ref["bounds"]):
                    self.assertAlmostEqual(ours, theirs, delta=2e-6)

    def test_published_obrien_fleming_boundaries(self):
        # Lan-DeMets O'Brien-Fleming, three equally spaced analyses, one-sided 0.025 (gsDesign, sfLDOF)
        designed = spec(design={"analyses": [1 / 3, 2 / 3, 1.0], "maximum": 300, "power": None}).data["designed"]
        for ours, published in zip(designed["bounds"], (3.7103, 2.5114, 1.9930)):
            self.assertAlmostEqual(ours, published, places=4)


class MaximumTest(unittest.TestCase):
    """With one analysis the design is a fixed-sample test: the maximum is the closed formula."""

    one = {"analyses": [1.0], "power": 0.9}

    def test_events_by_schoenfeld(self):
        designed = spec(design={**self.one, "assumptions": {"hazard_ratio": 0.7}}).data["designed"]
        events = (Z(0.975) + Z(0.9)) ** 2 / (0.25 * math.log(0.7) ** 2)
        self.assertEqual(designed["maximum"], math.ceil(events))
        self.assertEqual(designed["unit"], "events")

    def test_two_means(self):
        designed = spec("continuous", "greater", design={**self.one, "assumptions": {"effect": 0.4, "sigma": 1.2}}
                        ).data["designed"]
        n = 4 * (Z(0.975) + Z(0.9)) ** 2 * 1.2**2 / 0.4**2
        self.assertEqual(designed["maximum"], math.ceil(n - 1e-9))

    def test_two_proportions(self):
        designed = spec("binary", "greater", design={**self.one, "assumptions": {"p_control": 0.3,
                                                                                 "p_treatment": 0.45}}
                        ).data["designed"]
        p = (0.3 + 0.45) / 2
        h0 = math.sqrt(p * (1 - p) * 4)
        h1 = math.sqrt(0.45 * 0.55 / 0.5 + 0.3 * 0.7 / 0.5)
        n = ((Z(0.975) * h0 + Z(0.9) * h1) / 0.15) ** 2
        self.assertEqual(designed["maximum"], math.ceil(n))

    def test_group_sequential_inflation(self):
        one = spec(design=self.one).data["designed"]["maximum"]
        three = spec(design={"analyses": [1 / 3, 2 / 3, 1.0], "power": 0.9}).data["designed"]
        self.assertGreater(three["maximum"], one)
        self.assertAlmostEqual(three["power"], 0.9, places=6)
        self.assertLess(three["expected_fraction"], 1.0)


class CheckTest(unittest.TestCase):
    def test_check_finds_edits(self):
        s = spec(design={"futility": "obrien-fleming"})
        self.assertEqual(check_trial(s), [])
        edited = s.to_json()
        edited["designed"]["bounds"][1] += 0.01
        edited["designed"]["maximum"] += 1
        del edited["kind"]
        found = check_trial(TrialSpec(edited))
        self.assertTrue(any(p.startswith("bounds[1]") for p in found))
        self.assertTrue(any(p.startswith("maximum") for p in found))

    def test_plan_report(self):
        report = plan_report(spec(design={"futility": "obrien-fleming"}), generated_at=0)
        text = report.to_markdown()
        self.assertIn("Planned analyses", text)
        self.assertEqual(report.findings()[0].guarantee, "asymptotic")
        self.assertEqual(len(report.sections[0].tables[0].rows), 3)


if __name__ == "__main__":
    unittest.main()
