# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Novikov Laboratories LLC
# Commercial licenses for use outside the AGPL: see COMMERCIAL.md
"""Simulations of the whole pipeline: the type I error and the power are those of the design."""
import unittest

from sequentialist.clinical import simulate_trial

from .helpers import spec


class SimulationTest(unittest.TestCase):
    # a small binary design: 2 analyses, Pocock spending, 30 % against 60 %, power 0.9 (about 120 patients)
    design = {"analyses": [0.5, 1.0], "spending": "pocock", "assumptions": {"p_control": 0.3, "p_treatment": 0.6}}

    def test_type_one_error_and_power(self):
        s = spec("binary", "greater", design=self.design)
        designed = s.data["designed"]
        null = simulate_trial(s, under="null", runs=400, seed=1)
        lo, hi = null.rejected.interval
        self.assertLessEqual(lo, 0.025)
        self.assertLess(null.rejected.rate, 0.05)
        alternative = simulate_trial(s, under="alternative", runs=200, seed=2)
        lo, hi = alternative.rejected.interval
        self.assertTrue(lo <= designed["power"] <= hi, (lo, designed["power"], hi))
        self.assertGreater(alternative.stopped_early.rate, null.stopped_early.rate)

    def test_simulations_are_reproducible(self):
        s = spec(design={"analyses": [0.5, 1.0], "maximum": 60, "power": None})
        a = simulate_trial(s, under="alternative", runs=5, seed=3).to_json()
        b = simulate_trial(s, under="alternative", runs=5, seed=3).to_json()
        self.assertEqual(a, b)

    def test_a_simulation_needs_a_schedule(self):
        s = spec(design={"type": "anytime", "analyses": None, "power": None, "assumptions": {}})
        with self.assertRaises(ValueError):
            simulate_trial(s, runs=1)


if __name__ == "__main__":
    unittest.main()
