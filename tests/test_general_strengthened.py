"""Certification and budget checks for the general-graph solve pipeline."""

import copy
from itertools import combinations
import random
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from gurobipy import GRB

from our_formulation import compact_formulation as general_formulation
from our_formulation.validation import validate_schedule
from tests.oracle import exhaustive


class FakeCountModel:
    """Expose solver result attributes while preserving the real routing MIP."""

    def __init__(self, *, bound=0, incumbent=1, status=GRB.TIME_LIMIT):
        self.Params = SimpleNamespace()
        self.Status = status
        self.ObjBound = bound
        self.ObjVal = incumbent
        self.SolCount = int(incumbent is not None)
        self.IsMIP = True
        self.Runtime = self.Work = self.NodeCount = self.IterCount = 0
        self.MIPGap = 1
        self.NumVars = self.NumBinVars = self.NumConstrs = self.NumNZs = 0
        self.optimized = self.disposed = False

    def optimize(self):
        self.optimized = True

    def dispose(self):
        self.disposed = True


class GeneralStrengthenedTests(unittest.TestCase):
    def setUp(self):
        self.gates = [(0, 1), (1, 2), (0, 2)]
        self.dag = [(0, 1), (1, 2)]
        self.edges = [(0, 1), (1, 2), (2, 3)]
        self.n, self.horizon = 4, 4
        self.witness = {
            "placement": [[0, 1, 2, 3], [0, 1, 2, 3], [0, 1, 2, 3], [1, 0, 2, 3]],
            "gate_times": {0: 0, 1: 1, 2: 3}, "swaps": [[2, 0]],
        }
        validate_schedule(self.gates, self.dag, self.edges, self.n, self.horizon,
                          self.witness, 1)

    def heuristic(self):
        return {"schedule": copy.deepcopy(self.witness), "objective": 1,
                "seconds": 0, "generated": 0}

    def solve(self, **options):
        return general_formulation.solve(self.gates, self.dag, self.edges,
                                         self.n, self.horizon, threads=1,
                                         time_limit=2, **options)

    def assert_valid(self, result, gates, dag, edges, n, horizon, expected):
        self.assertIn(result["status"], ("optimal", "infeasible", "time_limit"))
        if result["status"] == "optimal":
            self.assertTrue(result["proven_optimal"])
            self.assertEqual(result["objective"], expected)
            self.assertEqual(result["bound"], expected)
            self.assertEqual(result["gap"], 0)
        elif result["status"] == "infeasible":
            self.assertIsNone(expected)
            self.assertIsNone(result["objective"])
            self.assertIsNone(result["schedule"])
        if expected is not None:
            self.assertLessEqual(result["bound"], expected)
        if result["schedule"] is not None:
            self.assertTrue(result["schedule_validated"])
            self.assertIsNotNone(expected)
            self.assertGreaterEqual(result["objective"], expected)
            validate_schedule(gates, dag, edges, n, horizon,
                              result["schedule"], result["objective"])

    def test_count_incumbent_is_not_used_as_a_lower_bound(self):
        count = FakeCountModel(bound=0, incumbent=1)
        with patch("our_formulation.general_beam_start.construct", side_effect=lambda *a, **k: self.heuristic()), \
             patch("our_formulation.general_subset_bound.build_bound", return_value=(count, {})):
            result = self.solve()
        self.assertTrue(count.optimized)
        self.assertTrue(count.disposed)
        self.assertEqual(result["movement_bound"]["objective"], 1)
        self.assertEqual(result["movement_bound"]["bound"], 0)
        # A false use of the incumbent would have skipped the actual routing MIP.
        self.assertIsNotNone(result["routing"])
        self.assertEqual(result["proof_source"], "routing_mip")
        self.assert_valid(result, self.gates, self.dag, self.edges, 4, 4, 1)

    def test_time_limited_count_dual_bound_can_certify_a_valid_schedule(self):
        count = FakeCountModel(bound=1, incumbent=2)
        with patch("our_formulation.general_beam_start.construct", side_effect=lambda *a, **k: self.heuristic()), \
             patch("our_formulation.general_subset_bound.build_bound", return_value=(count, {})):
            result = self.solve()
        self.assertIsNone(result["routing"])
        self.assertEqual(result["proof_source"], "subset_movement_bound")
        self.assertEqual(result["movement_bound"]["status"], "time_limit")
        self.assert_valid(result, self.gates, self.dag, self.edges, 4, 4, 1)

    def test_expiry_retains_a_valid_heuristic_without_claiming_optimality(self):
        clock = [0.0]

        def construct(*args, **kwargs):
            clock[0] = 3.0  # Total budget is two seconds; no sleep is required.
            return self.heuristic()

        with patch("our_formulation.compact_formulation.time.perf_counter", side_effect=lambda: clock[0]), \
             patch("our_formulation.general_beam_start.construct", side_effect=construct), \
             patch("our_formulation.general_subset_bound.build_bound",
                   side_effect=AssertionError("No budget remains for count preprocessing")):
            result = self.solve()
        self.assertEqual(result["status"], "time_limit")
        self.assertFalse(result["proven_optimal"])
        self.assertEqual(result["objective"], 1)
        self.assertEqual(result["bound"], 0)
        self.assertIsNone(result["routing"])
        self.assertIsNone(result["movement_bound"])
        self.assert_valid(result, self.gates, self.dag, self.edges, 4, 4, 1)

    def test_missing_heuristic_falls_back_to_the_exact_routing_model(self):
        count = FakeCountModel(bound=0, incumbent=None)
        with patch("our_formulation.general_beam_start.construct",
                   return_value={"schedule": None, "objective": None, "seconds": 0}), \
             patch("our_formulation.general_subset_bound.build_bound", return_value=(count, {})):
            result = self.solve()
        self.assertIsNotNone(result["routing"])
        self.assertEqual(result["proof_source"], "routing_mip")
        self.assert_valid(result, self.gates, self.dag, self.edges, 4, 4, 1)

    def test_plain_mode_bypasses_both_strengthening_helpers(self):
        with patch("our_formulation.general_beam_start.construct", side_effect=AssertionError("Unexpected heuristic")), \
             patch("our_formulation.general_subset_bound.build_bound", side_effect=AssertionError("Unexpected bound")):
            result = self.solve(strengthen=False)
        self.assertIsNone(result["heuristic"])
        self.assertIsNone(result["movement_bound"])
        self.assertFalse(result["strengthen"])
        self.assert_valid(result, self.gates, self.dag, self.edges, 4, 4, 1)

    def test_zero_swap_schedule_needs_only_nonnegativity(self):
        with patch("our_formulation.general_subset_bound.build_bound",
                   side_effect=AssertionError("Zero SWAPs already proves optimality")):
            result = general_formulation.solve([(0, 1), (2, 3)], [],
                                               [(0, 1), (1, 2), (2, 3), (3, 0)],
                                               4, 1, threads=1, time_limit=1)
        self.assertEqual(result["proof_source"], "nonnegative_cost")
        self.assertIsNone(result["routing"])
        self.assertIsNone(result["movement_bound"])
        self.assertEqual(result["objective"], 0)

    def test_certificate_contradictions_fail_closed(self):
        for count in (FakeCountModel(bound=2, incumbent=2),
                      FakeCountModel(bound=float("inf"), incumbent=None, status=GRB.INFEASIBLE)):
            with self.subTest(status=count.Status), \
                 patch("our_formulation.general_beam_start.construct", side_effect=lambda *a, **k: self.heuristic()), \
                 patch("our_formulation.general_subset_bound.build_bound", return_value=(count, {})):
                with self.assertRaisesRegex(RuntimeError, "contradicts"):
                    self.solve()
            self.assertTrue(count.disposed)

    def test_infeasible_count_relaxation_is_a_valid_infeasibility_proof(self):
        count = FakeCountModel(bound=float("inf"), incumbent=None, status=GRB.INFEASIBLE)
        with patch("our_formulation.general_beam_start.construct",
                   return_value={"schedule": None, "objective": None, "seconds": 0}), \
             patch("our_formulation.general_subset_bound.build_bound", return_value=(count, {})):
            result = general_formulation.solve([(0, 1)], [], [], 4, 3,
                                               time_limit=1, threads=1)
        self.assertEqual(result["status"], "infeasible")
        self.assertEqual(result["proof_source"], "subset_movement_bound")
        self.assertIsNone(result["routing"])
        self.assertIsNone(result["schedule"])

    def test_small_general_instances_match_independent_oracle(self):
        rng = random.Random(31281026)
        cases = [(self.gates, self.dag, self.edges, 4, 4)]
        for _ in range(20):
            n = rng.randint(2, 4)
            gates = [tuple(rng.sample(range(n), 2)) for _ in range(rng.randint(1, 5))]
            ids = rng.sample(range(len(gates)), len(gates))
            dag = [(ids[i], ids[j]) for i in range(len(gates))
                   for j in range(i + 1, len(gates)) if rng.random() < .35]
            edges = [pair for pair in combinations(range(n), 2) if rng.random() < .55]
            cases.append((gates, dag, edges, n, rng.randint(1, 6)))
        for case, (gates, dag, edges, n, horizon) in enumerate(cases):
            expected = exhaustive(gates, dag, edges, n, horizon)
            result = general_formulation.solve(gates, dag, edges, n, horizon,
                                               time_limit=2, threads=1)
            with self.subTest(case=case):
                self.assertIn(result["status"], ("optimal", "infeasible"))
                self.assert_valid(result, gates, dag, edges, n, horizon, expected)

    def test_integer_bound_rounding_and_invalid_parameters(self):
        self.assertEqual(general_formulation._integer_bound(None), 0)
        self.assertEqual(general_formulation._integer_bound(0), 0)
        self.assertEqual(general_formulation._integer_bound(.1), 1)
        self.assertEqual(general_formulation._integer_bound(1 + 1e-8), 1)
        self.assertEqual(general_formulation._integer_bound(1.01), 2)
        for options in (dict(time_limit=0), dict(time_limit=float("nan")),
                        dict(threads=0), dict(seed=-1)):
            with self.assertRaises(ValueError):
                general_formulation.solve(self.gates, self.dag, self.edges, 4, 4, **options)


if __name__ == "__main__":
    unittest.main()
