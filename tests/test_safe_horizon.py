"""Independent tests of automatic horizons and global optimality reporting."""

import copy
from heapq import heappop, heappush
from itertools import combinations, permutations
import random
import unittest
from unittest.mock import patch

import gurobipy as gp

from our_formulation import compact_formulation as formulation
from our_formulation.general_subset_bound import build_bound
from our_formulation.serial_heuristic import SafeHorizonUnavailable, construct
from our_formulation.validation import validate_schedule
from tests.test_general_strengthened import FakeCountModel


def unrestricted_oracle(gates, dag, edges, n):
    """Tiny full-state shortest-path search with no slot deadline.

    This test-only oracle does not call production routing, horizon selection,
    bounds, or MIP builders. Zero-cost gate transitions and unit-cost physical
    transpositions directly define the unrestricted minimum SWAP count.
    """
    distance = {(p, 0): 0 for p in permutations(range(n))}
    queue = [(0, p, 0) for p, _ in distance]
    hardware = {frozenset(edge) for edge in edges}
    finished = (1 << len(gates)) - 1
    while queue:
        cost, placement, done = heappop(queue)
        if cost != distance[placement, done]:
            continue
        if done == finished:
            return cost
        actions = []
        for g, (q, r) in enumerate(gates):
            if (not done >> g & 1 and all(done >> u & 1 for u, v in dag if v == g)
                    and frozenset((placement[q], placement[r])) in hardware):
                actions.append((cost, placement, done | (1 << g)))
        for i, j in edges:
            after = tuple(j if pos == i else i if pos == j else pos for pos in placement)
            actions.append((cost + 1, after, done))
        for new_cost, after, new_done in actions:
            if new_cost < distance.get((after, new_done), float("inf")):
                distance[after, new_done] = new_cost
                heappush(queue, (new_cost, after, new_done))
    return None


class SafeHorizonTests(unittest.TestCase):
    def setUp(self):
        self.gates = [(0, 1), (1, 2), (0, 2)]
        self.dag = [(0, 1), (1, 2)]
        self.edges = [(0, 1), (1, 2), (2, 3)]

    def test_serial_witness_absorbs_leading_swaps_and_is_reproducible(self):
        gates, dag = [(0, 3), (0, 1), (2, 3), (0, 2)], [(0, 2), (1, 2), (2, 3)]
        first = construct(gates, dag, self.edges, 4)
        second = construct(gates, dag, self.edges, 4)
        self.assertEqual(first["schedule"], second["schedule"])
        self.assertEqual(first["horizon"], len(gates) + first["objective"])
        self.assertEqual(first["schedule"]["gate_times"][0], 0)
        self.assertNotEqual(first["schedule"]["placement"][0], list(range(4)))
        self.assertTrue(all(0 < t < first["horizon"] - 1 for t, _ in first["schedule"]["swaps"]))
        validate_schedule(gates, dag, self.edges, 4, first["horizon"],
                          first["schedule"], first["objective"])

    def test_disconnected_hardware_embeds_interaction_components(self):
        # Identity placement puts (0,3) in distinct hardware components.
        gates, dag = [(0, 3), (1, 2), (0, 3)], [(0, 2)]
        edges = [(0, 1), (2, 3)]
        witness = construct(gates, dag, edges, 5)
        self.assertEqual(witness["objective"], 0)
        validate_schedule(gates, dag, edges, 5, witness["horizon"], witness["schedule"], 0)
        result = formulation.solve(gates, dag, edges, 5, threads=1, time_limit=2)
        self.assertTrue(result["globally_optimal"])

    def test_failure_to_construct_safe_witness_is_not_an_infeasibility_result(self):
        with patch.object(formulation, "build_model", side_effect=AssertionError("No safe horizon")):
            with self.assertRaisesRegex(SafeHorizonUnavailable, "Automatic horizon unavailable"):
                formulation.solve(self.gates, self.dag, [(0, 1), (2, 3)], 4)
            with self.assertRaisesRegex(SafeHorizonUnavailable, "nonempty"):
                formulation.solve([], [], self.edges, 4)

    def test_auto_horizon_matches_independent_unrestricted_oracle(self):
        rng = random.Random(20261104)
        cases = [(self.gates, self.dag, self.edges, 4)]
        for _ in range(16):
            n = rng.randint(2, 4)
            edges = [(q, q + 1) for q in range(n - 1)]
            edges += [e for e in combinations(range(n), 2) if e not in edges and rng.random() < .35]
            gates = [tuple(rng.sample(range(n), 2)) for _ in range(rng.randint(1, 5))]
            ids = rng.sample(range(len(gates)), len(gates))
            dag = [(ids[i], ids[j]) for i in range(len(gates))
                   for j in range(i + 1, len(gates)) if rng.random() < .4]
            cases.append((gates, dag, edges, n))
        for index, (gates, dag, edges, n) in enumerate(cases):
            expected = unrestricted_oracle(gates, dag, edges, n)
            for strengthen in (False, True):
                with self.subTest(case=index, strengthen=strengthen):
                    result = formulation.solve(gates, dag, edges, n, time_limit=2,
                                               threads=1, strengthen=strengthen)
                    self.assertEqual(result["status"], "optimal")
                    self.assertEqual(result["objective"], expected)
                    self.assertEqual(result["horizon"], len(gates) + result["heuristic_swap_count"])
                    self.assertEqual(result["optimality_scope"], "global")
                    self.assertTrue(result["globally_optimal"])
                    self.assertEqual(result["horizon_source"], "serial_heuristic")
                    validate_schedule(gates, dag, edges, n, result["horizon"],
                                      result["schedule"], result["objective"])

    def test_explicit_deadline_never_claims_unrestricted_optimality(self):
        for horizon, expected in ((3, "infeasible"), (4, "optimal")):
            result = formulation.solve(self.gates, self.dag, self.edges, 4, horizon,
                                       time_limit=2, threads=1)
            self.assertEqual(result["status"], expected)
            self.assertEqual(result["optimality_scope"], "fixed_horizon")
            self.assertEqual(result["horizon_source"], "explicit")
            self.assertFalse(result["globally_optimal"])
            self.assertIsNone(result["heuristic_swap_count"])

    def test_invalid_heuristic_cannot_supply_safe_horizon_or_certificate(self):
        witness = construct(self.gates, self.dag, self.edges, 4)
        invalid = copy.deepcopy(witness)
        invalid["schedule"]["gate_times"].pop(0)
        with patch("our_formulation.serial_heuristic.construct", return_value=invalid), \
             patch.object(formulation, "build_model", side_effect=AssertionError("Must validate first")):
            with self.assertRaisesRegex(ValueError, "Missing gate"):
                formulation.solve(self.gates, self.dag, self.edges, 4)
        invalid = copy.deepcopy(witness)
        invalid["horizon"] += 1
        with patch("our_formulation.serial_heuristic.construct", return_value=invalid):
            with self.assertRaisesRegex(ValueError, "Safe horizon disagrees"):
                formulation.solve(self.gates, self.dag, self.edges, 4)

    def test_failed_beam_search_retains_mandatory_witness(self):
        with patch("our_formulation.general_beam_start.construct",
                   return_value={"schedule": None, "objective": None, "seconds": 0}):
            result = formulation.solve(self.gates, self.dag, self.edges, 4, time_limit=2, threads=1)
        self.assertTrue(result["globally_optimal"])
        self.assertEqual(result["objective"], 1)
        self.assertIsNone(result["heuristic"]["objective"])

    def test_horizon_preprocessing_consumes_total_budget(self):
        witness = construct(self.gates, self.dag, self.edges, 4)
        clock = [0.0]

        def slow_construct(*args, **kwargs):
            self.assertEqual(kwargs["deadline"], 2)
            clock[0] = 3.0
            return witness

        with patch("our_formulation.compact_formulation.time.perf_counter", side_effect=lambda: clock[0]), \
             patch("our_formulation.serial_heuristic.construct", side_effect=slow_construct), \
             patch("our_formulation.general_beam_start.construct", side_effect=AssertionError("Budget spent")), \
             patch("our_formulation.general_subset_bound.build_bound", side_effect=AssertionError("Budget spent")):
            result = formulation.solve(self.gates, self.dag, self.edges, 4, time_limit=2, threads=1)
        self.assertEqual(result["status"], "time_limit")
        self.assertFalse(result["globally_optimal"])
        self.assertIsNone(result["routing"])
        self.assertEqual(result["production_seconds"], 3)
        self.assertTrue(result["schedule_validated"])

    def test_symmetry_defaults_and_overrides_reach_both_mips(self):
        with gp.Env(empty=True) as env:
            env.setParam("OutputFlag", 0)
            env.start()
            problem = formulation.build_model(self.gates, self.dag, self.edges, 4, 4, env=env)
            try:
                self.assertEqual(problem.model.Params.Symmetry, 2)
            finally:
                problem.model.dispose()
            bound, _ = build_bound(self.gates, self.dag, self.edges, 4, 4,
                                   size=2, deadline=0, env=env)
            try:
                self.assertEqual(bound.Params.Symmetry, 2)
            finally:
                bound.dispose()
        for symmetry in (-1, 0, 1, 2):
            count = FakeCountModel(bound=1, incumbent=1)
            with patch.object(formulation, "build_model", wraps=formulation.build_model) as build, \
                 patch("our_formulation.general_subset_bound.build_bound", return_value=(count, {})):
                result = formulation.solve(self.gates, self.dag, self.edges, 4,
                                           time_limit=2, threads=1, symmetry=symmetry)
            self.assertEqual(build.call_args.kwargs["symmetry"], symmetry)
            self.assertEqual(count.Params.Symmetry, symmetry)
            self.assertEqual(result["symmetry"], symmetry)
            self.assertTrue(result["globally_optimal"])
        for invalid in (-2, 3, True, .5):
            with self.assertRaisesRegex(ValueError, "symmetry"):
                formulation.solve(self.gates, self.dag, self.edges, 4, symmetry=invalid)


if __name__ == "__main__":
    unittest.main()
