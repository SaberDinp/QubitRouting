"""Independent validity checks for general-hardware subset movement bounds."""

from itertools import combinations, permutations, product
import random
import time
import unittest
from unittest.mock import patch

import gurobipy as gp
from gurobipy import GRB

import our_formulation.general_subset_bound as subset_bound
from our_formulation.compact_formulation import build_model, extract_schedule
from our_formulation.validation import validate_schedule
from tests.oracle import exhaustive


def alternating_witness(n, edges, seed):
    """Generate a valid route directly, independently of any routing solver."""
    rng = random.Random(seed)
    current = rng.sample(range(n), n)
    placements, gates, gate_times, swaps = [], [], {}, []
    for t in range(9):
        placements.append(current[:])
        e = rng.randrange(len(edges))
        i, j = edges[e]
        if t % 2 == 0:
            gate_times[len(gates)] = t
            gates.append((current.index(i), current.index(j)))
        else:
            swaps.append([t, e])
            current = [j if node == i else i if node == j else node for node in current]
    dag = [(g, g + 1) for g in range(len(gates) - 1)]
    schedule = dict(placement=placements, gate_times=gate_times, swaps=swaps)
    validate_schedule(gates, dag, edges, n, 9, schedule, len(swaps))
    return gates, dag, schedule


def movement_counts(schedule, edges, n):
    counts = [0] * n
    for t, e in schedule["swaps"]:
        i, j = edges[e]
        counts[schedule["placement"][t].index(i)] += 1
        counts[schedule["placement"][t].index(j)] += 1
    return counts


class GeneralSubsetBoundTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.env = gp.Env(empty=True)
        cls.env.setParam("OutputFlag", 0)
        cls.env.start()

    @classmethod
    def tearDownClass(cls):
        cls.env.dispose()

    def check_witness(self, gates, dag, edges, n, horizon, schedule, selections,
                      *, allow_one_anonymous=False):
        objective = len(schedule["swaps"])
        validate_schedule(gates, dag, edges, n, horizon, schedule, objective)
        counts = movement_counts(schedule, edges, n)
        pred = subset_bound.ancestors(gates, dag)
        checked = 0
        for selected in selections:
            k = len(selected)
            patterns = [(1,) * k, (2,) * k, tuple(0 if j == 0 else 2 if j == 1 else 1
                                                 for j in range(k))]
            for weights in patterns:
                costs, info = subset_bound.conditional_costs(
                    gates, pred, edges, n, selected, weights=weights,
                    allow_one_anonymous=allow_one_anonymous)
                initial = tuple(schedule["placement"][0][q] for q in selected)
                self.assertIn(initial, costs)
                actual = sum(weights[j] * counts[q] for j, q in enumerate(selected))
                self.assertLessEqual(costs[initial], actual)
                self.assertLess(len(selected), n)
                self.assertEqual(info["anonymous_qubits"], n - k)
                checked += 1
        return checked

    def test_weighted_inequalities_on_constructed_general_graph_routes(self):
        graphs = [
            (5, [(0, 1), (1, 2), (2, 3), (3, 4)]),
            (5, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]),
            (5, [(0, 1), (0, 2), (0, 3), (0, 4)]),
            (5, [(0, 1), (1, 2), (3, 4)]),
            (6, [(0, 1), (1, 2), (3, 4), (4, 5), (0, 3), (1, 4), (2, 5)]),
        ]
        checked = 0
        for index, (n, edges) in enumerate(graphs):
            gates, dag, schedule = alternating_witness(n, edges, 4182 + index)
            selections = list(combinations(range(n), 3))
            with self.subTest(edges=edges):
                checked += self.check_witness(gates, dag, edges, n, 9, schedule, selections)
        self.assertEqual(checked, 180)

    def test_weighted_inequalities_on_cycle_and_grid_routes(self):
        graphs = [
            [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)],
            [(0, 1), (1, 2), (3, 4), (4, 5), (0, 3), (1, 4), (2, 5)],
        ]
        # Build reproducible witnesses without private experiment results.
        for graph_index, edges in enumerate(graphs):
            for seed in (712, 917):
                gates, dag, schedule = alternating_witness(6, edges, seed)
                with self.subTest(graph=graph_index, seed=seed):
                    self.check_witness(gates, dag, edges, 6, 9, schedule,
                                       [(0, 1, 2), (1, 3, 5), (0, 2, 4, 5)])

    def test_inherited_precedence_through_deleted_gates(self):
        gates = [(0, 1), (0, 3), (1, 2), (0, 4), (0, 2)]
        dag = [(g, g + 1) for g in range(4)]
        edges = [(0, 1), (1, 2), (2, 3), (3, 4)]
        inherited, _ = subset_bound.conditional_costs(
            gates, subset_bound.ancestors(gates, dag), edges, 5, (0, 1, 2))
        retained = [gates[g] for g in (0, 2, 4)]
        direct, _ = subset_bound.conditional_costs(
            retained, subset_bound.ancestors(retained, [(0, 1), (1, 2)]),
            edges, 5, (0, 1, 2))
        self.assertEqual(inherited, direct)
        # Reverse labels, including dependencies, without changing the circuit.
        reverse = list(reversed(gates))
        reverse_dag = [(4 - g, 4 - h) for g, h in dag]
        reversed_costs, _ = subset_bound.conditional_costs(
            reverse, subset_bound.ancestors(reverse, reverse_dag), edges, 5, (0, 1, 2))
        self.assertEqual(inherited, reversed_costs)

    def test_disconnected_graph_omits_only_unreachable_initial_positions(self):
        edges = [(0, 1), (2, 3)]
        gates = [(0, 1)]
        costs, _ = subset_bound.conditional_costs(gates, [0], edges, 5, (0, 1, 2))
        expected = {pos for pos in permutations(range(5), 3)
                    if frozenset(pos[:2]) in {frozenset(e) for e in edges}}
        self.assertEqual(set(costs), expected)
        self.assertTrue(all(value == 0 for value in costs.values()))
        triangle = [(0, 1), (1, 2), (0, 2)]
        impossible, _ = subset_bound.conditional_costs(
            triangle, subset_bound.ancestors(triangle, [(0, 1), (1, 2)]),
            edges, 5, (0, 1, 2))
        self.assertEqual(impossible, {})
        model, _ = subset_bound.build_bound(triangle, [(0, 1), (1, 2)], edges,
                                             5, 8, env=self.env, size=3, weighted=True)
        try:
            model.Params.Threads = 1
            model.optimize()
            self.assertEqual(model.Status, GRB.INFEASIBLE)
        finally:
            model.dispose()

    def test_weight_scaling_and_zero_cost_moves(self):
        gates = [(0, 1), (1, 2), (0, 2)]
        pred = subset_bound.ancestors(gates, [(0, 1), (1, 2)])
        edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]
        tables = {}
        for weights in ((1, 1, 1), (2, 2, 2), (0, 0, 0), (0, 1, 0),
                        (10 ** 9, 10 ** 9, 10 ** 9)):
            tables[weights], _ = subset_bound.conditional_costs(
                gates, pred, edges, 5, (0, 1, 2), weights=weights)
        for pos, value in tables[(1, 1, 1)].items():
            self.assertEqual(tables[(2, 2, 2)][pos], 2 * value)
            self.assertEqual(tables[(10 ** 9, 10 ** 9, 10 ** 9)][pos], 10 ** 9 * value)
            self.assertEqual(tables[(0, 0, 0)][pos], 0)
            self.assertLessEqual(tables[(0, 1, 0)][pos], value)
        self.assertEqual(tables[(0, 1, 0)][0, 1, 2], 0)

    def test_all_51_pair_weight_patterns_on_a_valid_general_graph_witness(self):
        n = 6
        edges = [(0, 1), (1, 2), (3, 4), (4, 5), (0, 3), (1, 4), (2, 5)]
        gates, dag, schedule = alternating_witness(n, edges, 781026)
        counts = movement_counts(schedule, edges, n)
        selected = (0, 1, 2, 3, 4)
        # Independent definition: zero, one, or two coordinates differ from 1.
        patterns = {weights for weights in product((0, 1, 2), repeat=5)
                    if sum(w != 1 for w in weights) <= 2}
        self.assertEqual(len(patterns), 51)
        graph = subset_bound._projection_graph(
            gates, subset_bound.ancestors(gates, dag), edges, n, selected,
            allow_one_anonymous=True)
        initial = tuple(schedule["placement"][0][q] for q in selected)
        for weights in patterns:
            costs, _ = subset_bound._costs_on_graph(graph, weights)
            self.assertIn(initial, costs)
            self.assertLessEqual(costs[initial], sum(weights[j] * counts[q]
                                                    for j, q in enumerate(selected)))

        # Fix the complete feasible witness in the count MIP. Every generated
        # inequality (all six subsets and all 51 patterns) must accept it.
        model, info = subset_bound.build_bound(
            gates, dag, edges, n, 9, size=5, pair_weights=True,
            allow_one_anonymous=True, env=self.env)
        try:
            reports = info["projection_stats"]
            self.assertEqual(len(reports), 6 * 51)
            for selection in combinations(range(n), 5):
                self.assertEqual({tuple(row["weights"]) for row in reports
                                  if tuple(row["subset"]) == selection}, patterns)
            for q in range(n):
                for i in range(n):
                    variable = model.getVarByName(f"initial[{q},{i}]")
                    variable.LB = variable.UB = int(schedule["placement"][0][q] == i)
                variable = model.getVarByName(f"movement_count[{q}]")
                variable.LB = variable.UB = counts[q]
            variable = model.getVarByName("swap_count")
            variable.LB = variable.UB = len(schedule["swaps"])
            model.Params.Threads = 1
            model.optimize()
            self.assertEqual(model.Status, GRB.OPTIMAL)
        finally:
            model.dispose()

    def test_one_anonymous_qubit_requires_opt_in_and_full_subset_is_rejected(self):
        gates, dag = [(0, 1), (1, 2), (0, 2)], [(0, 1), (1, 2)]
        pred = subset_bound.ancestors(gates, dag)
        edges = [(0, 1), (1, 2), (2, 3), (3, 0)]
        with self.assertRaises(ValueError):
            subset_bound.conditional_costs(gates, pred, edges, 4, (0, 1, 2))
        costs, info = subset_bound.conditional_costs(
            gates, pred, edges, 4, (0, 1, 2), allow_one_anonymous=True)
        self.assertEqual(info["anonymous_qubits"], 1)
        self.assertEqual(info["placement_count"], 24)
        # Qubit 0 swaps with the anonymous token at site 3 before gate (0,2).
        self.assertEqual(costs[0, 1, 2], 1)
        for opt_in in (False, True):
            with self.assertRaises(ValueError):
                subset_bound.conditional_costs(gates, pred, edges, 4, (0, 1, 2, 3),
                                               allow_one_anonymous=opt_in)

    def test_default_projection_does_not_enumerate_complete_placements(self):
        real_permutations = subset_bound.permutations
        sizes = []

        def guarded(values, size):
            values = tuple(values)
            self.assertLessEqual(size, len(values) - 2)
            sizes.append(size)
            return real_permutations(values, size)

        with patch.object(subset_bound, "permutations", side_effect=guarded):
            model, info = subset_bound.build_bound(
                [(0, 1), (1, 2), (0, 2)], [(0, 1), (1, 2)],
                [(0, 1), (1, 2), (2, 3), (3, 4)], 5, 5, size=3,
                weighted=True, env=self.env)
            model.dispose()
        self.assertTrue(sizes)
        self.assertTrue(all(item["anonymous_qubits"] >= 2 for item in info["projection_stats"]))

    def test_expired_preprocessing_returns_a_valid_weak_count_master(self):
        gates, dag = [(0, 1), (1, 2), (0, 2)], [(0, 1), (1, 2)]
        edges = [(0, 1), (1, 2), (2, 3), (3, 4)]
        deadline = time.perf_counter() - 1
        with self.assertRaises(TimeoutError):
            subset_bound.conditional_costs(gates, subset_bound.ancestors(gates, dag),
                                           edges, 5, (0, 1, 2), deadline=deadline)
        model, info = subset_bound.build_bound(gates, dag, edges, 5, 5, env=self.env,
                                                 size=3, pair_weights=True, deadline=deadline)
        try:
            self.assertFalse(info["precompute_complete"])
            self.assertTrue(info["precompute_truncated"])
            self.assertEqual(info["completed_projection_tables"], 0)
            self.assertEqual(info["projection_stats"], [])
            self.assertIsNotNone(info["precompute_reason"])
            model.Params.Threads = 1
            model.optimize()
            self.assertEqual(model.Status, GRB.OPTIMAL)
            self.assertEqual(model.ObjVal, 0)
        finally:
            model.dispose()

    def test_partial_preprocessing_keeps_only_finished_weight_tables(self):
        gates, dag = [(0, 1), (1, 2), (0, 2)], [(0, 1), (1, 2)]
        edges = [(0, 1), (1, 2), (2, 3), (3, 4)]
        real_costs = subset_bound._costs_on_graph
        for finish in (0, 2):
            calls = 0

            def interrupt(graph, weights, **kwargs):
                nonlocal calls
                calls += 1
                if calls > finish:
                    raise TimeoutError("Deliberate test interruption")
                return real_costs(graph, weights, **kwargs)

            with patch.object(subset_bound, "_costs_on_graph", side_effect=interrupt):
                model, info = subset_bound.build_bound(
                    gates, dag, edges, 5, 5, env=self.env, size=3, pair_weights=True)
            try:
                self.assertFalse(info["precompute_complete"])
                self.assertTrue(info["precompute_truncated"])
                self.assertEqual(info["completed_projection_tables"], finish)
                self.assertEqual(len(info["projection_stats"]), finish)
                if finish:
                    self.assertEqual(tuple(info["projection_stats"][0]["weights"]), (1, 1, 1))
                model.Params.Threads = 1
                model.optimize()
                self.assertEqual(model.Status, GRB.OPTIMAL)
                self.assertLessEqual(model.ObjVal, 1)
            finally:
                model.dispose()

        # Interrupt during Dijkstra, rather than between tables. Partial
        # distances must raise and must never masquerade as completed costs.
        graph = subset_bound._projection_graph([], [], [(0, 1), (1, 2), (2, 3),
                                                        (3, 4), (4, 5)], 6, (0, 1, 2, 3))
        with patch.object(subset_bound, "_check_deadline",
                          side_effect=[None, None, TimeoutError("During Dijkstra")]):
            with self.assertRaises(TimeoutError):
                subset_bound._costs_on_graph(graph, (1, 1, 1, 1))

    def test_resource_limits_skip_projections_without_claiming_infeasibility(self):
        gates, dag = [(0, 1), (1, 2), (0, 2)], [(0, 1), (1, 2)]
        edges = [(0, 1), (1, 2), (2, 3), (3, 4)]
        for limits in (dict(state_limit=1), dict(arc_limit=1)):
            with self.subTest(limits=limits):
                model, info = subset_bound.build_bound(
                    gates, dag, edges, 5, 5, env=self.env, size=3,
                    pair_weights=True, **limits)
                try:
                    self.assertFalse(info["precompute_complete"])
                    self.assertFalse(info["precompute_truncated"])
                    self.assertEqual(info["completed_projection_tables"], 0)
                    self.assertEqual(len(info["skipped_projections"]), 10)
                    model.Params.Threads = 1
                    model.optimize()
                    self.assertEqual(model.Status, GRB.OPTIMAL)
                    self.assertEqual(model.ObjVal, 0)
                finally:
                    model.dispose()

    def test_count_bounds_against_small_independent_full_problem_oracle(self):
        rng = random.Random(67281026)
        cases = [
            ([(0, 1), (1, 2), (0, 2)], [(0, 1), (1, 2)],
             [(0, 1), (1, 2), (2, 3), (3, 4)], 5),
            ([(0, 1), (1, 2), (0, 2)] * 2, [(g, g + 1) for g in range(5)],
             [(0, 1), (1, 2), (2, 3), (3, 4)], 8),
        ]
        for _ in range(8):
            gates = [tuple(rng.sample(range(5), 2)) for _ in range(rng.randint(2, 4))]
            labels = rng.sample(range(len(gates)), len(gates))
            dag = [(labels[i], labels[j]) for i in range(len(gates))
                   for j in range(i + 1, len(gates)) if rng.random() < .35]
            edges = [pair for pair in combinations(range(5), 2) if rng.random() < .5]
            cases.append((gates, dag, edges, rng.randint(3, 6)))
        positive = 0
        for case, (gates, dag, edges, horizon) in enumerate(cases):
            expected = exhaustive(gates, dag, edges, 5, horizon)
            model, _ = subset_bound.build_bound(gates, dag, edges, 5, horizon,
                                                 size=3, weighted=True, env=self.env)
            try:
                model.Params.Threads = 1
                model.Params.TimeLimit = 5
                model.Params.MIPGap = model.Params.MIPGapAbs = 0
                model.optimize()
                with self.subTest(case=case):
                    self.assertIn(model.Status, (GRB.OPTIMAL, GRB.INFEASIBLE))
                    if expected is not None:
                        self.assertEqual(model.Status, GRB.OPTIMAL)
                        self.assertLessEqual(model.ObjVal, expected + 1e-5)
                        self.assertLessEqual(model.ObjBound, expected + 1e-5)
                        positive += expected > 0
            finally:
                model.dispose()
        self.assertGreaterEqual(positive, 2)

    def test_routing_integration_preserves_oracle_optimum(self):
        gates, dag = [(0, 1), (1, 2), (0, 2)], [(0, 1), (1, 2)]
        for edges in ([(0, 1), (1, 2), (2, 3), (3, 4)],
                      [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)],
                      [(0, 1), (0, 2), (0, 3), (0, 4)]):
            expected = exhaustive(gates, dag, edges, 5, 5)
            problem = build_model(gates, dag, edges, 5, 5, env=self.env)
            try:
                subset_bound.add_to_routing(problem, size=3, weighted=True)
                problem.model.Params.Threads = 1
                problem.model.Params.TimeLimit = 5
                problem.model.Params.MIPGap = problem.model.Params.MIPGapAbs = 0
                problem.model.optimize()
                self.assertEqual(problem.model.Status, GRB.OPTIMAL)
                self.assertEqual(round(problem.model.ObjVal), expected)
                schedule = extract_schedule(problem)
                validate_schedule(gates, dag, edges, 5, 5, schedule, expected)
            finally:
                problem.model.dispose()

    def test_preprocessing_limits_and_invalid_weights(self):
        gates = [(0, 1), (1, 2), (0, 2)]
        pred = subset_bound.ancestors(gates, [(0, 1), (1, 2)])
        edges = [(0, 1), (1, 2), (2, 3), (3, 4)]
        for options in (dict(weights=(1, 2)), dict(weights=(1, -1, 1)),
                        dict(weights=(1, .5, 1)), dict(ideal_limit=1), dict(state_limit=1)):
            with self.subTest(options=options), self.assertRaises(ValueError):
                subset_bound.conditional_costs(gates, pred, edges, 5, (0, 1, 2), **options)
        for min_size in (1, 4):
            with self.assertRaises(ValueError):
                subset_bound.build_bound(gates, [(0, 1), (1, 2)], edges, 5, 5,
                                           size=3, min_size=min_size, env=self.env)


if __name__ == "__main__":
    unittest.main()
