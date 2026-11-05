"""Independent finite-slot checks of the compact general-hardware MIP.

Run: .venv/bin/python -B -m unittest tests.test_general_formulation -v
The oracle enumerates permutations and legal actions; it does not reuse the
formulation's constraints, routing variables, or lower bounds.
"""

from itertools import combinations
import random
import unittest

import gurobipy as gp
from gurobipy import GRB

from our_formulation.compact_formulation import build_model, extract_schedule
from our_formulation.validation import validate_schedule
from tests.oracle import exhaustive


class GeneralFormulationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.env = gp.Env(empty=True)
        cls.env.setParam("OutputFlag", 0)
        cls.env.start()

    @classmethod
    def tearDownClass(cls):
        cls.env.dispose()

    def optimize(self, gates, dag, edges, n, horizon, *, coverage=False,
                 initial=None, final=None, gate_times=None, check_lp=False):
        problem = build_model(gates, dag, edges, n, horizon,
                              env=self.env, coverage=coverage)
        model = problem.model
        try:
            for t, placement in ((0, initial), (horizon - 1, final)):
                if placement is not None:
                    for q in range(n):
                        for i in range(n):
                            problem.w[t, q, i].LB = problem.w[t, q, i].UB = int(placement[q] == i)
            if gate_times is not None:
                for (t, g), variable in problem.a.items():
                    variable.LB = variable.UB = int(gate_times[g] == t)
            model.Params.Threads = 1
            model.Params.TimeLimit = 10
            model.Params.MIPGap = model.Params.MIPGapAbs = 0
            model.update()
            lp_value = None
            if check_lp:
                relaxation = model.relax()
                try:
                    relaxation.optimize()
                    self.assertIn(relaxation.Status, (GRB.OPTIMAL, GRB.INFEASIBLE))
                    if relaxation.SolCount:
                        lp_value = relaxation.ObjVal
                finally:
                    relaxation.dispose()
            model.optimize()
            self.assertIn(model.Status, (GRB.OPTIMAL, GRB.INFEASIBLE))
            if not model.SolCount:
                return None
            objective = round(model.ObjVal)
            self.assertAlmostEqual(model.ObjVal, objective, places=5)
            self.assertTrue(all(v.VType == GRB.BINARY for v in problem.w.values()))
            self.assertTrue(all(v.VType == GRB.BINARY for v in problem.a.values()))
            for family in (problem.f, problem.x):
                self.assertTrue(all(v.VType == GRB.CONTINUOUS for v in family.values()))
                self.assertTrue(all(abs(v.X - round(v.X)) < 1e-5 for v in family.values()))
            schedule = extract_schedule(problem)
            validate_schedule(gates, dag, edges, n, horizon, schedule, objective)
            if initial is not None:
                self.assertEqual(schedule["placement"][0], list(initial))
            if check_lp:
                self.assertIsNotNone(lp_value)
                self.assertLessEqual(lp_value, objective + 1e-5)
            return objective
        finally:
            model.dispose()

    def test_random_general_graphs_against_independent_oracle(self):
        rng = random.Random(28112026)
        for case in range(40):
            n = rng.randint(2, 4)
            gates = [tuple(rng.sample(range(n), 2))
                     for _ in range(rng.randint(1, 5))]
            edges = [edge for edge in combinations(range(n), 2)
                     if rng.random() < .55]
            order = rng.sample(range(len(gates)), len(gates))
            dag = [(order[i], order[j]) for i in range(len(gates))
                   for j in range(i + 1, len(gates)) if rng.random() < .3]
            horizon = rng.randint(1, 6)
            expected = exhaustive(gates, dag, edges, n, horizon)
            for coverage in (False, True):
                with self.subTest(case=case, n=n, T=horizon, coverage=coverage):
                    self.assertEqual(self.optimize(gates, dag, edges, n, horizon,
                                                   coverage=coverage, check_lp=True), expected)

    def test_triangle_and_branching_star_allow_stationary_gates(self):
        cases = [
            ([(0, 1), (1, 2), (0, 2)], [(0, 1), (1, 2), (0, 2)], 3),
            ([(0, 1), (0, 2), (0, 3)], [(0, 1), (0, 2), (0, 3)], 4),
        ]
        for gates, edges, n in cases:
            with self.subTest(edges=edges):
                self.assertEqual(self.optimize(gates, [(0, 1), (1, 2)], edges,
                                               n, 3, coverage=True, check_lp=True), 0)

    def test_nonpath_graphs_with_positive_swap_optimum(self):
        gates = list(combinations(range(4), 2))
        dag = [(g, g + 1) for g in range(len(gates) - 1)]
        graphs = [
            [(0, 1), (1, 2), (2, 3), (3, 0)],
            [(0, 1), (0, 2), (0, 3)],
        ]
        for edges in graphs:
            expected = exhaustive(gates, dag, edges, 4, 9)
            self.assertIsNotNone(expected)
            self.assertGreater(expected, 0)
            for coverage in (False, True):
                with self.subTest(edges=edges, coverage=coverage):
                    self.assertEqual(self.optimize(gates, dag, edges, 4, 9,
                                                   coverage=coverage, check_lp=True), expected)

    def test_disconnected_graphs(self):
        edges = [(0, 1), (2, 3)]
        self.assertEqual(self.optimize([(0, 1), (2, 3)], [], edges, 4, 1), 0)
        # One logical qubit cannot meet two partners in two-site components.
        self.assertIsNone(self.optimize([(0, 1), (0, 2)], [(0, 1)], edges,
                                        4, 5, coverage=True))
        self.assertIsNone(self.optimize([(0, 1)], [], [], 2, 3, coverage=True))
        # An isolated unused site does not prevent gates in the other component.
        self.assertEqual(self.optimize([(0, 1), (1, 2)], [(0, 1)],
                                       [(0, 1), (1, 2), (0, 2)], 4, 2), 0)

    def test_one_and_two_slots_parallelism_and_strict_dependencies(self):
        gates, edges = [(0, 1), (2, 3)], [(0, 1), (1, 2), (2, 3), (3, 0)]
        self.assertEqual(self.optimize(gates, [], edges, 4, 1), 0)
        # Gate ids deliberately run opposite to their execution order.
        self.assertIsNone(self.optimize(gates, [(1, 0)], edges, 4, 1))
        self.assertEqual(self.optimize(gates, [(1, 0)], edges, 4, 2), 0)
        self.assertIsNone(self.optimize([(0, 1), (0, 2)], [], edges, 4, 1))

    def test_fixed_layout_and_nonconsecutive_edge_swap(self):
        # The triangle's (0,2) edge exchanges endpoints of the numerical order
        # in one SWAP. This reverses three pair orders, so tests path assumptions.
        edges = [(0, 1), (1, 2), (0, 2)]
        gates, dag = [(0, 1), (0, 1)], [(0, 1)]
        self.assertEqual(self.optimize(gates, dag, edges, 3, 4,
                                       initial=(0, 1, 2), final=(2, 1, 0),
                                       gate_times={0: 0, 1: 2}, check_lp=True), 1)
        # A fixed nonadjacent first gate cannot use a first-slot SWAP.
        self.assertIsNone(self.optimize([(0, 2)], [], [(0, 1), (1, 2)], 3, 4,
                                        initial=(0, 1, 2)))

    def test_idle_slots_are_compressed_by_extractor(self):
        # Deliberately leave an internal idle slot in a feasible MIP solution.
        self.assertEqual(self.optimize([(0, 1), (0, 1)], [(0, 1)],
                                       [(0, 1)], 2, 4,
                                       gate_times={0: 0, 1: 3}), 0)

    def test_empty_circuit_remains_infeasible(self):
        # Preserve legacy_model.py's compulsory first-slot gate convention.
        self.assertIsNone(self.optimize([], [], [(0, 1)], 2, 3))

    def test_invalid_inputs_are_rejected(self):
        bad = [
            ([(0, 0)], [], [(0, 1)], 2, 3),
            ([(0, 2)], [], [(0, 1)], 2, 3),
            ([(0, 1)], [], [(0, 0)], 2, 3),
            ([(0, 1)], [], [(0, 2)], 2, 3),
            ([(0, 1)], [], [(0, 1), (1, 0)], 2, 3),
            ([(0, 1)], [(0, 1)], [(0, 1)], 2, 3),
            ([(0, 1)], [(0, 0)], [(0, 1)], 2, 3),
            ([(0, 1), (1, 0)], [(0, 1), (1, 0)], [(0, 1)], 2, 3),
            ([(0, 1)], [], [(0, 1)], 2, 0),
            ([(0, 1)], [], [(0, 1)], 0, 3),
        ]
        for args in bad:
            with self.subTest(args=args):
                with self.assertRaises(ValueError):
                    problem = build_model(*args, env=self.env)
                    problem.model.dispose()


if __name__ == "__main__":
    unittest.main()
