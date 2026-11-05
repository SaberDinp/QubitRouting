"""Independent exhaustive checks. Standard library only unless explicitly opted in.

    .venv/bin/python -m unittest luteberget.test_solver -v
    QMRS_TEST_GUROBI=1 .venv/bin/python -m unittest luteberget.test_solver -v
"""

from functools import lru_cache
from itertools import permutations
import json
import os
from pathlib import Path
import random
import tempfile
import unittest

from load_instance import load_instance
from luteberget import solve
from luteberget.validation import validate_solution


ROOT = Path(__file__).resolve().parents[1]


def exhaustive(gates, dag, edges, n, horizon, *, qmrs=True, sw=1, dw=0):
    """Enumerate full placements and ALL legal time-slot actions, including idle.

    No lazy assignment, paper bounds, Pareto pruning, or maximal-subset rule.
    Returns the minimum weighted objective for unit-duration operations.
    """
    full = (1 << len(gates)) - 1
    pred = [0] * len(gates)
    for u, v in dag:
        pred[v] |= 1 << u
    matchings = [((), frozenset())]
    for e, edge in enumerate(edges):
        matchings += [(es + (e,), used | set(edge)) for es, used in matchings
                      if not used.intersection(edge)]
    edge_set = {frozenset(e) for e in edges}

    @lru_cache(None)
    def actions(placement, done, swaps_allowed):
        positions = [placement.index(q) for q in range(n)]
        ready = [g for g, (p, q) in enumerate(gates)
                 if not (done >> g & 1) and pred[g] & done == pred[g]
                 and frozenset((positions[p], positions[q])) in edge_set]
        result = []
        for swaps, used in matchings if swaps_allowed else matchings[:1]:
            subsets = [(0, used)]
            for g in ready:
                p, q = gates[g]
                sites = {positions[p], positions[q]}
                subsets += [(mask | (1 << g), busy | sites) for mask, busy in subsets
                            if not busy.intersection(sites)]
            nxt = list(placement)
            for e in swaps:
                u, v = edges[e]
                nxt[u], nxt[v] = nxt[v], nxt[u]
            result.extend((tuple(nxt), done | mask, len(swaps), mask) for mask, _ in subsets)
        return result

    states = {(p, 0): 0 for p in permutations(range(n))}
    best = None
    for t in range(horizon):
        following = {}
        for (placement, done), cost in states.items():
            for nxt, new_done, extra, mask in actions(placement, done, not (qmrs and t == 0)):
                if qmrs and t == 0 and not mask:
                    continue
                value = cost + extra
                if new_done == full:
                    objective = sw * value + dw * (t + 1)
                    best = objective if best is None else min(best, objective)
                elif value < following.get((nxt, new_done), float("inf")):
                    following[nxt, new_done] = value
        states = following
    return best


class SolverTests(unittest.TestCase):
    def test_random_instances_against_exhaustive_schedules(self):
        rng = random.Random(250821718)
        for case in range(36):
            n = rng.choice((2, 3, 4))
            gates = [tuple(rng.sample(range(n), 2)) for _ in range(rng.randint(1, 5))]
            # Includes disjoint-pair dependencies and unordered gates on the
            # same qubit, unlike a conventional circuit DAG.
            dag = [(u, v) for u in range(len(gates)) for v in range(u + 1, len(gates))
                   if rng.random() < .3]
            edges = [(i, i + 1) for i in range(n - 1)]
            if n == 4 and case % 2:
                edges = [(0, 1), (0, 2), (0, 3)]
            if n > 2 and case % 5 == 0 and (0, n - 1) not in edges:
                edges.append((0, n - 1))
            horizon = rng.randint(2, 6)
            for mode, sw, dw in (("qmrs", 1, 0), ("paper", 0, 1), ("paper", 2, 1)):
                with self.subTest(case=case, mode=mode, weights=(sw, dw)):
                    expected = exhaustive(gates, dag, edges, n, horizon,
                                          qmrs=mode == "qmrs", sw=sw, dw=dw)
                    result = solve(gates, dag, edges, n, horizon=horizon, mode=mode,
                                   swap_weight=sw, depth_weight=dw, time_limit=10)
                    self.assertTrue(result["proven"], result)
                    self.assertEqual(result["objective"], expected)

    def test_requested_benchmark(self):
        gates, dag, edges, n, _ = load_instance(ROOT / "benchmark/qv6_depth6_01_opt7.txt")
        result = solve(gates, dag, edges, n, horizon=12, time_limit=30)
        self.assertEqual(result["status"], "optimal")
        self.assertEqual(result["objective"], 7)
        self.assertEqual(result["bound"], 7)
        self.assertLessEqual(result["makespan"], 12)

    def test_horizon_changes_feasibility(self):
        gates = [(0, 1), (1, 2), (0, 2)]
        dag = [(0, 1), (1, 2)]
        edges = [(0, 1), (1, 2)]
        impossible = solve(gates, dag, edges, 3, horizon=3)
        self.assertEqual(impossible["status"], "infeasible")
        possible = solve(gates, dag, edges, 3, horizon=4)
        self.assertEqual(possible["objective"], 1)

    def test_variable_durations_and_unlimited_horizon(self):
        gates = [(0, 1), (1, 2), (0, 2)]
        dag = [(0, 1), (1, 2)]
        edges = [(0, 1), (1, 2)]
        for sw, dw, expected in ((1, 0, 1), (0, 1, 14), (2, 3, 44)):
            result = solve(gates, dag, edges, 3, mode="paper", horizon=None,
                           gate_durations=(2, 3, 4), swap_duration=5,
                           swap_weight=sw, depth_weight=dw)
            self.assertEqual(result["status"], "optimal")
            self.assertEqual(result["objective"], expected)
            self.assertIsNone(result["schedule"])

    def test_unlimited_swap_objective(self):
        rng = random.Random(1729)
        for case in range(8):
            gates = [tuple(rng.sample(range(3), 2)) for _ in range(4)]
            dag = [(g, g + 1) for g in range(3)]
            edges = [(0, 1), (1, 2)]
            expected = exhaustive(gates, dag, edges, 3, 8, qmrs=False)
            result = solve(gates, dag, edges, 3, mode="paper", horizon=None)
            self.assertEqual(result["objective"], expected)
            self.assertEqual(result["status"], "optimal")

    def test_limits_are_not_optimality_claims(self):
        data = load_instance(ROOT / "benchmark/qv6_depth6_01_opt7.txt")
        for kwargs, status in (({"node_limit": 1}, "node_limit"),
                               ({"time_limit": 1e-12}, "time_limit")):
            result = solve(*data[:4], horizon=12, **kwargs)
            self.assertEqual(result["status"], status)
            self.assertFalse(result["proven"])
            self.assertLessEqual(result["bound"], 7)
            json.dumps(result, allow_nan=False)

    def test_disconnected_hardware_empty_circuit_and_validation(self):
        result = solve([(0, 1), (1, 2)], [(0, 1)], [(0, 1)], 3, mode="paper", horizon=None)
        self.assertEqual(result["status"], "infeasible")
        self.assertEqual(solve([], [], [], 2, mode="paper", horizon=None)["objective"], 0)
        self.assertEqual(solve([], [], [], 2, horizon=12)["status"], "infeasible")
        invalid = [dict(gate_dag=[(0, 1), (1, 0)]), dict(hw_edges=[(0, 1), (1, 0)]),
                   dict(n=0), dict(horizon=0), dict(time_limit=-1), dict(gate_duration=0)]
        for change in invalid:
            args = dict(gates=[(0, 1), (0, 1)], gate_dag=[], hw_edges=[(0, 1)], n=2, horizon=12)
            args.update(change)
            with self.assertRaises(ValueError):
                solve(**args)

    def test_reverse_swaps_through_unassigned_sites(self):
        # This witness assigns q2 for the first time at physical site 1 AFTER
        # swapping it with site 2. Its initial position must therefore be 2.
        from luteberget.solver import _Node, _Problem, _solution
        problem = _Problem([(0, 1), (1, 2)], [(0, 1)], [(0, 1), (1, 2)], 3,
                           [1, 1], 1, False, lambda: None)
        root = _Node(0, (-1, -1, -1), (0, 0, 0), (), 0)
        first = next(c for c in problem.children(root, True, 12) if c.assignment == (0, 1, -1))
        moved = next(c for c in problem.children(first, True, 12) if c.assignment == (0, 2, -1))
        last = next(c for c in problem.children(moved, True, 12) if c.done == 3)
        result = _solution(last, problem)
        self.assertEqual(result["initial_placement"], [0, 1, 2])
        self.assertTrue(validate_solution(problem.gates, problem.dag, problem.edges, 3,
                                          result, horizon=12, qmrs=True))

    def test_all_minimal_paths_include_longer_detour(self):
        from luteberget.solver import _Problem
        edges = [(0, 1), (1, 4), (0, 2), (2, 3), (3, 4)]
        problem = _Problem([], [], edges, 5, [], 1, False, lambda: None)
        self.assertEqual(set(problem.paths(0, 4)), {(0, 1, 4), (0, 2, 3, 4)})

    def test_layering_and_label_independence(self):
        from luteberget.__main__ import main
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "incorrect_label_opt999.txt"
            path.write_text("N = 2\noptimal_swaps = 999\nG =\n0 1\ngate_dag =\nhw_edges =\n0 1\n")
            output = Path(directory) / "results.jsonl"
            self.assertEqual(main([str(path), "--horizon", "12", "--output", str(output)]), 0)
            result = json.loads(output.read_text())
            self.assertEqual(result["objective"], 0)
            self.assertEqual(result["label_swaps"], 999)
        result = solve([(0, 1), (2, 3), (1, 2)], [(0, 2), (1, 2)],
                       [(0, 1), (1, 2), (2, 3)], 4, mode="paper", horizon=None,
                       layered=True)
        self.assertEqual(result["status"], "optimal")


@unittest.skipUnless(os.environ.get("QMRS_TEST_GUROBI") == "1", "optional Gurobi integration")
class LegacyIntegrationTests(unittest.TestCase):
    def test_benchmark_witness_satisfies_legacy_model(self):
        import gurobipy as gp
        from gurobipy import GRB
        from luteberget.compare import build_legacy_model, _witness
        path = ROOT / "benchmark/qv6_depth6_01_opt7.txt"
        gates, dag, edges, n, _ = load_instance(path)
        result = solve(gates, dag, edges, n, horizon=12)
        schedule = result["schedule"]
        with gp.Env(empty=True) as env:
            env.setParam("OutputFlag", 0)
            env.start()
            scope = build_legacy_model(path, 12, env)
            model = scope["model"]
            try:
                swaps = {tuple(s) for s in schedule["swaps"]}
                for (t, g), var in scope["a"].items():
                    var.LB = var.UB = int(schedule["gate_times"][g] == t)
                for (t, e), var in scope["s"].items():
                    var.LB = var.UB = int((t, e) in swaps)
                for (t, q, i), var in scope["w"].items():
                    var.LB = var.UB = int(schedule["placement"][t][q] == i)
                model.optimize()
                self.assertEqual(model.Status, GRB.OPTIMAL)
                self.assertEqual(model.ObjVal, 7)
                _witness(scope)
            finally:
                model.dispose()


if __name__ == "__main__":
    unittest.main()
