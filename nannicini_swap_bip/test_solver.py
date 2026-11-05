"""Behavioral checks against exhaustive placements and legal SWAP matchings."""

from collections import deque
import copy
from itertools import permutations
import math
from pathlib import Path
import random
import unittest

import gurobipy as gp
from gurobipy import GRB

from load_instance import load_instance
from .exact import exact_layered_optimum
from .layers import greedy_layers, path_order, sufficient_dummy_slots
from .solver import build_model, solve
from .validation import validate_schedule


def matchings(edges):
    result = [()]
    for edge in edges:
        result += [chosen + (edge,) for chosen in result
                   if not set(edge).intersection(v for e in chosen for v in e)]
    return result


def exhaustive_slots(gates, slots, edges, n):
    """Full time-slot DP with physical permutations and disjoint operations."""
    edge_set = {frozenset(e) for e in edges}
    choices = matchings(edges)
    states = {p: 0 for p in permutations(range(n))}  # physical -> logical
    for t, layer in enumerate(slots):
        following = {}
        for placement, cost in states.items():
            pos = {q: i for i, q in enumerate(placement)}
            if any(frozenset(pos[q] for q in gates[g]) not in edge_set for g in layer):
                continue
            active = {pos[q] for g in layer for q in gates[g]}
            for chosen in choices if t + 1 < len(slots) else [()]:
                if active.intersection(v for e in chosen for v in e):
                    continue
                nxt = list(placement)
                for i, j in chosen:
                    nxt[i], nxt[j] = nxt[j], nxt[i]
                nxt = tuple(nxt)
                following[nxt] = min(following.get(nxt, math.inf), cost + len(chosen))
        states = following
    return min(states.values(), default=None)


class LayeringTests(unittest.TestCase):
    def test_dag_and_shared_qubit_conflicts(self):
        gates = [(0, 1), (2, 3), (1, 2), (0, 3)]
        self.assertEqual(greedy_layers(gates, [(0, 1)]), ((0,), (1,), (2, 3)))
        # Backward gate IDs are legal; the DAG is authoritative.
        self.assertEqual(greedy_layers([(0, 1), (2, 3)], [(1, 0)]), ((1,), (0,)))
        with self.assertRaises(ValueError):
            greedy_layers(gates, [(0, 1), (1, 0)])

    def test_odd_even_bound_preserves_minimum_swap_count(self):
        for n in range(2, 8):
            edges = [(i, i + 1) for i in range(n - 1)]
            bound = sufficient_dummy_slots(edges, n)
            for p in permutations(range(n)):
                inversion_count = sum(p[i] > p[j] for i in range(n) for j in range(i + 1, n))
                values, count = list(p), 0
                for t in range(bound):
                    for i in range(t % 2, n - 1, 2):
                        if values[i] > values[i + 1]:
                            values[i], values[i + 1] = values[i + 1], values[i]
                            count += 1
                self.assertEqual(values, list(range(n)))
                self.assertEqual(count, inversion_count)
        self.assertEqual(path_order([(2, 0), (0, 3), (3, 1)], 4), (1, 3, 0, 2))
        self.assertEqual(sufficient_dummy_slots([(0, 1), (0, 2), (0, 3)], 4), 6)

    def test_five_rounds_do_not_allow_every_six_path_permutation(self):
        n = 6
        choices = matchings([(i, i + 1) for i in range(n - 1)])
        initial = tuple(range(n))
        distance, queue = {initial: 0}, deque([initial])
        while queue:
            p = queue.popleft()
            for chosen in choices:
                nxt = list(p)
                for i, j in chosen:
                    nxt[i], nxt[j] = nxt[j], nxt[i]
                nxt = tuple(nxt)
                if nxt not in distance:
                    distance[nxt] = distance[p] + 1
                    queue.append(nxt)
        self.assertEqual(len(distance), math.factorial(n))
        self.assertEqual(max(distance.values()), 6)
        self.assertEqual(distance[tuple(reversed(initial))], 6)


class ModelTests(unittest.TestCase):
    def fixed_transition(self, n, edges, final):
        with gp.Env(params={"OutputFlag": 0}) as env:
            p = build_model([(0, 1), (0, 1)], [(0, 1)], edges, n,
                            dummy_slots=0, env=env)
            try:
                for q in range(n):
                    p.model.addConstr(p.w[0, q, q] == 1)
                    p.model.addConstr(p.w[1, q, final[q]] == 1)
                p.model.Params.Threads = 1
                p.model.optimize()
                return p.model.Status, p.model.ObjVal if p.model.SolCount else None
            finally:
                p.model.dispose()

    def test_gate_freeze_idle_swaps_and_count(self):
        path = [(0, 1), (1, 2), (2, 3)]
        self.assertEqual(self.fixed_transition(4, path, [1, 0, 2, 3])[0], GRB.INFEASIBLE)
        self.assertEqual(self.fixed_transition(4, path, [0, 2, 1, 3])[0], GRB.INFEASIBLE)
        self.assertEqual(self.fixed_transition(4, path, [0, 1, 3, 2]), (GRB.OPTIMAL, 1))
        path6 = [(i, i + 1) for i in range(5)]
        self.assertEqual(self.fixed_transition(6, path6, [0, 1, 3, 2, 5, 4]), (GRB.OPTIMAL, 2))

    def test_directed_cycle_is_not_a_swap(self):
        edges = [(0, 1), (1, 2), (2, 3), (3, 4), (2, 4)]
        self.assertEqual(self.fixed_transition(5, edges, [0, 1, 3, 4, 2])[0], GRB.INFEASIBLE)

    def test_small_models_against_slot_enumeration(self):
        rng = random.Random(42)
        for n in range(2, 5):
            for trial in range(5):
                edges = [(i, i + 1) for i in range(n - 1)]
                if n > 2 and trial % 2:
                    edges.append((0, n - 1))
                gates = [tuple(rng.sample(range(n), 2)) for _ in range(4)]
                dag = [(i, j) for i in range(4) for j in range(i + 1, 4) if rng.random() < .3]
                for dummy in (0, 1, None):
                    with self.subTest(n=n, trial=trial, dummy=dummy):
                        r = solve(gates, dag, edges, n, dummy_slots=dummy, time_limit=15)
                        layers = r["layers"]
                        slots = [()] * r["total_slots"]
                        for layer, t in zip(layers, r["layer_times"]):
                            slots[t] = layer
                        expected = exhaustive_slots(gates, slots, edges, n)
                        self.assertIn(r["status"], ("OPTIMAL", "INFEASIBLE"))
                        self.assertEqual(r["objective"], expected)
                        if dummy is None:
                            self.assertEqual(expected, exact_layered_optimum(gates, layers, edges, n)["objective"])

    def test_empty_single_layer_and_infeasible_group(self):
        self.assertEqual(solve([], [], [], 1)["objective"], 0)
        self.assertEqual(solve([(0, 1)], [], [(0, 1)], 2)["objective"], 0)
        r = solve([(0, 1), (2, 3)], [], [(0, 1), (0, 2), (0, 3)], 4)
        self.assertEqual(r["status"], "INFEASIBLE")
        with self.assertRaises(ValueError):
            solve([(0, 1)], [], [], 2)
        with self.assertRaises(ValueError):
            solve([(0, 1)], [], [(0, 1), (1, 0)], 2)

    def test_compaction_and_validator(self):
        path = Path(__file__).resolve().parents[1] / "benchmark/instance_01_opt2.txt"
        gates, dag, edges, n, _ = load_instance(path)
        a = solve(gates, dag, edges, n, compact=True)
        b = solve(gates, dag, edges, n, compact=False)
        self.assertEqual(a["objective"], b["objective"])
        bad = copy.deepcopy(a["schedule"])
        bad[0]["placement"][0] = bad[0]["placement"][1]
        slots = [step["gates"] for step in a["schedule"]]
        with self.assertRaises(ValueError):
            validate_schedule(gates, dag, edges, n, a["layers"], slots, bad, a["objective"])


if __name__ == "__main__":
    unittest.main()
