"""Conditional proper-subset routing lower bounds for arbitrary hardware.

Selected logical qubits keep their physical positions.  Other qubits are
anonymous, not deleted physical vertices. A swap moves one or two selected
qubits and costs that many tracked-qubit movements.  Small subset shortest
paths supply conditional inequalities on total per-qubit movement counts.
The explicit one-anonymous option represents every physical placement while
omitting that qubit's gates. If it is unused, all gate progress is retained;
otherwise this is a proper circuit relaxation. The helper's default requires
at least two anonymous qubits; the solver explicitly opts in for small cases.
"""

from __future__ import annotations

from heapq import heappop, heappush
from itertools import combinations, permutations
from math import perm
import time

import gurobipy as gp
from gurobipy import GRB

from .compact_formulation import _input_and_windows


class ProjectionResourceLimit(ValueError):
    """A local projection is deliberately skipped, not declared infeasible."""


def _check_deadline(deadline):
    if deadline is not None and time.perf_counter() >= deadline:
        raise TimeoutError("Projection preprocessing time budget exhausted")


def ancestors(gates, dag, deadline=None):
    pred = [0] * len(gates)
    for g, h in dag:
        pred[h] |= 1 << g
    changed = True
    while changed:
        _check_deadline(deadline)
        changed = False
        for h in range(len(gates)):
            previous = pred[h]
            for g in range(len(gates)):
                if pred[h] >> g & 1:
                    pred[h] |= pred[g]
            changed |= previous != pred[h]
    return pred


def _projection_graph(gates, predecessor_masks, edges, n, subset, *, ideal_limit=20000,
                      allow_one_anonymous=False, state_limit=250000,
                      arc_limit=1000000, deadline=None):
    """Build once; arcs retain which tracked qubits move, independently of cost."""
    subset = tuple(subset)
    _check_deadline(deadline)
    if len(subset) >= n:
        raise ValueError("Do not project onto the complete set of logical qubits")
    if len(subset) >= n - 1 and not allow_one_anonymous:
        raise ValueError("Keep at least two anonymous qubits; do not enumerate full placements")
    index = {q: j for j, q in enumerate(subset)}
    retained = [g for g, pair in enumerate(gates) if all(q in index for q in pair)]
    local_gates = [(index[gates[g][0]], index[gates[g][1]]) for g in retained]
    pred = [sum(1 << j for j, g in enumerate(retained) if predecessor_masks[h] >> g & 1)
            for h in retained]
    full = (1 << len(retained)) - 1
    ideals, seen = [0], {0}
    for mask in ideals:
        _check_deadline(deadline)
        for g in range(len(retained)):
            if not (mask >> g & 1) and pred[g] & mask == pred[g]:
                nxt = mask | (1 << g)
                if nxt not in seen:
                    seen.add(nxt)
                    ideals.append(nxt)
                    if len(ideals) > ideal_limit:
                        raise ProjectionResourceLimit(f"Subset {subset} exceeds ideal preprocessing limit")
    if perm(n, len(subset)) > state_limit:
        raise ProjectionResourceLimit(f"Subset {subset} exceeds local placement limit {state_limit}")
    if perm(n, len(subset)) * len(ideals) > 4 * state_limit:
        raise ProjectionResourceLimit(f"Subset {subset} exceeds closure preprocessing limit {4 * state_limit}")
    placements = tuple(permutations(range(n), len(subset)))
    placement_index = {pos: i for i, pos in enumerate(placements)}
    neighbors = [set() for _ in range(n)]
    for i, j in edges:
        neighbors[i].add(j)
        neighbors[j].add(i)
    closures, moves = [], []
    for pos in placements:
        _check_deadline(deadline)
        adjacent = [pos[v] in neighbors[pos[u]] for u, v in local_gates]
        closed = {}
        for index_in_ideals, initial in enumerate(ideals):
            if index_in_ideals % 256 == 0:
                _check_deadline(deadline)
            mask = initial
            changed = True
            while changed:
                changed = False
                for g in range(len(retained)):
                    if adjacent[g] and not(mask >> g & 1) and pred[g] & mask == pred[g]:
                        mask |= 1 << g
                        changed = True
            closed[initial] = mask
        closures.append(closed)
        occupant = {node: q for q, node in enumerate(pos)}
        transitions = []
        for q, i in enumerate(pos):
            for j in neighbors[i]:
                r = occupant.get(j)
                if r is not None and r < q:
                    continue
                after = list(pos)
                after[q] = j
                if r is not None:
                    after[r] = i
                transitions.append((placement_index[tuple(after)], q, r))
        moves.append(transitions)
    states = []
    for i in range(len(placements)):
        _check_deadline(deadline)
        for mask in ideals:
            if closures[i][mask] == mask:
                states.append((i, mask))
                if len(states) > state_limit:
                    raise ProjectionResourceLimit(f"Subset {subset} exceeds local state limit {state_limit}")
    state_index = {state: i for i, state in enumerate(states)}
    reverse = [[] for _ in states]
    arc_count = 0
    for source, (i, mask) in enumerate(states):
        if source % 256 == 0:
            _check_deadline(deadline)
        for j, q, r in moves[i]:
            destination = state_index[j, closures[j][mask]]
            reverse[destination].append((source, q, r))
            arc_count += 1
            if arc_count > arc_limit:
                raise ProjectionResourceLimit(f"Subset {subset} exceeds local arc limit {arc_limit}")
    return {"placements": placements, "reverse": reverse,
            "initial": [state_index[i, closures[i][0]] for i in range(len(placements))],
            "terminal": [state_index[i, full] for i in range(len(placements))],
            "stats": {"subset": subset, "gates": len(retained), "ideals": len(ideals),
                      "placement_count": len(placements), "anonymous_qubits": n - len(subset),
                      "states": len(states), "arcs": arc_count}}


def _costs_on_graph(graph, weights, deadline=None):
    _check_deadline(deadline)
    weights = tuple(weights)
    if len(weights) != len(graph["stats"]["subset"]) or any(type(w) is not int or w < 0 for w in weights):
        raise ValueError("Weights must be one nonnegative integer per tracked qubit")
    reverse = graph["reverse"]
    infinity = float("inf")
    distance = [infinity] * len(reverse)
    queue = []
    for terminal in graph["terminal"]:
        distance[terminal] = 0
        heappush(queue, (0, terminal))
    processed = 0
    while queue:
        if processed % 256 == 0:
            _check_deadline(deadline)
        processed += 1
        cost, node = heappop(queue)
        if cost != distance[node]:
            continue
        for previous, q, r in reverse[node]:
            candidate = cost + weights[q] + (weights[r] if r is not None else 0)
            if candidate < distance[previous]:
                distance[previous] = candidate
                heappush(queue, (candidate, previous))
    result = {pos: distance[state] for pos, state in zip(graph["placements"], graph["initial"])}
    result = {pos: cost for pos, cost in result.items() if cost != infinity}
    return result, {**graph["stats"], "weights": weights,
                    "min_cost": min(result.values(), default=None),
                    "max_cost": max(result.values(), default=None)}


def conditional_costs(gates, predecessor_masks, edges, n, subset, *, ideal_limit=20000,
                      weights=None, allow_one_anonymous=False, state_limit=250000,
                      arc_limit=1000000, deadline=None):
    """Return exact no-deadline projection costs, in tracked-movement units."""
    subset = tuple(subset)
    weights = (1,) * len(subset) if weights is None else tuple(weights)
    graph = _projection_graph(gates, predecessor_masks, edges, n, subset,
                              ideal_limit=ideal_limit, allow_one_anonymous=allow_one_anonymous,
                              state_limit=state_limit, arc_limit=arc_limit, deadline=deadline)
    return _costs_on_graph(graph, weights, deadline=deadline)


def _prepare_tables(gates, dag, edges, n, size, min_size, weighted, allow_one_anonymous,
                    pair_weights=False, deadline=None, state_limit=250000,
                    arc_limit=1000000):
    if size >= n or (size >= n - 1 and not allow_one_anonymous) or size < 2:
        raise ValueError("Subset size must be at least two and strictly below n; n-1 requires explicit opt-in")
    if min_size is None:
        min_size = size
    if type(min_size) is not int or not 2 <= min_size <= size:
        raise ValueError("min_size must be an integer between two and size")
    start = time.perf_counter()
    try:
        _check_deadline(deadline)
        pred = ancestors(gates, dag, deadline=deadline)
    except TimeoutError as error:
        return [], time.perf_counter() - start, {
            "precompute_complete": False, "precompute_truncated": True,
            "precompute_reason": str(error), "skipped_projections": [],
            "completed_projection_tables": 0}
    data = []
    skipped, timed_out = [], False
    reason = None
    for k in range(min_size, size + 1):
        for subset in combinations(range(n), k):
            patterns = [(1,) * k]
            if weighted or pair_weights:
                patterns.extend(tuple(value if j == q else 1 for j in range(k))
                                for q in range(k) for value in (0, 2))
            if pair_weights:
                for q, r in combinations(range(k), 2):
                    for u, v in ((0, 0), (2, 2), (0, 2), (2, 0)):
                        patterns.append(tuple(u if j == q else v if j == r else 1
                                              for j in range(k)))
            tables, stats = {}, []
            try:
                graph = _projection_graph(gates, pred, edges, n, subset,
                                          allow_one_anonymous=allow_one_anonymous,
                                          deadline=deadline, state_limit=state_limit,
                                          arc_limit=arc_limit)
                for weights in patterns:
                    costs, report = _costs_on_graph(graph, weights, deadline=deadline)
                    tables[weights] = costs
                    stats.append(report)
            except ProjectionResourceLimit as error:
                skipped.append({"subset": subset, "reason": str(error)})
                reason = str(error)
            except TimeoutError as error:
                timed_out = True
                reason = str(error)
            # An incomplete Dijkstra result is never returned by _costs_on_graph.
            # The first completed all-ones table determines the exact support.
            if (1,) * k in tables:
                data.append((subset, tables, stats))
            if timed_out:
                break
        if timed_out:
            break
    precompute_seconds = time.perf_counter() - start
    return data, precompute_seconds, {
        "precompute_complete": not timed_out and not skipped,
        "precompute_truncated": timed_out,
        "precompute_reason": reason,
        "skipped_projections": skipped,
        "completed_projection_tables": sum(len(tables) for _, tables, _ in data)}


def _add_conditional_rows(m, w, c, n, data):
    all_weights = {}
    for subset, tables, stats in data:
        costs = tables[(1,) * len(subset)]
        weights = m.addVars(list(costs), name="initial_subset_" + "_".join(map(str, subset)))
        all_weights[subset] = weights
        m.addConstr(weights.sum() == 1)
        for j, q in enumerate(subset):
            for i in range(n):
                m.addConstr(gp.quicksum(v for pos, v in weights.items() if pos[j] == i) == w[q, i])
        for movement_weights, costs in tables.items():
            m.addConstr(gp.quicksum(movement_weights[j] * c[q] for j, q in enumerate(subset))
                        >= gp.quicksum(costs[pos] * var for pos, var in weights.items()))
    return all_weights


def build_bound(gates, dag, edges, n, horizon, *, size=3, env=None, min_size=None,
                weighted=False, allow_one_anonymous=False, pair_weights=False,
                deadline=None, state_limit=250000, arc_limit=1000000):
    gates, dag, edges, _, _ = _input_and_windows(gates, dag, edges, n, horizon)
    data, precompute_seconds, metadata = _prepare_tables(gates, dag, edges, n, size, min_size,
                                                        weighted, allow_one_anonymous, pair_weights,
                                                        deadline, state_limit, arc_limit)
    m = gp.Model("general_subset_movement_bound", env=env)
    w = m.addVars(n, n, vtype=GRB.BINARY, name="initial")
    c = m.addVars(n, vtype=GRB.INTEGER, name="movement_count")
    objective = m.addVar(vtype=GRB.INTEGER, name="swap_count")
    m.addConstr(2 * objective == c.sum())
    for q in range(n):
        m.addConstr(w.sum(q, "*") == 1)
        m.addConstr(w.sum("*", q) == 1)
        m.addConstr(c[q] <= min(max(0, horizon - 2), horizon - sum(q in gate for gate in gates)))
        m.addConstr(c[q] <= objective)
    _add_conditional_rows(m, w, c, n, data)
    m.setObjective(objective, GRB.MINIMIZE)
    m.update()
    return m, {**metadata, "precompute_seconds": precompute_seconds,
               "projection_stats": [item for _, _, stats in data for item in stats]}


def add_to_routing(problem, *, size=4, min_size=None, weighted=False,
                   allow_one_anonymous=False, pair_weights=False, deadline=None,
                   state_limit=250000, arc_limit=1000000):
    """Tie conditional subset inequalities to an existing general routing MIP."""
    p = problem
    start = time.perf_counter()
    data, precompute_seconds, metadata = _prepare_tables(p.gates, p.dag, p.edges, p.n, size,
                                                        min_size, weighted, allow_one_anonymous,
                                                        pair_weights, deadline, state_limit, arc_limit)
    m = p.model
    c = m.addVars(p.n, vtype=GRB.INTEGER, name="subset_movement_count")
    total = m.addVar(vtype=GRB.INTEGER, name="subset_swap_count")
    m.addConstr(2 * total == c.sum())
    for q in range(p.n):
        m.addConstr(c[q] == gp.quicksum(v for (t, r, i, j), v in p.f.items() if r == q))
        m.addConstr(c[q] <= total)
    w = {(q, i): p.w[0, q, i] for q in range(p.n) for i in range(p.n)}
    weights = _add_conditional_rows(m, w, c, p.n, data)
    m.update()
    return {**metadata, "movement": c, "total": total, "weights": weights,
            "precompute_seconds": precompute_seconds,
            "build_seconds": time.perf_counter() - start,
            "projection_stats": [item for _, _, stats in data for item in stats]}
