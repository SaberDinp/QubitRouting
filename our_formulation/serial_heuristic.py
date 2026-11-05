"""Deterministic feasible routing witness used to choose a safe time horizon.

This router does not read reference optima or solve an optimization model.
It places each logical interaction component within one hardware component,
then routes a topological gate order along shortest hardware paths. Leading
SWAPs are absorbed into the free initial placement.
"""

from collections import deque
import time

from .validation import validate_schedule


class SafeHorizonUnavailable(ValueError):
    """A heuristic witness was not constructed; this is not an MIP proof."""


def _components(n, pairs):
    neighbors = [set() for _ in range(n)]
    for i, j in pairs:
        neighbors[i].add(j)
        neighbors[j].add(i)
    unseen, result = set(range(n)), []
    while unseen:
        source = min(unseen)
        unseen.remove(source)
        reached, queue = [source], deque([source])
        while queue:
            for j in sorted(neighbors[queue.popleft()]):
                if j in unseen:
                    unseen.remove(j)
                    reached.append(j)
                    queue.append(j)
        result.append(tuple(sorted(reached)))
    return sorted(result, key=lambda component: (-len(component), component))


def _initial_placement(gates, edges, n, deadline):
    logical = _components(n, gates)
    hardware = _components(n, edges)
    remaining = [len(component) for component in hardware]
    assignment = [-1] * len(logical)

    def assign(k):
        if deadline is not None and time.perf_counter() >= deadline:
            raise SafeHorizonUnavailable("Automatic horizon unavailable: component placement exceeded the time budget")
        if k == len(logical):
            return True
        tried = set()
        for h, capacity in enumerate(remaining):
            if capacity < len(logical[k]) or capacity in tried:
                continue
            tried.add(capacity)
            remaining[h] -= len(logical[k])
            assignment[k] = h
            if assign(k + 1):
                return True
            remaining[h] += len(logical[k])
        return False

    if not assign(0):
        raise SafeHorizonUnavailable(
            "Automatic horizon unavailable: no placement of logical interaction components "
            "within hardware components was found; no routing MIP was solved")
    offsets = [0] * len(hardware)
    placement = [-1] * n
    for component, h in zip(logical, assignment):
        for q in component:
            placement[q] = hardware[h][offsets[h]]
            offsets[h] += 1
    return placement


def construct(gates, dag, edges, n, *, deadline=None):
    """Return a validated serial witness and T = gate count + SWAP count.

    With free initial and final placements, every feasible route using S SWAPs
    can be serialized into at most m + S operation slots: order parallel
    disjoint operations arbitrarily, absorb leading SWAPs into the initial
    placement, and discard trailing SWAPs. Our witness gives a feasible S_h,
    so a better route has S <= S_h and fits T = m + S_h. Closing the MIP gap
    at this T therefore certifies the unrestricted minimum SWAP count.
    """
    from .compact_formulation import _input_and_windows

    start = time.perf_counter()
    gates, dag, edges, _, _ = _input_and_windows(gates, dag, edges, n, 1)
    if not gates:
        raise SafeHorizonUnavailable(
            "Automatic horizon requires a nonempty circuit under the first-slot gate convention")
    placement = _initial_placement(gates, edges, n, deadline)
    neighbors = [[] for _ in range(n)]
    edge_ids = {}
    for e, (i, j) in enumerate(edges):
        neighbors[i].append(j)
        neighbors[j].append(i)
        edge_ids[frozenset((i, j))] = e
    for row in neighbors:
        row.sort()
    successors = [[] for _ in gates]
    indegree = [0] * len(gates)
    for g, h in dag:
        successors[g].append(h)
        indegree[h] += 1
    ready = deque(g for g, degree in enumerate(indegree) if degree == 0)
    placements, gate_times, swaps = [], {}, []
    while ready:
        if deadline is not None and time.perf_counter() >= deadline:
            raise SafeHorizonUnavailable("Automatic horizon unavailable: serial routing exceeded the time budget")
        g = ready.popleft()
        q, r = gates[g]
        source, target = placement[q], placement[r]
        # Compute a shortest path from the current physical source.
        predecessor = {source: None}
        queue = deque([source])
        while target not in predecessor and queue:
            i = queue.popleft()
            for j in neighbors[i]:
                if j not in predecessor:
                    predecessor[j] = i
                    queue.append(j)
        if target not in predecessor:
            raise SafeHorizonUnavailable("Automatic horizon unavailable: gate endpoints are disconnected")
        path, node = [target], target
        while predecessor[node] is not None:
            node = predecessor[node]
            path.append(node)
        path.reverse()
        for i, j in zip(path[:-2], path[1:-1]):
            # Before the first gate, this only changes the free initial map.
            if gate_times:
                swaps.append([len(placements), edge_ids[frozenset((i, j))]])
                placements.append(placement[:])
            placement = [j if pos == i else i if pos == j else pos for pos in placement]
        gate_times[g] = len(placements)
        placements.append(placement[:])
        for h in sorted(successors[g]):
            indegree[h] -= 1
            if indegree[h] == 0:
                ready.append(h)
    horizon = len(gates) + len(swaps)
    schedule = dict(placement=placements, gate_times=gate_times, swaps=swaps)
    validate_schedule(gates, dag, edges, n, horizon, schedule, len(swaps))
    return dict(schedule=schedule, objective=len(swaps), horizon=horizon,
                seconds=time.perf_counter() - start,
                algorithm="deterministic_serial_shortest_path")
