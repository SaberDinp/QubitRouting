"""Best-bound tree search with lazy assignment and Pareto dominance.

Implements Algorithms 1--3 and bounds (6)--(15) of Luteberget et al.
The optional deadline and no-swap-at-zero rule adapt the search to QMRS.
Only the supplied DAG imposes gate ordering; list order is not a constraint.
"""

from collections import deque
from dataclasses import dataclass
from functools import lru_cache
import heapq
import math
import time


class _LimitReached(Exception):
    pass


@dataclass(slots=True, eq=False)
class _Node:
    done: int
    assignment: tuple  # logical -> physical; -1 means not assigned yet
    depths: tuple  # earliest free time of each physical qubit
    releases: tuple  # extra DAG dependencies between disjoint logical pairs
    swaps: int
    parent: object = None
    action: object = None  # (gate index, physical u, physical v, start, end); -1=SWAP
    bound: float = 0
    active: bool = True


class _Problem:
    def __init__(self, gates, dag, edges, n, durations, swap_duration, layered,
                 check_time):
        if not isinstance(n, int) or n < 1:
            raise ValueError("N must be a positive integer")
        self.n = n
        self.gates = tuple(tuple(g) for g in gates)
        self.edges = tuple(tuple(sorted(e)) for e in edges)
        self.dag = tuple(dict.fromkeys(tuple(e) for e in dag))
        self.m = len(self.gates)
        self.full = (1 << self.m) - 1
        self.durations = tuple(durations)
        self.swap_duration = swap_duration
        self.layered = layered
        self.check_time = check_time
        if len(self.durations) != self.m or any(
                not isinstance(d, int) or d <= 0 for d in self.durations):
            raise ValueError("Provide one positive integer duration per gate")
        if not isinstance(swap_duration, int) or swap_duration <= 0:
            raise ValueError("swap_duration must be a positive integer")
        for name, pairs, size in (("G", self.gates, n),
                                  ("hw_edges", self.edges, n),
                                  ("gate_dag", self.dag, self.m)):
            if any(len(pair) != 2 or any(not isinstance(v, int) or v < 0 or v >= size
                                        for v in pair) or pair[0] == pair[1]
                   for pair in pairs):
                raise ValueError(f"Invalid pair in {name}")
        if len(set(self.edges)) != len(self.edges):
            raise ValueError("Duplicate undirected hardware edge")
        self.edge_set = set(self.edges)
        self.neighbors = [set() for _ in range(n)]
        for u, v in self.edges:
            self.neighbors[u].add(v)
            self.neighbors[v].add(u)
        self.pred = [0] * self.m
        self.successors = [[] for _ in self.gates]
        indegree = [0] * self.m
        for u, v in self.dag:
            self.pred[v] |= 1 << u
            self.successors[u].append(v)
            indegree[v] += 1
        queue = deque(g for g in range(self.m) if not indegree[g])
        order = []
        while queue:
            g = queue.popleft()
            order.append(g)
            for h in self.successors[g]:
                indegree[h] -= 1
                if not indegree[h]:
                    queue.append(h)
        if len(order) != self.m:
            raise ValueError("gate_dag must be acyclic")
        self.ancestors = [0] * self.m
        self.layers = [0] * self.m
        for g in order:
            for h in self.successors[g]:
                self.ancestors[h] |= self.ancestors[g] | (1 << g)
                self.layers[h] = max(self.layers[h], self.layers[g] + 1)
        self.delta = list(self.durations)  # Eq. (6), remaining critical path
        for g in reversed(order):
            self.delta[g] += max((self.delta[h] for h in self.successors[g]), default=0)
        self.qubit_gates = [[g for g, pair in enumerate(self.gates) if q in pair]
                            for q in range(n)]
        # The paper's state is sufficient for circuit dependencies on a shared
        # qubit. General input DAGs may also connect disjoint pairs: retain the
        # necessary release times, and include them in Pareto comparisons.
        extra = [(u, v) for u, v in self.dag
                 if not set(self.gates[u]).intersection(self.gates[v])]
        self.release_targets = tuple(sorted({v for _, v in extra}))
        self.release_index = {g: i for i, g in enumerate(self.release_targets)}
        self.release_successors = [[] for _ in self.gates]
        for u, v in extra:
            self.release_successors[u].append(self.release_index[v])
        self.distances = []
        for source in range(n):
            check_time()
            dist = [math.inf] * n
            dist[source] = 0
            queue = deque([source])
            while queue:
                u = queue.popleft()
                for v in self.neighbors[u]:
                    if dist[v] == math.inf:
                        dist[v] = dist[u] + 1
                        queue.append(v)
            self.distances.append(dist)
        # Instance-local caches: separate solves do not retain previous trees.
        self.remaining = lru_cache(maxsize=8192)(self._remaining)
        self.paths = lru_cache(maxsize=4096)(self._paths)

    def _remaining(self, done):
        unscheduled = tuple(g for g in range(self.m) if not (done >> g & 1))
        ready = tuple(g for g in unscheduled if self.pred[g] & done == self.pred[g])
        if self.layered and unscheduled:
            layer = min(self.layers[g] for g in unscheduled)
            ready = tuple(g for g in ready if self.layers[g] == layer)
        tails = tuple(max((self.delta[g] for g in gs if not (done >> g & 1)),
                          default=0) for gs in self.qubit_gates)
        # On circuit DAGs these sums are exactly lambda(q,m_q)-lambda(q,g).
        # Ancestors also make the bound valid when shared-qubit gates are not
        # totally ordered by the input DAG.
        before = {}
        for g in unscheduled:
            before[g] = tuple(sum(self.durations[h] for h in self.qubit_gates[q]
                                  if not (done >> h & 1) and self.ancestors[g] >> h & 1)
                              for q in self.gates[g])
        return unscheduled, ready, tails, before

    def _paths(self, source, target):
        """All minimal (chordless) paths, not just shortest paths: Section 3.3."""
        result = []
        stack = [(source,)]
        while stack:
            self.check_time()
            path = stack.pop()
            if path[-1] == target:
                result.append(path)
                continue
            for v in sorted(self.neighbors[path[-1]]):
                if v not in path and not any(v in self.neighbors[u] for u in path[:-1]):
                    stack.append(path + (v,))
        return tuple(result)

    def bounds(self, node, need_depth):
        unscheduled, _, tails, before = self.remaining(node.done)
        a, depths = node.assignment, node.depths
        extra_swaps = 0
        depth_bound = max(depths)
        if need_depth:
            # Eq. (10); max of tails equals delta(m_q) for circuit DAGs.
            depth_bound = max(depth_bound, max(
                ((depths[a[q]] if a[q] >= 0 else 0) + tails[q]
                 for q in range(self.n)), default=0))
            for g, release in zip(self.release_targets, node.releases):
                if not (node.done >> g & 1):
                    depth_bound = max(depth_bound, release + self.delta[g])
        for g in unscheduled:
            p, q = self.gates[g]
            u, v = a[p], a[q]
            if u < 0 or v < 0:
                continue
            distance = self.distances[u][v]
            if distance == math.inf:
                return math.inf, math.inf
            extra_swaps = max(extra_swaps, distance - 1)  # Eq. (14)
            if not need_depth:
                continue
            p_before, q_before = before[g]
            # A direct edge attains the minimum of H, with no movement.
            earliest = max(depths[u] + p_before, depths[v] + q_before)
            if distance > 1:
                earliest = math.inf
                for path in self.paths(u, v):
                    length = len(path)
                    for j in range(length - 1):  # Eq. (11), zero-based j
                        value = max(depths[u] + p_before + j * self.swap_duration,
                                    depths[v] + q_before + (length - j - 2) * self.swap_duration)
                        for k in range(1, j + 1):
                            value = max(value, depths[path[k]] + (j + 1 - k) * self.swap_duration)
                        for k in range(j + 1, length - 1):
                            value = max(value, depths[path[k]] + (k - j) * self.swap_duration)
                        earliest = min(earliest, value)
            depth_bound = max(depth_bound, self.delta[g] + earliest)
        return depth_bound, node.swaps + extra_swaps

    def children(self, node, no_swap_at_zero, horizon):
        _, ready, _, _ = self.remaining(node.done)
        a, depths = node.assignment, node.depths
        occupied = set(a) - {-1}
        for g in ready:
            p, q = self.gates[g]
            if a[p] >= 0 and a[q] >= 0:
                positions = [(a[p], a[q])] if a[q] in self.neighbors[a[p]] else []
            elif a[p] >= 0:
                positions = [(a[p], v) for v in sorted(self.neighbors[a[p]]) if v not in occupied]
            elif a[q] >= 0:
                positions = [(u, a[q]) for u in sorted(self.neighbors[a[q]]) if u not in occupied]
            else:
                # Both orientations are essential even though hardware is undirected.
                positions = [(u, v) for edge in self.edges if not occupied.intersection(edge)
                             for u, v in (edge, edge[::-1])]
            for u, v in positions:
                release = node.releases[self.release_index[g]] if g in self.release_index else 0
                start = max(depths[u], depths[v], release)
                end = start + self.durations[g]
                if horizon is not None and end > horizon:
                    continue
                assignment = list(a)
                assignment[p], assignment[q] = u, v
                new_depths = list(depths)
                new_depths[u] = new_depths[v] = end
                releases = list(node.releases)
                if g in self.release_index:
                    releases[self.release_index[g]] = 0
                for i in self.release_successors[g]:
                    releases[i] = max(releases[i], end)
                yield _Node(node.done | (1 << g), tuple(assignment), tuple(new_depths),
                            tuple(releases), node.swaps, node, (g, u, v, start, end))
        for u, v in self.edges:
            if u not in occupied and v not in occupied:
                continue  # Algorithm 2 skips SWAPs between two unassigned sites.
            start = max(depths[u], depths[v], int(no_swap_at_zero))
            end = start + self.swap_duration
            # QMRS has no swap variables at time T-1.
            if horizon is not None and (end > horizon or (no_swap_at_zero and end >= horizon)):
                continue
            assignment = tuple(v if pos == u else u if pos == v else pos for pos in a)
            new_depths = list(depths)
            new_depths[u] = new_depths[v] = end
            yield _Node(node.done, assignment, tuple(new_depths), node.releases,
                        node.swaps + 1, node, (-1, u, v, start, end))


def _dominates(a, b):
    return (a.swaps <= b.swaps
            and all(x <= y for x, y in zip(a.depths, b.depths))
            and all(x <= y for x, y in zip(a.releases, b.releases)))


def _solution(node, problem):
    actions = []
    current = node
    while current.parent is not None:
        actions.append(current.action)
        current = current.parent
    actions.reverse()
    final = list(node.assignment)
    unused = iter(sorted(set(range(problem.n)) - set(final)))
    for q, pos in enumerate(final):
        if pos < 0:
            final[q] = next(unused)
    initial = final[:]
    # Reverse swaps through empty placeholders as well. A qubit first assigned
    # after a SWAP need not have started at its first gate's physical location.
    for g, u, v, _, _ in reversed(actions):
        if g == -1:
            initial = [v if pos == u else u if pos == v else pos for pos in initial]
    operations = []
    for g, u, v, start, end in actions:
        op = {"kind": "swap" if g == -1 else "gate", "physical": [u, v],
              "start": start, "duration": end - start}
        if g >= 0:
            op["gate"] = g
            op["logical"] = list(problem.gates[g])
        operations.append(op)
    operations.sort(key=lambda op: (op["start"], op["kind"], op.get("gate", -1), op["physical"]))
    return {"initial_placement": initial, "final_placement": final,
            "operations": operations, "swap_count": node.swaps,
            "makespan": max(node.depths)}


def solve(gates, gate_dag, hw_edges, n, *, horizon=None, mode="qmrs",
          gate_duration=1, gate_durations=None, swap_duration=1,
          swap_weight=1, depth_weight=0, layered=False,
          time_limit=60.0, node_limit=None):
    """Return a JSON-serializable result, including an independently checked witness.

    ``qmrs`` requires an explicit horizon and matches legacy_model.py: unit durations,
    a finite T, no swaps at t=0,
    and swap minimization. ``paper`` allows no deadline (horizon=None), other
    positive integer durations, and nonnegative weighted depth/SWAP objectives.
    A limit result is never reported as optimal unless the bound proves it.
    The benchmark's optimal_swaps annotation is deliberately not an argument.
    """
    started = time.perf_counter()
    if mode not in ("qmrs", "paper"):
        raise ValueError("mode must be 'qmrs' or 'paper'")
    if horizon is not None and (not isinstance(horizon, int) or horizon < 1):
        raise ValueError("horizon must be a positive integer or None")
    if (not math.isfinite(swap_weight) or not math.isfinite(depth_weight)
            or min(swap_weight, depth_weight) < 0 or swap_weight + depth_weight <= 0):
        raise ValueError("Objective weights must be finite, nonnegative, and not both zero")
    if time_limit is not None and (not math.isfinite(time_limit) or time_limit <= 0):
        raise ValueError("time_limit must be positive and finite, or None")
    if node_limit is not None and (not isinstance(node_limit, int) or node_limit < 1):
        raise ValueError("node_limit must be a positive integer or None")
    gates = tuple(gates)
    gate_dag = tuple(gate_dag)
    hw_edges = tuple(hw_edges)
    durations = tuple(gate_durations) if gate_durations is not None else (gate_duration,) * len(gates)
    if mode == "qmrs" and (horizon is None or swap_duration != 1 or any(d != 1 for d in durations)
                           or swap_weight != 1 or depth_weight != 0 or layered):
        raise ValueError("qmrs mode requires a horizon, unit durations, swap objective, and no layering")
    deadline = math.inf if time_limit is None else started + time_limit

    def check_time():
        if time.perf_counter() >= deadline:
            raise _LimitReached

    stats = dict(expanded=0, generated=0, inserted=0, dominated=0,
                 bound_pruned=0, deadline_pruned=0, peak_open=0, pareto_states=0)
    incumbent = None
    incumbent_value = math.inf
    bound = 0
    status = "time_limit"
    problem = None
    fronts = {}
    need_depth = horizon is not None or depth_weight > 0
    try:
        problem = _Problem(gates, gate_dag, hw_edges, n, durations, swap_duration,
                           layered, check_time)
        root = _Node(0, (-1,) * n, (0,) * n, (0,) * len(problem.release_targets), 0)
        queue = []
        sequence = 0

        def insert(node):
            nonlocal sequence, incumbent, incumbent_value
            stats["generated"] += 1
            check_time()
            hd, hs = problem.bounds(node, need_depth)
            if math.isinf(hs) or (horizon is not None and hd > horizon):
                stats["deadline_pruned"] += 1
                return
            node.bound = max(node.parent.bound if node.parent else 0,
                             depth_weight * hd + swap_weight * hs)
            if node.bound >= incumbent_value:
                stats["bound_pruned"] += 1
                return
            key = (node.done, node.assignment)
            front = fronts.get(key, [])
            if any(_dominates(other, node) for other in front):
                stats["dominated"] += 1
                return
            kept = []
            for other in front:
                if _dominates(node, other):
                    other.active = False
                else:
                    kept.append(other)
            kept.append(node)
            fronts[key] = kept
            stats["inserted"] += 1
            if node.done == problem.full:
                # legacy_model.py requires at least one gate at t=0; the nonempty search
                # has this automatically because its first action is a gate.
                if mode == "qmrs" and not problem.m:
                    return
                incumbent = node
                incumbent_value = swap_weight * node.swaps + depth_weight * max(node.depths)
                return
            sequence += 1
            heapq.heappush(queue, (node.bound, -node.done.bit_count(), max(node.depths),
                                   sequence, node))
            stats["peak_open"] = max(stats["peak_open"], len(queue))

        insert(root)
        while queue:
            while queue and not queue[0][-1].active:
                heapq.heappop(queue)
            if not queue:
                break
            bound = min(queue[0][0], incumbent_value)
            if bound >= incumbent_value:
                break
            check_time()
            if node_limit is not None and stats["expanded"] >= node_limit:
                status = "node_limit"
                break
            *_, node = heapq.heappop(queue)
            stats["expanded"] += 1
            for child in problem.children(node, mode == "qmrs", horizon):
                insert(child)
        else:
            status = "optimal" if incumbent is not None else "infeasible"
        if status != "node_limit":
            status = "optimal" if incumbent is not None else "infeasible"
        if status == "optimal":
            bound = incumbent_value
        elif status == "infeasible":
            bound = None
    except _LimitReached:
        # If interrupted during expansion, `bound` still includes the popped
        # parent and therefore covers all its not-yet-generated descendants.
        status = "time_limit"
    stats["pareto_states"] = len(fronts)
    search_seconds = time.perf_counter() - started
    solution = _solution(incumbent, problem) if incumbent is not None else None
    schedule = None
    if solution is not None:
        from .validation import validate_solution, unit_schedule
        validate_solution(gates, gate_dag, hw_edges, n, solution, durations,
                          swap_duration, horizon=horizon, qmrs=(mode == "qmrs"))
        if swap_duration == 1 and all(d == 1 for d in durations):
            schedule = unit_schedule(solution, hw_edges, horizon or solution["makespan"])
    value = None if solution is None else incumbent_value
    if value is not None and bound is not None and bound >= value:
        status, bound = "optimal", value
    return {"method": "luteberget", "status": status,
            "proven": status in ("optimal", "infeasible"), "mode": mode,
            "horizon": horizon, "layered": layered, "objective": value,
            "swap_count": None if solution is None else solution["swap_count"],
            "makespan": None if solution is None else solution["makespan"],
            "bound": bound, "gap": None if value is None or bound is None else
            (value - bound) / max(abs(value), 1e-10),
            "swap_weight": swap_weight, "depth_weight": depth_weight,
            "gate_durations": list(durations), "swap_duration": swap_duration,
            "time_limit": time_limit, "node_limit": node_limit,
            "seconds": search_seconds, "total_seconds": time.perf_counter() - started,
            "stats": stats, "solution": solution, "schedule": schedule}
