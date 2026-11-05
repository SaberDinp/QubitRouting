"""Bounded beam heuristic for arbitrary undirected routing hardware.

Sampled placements and a capped frontier supply feasible MIP starts only.
No lower bound or optimality claim comes from this heuristic, and it never
enumerates the complete placement/state graph. Returned schedules are replayed
by the independent validator used for the exact formulations.
"""

from collections import deque
import math
import random
import time

from .compact_formulation import _input_and_windows
from .validation import validate_schedule


def _hardware_data(edges, n, rng, matching_limit):
    neighbors = [[] for _ in range(n)]
    for e, (i, j) in enumerate(edges):
        neighbors[i].append((j, e))
        neighbors[j].append((i, e))
    distance = [[float("inf")] * n for _ in range(n)]
    for source in range(n):
        distance[source][source] = 0
        queue = deque([source])
        while queue:
            i = queue.popleft()
            for j, _ in neighbors[i]:
                if distance[source][j] == float("inf"):
                    distance[source][j] = distance[source][i] + 1
                    queue.append(j)

    # Always include the empty matching and individual hardware edges.
    # Small sparse graphs are enumerated exactly up to the cap; dense graphs
    # additionally get sampled maximal/partial matchings for more variety.
    limit = max(matching_limit, len(edges) + 1)
    found = {(): 0}
    for e, (i, j) in enumerate(edges):
        found[(e,)] = (1 << i) | (1 << j)
    enumeration_limit = max(len(found), limit // 2)
    capped = False

    def enumerate_matchings(available, chosen, occupied):
        nonlocal capped
        if len(found) >= enumeration_limit:
            capped = True
            return
        if not available:
            found[tuple(sorted(chosen))] = occupied
            return
        bit = available & -available
        i = bit.bit_length() - 1
        rest = available ^ bit
        enumerate_matchings(rest, chosen, occupied)
        for j, e in neighbors[i]:
            jbit = 1 << j
            if rest & jbit:
                enumerate_matchings(rest ^ jbit, chosen + (e,), occupied | bit | jbit)

    enumerate_matchings((1 << n) - 1, (), 0)
    if capped:
        ids = list(range(len(edges)))
        for attempt in range(8 * limit):
            if len(found) >= limit:
                break
            rng.shuffle(ids)
            occupied, chosen = 0, []
            for e in ids:
                i, j = edges[e]
                bits = (1 << i) | (1 << j)
                if not occupied & bits and (attempt % 2 == 0 or rng.random() < .75):
                    occupied |= bits
                    chosen.append(e)
            found[tuple(sorted(chosen))] = occupied
    matchings = sorted(found.items(), key=lambda item: (len(item[0]), item[0]))
    return distance, matchings, capped


def construct(gates, dag, edges, n, horizon, *, width=100, seed=0,
              weights=(1.0, 2.0, 0.5, 0.125), time_limit=1.0, matching_limit=256):
    start = time.perf_counter()
    if type(width) is not int or width < 1:
        raise ValueError("width must be a positive integer")
    if type(matching_limit) is not int or matching_limit < 1:
        raise ValueError("matching_limit must be a positive integer")
    if not math.isfinite(time_limit) or time_limit <= 0:
        raise ValueError("time_limit must be positive and finite")
    gates, dag, canonical_edges, _, _ = _input_and_windows(gates, dag, edges, n, horizon)
    # Canonicalization preserves edge indices, used by the schedule validator.
    edges = canonical_edges
    rng = random.Random(seed)
    distance, matchings, matchings_capped = _hardware_data(edges, n, rng, matching_limit)
    pred = [0] * len(gates)
    predecessors = [[] for _ in gates]
    successors = [[] for _ in gates]
    indegree = [0] * len(gates)
    for g, h in dag:
        pred[h] |= 1 << g
        predecessors[h].append(g)
        successors[g].append(h)
        indegree[h] += 1
    queue = deque(g for g in range(len(gates)) if not indegree[g])
    topological = []
    while queue:
        g = queue.popleft()
        topological.append(g)
        for h in successors[g]:
            indegree[h] -= 1
            if not indegree[h]:
                queue.append(h)
    all_done = (1 << len(gates)) - 1
    roots = [g for g in range(len(gates)) if not pred[g]]
    ready_cache = {}

    def ready(done):
        if done not in ready_cache:
            ready_cache[done] = [g for g in range(len(gates))
                                 if not done >> g & 1 and pred[g] & done == pred[g]]
        return ready_cache[done]

    # Embed disjoint source gates in sampled physical matchings. Unlike a
    # chain-block construction, this uses only edges actually in the hardware.
    initial = set()
    preferred = sorted(matchings, key=lambda item: -len(item[0]))
    largest = len(preferred[0][0]) if preferred else 0
    preferred = [item for item in preferred if len(item[0]) == largest]
    for _ in range(min(300, max(width, 64))):
        p = list(range(n))
        rng.shuffle(p)
        initial.add(tuple(p))
        root_order = roots[:]
        rng.shuffle(root_order)
        blocks, used = [], set()
        for g in root_order:
            if not used.intersection(gates[g]):
                blocks.append(gates[g])
                used.update(gates[g])
        chosen, _ = rng.choice(preferred)
        chosen = list(chosen)
        rng.shuffle(chosen)
        p = [-1] * n
        assigned = set()
        for (q, r), e in zip(blocks, chosen):
            i, j = edges[e]
            if rng.randrange(2):
                i, j = j, i
            p[i], p[j] = q, r
            assigned.update((q, r))
        unused = [q for q in range(n) if q not in assigned]
        rng.shuffle(unused)
        for i in range(n):
            if p[i] == -1:
                p[i] = unused.pop()
        initial.add(tuple(p))

    best = None
    generated = 0
    inspected = 0
    timed_out = False
    for weight in weights:
        beam = [(p, 0, 0, []) for p in sorted(initial)]
        for t in range(horizon):
            candidates = {}
            for p, done, cost, hist in beam:
                inspected += 1
                if inspected % 16 == 0 and time.perf_counter() - start >= time_limit:
                    timed_out = True
                    break
                pos = [0] * n
                for i, q in enumerate(p):
                    pos[q] = i
                executable = [g for g in ready(done)
                              if distance[pos[gates[g][0]]][pos[gates[g][1]]] == 1]
                options = matchings if 0 < t < horizon - 1 else matchings[:1]
                for chosen, occupied in options:
                    new_cost = cost + len(chosen)
                    if best is not None and new_cost >= best[0]:
                        continue
                    gate_mask = 0
                    used_nodes = occupied
                    for g in executable:
                        q, r = gates[g]
                        bits = (1 << pos[q]) | (1 << pos[r])
                        if not bits & used_nodes:
                            used_nodes |= bits
                            gate_mask |= 1 << g
                    if not chosen and not gate_mask:
                        continue
                    new_done = done | gate_mask
                    next_p = list(p)
                    for e in chosen:
                        i, j = edges[e]
                        next_p[i], next_p[j] = next_p[j], next_p[i]
                    next_p = tuple(next_p)
                    key = (next_p, new_done)
                    if key in candidates and candidates[key][2] <= new_cost:
                        continue
                    new_hist = hist + [(p, gate_mask, chosen)]
                    generated += 1
                    if new_done == all_done:
                        best = (new_cost, new_hist, next_p)
                        continue
                    next_pos = [0] * n
                    for i, q in enumerate(next_p):
                        next_pos[q] = i
                    finish = [0] * len(gates)
                    for g in topological:
                        if new_done >> g & 1:
                            continue
                        if pred[g] & new_done == pred[g]:
                            q, r = gates[g]
                            gap = distance[next_pos[q]][next_pos[r]]
                            finish[g] = float("inf") if gap == float("inf") else 1 + gap // 2
                        else:
                            finish[g] = 1 + max((finish[h] for h in predecessors[g]), default=0)
                    critical = max(finish, default=0)
                    if critical > horizon - t - 1:
                        continue
                    separation = sum(max(0, distance[next_pos[gates[g][0]]][next_pos[gates[g][1]]] - 1)
                                     for g in ready(new_done))
                    score = (new_cost + weight * (len(gates) - new_done.bit_count())
                             + .2 * separation + .3 * critical + rng.random() * .001)
                    candidates[key] = (next_p, new_done, new_cost, new_hist, score)
            if timed_out:
                break
            ordered = sorted(candidates.values(), key=lambda v: (v[4], v[2], v[0], v[1]))[:width]
            beam = [v[:4] for v in ordered]
            if time.perf_counter() - start >= time_limit:
                timed_out = True
                break
            if not beam:
                break
        if timed_out or (best is not None and best[0] == 0):
            break

    result = dict(seconds=time.perf_counter() - start, width=width, seed=seed,
                  generated=generated, timed_out=timed_out, schedule=None, objective=None,
                  sampled_initial_placements=len(initial), hardware_matchings=len(matchings),
                  hardware_matchings_capped=matchings_capped)
    if best is not None:
        cost, history, final = best
        placement, gate_times, swaps = [], {}, []
        for t in range(horizon):
            p, mask, chosen = history[t] if t < len(history) else (final, 0, ())
            placement.append([p.index(q) for q in range(n)])
            gate_times.update({g: t for g in range(len(gates)) if mask >> g & 1})
            swaps.extend([t, e] for e in chosen)
        schedule = dict(placement=placement, gate_times=gate_times, swaps=swaps)
        validate_schedule(gates, dag, edges, n, horizon, schedule, cost)
        result.update(schedule=schedule, objective=cost)
    result["seconds"] = time.perf_counter() - start
    return result


def set_general_start(problem, schedule):
    """Populate the four core variable families of the general routing MIP."""
    if schedule is None:
        return
    placement, gate_times = schedule["placement"], schedule["gate_times"]
    chosen = {tuple(item) for item in schedule["swaps"]}
    edges = {frozenset(pair): e for e, pair in enumerate(problem.edges)}
    for (t, q, i), variable in problem.w.items():
        variable.Start = int(placement[t][q] == i)
    for (t, g), variable in problem.a.items():
        variable.Start = int(gate_times[g] == t)
    for (t, q, i, j), variable in problem.f.items():
        variable.Start = int(placement[t][q] == i and (t, edges[frozenset((i, j))]) in chosen)
    for (t, g, i, j), variable in problem.x.items():
        q, r = problem.gates[g]
        variable.Start = int(gate_times[g] == t and placement[t][q] == i and placement[t][r] == j)
