"""Independent factorial-size oracle; never used to build or warm-start the BIP."""

import heapq
from itertools import permutations
import math
import time


def exact_layered_optimum(gates, layers, edges, n):
    """Minimum swaps with fixed simultaneous groups and unlimited routing time.

    Enumerate full placements. Between layers take the shortest-path closure
    under single physical-edge swaps (unit cost), then retain placements where
    every gate in the next layer is adjacent. Disjoint swaps can be serialized;
    there is no depth objective. No label, Gurobi result, or lower bound is read.
    """
    if n > 7:
        raise ValueError("The independent exhaustive oracle is limited to N <= 7")
    start = time.perf_counter()
    placements = list(permutations(range(n)))  # physical -> logical
    ids = {p: i for i, p in enumerate(placements)}
    adjacent = [{frozenset((p[i], p[j])) for i, j in edges} for p in placements]
    neighbors = []
    for p in placements:
        row = []
        for i, j in edges:
            other = list(p)
            other[i], other[j] = other[j], other[i]
            row.append(ids[tuple(other)])
        neighbors.append(row)
    costs = [0] * len(placements)
    for k, layer in enumerate(layers):
        if k:
            queue = [(cost, i) for i, cost in enumerate(costs) if cost < math.inf]
            heapq.heapify(queue)
            while queue:
                cost, i = heapq.heappop(queue)
                if cost != costs[i]:
                    continue
                for j in neighbors[i]:
                    if cost + 1 < costs[j]:
                        costs[j] = cost + 1
                        heapq.heappush(queue, (cost + 1, j))
        pairs = {frozenset(gates[g]) for g in layer}
        costs = [cost if pairs <= adjacent[i] else math.inf for i, cost in enumerate(costs)]
    best = min(costs)
    return {"objective": int(best) if best < math.inf else None,
            "status": "OPTIMAL" if best < math.inf else "INFEASIBLE",
            "seconds": time.perf_counter() - start, "placements": len(placements)}
