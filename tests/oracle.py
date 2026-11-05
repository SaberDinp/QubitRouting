"""Small exhaustive schedule oracle, independent of the production models."""

from itertools import permutations


def exhaustive(gates, dag, edges, n, horizon):
    """Tiny independent slot enumeration used only for correctness tests.

    Placements map logical qubits to physical sites.
    Idle slots are allowed: removing them must preserve the optimum.
    """
    frontier = {(p, frozenset()): 0 for p in permutations(range(n))}
    hardware = {frozenset(edge) for edge in edges}
    best = None
    for t in range(horizon):
        following = {}
        for (placement, done), cost in frontier.items():
            ready = [g for g, (p, q) in enumerate(gates)
                     if g not in done and all(u in done for u, v in dag if v == g)
                     and frozenset((placement[p], placement[q])) in hardware]
            for gate_bits in range(1 << len(ready)):
                selected = [g for k, g in enumerate(ready) if gate_bits >> k & 1]
                if t == 0 and not selected:
                    continue
                sites = [placement[q] for g in selected for q in gates[g]]
                if len(sites) != len(set(sites)):
                    continue
                swap_sets = range(1 << len(edges)) if 0 < t < horizon - 1 else (0,)
                for swap_bits in swap_sets:
                    chosen = [e for e in range(len(edges)) if swap_bits >> e & 1]
                    used = sites + [i for e in chosen for i in edges[e]]
                    if len(used) != len(set(used)):
                        continue
                    after = list(placement)
                    for e in chosen:
                        i, j = edges[e]
                        after = [j if pos == i else i if pos == j else pos for pos in after]
                    new_done = done | frozenset(selected)
                    new_cost = cost + len(chosen)
                    if len(new_done) == len(gates):
                        best = new_cost if best is None else min(best, new_cost)
                    key = (tuple(after), new_done)
                    following[key] = min(following.get(key, float("inf")), new_cost)
        frontier = following
    return best
