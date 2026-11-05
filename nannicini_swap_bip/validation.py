"""Witness replay independent of Gurobi variables and flow equations."""


def validate_schedule(gates, dag, edges, n, layers, slots, schedule, objective):
    def check(condition, message):
        if not condition:
            raise ValueError(message)

    edge_set = {tuple(sorted(e)) for e in edges}
    check(len(schedule) == len(slots), "Wrong number of slots")
    check([tuple(s) for s in slots if s] == [tuple(g) for g in layers], "Wrong layers")
    gate_times, swap_count = {}, 0
    for t, step in enumerate(schedule):
        position = step["placement"]  # logical -> physical
        check(step["t"] == t, "Wrong slot index")
        check(sorted(position) == list(range(n)), "Placement is not a bijection")
        check(tuple(step["gates"]) == tuple(slots[t]), "Gate assigned to wrong layer")
        used = set()
        for g in step["gates"]:
            check(g not in gate_times, "Repeated gate")
            gate_times[g] = t
            p, q = gates[g]
            physical = {position[p], position[q]}
            check(tuple(sorted(physical)) in edge_set, "Gate is not adjacent")
            check(not used.intersection(physical), "Overlapping gates")
            used.update(physical)
        after = position.copy()
        inverse = {i: q for q, i in enumerate(position)}
        for i, j in step["swaps"]:
            check(tuple(sorted((i, j))) in edge_set, "Non-edge swap")
            check(not used.intersection((i, j)), "Gate/SWAP or SWAP/SWAP overlap")
            used.update((i, j))
            after[inverse[i]], after[inverse[j]] = j, i
            swap_count += 1
        if t + 1 < len(schedule):
            check(after == schedule[t + 1]["placement"], "Incorrect placement transition")
        else:
            check(not step["swaps"], "Swap in final slot")
    check(set(gate_times) == set(range(len(gates))), "Missing gates")
    check(all(gate_times[u] < gate_times[v] for u, v in dag), "Precedence violation")
    check(swap_count == objective, "Incorrect SWAP count")
    return True
