"""Independent schedule replay and Gurobi result serialization.

These helpers have no dependency on a routing formulation.
"""

import math

from gurobipy import GRB


def validate_schedule(gates, dag, edges, n, horizon, schedule, objective):
    """Replay every operation independently of graph generation and MIP rows."""
    placement = schedule["placement"]
    times = {int(g): t for g, t in schedule["gate_times"].items()}
    swaps = schedule["swaps"]
    if len(placement) != horizon or any(sorted(p) != list(range(n)) for p in placement):
        raise ValueError("Invalid placement")
    if set(times) != set(range(len(gates))) or any(type(t) is not int or not 0 <= t < horizon for t in times.values()):
        raise ValueError("Missing gate or invalid gate time")
    if not any(t == 0 for t in times.values()) or any(times[u] >= times[v] for u, v in dag):
        raise ValueError("First-slot gate rule or strict precedence violated")
    if len({tuple(s) for s in swaps}) != len(swaps):
        raise ValueError("Duplicate swap")
    if any(type(t) is not int or not 0 < t < horizon - 1 or type(e) is not int
           or not 0 <= e < len(edges) for t, e in swaps):
        raise ValueError("Invalid swap time or edge")
    hardware = {frozenset(e) for e in edges}
    idle_seen = False
    for t in range(horizon):
        used = set()
        for g, when in times.items():
            if when == t:
                p, q = gates[g]
                pair = {placement[t][p], placement[t][q]}
                if frozenset(pair) not in hardware or used & pair:
                    raise ValueError("Nonadjacent gate or gate conflict")
                used.update(pair)
        after = placement[t][:]
        for when, e in swaps:
            if when == t:
                i, j = edges[e]
                if used.intersection((i, j)):
                    raise ValueError("Swap resource conflict")
                used.update((i, j))
                after = [j if pos == i else i if pos == j else pos for pos in after]
        if t < horizon - 1 and after != placement[t + 1]:
            raise ValueError("Invalid placement transition")
        if used and idle_seen:
            raise ValueError("Activity is not a prefix")
        idle_seen |= not used
    if abs(len(swaps) - objective) > 1e-5:
        raise ValueError("Reported objective differs from swap count")
    return True


def finite(value):
    return float(value) if math.isfinite(value) and abs(value) < 1e90 else None


def model_result(model):
    statuses = {GRB.OPTIMAL: "optimal", GRB.INFEASIBLE: "infeasible",
                GRB.TIME_LIMIT: "time_limit", GRB.INTERRUPTED: "interrupted",
                GRB.INF_OR_UNBD: "infeasible_or_unbounded"}
    return {
        "status": statuses.get(model.Status, f"gurobi_status_{model.Status}"),
        "objective": finite(model.ObjVal) if model.SolCount else None,
        "bound": finite(model.ObjBound) if model.IsMIP else
                 finite(model.ObjVal) if model.Status == GRB.OPTIMAL else None,
        "solve_seconds": model.Runtime, "work": model.Work,
        "nodes": model.NodeCount if model.IsMIP else 0,
        "gap": finite(model.MIPGap) if model.IsMIP and model.SolCount else None,
        "variables": model.NumVars, "binaries": model.NumBinVars,
        "constraints": model.NumConstrs, "nonzeros": model.NumNZs,
    }
