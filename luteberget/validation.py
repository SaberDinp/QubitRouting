"""Replay witnesses independently of search nodes, bounds, and dominance."""


def validate_solution(gates, dag, edges, n, solution, durations=None,
                      swap_duration=1, *, horizon=None, qmrs=False):
    durations = tuple(durations) if durations is not None else (1,) * len(gates)
    edges = {frozenset(e) for e in edges}
    placement = list(solution["initial_placement"])
    if sorted(placement) != list(range(n)):
        raise ValueError("Initial placement is not bijective")
    free = [0] * n
    gate_times = {}
    gate_ends = {}
    swap_count = 0
    operations = sorted(solution["operations"], key=lambda op: op["start"])
    for op in operations:
        u, v = op["physical"]
        start, duration = op["start"], op["duration"]
        if not isinstance(start, int) or start < 0 or frozenset((u, v)) not in edges:
            raise ValueError("Invalid start time or hardware edge")
        if start < max(free[u], free[v]):
            raise ValueError("Overlapping operations on a physical qubit")
        end = start + duration
        if horizon is not None and end > horizon:
            raise ValueError("Operation exceeds horizon")
        free[u] = free[v] = end
        if op["kind"] == "gate":
            g = op["gate"]
            if g not in range(len(gates)) or g in gate_times or duration != durations[g]:
                raise ValueError("Duplicate/invalid gate or duration")
            p, q = gates[g]
            if [placement[p], placement[q]] != [u, v]:
                raise ValueError("Gate executed on the wrong logical qubits")
            gate_times[g], gate_ends[g] = start, end
        elif op["kind"] == "swap":
            if duration != swap_duration:
                raise ValueError("Incorrect SWAP duration")
            if qmrs and (start == 0 or start >= horizon - 1):
                raise ValueError("QMRS forbids swaps in the first and last slots")
            # Future operations on either site cannot start before completion,
            # so updating immediately also works for unequal durations.
            placement = [v if pos == u else u if pos == v else pos for pos in placement]
            swap_count += 1
        else:
            raise ValueError("Unknown operation kind")
    if set(gate_times) != set(range(len(gates))):
        raise ValueError("Missing circuit gate")
    if any(gate_ends[u] > gate_times[v] for u, v in dag):
        raise ValueError("DAG precedence violated")
    if qmrs and 0 not in gate_times.values():
        raise ValueError("QMRS requires a gate at time zero")
    if swap_count != solution["swap_count"] or max(free) != solution["makespan"]:
        raise ValueError("Incorrect reported cost or makespan")
    if placement != solution["final_placement"]:
        raise ValueError("Incorrect final placement")
    return True


def unit_schedule(solution, edges, horizon):
    """Export the shared placement/gate_times/swaps schedule schema."""
    ids = {frozenset(e): i for i, e in enumerate(edges)}
    gate_times = {op["gate"]: op["start"] for op in solution["operations"] if op["kind"] == "gate"}
    swaps = [[op["start"], ids[frozenset(op["physical"])]]
             for op in solution["operations"] if op["kind"] == "swap"]
    placement = []
    current = list(solution["initial_placement"])
    for t in range(horizon):
        placement.append(current[:])
        for start, e in swaps:
            if start == t:
                u, v = edges[e]
                current = [v if pos == u else u if pos == v else pos for pos in current]
    return {"placement": placement, "gate_times": gate_times, "swaps": swaps}
