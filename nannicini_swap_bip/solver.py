"""Direct binary w/y/x implementation of the handwritten flow formulation."""

from dataclasses import dataclass
import math
import time

import gurobipy as gp
from gurobipy import GRB

from .layers import (greedy_layers, insert_dummy_slots, sufficient_dummy_slots,
                     validate_input)
from .validation import validate_schedule


@dataclass
class RoutingModel:
    model: gp.Model
    w: gp.tupledict
    y: gp.tupledict
    x: gp.tupledict
    z: gp.tupledict
    gates: tuple
    dag: tuple
    edges: tuple
    n: int
    layers: tuple
    slots: tuple
    layer_times: tuple
    dummy_slots: int
    sufficient_slots: int
    compact: bool


def build_model(gates, dag, edges, n, *, dummy_slots=None, compact=True, env=None):
    gates, dag, edges = validate_input(gates, dag, edges, n)
    layers = greedy_layers(gates, dag)
    sufficient = sufficient_dummy_slots(edges, n)
    if dummy_slots is None:
        dummy_slots = sufficient
    if type(dummy_slots) is not int or dummy_slots < 0:
        raise ValueError("dummy_slots must be a nonnegative integer or None")
    slots, layer_times, blocks = insert_dummy_slots(layers, dummy_slots)
    T = len(slots)
    arcs = tuple(arc for i, j in edges for arc in ((i, j), (j, i)))
    flow_arcs = arcs + tuple((i, i) for i in range(n))
    neighbors = [tuple(j for u, j in arcs if u == i) for i in range(n)]
    active = [set(q for g in slot for q in gates[g]) for slot in slots]
    model = gp.Model("nannicini_layered_swap_count", env=env)
    w = model.addVars(T, n, n, vtype=GRB.BINARY, name="w")
    y = model.addVars([(t, g, i, j) for t, slot in enumerate(slots)
                       for g in slot for i, j in arcs], vtype=GRB.BINARY, name="y")
    x = model.addVars([(t, q, i, j) for t in range(T - 1) for q in range(n)
                       for i, j in flow_arcs], vtype=GRB.BINARY, name="x")

    for t in range(T):
        for q in range(n):
            model.addConstr(w.sum(t, q, "*") == 1, name=f"qubit[{t},{q}]")
        for i in range(n):
            model.addConstr(w.sum(t, "*", i) == 1, name=f"node[{t},{i}]")
        for g in slots[t]:
            p, q = gates[g]
            model.addConstr(y.sum(t, g, "*", "*") == 1, name=f"gate[{t},{g}]")
            for i, j in arcs:
                yg = y[t, g, i, j]
                # Handwritten strengthened McCormick; final slot has no x.
                p_at_i = x[t, p, i, i] if t < T - 1 else w[t, p, i]
                q_at_j = x[t, q, j, j] if t < T - 1 else w[t, q, j]
                model.addConstr(yg <= p_at_i, name=f"gate_tail[{t},{g},{i},{j}]")
                model.addConstr(yg <= q_at_j, name=f"gate_head[{t},{g},{i},{j}]")
                model.addConstr(yg >= w[t, p, i] + w[t, q, j] - 1,
                                name=f"gate_product[{t},{g},{i},{j}]")

    for t in range(T - 1):
        for q in range(n):
            for i in range(n):
                model.addConstr(x[t, q, i, i] + gp.quicksum(
                    x[t, q, i, j] for j in neighbors[i]) == w[t, q, i],
                    name=f"outflow[{t},{q},{i}]")
                model.addConstr(x[t, q, i, i] + gp.quicksum(
                    x[t, q, j, i] for j in neighbors[i]) == w[t + 1, q, i],
                    name=f"inflow[{t},{q},{i}]")
                if q in active[t]:
                    model.addConstr(x[t, q, i, i] == w[t, q, i],
                                    name=f"gate_stays[{t},{q},{i}]")
        idle = tuple(q for q in range(n) if q not in active[t])
        # One equality per undirected edge suffices (the reverse duplicates it).
        for i, j in edges:
            model.addConstr(gp.quicksum(x[t, q, i, j] for q in idle) ==
                            gp.quicksum(x[t, q, j, i] for q in idle),
                            name=f"reciprocal[{t},{i},{j}]")

    # Paper Eq. (4): remove gaps only WITHIN each consecutive dummy block.
    dummy_times = [t for block in blocks for t in block] if compact else []
    z = model.addVars(dummy_times, vtype=GRB.BINARY, name="dummy_active")
    for t in dummy_times:
        movement = gp.quicksum(x[t, q, i, j] for q in range(n) for i, j in arcs)
        model.addConstr(movement <= 2 * (n // 2) * z[t], name=f"active_upper[{t}]")
        model.addConstr(movement >= 2 * z[t], name=f"active_lower[{t}]")
    if compact:
        for block in blocks:
            for t, u in zip(block, block[1:]):
                model.addConstr(z[t] >= z[u], name=f"compact[{t},{u}]")

    model.setObjective(0.5 * gp.quicksum(x[t, q, i, j] for t in range(T - 1)
                                       for q in range(n) for i, j in arcs), GRB.MINIMIZE)
    model.update()
    return RoutingModel(model, w, y, x, z, gates, dag, edges, n, layers, slots,
                        layer_times, dummy_slots, sufficient, compact)


def extract_schedule(problem):
    w, x = problem.w, problem.x
    steps = []
    for t, gates in enumerate(problem.slots):
        placement = [next(i for i in range(problem.n) if w[t, q, i].X > 0.5)
                     for q in range(problem.n)]
        swaps = []
        if t < len(problem.slots) - 1:
            swaps = [(i, j) for i, j in problem.edges
                     if sum(x[t, q, i, j].X for q in range(problem.n)) > 0.5]
        steps.append({"t": t, "gates": list(gates), "swaps": swaps,
                      "placement": placement})
    return steps


def solve(gates, dag, edges, n, *, dummy_slots=None, compact=True,
          time_limit=60, threads=1, seed=0, log_file=None, output=False,
          model_file=None):
    """Solve without using instance labels, external bounds, or MIP starts.

    `time_limit` applies to Gurobi optimization, not model construction.
    An explicitly smaller dummy_slots value restricts the routing horizon.
    """
    if time_limit is not None and (not math.isfinite(time_limit) or time_limit <= 0):
        raise ValueError("time_limit must be positive and finite, or None")
    if type(threads) is not int or threads < 1 or type(seed) is not int or seed < 0:
        raise ValueError("threads must be positive and seed nonnegative integers")
    start = time.perf_counter()
    with gp.Env(empty=True) as env:
        env.setParam("OutputFlag", 1 if output or log_file else 0)
        env.setParam("LogToConsole", int(output))
        env.start()
        problem = build_model(gates, dag, edges, n, dummy_slots=dummy_slots,
                              compact=compact, env=env)
        model = problem.model
        try:
            model.Params.Threads = threads
            model.Params.Seed = seed
            model.Params.MIPGap = 0
            model.Params.MIPGapAbs = 0
            if time_limit is not None:
                model.Params.TimeLimit = time_limit
            if log_file:
                model.Params.LogFile = str(log_file)
            if model_file:
                model.write(str(model_file))
            build_seconds = time.perf_counter() - start
            model.optimize()
            statuses = {getattr(GRB, name): name for name in (
                "LOADED", "OPTIMAL", "INFEASIBLE", "INF_OR_UNBD", "UNBOUNDED",
                "CUTOFF", "ITERATION_LIMIT", "NODE_LIMIT", "TIME_LIMIT",
                "SOLUTION_LIMIT", "INTERRUPTED", "NUMERIC", "SUBOPTIMAL",
                "INPROGRESS", "USER_OBJ_LIMIT", "WORK_LIMIT", "MEM_LIMIT")}
            objective = int(round(model.ObjVal)) if model.SolCount else None
            schedule = extract_schedule(problem) if model.SolCount else None
            if schedule is not None:
                validate_schedule(problem.gates, problem.dag, problem.edges, n,
                                  problem.layers, problem.slots, schedule, objective)
            bound = model.ObjBound
            result = {
                "status": statuses.get(model.Status, str(model.Status)),
                "objective": objective,
                "bound": bound if math.isfinite(bound) and abs(bound) < GRB.INFINITY else None,
                "gap": model.MIPGap if model.SolCount and model.IsMIP else None,
                "proven_optimal": model.Status == GRB.OPTIMAL,
                "build_seconds": build_seconds,
                "solve_seconds": model.Runtime,
                "total_seconds": time.perf_counter() - start,
                "work": model.Work, "nodes": model.NodeCount,
                "variables": model.NumVars, "binaries": model.NumBinVars,
                "constraints": model.NumConstrs, "nonzeros": model.NumNZs,
                "gurobi_version": list(gp.gurobi.version()),
                "threads": threads, "seed": seed, "time_limit": time_limit,
                "compact_dummy_blocks": compact,
                "n": n, "gate_count": len(problem.gates),
                "layers": problem.layers, "layer_times": problem.layer_times,
                "dummy_slots_per_gap": problem.dummy_slots,
                "sufficient_dummy_slots_per_gap": problem.sufficient_slots,
                "unrestricted_rerouting_guaranteed": (
                    problem.dummy_slots >= problem.sufficient_slots or len(problem.layers) <= 1),
                "total_slots": len(problem.slots),
                "schedule": schedule, "schedule_validated": schedule is not None,
            }
            return result
        finally:
            model.dispose()
