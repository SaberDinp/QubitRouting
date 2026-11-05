"""Compact routing MIP for arbitrary undirected hardware graphs.

Four core variable families: binary placements w and gate times a; continuous
directed movements f and gate allocations x. Reciprocal movement forces SWAPs.
A small subset-movement MIP supplies a stronger certified lower bound by default.
The original DAG and deadline are preserved; input gates are not fixed to layers.
See our_formulation/COMPACT_FORMULATION.md for the equations and correctness argument.
"""

from __future__ import annotations

import argparse
from collections import deque
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
import time

import gurobipy as gp
from gurobipy import GRB

from load_instance import load_instance
from .validation import model_result, validate_schedule


@dataclass
class RoutingModel:
    model: gp.Model
    w: gp.tupledict
    a: gp.tupledict
    f: gp.tupledict
    x: gp.tupledict
    gates: tuple
    dag: tuple
    edges: tuple
    n: int
    horizon: int


def _input_and_windows(gates, dag, edges, n, horizon):
    if type(n) is not int or n < 1:
        raise ValueError("n must be a positive integer")
    if type(horizon) is not int or horizon < 1:
        raise ValueError("horizon must be a positive integer")
    gates, dag, edges = tuple(map(tuple, gates)), tuple(map(tuple, dag)), tuple(map(tuple, edges))
    for label, pairs, size in (("gate", gates, n), ("dependency", dag, len(gates)),
                               ("hardware edge", edges, n)):
        if any(len(p) != 2 or any(type(i) is not int or not 0 <= i < size for i in p)
               or p[0] == p[1] for p in pairs):
            raise ValueError(f"Invalid {label}")
    edges = tuple(tuple(sorted(e)) for e in edges)
    if len(set(edges)) != len(edges) or len(set(dag)) != len(dag):
        raise ValueError("Duplicate hardware edge or dependency")
    successors = [[] for _ in gates]
    indegree = [0] * len(gates)
    for g, h in dag:
        successors[g].append(h)
        indegree[h] += 1
    queue = deque(g for g in range(len(gates)) if indegree[g] == 0)
    order, earliest = [], [0] * len(gates)
    while queue:
        g = queue.popleft()
        order.append(g)
        for h in successors[g]:
            earliest[h] = max(earliest[h], earliest[g] + 1)
            indegree[h] -= 1
            if indegree[h] == 0:
                queue.append(h)
    if len(order) != len(gates):
        raise ValueError("gate_dag must be acyclic")
    latest = [horizon - 1] * len(gates)
    for g in reversed(order):
        for h in successors[g]:
            latest[g] = min(latest[g], latest[h] - 1)
    return gates, dag, edges, earliest, latest


def _validate_symmetry(symmetry):
    if type(symmetry) is not int or symmetry not in (-1, 0, 1, 2):
        raise ValueError("symmetry must be -1 (automatic), 0 (off), 1, or 2 (aggressive)")


def build_model(gates, dag, edges, n, horizon, *, env=None, coverage=False, symmetry=2):
    """Build the MIP without solving it or reading any known optimum.

    Logical and physical qubit counts are both n, as in the legacy model. Disconnected
    graphs are accepted; infeasibility is left to the model. First-slot SWAPs
    are prohibited and at least one gate must execute in that slot. Globally
    idle slots are allowed and removed when extracting a schedule.
    """
    _validate_symmetry(symmetry)
    gates, dag, edges, earliest, latest = _input_and_windows(gates, dag, edges, n, horizon)
    arcs = tuple((i, j) for u, v in edges for i, j in ((u, v), (v, u)))
    neighbors = [set() for _ in range(n)]
    for i, j in arcs:
        neighbors[i].add(j)
    gate_slots = [(t, g) for g in range(len(gates))
                  for t in range(earliest[g], latest[g] + 1)]
    times_of = [[t for t in range(earliest[g], latest[g] + 1)] for g in range(len(gates))]
    gates_of = [[g for g, pair in enumerate(gates) if q in pair] for q in range(n)]
    movement_times = range(1, horizon - 1)

    model = gp.Model("compact_reciprocal_routing", env=env)
    model.Params.Symmetry = symmetry
    w = model.addVars(horizon, n, n, vtype=GRB.BINARY, name="w")
    a = model.addVars(gate_slots, vtype=GRB.BINARY, name="a")
    f = model.addVars([(t, q, i, j) for t in movement_times for q in range(n)
                      for i, j in arcs], ub=1, name="f")
    x = model.addVars([(t, g, i, j) for t, g in gate_slots for i, j in arcs], ub=1, name="x")

    # 1. Initial placement is a bijection. Conservation and reciprocity imply
    #    the same row/column sums at every later time, even in the relaxation.
    for q in range(n):
        model.addConstr(gp.quicksum(w[0, q, i] for i in range(n)) == 1,
                        name=f"initial_logical[{q}]")
    for i in range(n):
        model.addConstr(gp.quicksum(w[0, q, i] for q in range(n)) == 1,
                        name=f"initial_physical[{i}]")

    # 2. Every gate executes once, on a real directed hardware edge.
    for g in range(len(gates)):
        model.addConstr(gp.quicksum(a[t, g] for t in times_of[g]) == 1,
                        name=f"gate_once[{g}]")
    for t, g in gate_slots:
        model.addConstr(gp.quicksum(x[t, g, i, j] for i, j in arcs) == a[t, g],
                        name=f"gate_allocation[{t},{g}]")
    for g, h in dag:
        for t in range(earliest[h], latest[g] + 1):
            model.addConstr(gp.quicksum(a[tau, h] for tau in times_of[h] if tau <= t)
                            <= gp.quicksum(a[tau, g] for tau in times_of[g] if tau < t),
                            name=f"precedence[{t},{g},{h}]")
    model.addConstr(gp.quicksum(a[0, g] for g in range(len(gates)) if (0, g) in a) >= 1,
                    name="first_gate")

    # 3. A qubit can take part in one gate OR move out of its current node.
    #    These aggregated rows also imply all physical/logical resource limits.
    for t in range(horizon):
        for q in range(n):
            for i in range(n):
                gate_use = gp.quicksum(
                    x[t, g, i, j] if gates[g][0] == q else x[t, g, j, i]
                    for g in gates_of[q] if (t, g) in a for j in neighbors[i])
                outgoing = gp.quicksum(f[t, q, i, j] for j in neighbors[i]) if t in movement_times else 0
                model.addConstr(gate_use + outgoing <= w[t, q, i],
                                name=f"capacity[{t},{q},{i}]")
                # 4. Qubit flow between consecutive placements.
                if t < horizon - 1:
                    incoming = gp.quicksum(f[t, q, j, i] for j in neighbors[i]) if t in movement_times else 0
                    model.addConstr(w[t + 1, q, i] == w[t, q, i] - outgoing + incoming,
                                    name=f"conservation[{t},{q},{i}]")

    # 5. Opposite movements along each edge have equal totals: only SWAPs.
    for t in movement_times:
        for i, j in edges:
            model.addConstr(gp.quicksum(f[t, q, i, j] for q in range(n))
                            == gp.quicksum(f[t, q, j, i] for q in range(n)),
                            name=f"reciprocal[{t},{i},{j}]")

    # Optional single general-graph valid inequality, kept separate from core.
    # A swap creates at most |N(i) symmetric_difference N(j)| - 2 new pairs.
    if coverage:
        distinct_pairs = len({frozenset(g) for g in gates})
        model.addConstr(gp.quicksum(
            (len(neighbors[i] ^ neighbors[j]) - 2) * f[t, q, i, j]
            for t in movement_times for i, j in edges for q in range(n))
            >= distinct_pairs - len(edges), name="interaction_coverage")

    # Count one direction per undirected edge: one unit for each SWAP.
    model.setObjective(gp.quicksum(f[t, q, i, j]
                                  for t in movement_times for i, j in edges for q in range(n)),
                       GRB.MINIMIZE)
    model.update()
    return RoutingModel(model, w, a, f, x, gates, dag, edges, n, horizon)


def extract_schedule(problem):
    """Extract and independently replay a schedule, removing globally idle slots."""
    p = problem
    if not p.model.SolCount:
        raise ValueError("No incumbent to extract")
    for family in (p.w, p.a, p.f, p.x):
        if any(abs(v.X - round(v.X)) > 1e-5 for v in family.values()):
            raise ValueError("Incumbent has a fractional placement or operation")
    placements = [[next(i for i in range(p.n) if p.w[t, q, i].X > .5)
                   for q in range(p.n)] for t in range(p.horizon)]
    times = {g: t for (t, g), v in p.a.items() if v.X > .5}
    swaps = [(t, e) for t in range(1, p.horizon - 1) for e, (i, j) in enumerate(p.edges)
             if sum(p.f[t, q, i, j].X for q in range(p.n)) > .5]
    active = sorted(set(times.values()) | {t for t, _ in swaps})
    mapping = {t: k for k, t in enumerate(active)}
    compressed = [placements[t] for t in active]
    # Replay the last active slot to obtain the trailing placement (also valid
    # for a time-limited incumbent containing a redundant final SWAP).
    last = compressed[-1][:]
    for t, e in swaps:
        if t == active[-1]:
            i, j = p.edges[e]
            last = [j if v == i else i if v == j else v for v in last]
    compressed.extend(last[:] for _ in range(p.horizon - len(compressed)))
    schedule = {"placement": compressed, "gate_times": {g: mapping[t] for g, t in times.items()},
                "swaps": [[mapping[t], e] for t, e in swaps]}
    validate_schedule(p.gates, p.dag, p.edges, p.n, p.horizon, schedule, p.model.ObjVal)
    return schedule


def _integer_bound(value):
    """Round a certified lower bound using the integral SWAP objective."""
    return max(0, math.ceil(value - 1e-6)) if value is not None else 0


def solve(gates, dag, edges, n, horizon=None, *, time_limit=10, threads=4, seed=0,
          coverage=False, strengthen=True, diagnose_lp=False, output=False, log_file=None,
          symmetry=2):
    """Cold solve, including all bound preprocessing and heuristic time.

    The default horizon is the gate count plus SWAPs in an independently
    validated serial heuristic route. This retains a global minimum-SWAP
    solution. An explicit horizon instead proves only a fixed-horizon optimum.
    All preprocessing, including the serial router, consumes the time budget.
    The auxiliary MIP is a relaxation: only its dual bound can certify a
    routing lower bound. A feasible routing schedule is independently replayed.
    If these bounds meet, the time-indexed MIP need not be optimized.
    Otherwise it receives the certified bound and the feasible MIP start.
    Optional LP diagnostics run afterward, outside the production budget.
    """
    if not math.isfinite(time_limit) or time_limit <= 0:
        raise ValueError("time_limit must be positive and finite")
    if type(threads) is not int or threads < 1 or type(seed) is not int or seed < 0:
        raise ValueError("threads must be positive and seed nonnegative integers")
    _validate_symmetry(symmetry)
    start = time.perf_counter()
    automatic_horizon = horizon is None
    safe_witness = None
    gates, dag, edges = tuple(map(tuple, gates)), tuple(map(tuple, dag)), tuple(map(tuple, edges))
    if automatic_horizon:
        from .serial_heuristic import construct as serial_construct

        safe_witness = serial_construct(gates, dag, edges, n, deadline=start + time_limit)
        horizon = safe_witness["horizon"]
        if horizon != len(gates) + safe_witness["objective"]:
            raise ValueError("Safe horizon disagrees with heuristic gate and SWAP counts")
        # Validate at the consumer boundary too: no unchecked witness can
        # justify the horizon, supply an incumbent, or establish optimality.
        validate_schedule(gates, dag, edges, n, horizon,
                          safe_witness["schedule"], safe_witness["objective"])
    from .general_beam_start import set_general_start

    with gp.Env(empty=True) as env:
        env.setParam("OutputFlag", int(output or log_file is not None))
        env.setParam("LogToConsole", int(output))
        env.start()
        p = build_model(gates, dag, edges, n, horizon, env=env, coverage=coverage,
                        symmetry=symmetry)
        m = p.model
        bound_lp = None
        try:
            m.Params.Threads = threads
            m.Params.Seed = seed
            m.Params.Symmetry = symmetry
            m.Params.MIPGap = m.Params.MIPGapAbs = 0
            if log_file:
                m.Params.LogFile = str(log_file)
            build_seconds = time.perf_counter() - start
            schedule, objective, lower_bound = None, None, 0
            if safe_witness is not None:
                schedule, objective = safe_witness["schedule"], safe_witness["objective"]
            heuristic, movement_bound = None, None
            infeasible = False
            if strengthen and p.gates and time.perf_counter() < start + time_limit:
                from .general_beam_start import construct
                from .general_subset_bound import build_bound

                heuristic = construct(p.gates, p.dag, p.edges, n, horizon,
                                      width=256, seed=seed,
                                      time_limit=max(.001, min(2.0, .15 * time_limit,
                                                              start + time_limit - time.perf_counter())))
                candidate, candidate_objective = heuristic.pop("schedule"), heuristic["objective"]
                if candidate is not None:
                    validate_schedule(p.gates, p.dag, p.edges, n, horizon, candidate, candidate_objective)
                    if objective is None or candidate_objective < objective:
                        schedule, objective = candidate, candidate_objective
                # Zero is already a universal lower bound; no preprocessing
                # is needed when the heuristic finds a zero-SWAP schedule.
                if objective != 0 and n >= 3 and time.perf_counter() < start + .55 * time_limit:
                    bound_start = time.perf_counter()
                    size = min(5, n - 1) if n <= 6 else 4
                    b, info = build_bound(
                        p.gates, p.dag, p.edges, n, horizon, size=size, env=env,
                        symmetry=symmetry,
                        pair_weights=True, allow_one_anonymous=size == n - 1,
                        deadline=min(bound_start + 10.0, start + .55 * time_limit))
                    try:
                        b.Params.Threads = threads
                        b.Params.Seed = seed
                        b.Params.Symmetry = symmetry
                        b.Params.MIPGap = b.Params.MIPGapAbs = 0
                        if log_file:
                            path = Path(log_file)
                            b.Params.LogFile = str(path.with_name(path.stem + "_bound" + path.suffix))
                        b.Params.TimeLimit = max(.001, min(bound_start + 15.0,
                                                         start + .75 * time_limit)
                                                - time.perf_counter())
                        b.optimize()
                        movement_bound = {**info, **model_result(b), "subset_size": size,
                                          "total_seconds": time.perf_counter() - bound_start,
                                          "iterations": b.IterCount}
                        lower_bound = _integer_bound(movement_bound["bound"])
                        infeasible = b.Status == GRB.INFEASIBLE
                        if diagnose_lp:
                            bound_lp = b.relax()
                    finally:
                        b.dispose()

            if schedule is not None and (infeasible or lower_bound > objective):
                raise RuntimeError("Certified lower bound contradicts the validated schedule")
            certified = schedule is not None and lower_bound == objective
            routing = None
            # Keep the full routing model lean: one cost variable carries the
            # computed global bound; subset tables stay in the auxiliary MIP.
            if lower_bound:
                total = m.addVar(lb=lower_bound, vtype=GRB.INTEGER, name="certified_swap_count")
                m.addConstr(total == m.getObjective(), name="certified_swap_total")
                m.setObjective(total, GRB.MINIMIZE)
                if objective is not None:
                    total.Start = objective
                m.update()
            if not certified and not infeasible and time.perf_counter() < start + time_limit:
                if schedule is not None:
                    set_general_start(p, schedule)
                m.Params.TimeLimit = max(.001, start + time_limit - time.perf_counter())
                m.optimize()
                routing = model_result(m)
                if m.SolCount and (objective is None or m.ObjVal < objective + 1e-6):
                    schedule = extract_schedule(p)
                    objective = len(schedule["swaps"])
                lower_bound = max(lower_bound, _integer_bound(routing["bound"]))
                infeasible = m.Status == GRB.INFEASIBLE
                certified = schedule is not None and lower_bound == objective
            if schedule is not None:
                if infeasible or lower_bound > objective:
                    raise RuntimeError("Routing proof contradicts the validated schedule")
                validate_schedule(p.gates, p.dag, p.edges, n, horizon, schedule, objective)
            result = {
                "status": "optimal" if certified else "infeasible" if infeasible
                          else routing["status"] if routing else "time_limit",
                "objective": objective, "bound": None if infeasible else lower_bound,
                "gap": (objective - lower_bound) / max(1, objective) if objective is not None else None,
                "nodes": (routing["nodes"] if routing else 0)
                         + (movement_bound["nodes"] if movement_bound else 0),
                "solve_seconds": (routing["solve_seconds"] if routing else 0)
                                 + (movement_bound["solve_seconds"] if movement_bound else 0),
                "work": (routing["work"] if routing else 0)
                        + (movement_bound["work"] if movement_bound else 0),
                "variables": m.NumVars, "binaries": m.NumBinVars,
                "constraints": m.NumConstrs, "nonzeros": m.NumNZs,
                "routing": routing, "movement_bound": movement_bound, "heuristic": heuristic,
                "proof_source": "subset_movement_bound" if certified and routing is None and movement_bound
                                else "nonnegative_cost" if certified and routing is None
                                else "routing_mip" if certified or infeasible else None,
            }
            if infeasible and movement_bound and movement_bound["status"] == "infeasible":
                result["proof_source"] = "subset_movement_bound"
            iterations = (m.IterCount if routing else 0) + (movement_bound["iterations"] if movement_bound else 0)
            result.update(build_seconds=build_seconds, production_seconds=time.perf_counter() - start,
                          schedule=schedule, schedule_validated=schedule is not None,
                          horizon=horizon, n=n, gate_count=len(p.gates), hardware_edges=p.edges,
                          horizon_source="serial_heuristic" if automatic_horizon else "explicit",
                          optimality_scope="global" if automatic_horizon else "fixed_horizon",
                          heuristic_swap_count=safe_witness["objective"] if safe_witness else None,
                          heuristic_seconds=safe_witness["seconds"] if safe_witness else None,
                          safe_horizon={key: value for key, value in safe_witness.items() if key != "schedule"}
                                       if safe_witness else None,
                          coverage=coverage, strengthen=strengthen, threads=threads, seed=seed, time_limit=time_limit,
                          symmetry=symmetry,
                          proven_optimal=certified,
                          globally_optimal=certified and automatic_horizon,
                          iterations=iterations, iterations_per_node=iterations / max(1, result["nodes"]),
                          external_start=False, external_bound=False,
                          source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                          helper_sha256={name: hashlib.sha256((Path(__file__).parent / name).read_bytes()).hexdigest()
                                         for name in (("general_subset_bound.py", "general_beam_start.py")
                                                      if strengthen else ())
                                         + (("serial_heuristic.py",) if automatic_horizon else ())},
                          gurobi_version=list(gp.gurobi.version()))
            if diagnose_lp:
                diagnostic_start = time.perf_counter()
                lp = m.relax()
                try:
                    # Report the raw core LP, before the independently proved
                    # scalar bound. Otherwise a displayed LP=7 would just
                    # reflect the imported lower bound of seven.
                    cost = lp.getVarByName("certified_swap_count")
                    if cost is not None:
                        cost.LB = 0
                    lp.Params.OutputFlag = 0
                    lp.Params.Threads = threads
                    lp.Params.TimeLimit = time_limit
                    lp.optimize()
                    result["lp"] = model_result(lp)
                    result["lp_diagnostic_scope"] = "core_before_computed_bound"
                    if bound_lp is not None:
                        bound_lp.Params.OutputFlag = 0
                        bound_lp.Params.Threads = threads
                        bound_lp.Params.TimeLimit = time_limit
                        bound_lp.optimize()
                        movement_bound["lp"] = model_result(bound_lp)
                finally:
                    lp.dispose()
                result["lp_diagnostic_seconds"] = time.perf_counter() - diagnostic_start
            return result
        finally:
            if bound_lp is not None:
                bound_lp.dispose()
            m.dispose()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("instances", nargs="+", type=Path)
    parser.add_argument("--horizon", type=int, default=None,
                        help="explicit slot deadline; default: gates + heuristic SWAPs for global optimality")
    parser.add_argument("--limit", type=float, default=10)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--symmetry", type=int, choices=(-1, 0, 1, 2), default=2,
                        help="Gurobi symmetry detection: -1 automatic, 0 off, 1 conservative, 2 aggressive (default)")
    parser.add_argument("--coverage", action="store_true", help="add the single interaction-coverage inequality")
    parser.add_argument("--plain", action="store_true",
                        help="disable subset bounds and beam search; retain the automatic-horizon serial witness")
    parser.add_argument("--diagnose-lp", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.output and args.output.exists():
        parser.error("Choose a fresh output filename")
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
    for path in args.instances:
        gates, dag, edges, n, _ = load_instance(path)  # Reference optimum is deliberately ignored.
        from .serial_heuristic import SafeHorizonUnavailable

        try:
            result = solve(gates, dag, edges, n, args.horizon, time_limit=args.limit,
                           threads=args.threads, seed=args.seed, coverage=args.coverage,
                           symmetry=args.symmetry,
                           strengthen=not args.plain,
                           diagnose_lp=args.diagnose_lp)
        except SafeHorizonUnavailable as exc:
            parser.error(f"{path}: {exc}")
        result.update(instance=str(path), instance_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        if args.output:
            with args.output.open("a") as stream:
                stream.write(json.dumps(result, allow_nan=False) + "\n")
        print(f"{path.stem}: {result['status']} ({result['optimality_scope']}), "
              f"T={result['horizon']}, symmetry={result['symmetry']}, swaps={result['objective']}, "
              f"bound={result['bound']}, nodes={result['nodes']:g}, "
              f"total={result['production_seconds']:.3f}s", flush=True)


if __name__ == "__main__":
    main()
