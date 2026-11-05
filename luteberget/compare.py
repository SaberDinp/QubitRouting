"""Run our_formulation/legacy_model.py with controlled data and solver settings.

No snapshot, added cuts, bound from the other method, or warm start is used.
Gurobi is imported only when this optional comparison is requested.
"""

import ast
import hashlib
import math
from pathlib import Path
import time

from .validation import unit_schedule, validate_solution


LEGACY_MODEL = Path(__file__).resolve().parents[1] / "our_formulation" / "legacy_model.py"


def build_legacy_model(instance, horizon, env, source=LEGACY_MODEL):
    """Replace T, the input filename and Gurobi env; defer optimize()."""
    source = Path(source)
    code = source.read_text(encoding="utf-8")
    tree = ast.parse(code, filename=str(source))
    counts = dict(horizon=0, instance=0, environment=0, optimize=0)
    body = []
    for node in tree.body:
        if isinstance(node, ast.Assign):
            names = [t.id for t in node.targets if isinstance(t, ast.Name)]
            if names == ["T"]:
                node.value = ast.Constant(horizon)
                counts["horizon"] += 1
            elif (names == ["model"] and isinstance(node.value, ast.Call)
                  and isinstance(node.value.func, ast.Attribute)
                  and node.value.func.attr == "Model"):
                if any(kw.arg == "env" for kw in node.value.keywords):
                    raise ValueError("legacy_model.py already specifies an environment; update the comparison adapter")
                node.value.keywords.append(ast.keyword(arg="env", value=ast.Name(id="_env", ctx=ast.Load())))
                counts["environment"] += 1
            elif (isinstance(node.value, ast.Call) and isinstance(node.value.func, ast.Attribute)
                  and node.value.func.attr == "load_instance"):
                node.value.args = [ast.Constant(str(Path(instance).resolve()))]
                counts["instance"] += 1
        # The legacy script can wrap optimization in its optional relaxation
        # branch. Defer that whole solve/report block as well as a bare call.
        optimize_calls = [
            child for child in ast.walk(node)
            if isinstance(child, ast.Call)
            and isinstance(child.func, ast.Attribute)
            and isinstance(child.func.value, ast.Name)
            and child.func.value.id == "model"
            and child.func.attr == "optimize"
        ]
        if optimize_calls:
            counts["optimize"] += len(optimize_calls)
            continue
        body.append(node)
    if any(count != 1 for count in counts.values()):
        raise ValueError(f"legacy_model.py structure changed; update comparison adapter (matches: {counts})")
    tree.body = body
    ast.fix_missing_locations(tree)
    scope = {"_env": env, "__name__": "_qmrs_comparison_", "__file__": str(source)}
    exec(compile(tree, str(source), "exec"), scope)
    scope["model"].update()
    scope["source_sha256"] = hashlib.sha256(code.encode()).hexdigest()
    return scope


def _witness(scope):
    n, gates, edges, horizon = (scope[k] for k in ("N", "G", "hw_edges", "T"))
    for name in ("a", "s", "w"):
        if any(abs(v.X - round(v.X)) > 1e-5 for v in scope[name].values()):
            raise ValueError("legacy_model.py returned a nonintegral witness")
    placements = []
    for t in range(horizon):
        row = []
        for q in range(n):
            positions = [i for i in range(n) if scope["w"][t, q, i].X > .5]
            if len(positions) != 1:
                raise ValueError("legacy_model.py returned an invalid placement")
            row.append(positions[0])
        placements.append(row)
    operations = []
    for g, (p, q) in enumerate(gates):
        times = [t for t in range(horizon) if scope["a"][t, g].X > .5]
        if len(times) != 1:
            raise ValueError("legacy_model.py did not schedule each gate exactly once")
        t = times[0]
        operations.append({"kind": "gate", "gate": g, "logical": [p, q],
                           "physical": [placements[t][p], placements[t][q]],
                           "start": t, "duration": 1})
    swaps = []
    for (t, e), var in scope["s"].items():
        if var.X > .5:
            swaps.append([t, e])
            operations.append({"kind": "swap", "physical": list(edges[e]),
                               "start": t, "duration": 1})
    solution = {"initial_placement": placements[0], "final_placement": placements[-1],
                "operations": sorted(operations, key=lambda op: op["start"]),
                "swap_count": len(swaps),
                "makespan": max((op["start"] + 1 for op in operations), default=0)}
    validate_solution(gates, scope["gate_dag"], edges, n, solution,
                      horizon=horizon, qmrs=True)
    schedule = unit_schedule(solution, edges, horizon)
    if schedule["placement"] != placements or abs(len(swaps) - scope["model"].ObjVal) > 1e-5:
        raise ValueError("legacy_model.py witness transitions or objective are inconsistent")
    return solution, schedule


def run_legacy(instance, horizon, *, time_limit=60.0, threads=1, seed=0,
             log_file=None):
    import gurobipy as gp
    from gurobipy import GRB

    started = time.perf_counter()
    with gp.Env(empty=True) as env:
        env.setParam("OutputFlag", 0)
        env.start()
        scope = build_legacy_model(instance, horizon, env)
        model = scope["model"]
        try:
            model.Params.Threads = threads
            model.Params.Seed = seed
            model.Params.TimeLimit = time_limit
            model.Params.MIPGap = 0
            if log_file is not None:
                log_file = Path(log_file)
                log_file.parent.mkdir(parents=True, exist_ok=True)
                model.Params.LogFile = str(log_file)
                model.Params.LogToConsole = 0
                model.Params.OutputFlag = 1
            build_seconds = time.perf_counter() - started
            model.optimize()
            statuses = {GRB.OPTIMAL: "optimal", GRB.INFEASIBLE: "infeasible",
                        GRB.TIME_LIMIT: "time_limit", GRB.INTERRUPTED: "interrupted",
                        GRB.INF_OR_UNBD: "infeasible_or_unbounded"}
            status = statuses.get(model.Status, f"gurobi_status_{model.Status}")
            solution, schedule = _witness(scope) if model.SolCount else (None, None)
            value = model.ObjVal if model.SolCount else None
            bound = model.ObjBound
            bound = bound if math.isfinite(bound) and abs(bound) < 1e90 else None
            return {"method": "legacy_model", "status": status, "gurobi_status": model.Status,
                    "proven": status in ("optimal", "infeasible"), "horizon": horizon,
                    "mode": "qmrs", "objective": value, "bound": bound,
                    "gap": model.MIPGap if model.SolCount else None,
                    "swap_count": None if solution is None else solution["swap_count"],
                    "makespan": None if solution is None else solution["makespan"],
                    "build_seconds": build_seconds, "seconds": model.Runtime,
                    "total_seconds": time.perf_counter() - started,
                    "nodes": model.NodeCount, "work": model.Work,
                    "threads": threads, "seed": seed, "time_limit": time_limit,
                    "gurobi_version": list(gp.gurobi.version()),
                    "source_sha256": scope["source_sha256"],
                    "variables": model.NumVars, "constraints": model.NumConstrs,
                    "solution": solution, "schedule": schedule}
        finally:
            model.dispose()
