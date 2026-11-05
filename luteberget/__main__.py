"""Run python -m luteberget --help from the project directory."""

import argparse
import glob
import hashlib
import json
from pathlib import Path
import platform
import sys
import time

from .solver import solve


ROOT = Path(__file__).resolve().parents[1]


def _instances(patterns):
    paths = set()
    for pattern in patterns:
        matches = []
        for candidate in (Path(pattern), ROOT / pattern, ROOT / "benchmark" / pattern):
            matches = [Path(p).resolve() for p in glob.glob(str(candidate)) if Path(p).is_file()]
            if matches:
                break
        if not matches:
            raise ValueError(f"No instance matches {pattern!r}")
        paths.update(matches)
    return sorted(paths)


def _print_result(path, result):
    value = result["objective"]
    bound = result.get("bound")
    print(f"{path.name}  {result['method']:10} {result['status']:12} "
          f"horizon={result['horizon'] if result['horizon'] is not None else 'unlimited'}  "
          f"objective={value if value is not None else '—'}  "
          f"bound={bound if bound is not None else '—'}  "
          f"swaps={result.get('swap_count')}  makespan={result.get('makespan')}  "
          f"search={result['seconds']:.3f}s  total={result['total_seconds']:.3f}s", flush=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Luteberget et al. exact lazy-assignment branch and bound")
    parser.add_argument("instances", nargs="*", help="Instance paths or quoted globs (default: qv6_depth6_01_opt7.txt)")
    parser.add_argument("--mode", choices=("qmrs", "paper"), default="qmrs",
                        help="qmrs matches legacy_model.py; paper allows durations and other objectives")
    parser.add_argument("--horizon", type=int,
                        help="Number of time steps; required in qmrs mode, optional in paper mode")
    parser.add_argument("--time-limit", type=float, default=60, help="Seconds per search (default: 60)")
    parser.add_argument("--node-limit", type=int, help="Maximum expanded Luteberget nodes")
    parser.add_argument("--objective", choices=("swaps", "depth", "weighted"), default="swaps")
    parser.add_argument("--swap-weight", type=float, default=1, help="Used with --objective weighted")
    parser.add_argument("--depth-weight", type=float, default=1, help="Used with --objective weighted")
    parser.add_argument("--gate-duration", type=int, default=1)
    parser.add_argument("--swap-duration", type=int, default=1)
    parser.add_argument("--layered", action="store_true", help="Paper Section 3.5 expansion restriction")
    parser.add_argument("--compare-legacy", "--compare-main", dest="compare_legacy",
                        action="store_true", help="Also solve legacy_model.py with Gurobi")
    parser.add_argument("--threads", type=int, default=1, help="Gurobi comparison threads (default: 1)")
    parser.add_argument("--seed", type=int, default=0, help="Gurobi comparison seed")
    parser.add_argument("--output", type=Path, help="Append metrics and validated witnesses to JSONL")
    parser.add_argument("--show-schedule", action="store_true")
    args = parser.parse_args(argv)
    if args.mode == "qmrs" and args.horizon is None:
        parser.error("--horizon is required in qmrs mode; specify the number of time steps")
    horizon = args.horizon
    sw, dw = {"swaps": (1, 0), "depth": (0, 1),
              "weighted": (args.swap_weight, args.depth_weight)}[args.objective]
    if args.compare_legacy and args.mode != "qmrs":
        parser.error("--compare-legacy requires --mode qmrs so the two feasible sets and objectives match")
    if args.threads < 1:
        parser.error("--threads must be positive")
    try:
        paths = _instances(args.instances or ["benchmark/qv6_depth6_01_opt7.txt"])
    except ValueError as exc:
        parser.error(str(exc))
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
    # Loading is shared, but each method solves from scratch, with no transferred
    # witness, bound, or use of optimal_swaps as optimization input.
    from load_instance import load_instance
    source_hash = hashlib.sha256(Path(__file__).with_name("solver.py").read_bytes()).hexdigest()
    had_error = False
    for path in paths:
        started = time.perf_counter()
        try:
            gates, dag, edges, n, label = load_instance(path)
            load_seconds = time.perf_counter() - started
            result = solve(gates, dag, edges, n, horizon=horizon, mode=args.mode,
                           gate_duration=args.gate_duration, swap_duration=args.swap_duration,
                           swap_weight=sw, depth_weight=dw, layered=args.layered,
                           time_limit=args.time_limit, node_limit=args.node_limit)
        except (OSError, ValueError) as exc:
            parser.error(f"{path}: {exc}")
        result["source_sha256"] = source_hash

        def record(row):
            row.update(instance=path.name, instance_path=str(path),
                       instance_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                       label_swaps=label, load_seconds=load_seconds,
                       python_version=platform.python_version(), platform=platform.platform())
            if args.output:
                with args.output.open("a", encoding="utf-8") as handle:
                    handle.write(json.dumps(row, allow_nan=False) + "\n")
            if row["status"] != "error":
                _print_result(path, row)

        record(result)
        if args.show_schedule and result["solution"]:
            print("Initial logical -> physical:", result["solution"]["initial_placement"])
            for op in result["solution"]["operations"]:
                name = f"gate {op['gate']} {op['logical']}" if op["kind"] == "gate" else "SWAP"
                print(f"  t={op['start']:3}  {name:20} physical={op['physical']} duration={op['duration']}")
        if args.compare_legacy:
            try:
                from .compare import run_legacy
                log_file = (args.output.parent / (args.output.stem + "_logs") /
                            f"{path.stem}__legacy_model__seed{args.seed}.log") if args.output else None
                comparison = run_legacy(path, horizon, time_limit=args.time_limit,
                                      threads=args.threads, seed=args.seed, log_file=log_file)
                record(comparison)
            except Exception as exc:
                # Keep the already-computed independent B&B result if Gurobi
                # is missing, the license fails, or legacy_model.py changes shape.
                had_error = True
                record({"method": "legacy_model", "status": "error", "error": str(exc)})
                print(f"legacy_model.py comparison failed: {exc}", file=sys.stderr)
    return 1 if had_error else 0


if __name__ == "__main__":
    raise SystemExit(main())
