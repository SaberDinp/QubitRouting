"""Run the supplied benchmark or selected instance files reproducibly."""

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import sys
import time

from load_instance import load_instance
from .exact import exact_layered_optimum
from .solver import solve

ROOT = Path(__file__).resolve().parents[1]


def write_summary(records, output):
    lines = ["# Nannicini SWAP-count BIP benchmark", "",
             "Times are seconds. Gurobi receives no benchmark-label bound or external MIP start.",
             "The oracle is timed separately and runs after optimization.", "",
             "| Instance | Layers | Dummy/gap | Slots | Status | SWAPs | Bound | Gurobi s | Build s | Total s | Oracle |",
             "|---|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|"]
    for r in records:
        show = lambda x: "—" if x is None else str(x)
        lines.append(f"| {Path(r['instance']).stem} | {len(r['layers'])} | "
                     f"{r['dummy_slots_per_gap']} | {r['total_slots']} | {r['status']} | "
                     f"{show(r['objective'])} | {show(r['bound'])} | "
                     f"{r['solve_seconds']:.4f} | {r['build_seconds']:.4f} | "
                     f"{r['total_seconds']:.4f} | {show(r.get('oracle', {}).get('objective'))} |")
    if records:
        r = records[0]
        lines += ["", f"Gurobi {'.'.join(map(str, r['gurobi_version']))}; "
                  f"threads={r['threads']}; seed={r['seed']}; MIPGap=MIPGapAbs=0; "
                  f"time limit={r['time_limit']} seconds per model.", "",
                  f"Environment: {r['platform']}; Python {r['python_version']}.", "",
                  "Total includes environment startup, model construction, optimization, and witness validation; "
                  "it excludes input loading, oracle verification, and report writing. "
                  "The JSONL also records full benchmark wall time per instance."]
    output.with_suffix(".md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("instances", nargs="*", help="Paths or filenames under benchmark/; default all")
    parser.add_argument("--time-limit", type=float, default=60)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--dummy-slots", type=int, default=None,
                        help="Override the safe automatic bound; smaller values restrict rerouting")
    parser.add_argument("--no-compact", action="store_true")
    parser.add_argument("--verify-exact", action="store_true", help="Independent oracle after Gurobi; N <= 7")
    parser.add_argument("--console-log", action="store_true")
    parser.add_argument("--write-models", action="store_true")
    parser.add_argument("--output", type=Path, default=Path(__file__).parent / "results" / "benchmark.jsonl")
    args = parser.parse_args(argv)
    if args.output.exists():
        parser.error(f"Output already exists: {args.output}; choose a fresh --output path")
    paths = [Path(p) if Path(p).is_file() else ROOT / "benchmark" / p for p in args.instances]
    paths = paths or sorted((ROOT / "benchmark").glob("*.txt"))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    artifacts = args.output.parent / (args.output.stem + "_artifacts")
    artifacts.mkdir(parents=True, exist_ok=True)
    source_hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                     for p in sorted(Path(__file__).parent.glob("*.py"))}
    records = []
    with args.output.open("x", encoding="utf-8") as stream:
        for index, path in enumerate(paths):
            start = time.perf_counter()
            gates, dag, edges, n, label = load_instance(path)
            prefix = artifacts / f"{index:02d}_{path.stem}"
            print(f"Solving {path.name}: N={n}, gates={len(gates)}", flush=True)
            result = solve(gates, dag, edges, n, dummy_slots=args.dummy_slots,
                           compact=not args.no_compact, time_limit=args.time_limit,
                           threads=args.threads, seed=args.seed,
                           log_file=prefix.with_suffix(".log"), output=args.console_log,
                           model_file=prefix.with_suffix(".lp") if args.write_models else None)
            result.update({"instance": str(path.resolve()), "benchmark_label": label,
                           "instance_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                           "source_sha256": source_hashes,
                           "timestamp_utc": datetime.now(timezone.utc).isoformat(),
                           "platform": platform.platform(), "machine": platform.machine(),
                           "processor": platform.processor(), "python_version": platform.python_version(),
                           "command": sys.argv, "gates": gates, "dag": dag, "edges": edges})
            if args.verify_exact:
                oracle = exact_layered_optimum(gates, result["layers"], edges, n)
                result["oracle"] = oracle
                if result["unrestricted_rerouting_guaranteed"] and result["proven_optimal"]:
                    if result["objective"] != oracle["objective"]:
                        raise RuntimeError(f"Gurobi/oracle disagreement for {path}")
                if result["status"] == "INFEASIBLE" and result["unrestricted_rerouting_guaranteed"]:
                    if oracle["status"] != "INFEASIBLE":
                        raise RuntimeError(f"False infeasibility for {path}")
            result["benchmark_seconds"] = time.perf_counter() - start
            stream.write(json.dumps(result, allow_nan=False) + "\n")
            stream.flush()
            records.append(result)
            write_summary(records, args.output)
            print(f"  {result['status']}: swaps={result['objective']}, bound={result['bound']}, "
                  f"solve={result['solve_seconds']:.3f}s, total={result['total_seconds']:.3f}s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
