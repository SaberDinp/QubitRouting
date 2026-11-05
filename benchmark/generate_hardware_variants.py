"""Copy the six-qubit quantum-volume circuits onto cycle, star, and grid hardware.

Only the hardware graph changes. Gate order and every supplied dependency are
preserved, and path-specific optimum labels are discarded. Run with --check
to verify the checked-in files without changing anything.
"""

import argparse
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from load_instance import load_instance


SOURCES = {
    6: (
        "qv6_depth6_01_opt7.txt",
        "qv6_depth6_02_opt8.txt",
        "qv6_depth6_03_opt11.txt",
        "qv6_depth6_04_opt13.txt",
    ),
}


def hardware_graphs(n):
    """Return undirected edges, each represented once, in a fixed order."""
    if n != 6:
        raise ValueError("Only the existing six-qubit circuit suite is supported")
    cycle = [(i, i + 1) for i in range(n - 1)] + [(0, n - 1)]
    columns = n // 2
    grid = ([(row * columns + col, row * columns + col + 1)
             for row in range(2) for col in range(columns - 1)]
            + [(col, col + columns) for col in range(columns)])
    return {
        "cycle6": cycle,
        "star6": [(0, i) for i in range(1, n)],
        "grid2x3": grid,
    }


def expected_instances():
    """Yield filenames and deterministic text for the 12 hardware variants."""
    for n, source_names in SOURCES.items():
        for index, source_name in enumerate(source_names, 1):
            source = ROOT / "benchmark" / source_name
            gates, dag, _, source_n, _ = load_instance(source)
            if source_n != n:
                raise ValueError(f"Unexpected qubit count in {source}")
            for topology, edges in hardware_graphs(n).items():
                lines = [
                    f"# Source circuit: benchmark/{source_name}",
                    "# Gate pairs and gate_dag are copied unchanged from the source.",
                    f"# Hardware topology: {topology}; undirected edges listed once.",
                    "# Generator: benchmark/generate_hardware_variants.py",
                    "# optimal_swaps = -1 means unknown; source path optima do not apply.",
                    f"N = {n}",
                    "optimal_swaps = -1",
                ]
                for heading, pairs in (("G", gates), ("gate_dag", dag),
                                       ("hw_edges", edges)):
                    lines.extend(["", f"{heading} ="])
                    lines.extend(f"{u} {v}" for u, v in pairs)
                filename = f"qv{n}_depth{n}_{index:02}_{topology}.txt"
                yield filename, "\n".join(lines) + "\n"


def write_instances(directory, *, check=False):
    """Validate everything first; never overwrite differing existing files."""
    directory = Path(directory)
    pending = [(directory / name, content)
               for name, content in expected_instances()]
    for path, content in pending:
        if path.exists():
            if path.read_text(encoding="utf-8") != content:
                raise FileExistsError(f"Different existing instance: {path}")
        elif check:
            raise FileNotFoundError(f"Missing generated instance: {path}")
    if not check:
        directory.mkdir(parents=True, exist_ok=True)
        for path, content in pending:
            if not path.exists():
                path.write_text(content, encoding="utf-8")
    return [path for path, _ in pending]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path,
                        default=ROOT / "benchmark")
    parser.add_argument("--check", action="store_true",
                        help="Verify every generated instance; do not write files")
    args = parser.parse_args()
    paths = write_instances(args.output_dir, check=args.check)
    action = "Verified" if args.check else "Generated or retained"
    print(f"{action} {len(paths)} instances in {args.output_dir}")


if __name__ == "__main__":
    main()
