"""Generate reproducible quantum-volume routing circuits on path hardware.

The file format stores the routing structure only, not Haar-random SU(4)
matrices or measurement results. Each layer uses a uniformly shuffled pairing.
Dependencies link consecutive gates on each logical qubit, without barriers
between layers. The defaults create new six-qubit circuits in outputs/ rather
than replacing the older, optimum-labelled QV6 benchmark circuits.
"""

import argparse
import hashlib
import json
from pathlib import Path
import platform
import random


def generate(n=6, depth=6, seed=601):
    if type(n) is not int or n < 2 or n % 2:
        raise ValueError("n must be an even integer >= 2")
    if type(depth) is not int or depth < 1 or type(seed) is not int or seed < 0:
        raise ValueError("depth must be positive and seed a nonnegative integer")
    rng = random.Random(seed)
    gates, dag, layers = [], [], []
    previous = [None] * n
    for _ in range(depth):
        permutation = list(range(n))
        rng.shuffle(permutation)
        pairs = sorted(tuple(sorted(permutation[i:i + 2])) for i in range(0, n, 2))
        layer = []
        for p, q in pairs:
            g = len(gates)
            gates.append((p, q))
            dag.extend((parent, g) for parent in sorted(
                {previous[p], previous[q]} - {None}))
            previous[p] = previous[q] = g
            layer.append(g)
        layers.append(layer)
    return gates, dag, [(i, i + 1) for i in range(n - 1)], layers


def instance_text(n, depth, seed):
    gates, dag, edges, layers = generate(n, depth, seed)
    lines = [f"# QV routing structure: qubits={n}, depth={depth}, seed={seed}",
             "# Generator: benchmark/generate_quantum_volume.py; random.Random(seed).shuffle per layer.",
             "# optimal_swaps = -1 denotes an unknown optimum, not a lower bound.",
             "# No unitary matrices or fixed-layer scheduling barriers are encoded.",
             f"N = {n}", "optimal_swaps = -1"]
    for title, pairs in (("G", gates), ("gate_dag", dag), ("hw_edges", edges)):
        lines += ["", f"{title} ="] + [f"{u} {v}" for u, v in pairs]
    return "\n".join(lines) + "\n", layers


def write_instances(directory, n=6, depth=6, seeds=(601, 602, 603, 604)):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    entries = []
    files = []
    for index, seed in enumerate(seeds, 1):
        content, layers = instance_text(n, depth, seed)
        path = directory / f"qv{n}_depth{depth}_{index:02}.txt"
        if path.exists() and path.read_text(encoding="utf-8") != content:
            raise FileExistsError(f"Refusing to overwrite different existing instance: {path}")
        files.append((path, content))
        entries.append({"file": path.name, "seed": seed, "n": n, "depth": depth,
                        "gate_count": n * depth // 2, "layers": layers,
                        "optimal_swaps": None, "sha256": hashlib.sha256(content.encode()).hexdigest()})
    manifest = {"generator": "benchmark/generate_quantum_volume.py", "python_version": platform.python_version(),
                "generator_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                "hardware": "numbered path", "unknown_optimum_sentinel": -1,
                "instances": entries}
    manifest_path = directory / f"qv{n}_depth{depth}_manifest.json"
    if manifest_path.exists():
        previous = json.loads(manifest_path.read_text())
        if previous["instances"] != entries:
            raise FileExistsError(f"Refusing to overwrite different manifest: {manifest_path}")
    for path, content in files:
        if not path.exists():
            path.write_text(content, encoding="utf-8")
    if not manifest_path.exists():
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return [p for p, _ in files]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--qubits", type=int, default=6)
    parser.add_argument("--depth", type=int, default=6)
    parser.add_argument("--seeds", type=int, nargs="+", default=[601, 602, 603, 604])
    parser.add_argument("--output-dir", type=Path,
                        default=Path(__file__).resolve().parents[1] / "outputs" / "quantum_volume")
    args = parser.parse_args()
    for path in write_instances(args.output_dir, args.qubits, args.depth, args.seeds):
        print(path)


if __name__ == "__main__":
    main()
