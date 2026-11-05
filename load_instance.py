from pathlib import Path
import sys

def load_instance(filename):
    lines = Path(filename).read_text(encoding="utf-8").splitlines()

    N = None
    optimal_swaps = None
    sections = {"G": [], "gate_dag": [], "hw_edges": []}
    current = None

    for raw in lines:
        line = raw.strip()
        if not line or line.startswith("#"):
            continue

        if line.startswith("N ="):
            N = int(line.split("=", 1)[1].strip())
            current = None
        elif line.startswith("optimal_swaps ="):
            optimal_swaps = int(line.split("=", 1)[1].strip())
            current = None
        elif line in ("G =", "gate_dag =", "hw_edges ="):
            current = line[:-2].strip()
        else:
            if current is None:
                raise ValueError(f"Unexpected line: {line}")
            u, v = map(int, line.split())
            sections[current].append((u, v))

    if N is None or optimal_swaps is None:
        raise ValueError("Missing N or optimal_swaps.")

    G = sections["G"]
    gate_dag = sections["gate_dag"]
    hw_edges = sections["hw_edges"]

    return G, gate_dag, hw_edges, N, optimal_swaps


def print_as_python(filename):
    G, gate_dag, hw_edges, N, optimal_swaps = load_instance(filename)

    print("G = [")
    for e in G:
        print(f"    {e},")
    print("]\n")

    print("gate_dag = [")
    for e in gate_dag:
        print(f"    {e},")
    print("]\n")

    print(f"hw_edges = {hw_edges}\n")
    print(f"N = {N}")
    print(f"optimal_swaps = {optimal_swaps}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("Usage: python load_instance.py path/to/instance.txt")

    # Ordinary Python variables available here:
    G, gate_dag, hw_edges, N, optimal_swaps = load_instance(sys.argv[1])

    # Print them in Python-assignment form.
    print_as_python(sys.argv[1])
