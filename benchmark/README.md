# QubitRouting benchmark instances

The 22 circuit files in this directory use at most six qubits. Each file
contains the number of logical/physical qubits, a reference SWAP optimum when
known, two-qubit gate pairs, a gate-dependency DAG, and hardware edges. All
instances live directly in `benchmark/`.

| Circuit family | Hardware | Files | Two-qubit gates per circuit |
|---|---|---:|---:|
| `original_4q_opt3.txt`, `instance_*.txt` | Paths, 4–6 qubits | 6 | Varies |
| `qv6_depth6_*_opt*.txt` | Six-node path | 4 | 18 |
| `qv6_depth6_*_cycle6.txt` | Six-node cycle | 4 | 18 |
| `qv6_depth6_*_star6.txt` | Six-node star | 4 | 18 |
| `qv6_depth6_*_grid2x3.txt` | 2 × 3 grid | 4 | 18 |

`optimal_swaps = -1` means **unknown**. Known values are reference answers,
not optimization constraints or stopping targets. Hardware variants discard
the source path optimum because a different graph can change the answer.

Gates may be rescheduled subject to the DAG. Input layers do not add barriers.
Files do not prescribe a scheduling horizon. For an unrestricted minimum-SWAP
claim, our formulation uses the number of gates plus a validated heuristic
route's SWAP count as its default horizon; an arbitrarily smaller horizon can
exclude an optimal route. See the [repository README](../README.md) and
[formulation](../our_formulation/COMPACT_FORMULATION.md) for solve commands and
precise scheduling conventions.

## Reproduce hardware-topology variants

[`generate_hardware_variants.py`](generate_hardware_variants.py) copies each
of the four existing QV6 circuits onto three hardware graphs. It preserves
every gate pair and dependency, changing only the hardware and optimum label:

- `cycle6` adds `(0, 5)` to the numbered six-node path.
- `star6` connects node 0 to each other node.
- `grid2x3` numbers nodes by row and joins horizontal and vertical neighbors.

Each output file identifies its source circuit and topology in comments.
From the repository root:

```bash
# Create missing variants in benchmark/; keep identical existing files.
python benchmark/generate_hardware_variants.py

# Verify all 12 checked-in variants without changing files.
python benchmark/generate_hardware_variants.py --check
```

The script validates every existing output before creating missing files and
refuses to overwrite different content. It uses only the Python standard
library and the repository's instance loader; it does not solve the instances.

## Generate additional quantum-volume routing circuits

Quantum-volume (QV) routing structures consist of layers of randomly paired
qubits. These files encode the pair interactions and their dependencies;
they omit the random two-qubit unitary matrices and measurements of a full
quantum-volume circuit.

[`generate_quantum_volume.py`](generate_quantum_volume.py) creates new circuits
on a numbered path. Each layer shuffles the qubits with Python's
`random.Random(seed)` and pairs consecutive entries. A gate depends on the
previous gate on each of its two qubits. Repeated pairs are retained.

```bash
# Defaults: 6 qubits, depth 6, seeds 601–604, and outputs/quantum_volume/.
python benchmark/generate_quantum_volume.py

# An explicit destination for a new set of circuits.
python benchmark/generate_quantum_volume.py --qubits 6 --depth 6 \
  --seeds 611 612 --output-dir outputs/quantum_volume_custom
```

The default output directory is ignored by Git. New files use unknown optimum
labels and filenames without an `_opt...` suffix. A manifest records the seeds,
layer membership, generator and instance hashes, and Python version. Different
existing files are never overwritten.

The older, optimum-labelled QV6 path benchmarks are retained unchanged. Their
original generation seeds were not recorded, so this generator does not claim
to recreate them. The hardware-variant generator reproduces their derived
cycle, star, and grid instances exactly.
