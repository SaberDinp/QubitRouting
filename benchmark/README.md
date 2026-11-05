# QubitRouting benchmark instances

Each text file contains the number of logical/physical qubits, a reference SWAP
optimum when known, two-qubit gate pairs, a gate-dependency DAG, and hardware
edges. The original files in this directory use path hardware; the
[`general/`](general/README.md) suite adds cycles, grids, stars, and a branching
tree. Gates may be rescheduled subject to the DAG; grouping gates into input
layers does not impose a barrier between layers.

`optimal_swaps = -1` means **unknown**, not a negative optimum or a lower bound.
Do not use this field as an optimization constraint or stopping target. A known
SWAP optimum also does not guarantee feasibility at a particular time horizon.

## Eight-qubit quantum-volume routing instances

Four reproducible instances were added:

| File | Seed | Qubits | Input depth | Two-qubit gates |
|---|---:|---:|---:|---:|
| `qv8_depth8_01.txt` | 801 | 8 | 8 | 32 |
| `qv8_depth8_02.txt` | 802 | 8 | 8 | 32 |
| `qv8_depth8_03.txt` | 803 | 8 | 8 | 32 |
| `qv8_depth8_04.txt` | 804 | 8 | 8 | 32 |

Like the existing six-qubit/depth-six data, these are square QV **routing
structures**. Each input layer independently shuffles the qubits and pairs
consecutive entries, giving four disjoint two-qubit interactions. Each gate
depends on the immediately previous gate on each of its qubits, with duplicate
dependencies removed. Repeated pairs in different layers are retained. All
four seeds are included without selection by difficulty or solver outcome.

Actual quantum-volume circuits also specify random two-qubit SU(4) matrices;
those matrices and measurements are outside this repository's routing format.
See the [IBM definition](https://quantum.cloud.ibm.com/docs/en/api/qiskit/qiskit.circuit.library.quantum_volume).
The generator uses Python's `random.Random`, so these seed numbers do not claim
to reproduce Qiskit's complete seeded circuits.

```bash
# Reproduce the checked-in instances and their manifest; identical files are kept.
python benchmark/generate_qv.py

# Optional depth-six controls in a separate directory.
python benchmark/generate_qv.py --qubits 8 --depth 6 \
  --output-dir /tmp/qv8_depth6
```

[The manifest](qv8_depth8_manifest.json) records seeds, layer membership, source
and instance hashes, and the generator's Python version. The generator refuses
to overwrite different existing files. Exact optima remain unknown until
proved; there is deliberately no `_opt...` suffix on these filenames.

The files do not set a deadline. Choose and report the horizon when using
our formulation or Luteberget's `qmrs` mode; a longer circuit or a different
hardware topology can require more slots. See the
[repository README](../README.md) for solve commands and the
[compact formulation](../our_formulation/COMPACT_FORMULATION.md) for the model's
precise scheduling conventions.
