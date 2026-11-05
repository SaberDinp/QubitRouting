# Nannicini layered BIP with a SWAP-count objective

This package adapts
[Nannicini et al.'s network-flow BIP](https://arxiv.org/abs/2106.06446)
in Gurobi. It minimizes the number of SWAPs and forbids a SWAP on either qubit of a gate in the same slot. Gates are greedily
grouped into fixed simultaneous layers, respecting the supplied DAG and breaking
ties by gate index. The earlier free-scheduling model is retained as
[`legacy_model.py`](../our_formulation/legacy_model.py).

The model uses binary `w`, `y`, and `x` variables and the objective
`0.5 * sum(x[t,q,i,j] for i != j)`. Gate qubits have only self-loop movement;
idle-qubit flow must be reciprocal on each hardware edge. The full equations,
correctness argument, and adaptations from the paper are in
[FORMULATION.md](FORMULATION.md). The final slot uses ordinary McCormick upper
bounds because outgoing movement variables do not exist there.

On path hardware with `N >= 3`, the default inserts **N dummy slots per gap**.
This allows every placement permutation and a minimum-swap route between any two placements. Five slots do not provide
that guarantee on a six-node path once swaps involving gate qubits are forbidden.
For general connected hardware the conservative bound is `N*(N-1)//2`.
Initial and final placement are free. Its layer-and-gap horizon differs from
our free-scheduling model's automatic `T = number of gates + heuristic SWAP
count` cap; both horizon constructions preserve their respective SWAP optima,
but fixed layers can change the underlying routing problem.

## Run

From the repository root, with your Python environment and Gurobi license configured:

```bash
# A small instance with independent verification after the solve.
python -m nannicini_swap_bip benchmark/instance_01_opt2.txt \
  --time-limit 30 --threads 1 --seed 0 --verify-exact \
  --output nannicini_swap_bip/results/my_benchmark.jsonl

# A single instance; optionally retain the LP model and show Gurobi's log.
python -m nannicini_swap_bip qv6_depth6_01_opt7.txt \
  --time-limit 60 --verify-exact --write-models --console-log \
  --output nannicini_swap_bip/results/my_single_run.jsonl

# Behavioral tests, including an independent finite-slot enumeration.
python -m unittest nannicini_swap_bip.test_solver -v
```

`run.py` is also executable directly from an IDE or another working directory.
Use a fresh output filename: existing reports are protected from overwriting.
The CLI writes JSONL results and complete schedules, a Markdown timing table,
and one Gurobi log per instance. `--write-models` additionally saves LP files.
`--no-compact` disables the paper's within-gap idle-slot symmetry rule.
`--dummy-slots D` explicitly restricts the horizon if `D` is below the safe bound;
the output flags that loss of a general guarantee. This option is unnecessary
for ordinary benchmark runs.

```python
from load_instance import load_instance
from nannicini_swap_bip import solve, build_model

gates, dag, edges, n, reference_label = load_instance(
    "benchmark/qv6_depth6_01_opt7.txt"
)
result = solve(gates, dag, edges, n, time_limit=60, threads=1, seed=0)
print(result["status"], result["objective"], result["solve_seconds"])

# To inspect or customize the actual Gurobi model:
problem = build_model(gates, dag, edges, n)
try:
    problem.model.write("nannicini_swap_bip/model.lp")
finally:
    problem.model.dispose()
```

`schedule[t]["placement"][q]` is the physical location of logical qubit `q`
before slot `t`; swaps are pairs of physical nodes. Empty slots are retained
so every transition can be replayed. Reported solve time is Gurobi `Runtime`.
Total time includes environment startup, construction, optimization, and witness
validation. Exact-oracle verification is recorded separately, after optimization.

## Interpretation and validation

The reported optimum is for the fixed greedy groups and the chosen dummy-slot
allowance. These groups restrict the original DAG-scheduling problem; their
optimum may exceed that of a model with free gate scheduling. Reference
`optimal_swaps` labels are recorded as metadata and never supplied to Gurobi
as bounds or stopping targets.

Every returned incumbent is independently replayed for placement bijection,
gate adjacency, dependencies, resource conflicts, movement, and cost. The
optional `--verify-exact` oracle enumerates full placements for at most seven
qubits, runs after optimization, and is timed separately. All bundled circuits
have at most six qubits. With no instance argument, the CLI selects all text
files directly under `benchmark/`, including the hardware-topology variants.

Generated reports, logs, and witnesses are ignored by Git. Use a fresh output
path to record results in your own environment. A limited solve can return an
incumbent and a lower bound without proving optimality.

## Files

- `layers.py`: input validation, deterministic greedy groups, dummy-slot bounds.
- `solver.py`: binary flow model, solve metrics, witness extraction.
- `validation.py`: independent replay of returned schedules.
- `exact.py`: independent small-instance optimum check.
- `__main__.py` / `run.py`: benchmark and single-instance entry points.
- `test_solver.py`: mathematical and behavioral regression checks.
- `FORMULATION.md`: equations, assumptions, and proofs.
- `results/`: local generated reports, solver logs, and complete witnesses (ignored by Git).
