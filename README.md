# QubitRouting

Exact methods for **placing logical qubits on hardware, inserting SWAPs, and
scheduling two-qubit gates**. The aim is to execute every gate on a hardware
edge while respecting the circuit's dependencies and minimizing routing cost.

The repository contains our compact mixed-integer formulation, an independent
implementation of Luteberget et al.'s branch-and-bound method, and a SWAP-count
adaptation of Nannicini et al.'s network-flow formulation. Small routing
instances are included so each method can be run from the command line.

## The methods in a nutshell

| Method | Main idea | Scheduling and objective | Requirements |
|---|---|---|---|
| [Our formulation](our_formulation/COMPACT_FORMULATION.md) | Binary qubit placements and gate times, with continuous reciprocal movement flows representing SWAPs. A heuristic supplies a safe time horizon and a feasible schedule; a separate small-subset model supplies certified lower bounds. | Minimize SWAPs with free gate scheduling subject to the supplied dependency DAG. The default horizon preserves the unrestricted SWAP optimum. | Python and Gurobi |
| [Luteberget](luteberget/README.md) | Search partial gate/SWAP schedules, assign qubits lazily, and prune using admissible bounds and dominance. | Default `qmrs` mode matches our deadline and operation conventions. `paper` mode also supports duration and weighted objectives. | Python; Gurobi only for optional legacy comparison |
| [Nannicini adaptation](nannicini_swap_bip/README.md) | A binary network-flow model tracks qubits through fixed gate layers and intervening routing slots. | Minimize SWAPs for deterministic greedy gate groups; enough dummy slots allow rerouting between groups. | Python and Gurobi |

The Luteberget implementation follows Section 3 of
[*An Exact Branch and Bound Algorithm for the generalized Qubit Mapping Problem*](https://arxiv.org/abs/2508.21718v1).
It is a reimplementation, not the authors' original code. The Nannicini package
adapts [*Optimal qubit assignment and routing via integer programming*](https://arxiv.org/abs/2106.06446).
Unlike the source paper's gate/SWAP merging convention, our adaptation forbids
swapping a gate participant in the same slot and optimizes only SWAP count.

**Compare like problems.** Our formulation and Luteberget's `qmrs` mode allow
gates to move across input layers as long as the DAG is respected. The Nannicini
adaptation fixes simultaneous layers before solving and has a different total
horizon. Its optimum is for those fixed groups, so a larger SWAP count does not
by itself indicate a worse solver. For a matching deadline, give Luteberget's
`qmrs` mode the `horizon` reported by our solver. With unit durations and SWAP
minimization, its unrestricted `paper` mode also targets the unrestricted SWAP
optimum. A smaller explicit deadline can change that optimum.

## Setup

Use **Python 3.10 or newer**. From a terminal:

```bash
git clone https://github.com/SaberDinp/QubitRouting.git
cd QubitRouting
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

On Windows, activate with `.venv\Scripts\activate` instead. The mixed-integer
models require a working Gurobi license; installing `gurobipy` alone does not
provide an unrestricted license. Luteberget's standalone search uses only the
Python standard library and can run without installing the Gurobi dependency.

Run all commands below from the repository root with the environment activated.

## Run an instance

Start with the included four-qubit instance:

```bash
# Our formulation: automatic safe horizon, 10-second total solve budget.
python -m our_formulation benchmark/instance_01_opt2.txt \
  --limit 10 --symmetry 2 --output results/ours.jsonl

# Luteberget: SWAP minimization without a deadline.
python -m luteberget benchmark/instance_01_opt2.txt \
  --mode paper --time-limit 10 --output results/luteberget.jsonl

# Nannicini adaptation: fixed layers, automatic dummy-slot allowance.
python -m nannicini_swap_bip benchmark/instance_01_opt2.txt \
  --time-limit 10 --output results/nannicini.jsonl
```

Each command prints a result summary and writes JSON Lines (one JSON object per
result) at the requested path. Our and Nannicini's commands require a fresh
output filename; Luteberget appends to an existing file. Nannicini also writes a
Markdown summary and solver logs next to its report. Generated results are
ignored by Git.

To try a larger circuit or a different hardware graph:

```bash
python -m our_formulation benchmark/qv6_depth6_01_opt7.txt \
  --limit 30
python -m our_formulation benchmark/qv6_depth6_01_grid2x3.txt \
  --limit 30
```

Our default scheduling horizon is **`T = number of gates + heuristic SWAP count`**.
The heuristic first constructs and independently validates a feasible route.
Any route with no more SWAPs can be serialized within this many unit time slots,
so the cap cannot exclude a better SWAP solution. A closed optimality
certificate therefore establishes the unrestricted minimum SWAP count under
our operation conventions. `--horizon T` instead requests a deadline problem;
its certificate is reported as applying to that fixed horizon.

`--symmetry 2` is the default and selects Gurobi's maximum, aggressive symmetry
detection. Use `--symmetry -1` to compare with Gurobi's automatic setting;
symmetry handling changes search behavior, not the routing problem.
See [Gurobi's Symmetry parameter](https://docs.gurobi.com/projects/optimizer/en/current/reference/parameters.html#parameter.Symmetry).

Useful options for our formulation are `--coverage` (an additional interaction
coverage inequality), `--plain` (disable the optional beam heuristic and subset
bound, while retaining safe horizon construction), and
`--diagnose-lp` (solve diagnostic LP relaxations after the production run).
Diagnostics take additional time outside `--limit`. All entry points support
`--help`; see the method READMEs for the Python APIs and additional options.

The methods use different time-limit scopes. Our `--limit` includes model
construction, preprocessing, the heuristic, optimization, and replay, with
possible small solver overshoot. Luteberget limits its preprocessing/search;
Nannicini limits Gurobi optimization. Use the recorded total times for
end-to-end comparisons.

## Input format

An instance describes a simple undirected hardware graph and a circuit of
two-qubit interactions. There are `N` logical qubits and `N` physical nodes,
both numbered `0` through `N-1`. For example:

```text
N = 4
optimal_swaps = -1

G =
0 1
1 2
2 3
0 3
0 2

gate_dag =
0 1
1 2
2 3
3 4

hw_edges =
0 1
1 2
2 3
```

- `G`: one logical-qubit pair per gate. Gate IDs are zero-based line indices;
  repeated pairs are distinct gates.
- `gate_dag`: pairs of **gate IDs**. `g h` requires gate `g` to complete before
  gate `h` starts. List order adds no dependencies of its own.
- `hw_edges`: physical-node pairs defining allowed interactions and SWAPs.
  List each undirected edge once.
- `optimal_swaps`: required reference metadata; use `-1` if unknown. Solvers
  do not use this value as a bound, target, or initial solution.

Blank lines and `#` comment lines are accepted. The file has no scheduling
horizon; our solver computes a safe one unless an explicit deadline is supplied.
The [benchmark guide](benchmark/README.md) describes the included circuits of at
most six qubits, `generate_quantum_volume.py` for new circuit structures, and
`generate_hardware_variants.py` for alternative hardware topologies. These are routing structures,
not complete quantum circuits with gate matrices or measurement instructions.

Our formulation and the legacy model use unit-duration gates and SWAPs, free
initial placement, disjoint operations within each slot, at least one gate in
slot zero, and no SWAPs in the first or final slot. Our formulation supports
arbitrary simple undirected graphs, including disconnected hardware; a circuit
can be infeasible for the chosen graph or deadline. The Nannicini package
requires connected hardware and uses its own fixed-layer time convention.

## Read the output

A feasible schedule gives an **upper bound** on the minimum SWAP count. A
certified `bound` gives a **lower bound**. Optimality requires a proof that they
meet, or the solver's corresponding optimal status. A time-limited run with no
schedule does not prove infeasibility. An `optimal_swaps` filename or metadata
label does not guarantee feasibility at your chosen deadline.

Our JSON reports include `status`, `objective`, `bound`, `proven_optimal`,
`globally_optimal`, `optimality_scope`, `horizon`, `heuristic_swap_count`,
`production_seconds`, and the independently replayed `schedule`.
`proven_optimal` means the reported problem was solved to optimality;
`globally_optimal` additionally requires the automatic safe horizon.
An explicit deadline is conservatively marked `fixed_horizon` even if it
happens to be large enough. In our
schedules and Luteberget's `qmrs` schedules, `placement[t][q]` gives the physical
node of logical qubit `q` before slot `t`, `gate_times` maps gate IDs to slots, and each SWAP is
`[slot, edge_index]`. Nannicini uses a list of per-slot records with physical-node
SWAP pairs; its README documents that schema. Exact algorithms can still take
exponential time, so small examples are the best starting point.

## Repository layout

```text
our_formulation/        Compact formulation, bound/heuristic helpers, math notes
luteberget/             Independent branch-and-bound implementation and guide
nannicini_swap_bip/     Fixed-layer flow BIP, validation, and formulation notes
benchmark/             Input instances and reproducible instance generators
tests/                 Regression checks for our formulation and helpers
load_instance.py       Shared instance reader
legacy_model.py        Earlier time-indexed Gurobi model, retained for reference
requirements.txt       Gurobi Python dependency
```

`legacy_model.py` retains its source-configured instance and horizon; run it
with `python legacy_model.py` or use Luteberget's optional comparison adapter.
Use the package CLIs for new runs. To compare with it through the adapter:

```bash
python -m luteberget benchmark/instance_01_opt2.txt \
  --horizon 12 --time-limit 10 --compare-legacy
```

Run the regression checks from the repository root:

```bash
python -m unittest discover -v
```

The suite includes Gurobi-backed checks, which need a working license.

Development experiments, presentation sources/exports, and generated solver
artifacts remain local and are excluded through `.gitignore`. They are not
required to run the published methods.
