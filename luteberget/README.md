# Luteberget branch and bound

Independent Python implementation of Section 3 of
[Luteberget et al., *An Exact Branch and Bound Algorithm for the generalized
Qubit Mapping Problem*, arXiv:2508.21718v1](https://arxiv.org/abs/2508.21718v1).
This is a reimplementation from the paper, not the authors' original software.
It uses the repository's instance format and `load_instance.py`.
The search needs only Python 3.10 or newer; the optional comparison needs
`gurobipy` and a working Gurobi license.

## Run an instance

From the repository root, with your Python environment activated:

```bash
python -m luteberget benchmark/qv6_depth6_01_opt7.txt --horizon 12
```

The number of time steps must be supplied with `--horizon T` in the default
`qmrs` mode; there is no implicit 12-step horizon. Omitting the instance path
selects the same instance. Defaults are swap minimization and a 60-second search
limit. To print the complete placement and schedule:

```bash
python -m luteberget benchmark/qv6_depth6_01_opt7.txt --horizon 12 --show-schedule
```

The alternative entry point `luteberget/run.py` works in an IDE or from another
working directory. Instance paths are resolved against the working directory,
then the project directory, then `benchmark/`.

## Compare with the legacy model

```bash
python -m luteberget benchmark/qv6_depth6_01_opt7.txt \
  --horizon 12 --time-limit 30 --compare-legacy \
  --output luteberget/results/my_comparison.jsonl
```

`--compare-legacy` builds the retained `legacy_model.py` with the chosen
instance and horizon. It sets the Gurobi limit, thread count, seed, and zero relative MIP gap.
The older `--compare-main` spelling remains an alias. Both searches start
independently; neither gets the other's solution or bound. Default Gurobi settings are `--threads 1 --seed 0`.
Changes to the structure of `legacy_model.py` may require updating the adapter in
`compare.py`; it checks that each intended substitution occurs exactly once.

The default `--mode qmrs` matches the constraints in `legacy_model.py`:

- Free initial bijective placement; undirected hardware edges.
- The supplied `gate_dag`, including strict precedence between time slots.
- Unit gate and SWAP durations, with disjoint operations allowed in parallel.
- All gates complete within `T` slots; no swaps at slots `0` or `T-1`.
- At least one gate at slot zero and minimum total SWAP count.

The paper has no fixed deadline by default. Its unrestricted optimum can be
infeasible at your chosen `T`, which is why comparison is restricted to `qmrs`
mode. The benchmark's `optimal_swaps` field is recorded as metadata only; it
never supplies a bound, target, incumbent, or stopping condition.

Run other instances, or a batch, with the same command:

```bash
python -m luteberget benchmark/instance_04_opt4.txt --horizon 12
python -m luteberget benchmark/qv6_depth6_04_opt13.txt --horizon 14
python -m luteberget 'benchmark/*.txt' --horizon 14 --time-limit 60 \
  --compare-legacy --output luteberget/results/all_comparisons.jsonl
```

Output files **append** one JSON object per method and instance. Each successful
run records status, objective, lower bound, gap, SWAP count, makespan, timings,
search statistics, settings, source/instance hashes, and the independently
validated witness. `schedule` in `qmrs` mode uses the shared routing schema:
`placement[t][q]` is the physical position of logical qubit `q`, `gate_times`
maps zero-based gate indices to start slots, and `swaps` contains `[slot, edge_index]`.
`solution` also contains individual operations, including their durations.
Gurobi logs go beside the JSONL file in a corresponding `_logs` folder.

`seconds` measures B&B preprocessing and search, or Gurobi's optimization time.
`total_seconds` also includes reconstruction/validation, and for Gurobi includes
environment startup and model construction. Input loading is recorded separately
as `load_seconds`. The time limit applies to B&B preprocessing/search or to the
Gurobi optimization, respectively; reconstruction and model construction may
add wall time. Prefer total times when comparing end-to-end performance.

`optimal` and `infeasible` are proof statuses. `time_limit` and `node_limit`
retain any incumbent and a valid lower bound, without claiming optimality.
An absent incumbent on a limited run does not prove infeasibility.

## Paper objectives and durations

```bash
# Swap minimization without a deadline
python -m luteberget benchmark/qv6_depth6_01_opt7.txt --mode paper

# Duration minimization, with SWAP taking three times as long as a gate
python -m luteberget benchmark/qv6_depth6_01_opt7.txt \
  --mode paper --objective depth --gate-duration 1 --swap-duration 3

# Weighted objective: 2 * SWAPs + 1 * makespan
python -m luteberget benchmark/instance_04_opt4.txt \
  --mode paper --objective weighted --swap-weight 2 --depth-weight 1

# Optional Section 3.5 restriction on expansion order
python -m luteberget benchmark/instance_04_opt4.txt --mode paper --layered
```

In paper mode, `--horizon T` adds an optional duration deadline. `--layered`
requires scheduling all lower-layer gates before expanding a higher-layer gate,
as in Section 3.5; it does not force all gates in a layer to start simultaneously.
The input files contain no gate durations, so CLI durations are uniform. The
Python API also accepts one positive integer duration per gate. Integer time
units cover the paper's Table 1 durations; zero and fractional durations are not
supported. Single-qubit gates are not represented by this repository's format;
the caller must fold their durations into two-qubit gates if desired.

```python
from load_instance import load_instance
from luteberget import solve

gates, dag, edges, n, _ = load_instance("benchmark/qv6_depth6_01_opt7.txt")
result = solve(gates, dag, edges, n, horizon=12, time_limit=60)
print(result["status"], result["objective"], result["bound"])

# In paper mode, horizon=None (the API default) means no deadline.
result = solve(gates, dag, edges, n, mode="paper", horizon=None,
               gate_durations=[4] * len(gates), swap_duration=15,
               swap_weight=0, depth_weight=1, time_limit=60)
```

## Implementation correspondence

| Paper | Code in `solver.py` |
|---|---|
| Algorithm 1: best-bound tree search | `solve`, heap ordered by the admissible objective bound |
| Algorithm 2: lazy assignment and gate/SWAP expansion | `_Problem.children`; both orientations of an undirected edge are considered |
| Algorithm 3: Pareto insertion | `insert`, `_dominates`; same completed gates and partial assignment, componentwise depth and SWAP-count dominance |
| Equation (6): critical-path tails | `_Problem.delta` |
| Equations (10)–(13): depth bounds | `_Problem.bounds`, including all minimal/chordless hardware paths, not just shortest paths |
| Equation (14): SWAP bound | Current SWAPs plus maximum remaining assigned-pair distance minus one |
| Equation (15): weighted objective | `depth_weight * hD + swap_weight * hS` |
| Section 3.5: layering | `_Problem.layers` and `layered=True` |

The heap uses completed-gate count and depth to break equal-bound ties.
Dominated heap entries are removed lazily. An incumbent prunes nodes that cannot
improve it. Bound evaluation is skipped for the depth component when it has
neither weight nor a deadline. These choices preserve exactness; this code does
not use a beam width or another heuristic cutoff.

Two explicit adaptations handle this repository. First, QMRS mode adds the
deadline and endpoint-slot rules above. Second, the input DAG is authoritative:
list order does not add dependencies. If it has dependencies between disjoint
logical pairs, the state retains the corresponding release times and includes
them in Pareto dominance. For unordered gates on a shared qubit, the depth bound
uses DAG critical-path tails and mandatory ancestor durations. These reduce to
the paper's formulas on conventional circuit DAGs. The state space is still
exponential; large circuits may need a time or node limit.

`validation.py` replays every returned witness independently of the search,
checking placement, adjacency, resources, precedence, durations, deadline, and
reported cost. Initial placement is reconstructed by reversing every SWAP,
including those through positions unassigned at that point in the search.

## Verification

```bash
python -m unittest luteberget.test_solver -v
QMRS_TEST_GUROBI=1 python -m unittest luteberget.test_solver -v
```

Tests compare against an independent exhaustive time-slot scheduler on 36 small
generated instances with three objective/mode combinations, plus unrestricted
routing, unequal durations, disconnected graphs, limits, lazy reconstruction,
longer minimal paths, and independence from instance labels. The optional
integration test fixes the example instance's B&B witness into the retained
`legacy_model.py` variables and checks all original Gurobi constraints.

Generated JSONL reports and logs are local artifacts and are ignored by Git.
Use the commands above to produce fresh results for your own environment.
A benchmark's reference SWAP count does not imply feasibility at every
horizon; report both the horizon and the mode with any comparison.
