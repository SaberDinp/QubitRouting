# General-hardware QV routing benchmarks

These 24 instances reuse all four existing six-qubit QV circuits and all four
existing eight-qubit QV circuits. Each circuit's gate pairs and gate-dependency
DAG are preserved exactly; only its undirected hardware graph changes.

| Circuit suite | Hardware variants | Instances per topology |
|---|---|---:|
| QV6, depth 6, 18 gates | `cycle6`, `star6`, `grid2x3` | 4 |
| QV8, depth 8, 32 gates | `cycle8`, `grid2x4`, `tree8` | 4 |

A cycle adds the edge `(0, n-1)` to the numbered path. The six-node star has
center 0. Grid nodes are numbered in row order, with horizontal edges in each
row and vertical edges between corresponding columns. The eight-node branching
tree has edges `(0,1), (0,2), (1,3), (1,4), (2,5), (2,6), (6,7)`.

Each file identifies its source circuit and hardware topology in comments.
All files use `optimal_swaps = -1`, meaning **unknown**, and their filenames have
no optimum suffix. The original path optimum is not a reference value for a
different graph. This field must never be used as a bound or stopping target.

The input format does not specify a scheduling horizon. A feasible horizon on
a path need not be feasible on another topology; in particular, every operation
on a star uses its center. Any experiment must report its horizon separately.
Input layers impose no additional barriers beyond the supplied DAG.

From the repository root:

```bash
# Generate missing files; refuse to overwrite differing existing files.
python benchmark/generate_general.py

# Verify reproduction without writing anything.
python benchmark/generate_general.py --check
```

The generator depends only on the Python standard library and the repository's
instance loader. It validates all existing output files before creating any
missing ones. This dataset construction does not solve the routing problems.
