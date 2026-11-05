# Our compact routing MIP for arbitrary hardware topologies

The model is implemented in
[`compact_formulation.py`](compact_formulation.py), with a runnable entry point
at `python -m our_formulation`. It accepts **any simple undirected hardware
graph**, including cycles, branching trees, grids, and disconnected graphs.
Initial placement and gate times are decisions. The original gate DAG is
preserved; input layers do not become scheduling barriers. By default, a
validated heuristic supplies a horizon large enough to preserve the
unrestricted minimum SWAP count. An explicit deadline remains optional.

The core idea is to **place qubits, execute gates on edges, and move qubits in
opposite pairs along edges**. Opposite movements are precisely SWAPs. Four
variable families define the routing MIP. Its LP relaxation can be weak, so
the default solver also computes a lower bound on the movement required by
small groups of logical qubits, using a separate auxiliary MIP.

The earlier model is retained as [`legacy_model.py`](../legacy_model.py).
No known-optimum metadata or saved solutions are used by this solver. For
setup and an overview of the other methods, see the [repository README](../README.md).

## A safe total time cap

Let m be the number of gates, and let `S_h` be the number of SWAPs in a feasible
heuristic route for this instance. The default solve uses

$$
\boxed{T=m+S_h.}
$$

Before building either MIP, a deterministic serial heuristic constructs a
route, absorbs any leading SWAPs into the free initial placement, and replays
the resulting witness independently. Its SWAP count determines T; reference
optimum labels never enter this calculation. The optional beam heuristic can
then improve the incumbent within this horizon. `--plain` disables the beam
heuristic and subset bound, but retains this safe-cap construction.

**Why this cap preserves the optimum.** The heuristic establishes a feasible
SWAP count `S_h`, so a minimum-SWAP route uses `S* <= S_h` SWAPs. For any route
with at most `S_h` SWAPs, serialize its operations while preserving their
original order between slots. Operations within a slot are disjoint, so they
can execute in any order without changing placements used by gates or
violating the DAG. Delete idle slots. Absorb all SWAPs before the first gate
into the freely chosen initial placement, and remove SWAPs after the last
gate because final placement is unrestricted. The resulting schedule begins
and ends with a gate and uses at most `m + S_h` slots. Pad its end with idle
slots if needed. It satisfies the model's first-slot gate requirement and
prohibition on first- and final-slot SWAPs. Consequently, at least one
unrestricted SWAP-optimal route is feasible in this finite horizon. Every
model solution is also an unrestricted route, so the two optimum values are
equal.

This argument assumes a nonempty circuit, unit-duration gates and SWAPs,
free initial and final placement, and the operation rules specified below.
Disconnected hardware does not change the proof once a feasible witness is
available. If the heuristic cannot produce a witness, automatic horizon
selection fails without claiming that the circuit is infeasible. An explicit
`--horizon T` can still be used to study that deadline problem.

A safe cap does not itself prove optimality: a validated incumbent must meet
a certified lower bound. With the automatic cap, that certificate proves the
unrestricted minimum SWAP count. With an explicit horizon, the report
conservatively certifies only that deadline, even if it happens to be large
enough. Earlier runs using a fixed `T=12` alone therefore establish optimality
only within that deadline.

## Sets and variables

There are n logical qubits and n physical nodes, as in [`legacy_model.py`](../legacy_model.py). Hardware edges
form E; A contains both directions of each edge. Time slots are `0,...,T-1`.
Gate g acts on logical qubits `(p_g,q_g)`. A dependency `(g,h)` requires g to
execute strictly before h.

| Variable | Meaning | Domain |
|---|---|---|
| `w[t,q,i]` | Logical qubit q occupies physical node i at the start of slot t | Binary |
| `a[t,g]` | Gate g executes in slot t | Binary |
| `f[t,q,i,j]` | Qubit q moves from node i to neighboring node j during slot t | Continuous, [0,1] |
| `x[t,g,i,j]` | Gate g executes with its first qubit at i and second at j | Continuous, [0,1] |

Movement variables exist only for hardware arcs and slots `1,...,T-2`.
They are zero elsewhere: first-slot SWAPs are prohibited, and final-slot SWAPs
cannot help execute any remaining gate. Gate-time variables outside the DAG's
earliest/latest windows are omitted. These windows only use dependency paths.

## 1. Place the qubits

Initial placement is a bijection:

$$
\sum_i w_{0qi}=1 \quad(\forall q), \qquad
\sum_q w_{0qi}=1 \quad(\forall i).
$$

Later placement sums follow from conservation and reciprocal movement below,
so repeating these assignment equations at every time is unnecessary.
The placement variables themselves remain binary at every time.

## 2. Execute every gate once, in dependency order

$$
\sum_t a_{tg}=1 \quad(\forall g), \qquad
\sum_{(i,j)\in A}x_{tgij}=a_{tg} \quad(\forall t,g).
$$

For every dependency `(g,h)` and every time t:

$$
\sum_{\tau=0}^{t}a_{\tau h}
\le \sum_{\tau=0}^{t-1}a_{\tau g}.
$$

At least one gate executes in slot zero, matching the original convention:

$$
\sum_g a_{0g}\ge1.
$$

An empty circuit is therefore infeasible under this convention.

## 3. A qubit either participates in a gate or moves

Define expressions, not additional variables:

$$
\begin{aligned}
C_{tqi}&=
\sum_{g:p_g=q}\sum_{j:\{i,j\}\in E}x_{tgij}
+\sum_{g:q_g=q}\sum_{j:\{i,j\}\in E}x_{tgji},\\
F^{\mathrm{out}}_{tqi}&=\sum_{j:\{i,j\}\in E}f_{tqij},\\
F^{\mathrm{in}}_{tqi}&=\sum_{j:\{i,j\}\in E}f_{tqji}.
\end{aligned}
$$

The single resource-capacity family is

$$
C_{tqi}+F^{\mathrm{out}}_{tqi}\le w_{tqi}.
$$

This links gate allocations to the actual placements. It also prevents two
gates from sharing a qubit, and prevents a gate participant from being swapped
in that slot. Summing over qubits gives the physical-node resource limits;
separate gate/SWAP conflict inequalities are unnecessary.

## 4. Conserve qubits and require reciprocal movement

For every qubit, node, and slot before the last:

$$
w_{t+1,q,i}=w_{tqi}-F^{\mathrm{out}}_{tqi}+F^{\mathrm{in}}_{tqi}.
$$

For every hardware edge and movement slot:

$$
\sum_q f_{tqij}=\sum_q f_{tqji} \quad(\forall\{i,j\}\in E).
$$

The first equation updates locations. The second forces the same amount of
movement in both directions along an edge. Together with binary placements,
these equations permit exactly disjoint SWAPs, rather than arbitrary token
rotations around a cycle.

## 5. Minimize SWAPs

Choose one orientation `(i,j)` for each undirected edge, denoted `E+`:

$$
\min\ \sum_{t=1}^{T-2}\sum_{(i,j)\in E^+}\sum_q f_{tqij}.
$$

Each SWAP contributes exactly one unit in the selected direction.
Equivalently, minimize half the movement over all directed arcs.

## Why continuous movement variables are sufficient

A binary placement gives each qubit one origin. Capacity forces its outgoing
flow from every other node to zero. Conservation and the next binary placement
then force it either to stay put or move with value one to a single neighbor.
Reciprocity forces that neighbor's occupant to return along the same edge.
Thus movements form disjoint SWAPs automatically.

Similarly, binary placement and a binary gate time leave only one possible
oriented hardware edge for an executed gate, so its allocation is integral.
This argument depends on the binary placements; it does **not** imply that the
LP relaxation is integral.

The model allows idle slots. The extractor removes globally idle slots and
pads the end, preserving strict dependencies, cost, and the deadline. This
replaces the original activity-prefix variables without changing the optimum.

## Optional coverage inequality

In addition to the default movement bound described below, `--coverage` adds one valid
inequality. Let D be the number of distinct unordered logical interaction
pairs and let `N(i)` denote the hardware neighbors of node i:

$$
c_{ij}=|N(i)|+|N(j)|-2-2|N(i)\cap N(j)|,
$$

$$
D\le |E|+\sum_{t=1}^{T-2}\sum_{(i,j)\in E^+}
c_{ij}\sum_q f_{tqij}.
$$

Initially there are exactly `|E|` adjacent logical pairs. Swapping on edge
`{i,j}` creates at most `c[i,j]` new pairs; common neighbors create none.
Every required interaction must appear among those pairs. Parallel disjoint
SWAPs can be serialized for this counting argument. This bound is valid on
any hardware graph, but it is often weak and can be vacuous on dense graphs.

## Core size and limitations

With m gates, the model has
`O(T*n*n + T*n*|E| + T*m*(1+|E|))` variables and
`O(T*n*n + T*|E| + T*m + T*|DAG|)` constraints. It does not enumerate
placements or complete routing states. Disconnected graphs are accepted;
automatic horizon selection first needs a feasible heuristic witness. With
an explicit deadline, Gurobi can also be used to test feasibility.

The formulation is exact for the stated horizon and conventions, and the
default horizon preserves the unrestricted SWAP optimum as proved above.
Exactness is not a runtime guarantee: fractional mixtures of static placements can
still give a zero-SWAP LP bound. The optional coverage inequality helps some
bounds but does not eliminate that weakness. This is why a stronger movement
bound is computed separately, rather than relying on a smaller core alone.

## Stronger bound: minimum movement of small qubit groups

The bound is implemented in
[`general_subset_bound.py`](general_subset_bound.py).

Let `c[q]` be the total number of SWAPs involving logical qubit q, and let C
be the total SWAP count. Every routing schedule within T slots satisfies

$$
2C=\sum_q c_q,\qquad 0\le c_q\le C,
\qquad c_q\le\min\{\max(0,T-2),\ T-d_q\},
$$

where `d[q]` counts gates involving q. A qubit cannot perform a gate and move
in the same slot. The count variables and C are nonnegative integers.
With the automatic cap, these T-dependent limits preserve a global optimum
by the serialization argument above.

Choose a small **proper** subset S of logical qubits. Keep their physical
positions and their mutual gates; the other logical qubits become anonymous.
**All physical nodes and hardware edges remain.** Retained gates inherit
precedence even when their dependency path passes through a removed gate.
Drop the deadline, allowing adjacent ready gates to execute at zero movement
cost, sequentially if necessary.

For nonnegative weights alpha on S and an initial injective placement p,
compute `L[S,alpha,p]`: the minimum weighted movement needed to finish those
retained gates. A SWAP involving one tracked qubit q costs `alpha[q]`; a SWAP
involving two tracked qubits q,r costs `alpha[q]+alpha[r]`. This computation
uses shortest paths in the **projected** problem: a graph for the selected
qubits and their retained gates. The graph is reused for every weight vector.

Every complete routing schedule projects to a feasible walk of this kind.
Consequently, it satisfies

$$
\sum_{q\in S}\alpha_q c_q\ge L_{S,\alpha,p(S)}.
$$

The lower bound depends on the chosen initial placement. We express that
dependence linearly using a continuous nonnegative variable `lambda[S,p]`
for each reachable projected initial placement, and a shared binary initial
assignment `b[q,i]`:

$$
\sum_p\lambda_{S,p}=1,\qquad
\sum_{p:p(q)=i}\lambda_{S,p}=b_{qi}\quad(\forall q\in S,i).
$$

The b variables satisfy the usual assignment row and column sums. When b is
binary, these equations force all lambda mass onto the selected positions
of S. All weight vectors for a subset share the same lambda variables.
The movement inequality becomes

$$
\boxed{\sum_{q\in S}\alpha_q c_q
\ge\sum_p L_{S,\alpha,p}\lambda_{S,p}.}
$$

**Minimize C subject to these inequalities and the count/assignment equations.**
This is the auxiliary MIP. It is a relaxation of the routing problem: its
count vectors need not describe any complete route. Only its certified
**dual bound (`ObjBound`)** is a routing lower bound. Its incumbent objective
must never be mistaken for a routing bound or a feasible route.

### Fixed choices and size

For three through six qubits the solver uses subsets of size `n-1`, capped
at five; for larger instances it uses size four. No subset model is built
for fewer than three qubits. The candidate weights are the all-ones vector
and vectors made by changing one or two entries to zero or two. Size five
gives 51 vectors; a preprocessing timeout may keep fewer completed tables.
This is a fixed, instance-independent rule, not a choice based on an optimum.

This preprocessing is not free or purely polynomial in circuit size:
projected gate-progress states can grow exponentially. Five tracked qubits on
six nodes have 720 placements, so these projections represent every physical
placement locally while omitting one qubit's gates. If that omitted qubit is
unused, such a projection can retain all gate progress as well as complete
placement information.

Preprocessing has time, state, and arc caps. Oversized projections are
skipped; a timeout keeps only fully computed shortest-path tables. An
unfinished Dijkstra search cannot supply valid minimum-cost tables.
Omitting projections only weakens the bound. The JSON report records
whether preparation was complete and which projections were skipped.

### How the full solver uses the bound

1. A serial heuristic constructs and independently validates a feasible route,
   setting `T=m+S_h` unless an explicit horizon was requested. A bounded beam
   heuristic can improve the candidate; all returned schedules are replayed.
2. The auxiliary MIP computes a lower bound. If its integer-rounded dual
   bound equals the validated schedule cost, optimality is proved immediately.
3. Otherwise the core routing MIP receives that schedule as a MIP start and
   the proved bound as a lower bound on one integer total-SWAP variable.
   It uses the remaining time budget.

The lambda variables stay in the auxiliary MIP, keeping the core small.
The bound method applies to cycles, stars, grids, and arbitrary undirected
hardware; it uses no path order or path-specific distance identity.
The heuristic supplies only upper bounds. All its time and all lower-bound
preprocessing are included in the reported production time.

Both MIPs use Gurobi `Symmetry=2` by default, its maximum setting for aggressive
symmetry detection. `--symmetry -1` restores Gurobi's automatic setting;
`0` disables detection and `1` selects the conservative setting. These settings
affect solver search, not the feasible routing schedules or optimality proof.
See the [Gurobi parameter reference](https://docs.gurobi.com/projects/optimizer/en/current/reference/parameters.html#parameter.Symmetry).

## Run and inspect results

From the repository root, after activating an environment with Gurobi:

```bash
# Safe automatic horizon and maximum symmetry detection (both defaults).
python -m our_formulation benchmark/qv6_depth6_01_opt7.txt \
  --limit 10 --symmetry 2 --output results/ours_qv6.jsonl

# The same formulation on a grid, plus optional bound diagnostics.
python -m our_formulation benchmark/qv6_depth6_01_grid2x3.txt \
  --limit 10 --coverage --diagnose-lp

# Core routing MIP with the mandatory safe-cap heuristic only.
python -m our_formulation benchmark/qv6_depth6_01_opt7.txt \
  --limit 10 --plain --coverage

# An explicitly restricted deadline; certificate scope is fixed_horizon.
python -m our_formulation benchmark/instance_01_opt2.txt \
  --horizon 12 --limit 10
```

The `--limit` budget includes construction, preprocessing, the heuristic,
both MIPs, and schedule replay. It is a soft wall limit: solver operations can
overshoot slightly. Timeouts do not certify optimality unless a valid lower
bound matches a validated schedule. `--diagnose-lp` runs after the production
solve and can add time outside the budget. All CLI options are listed by
`python -m our_formulation --help`.

The JSON report stores schedules, bounds, branch-and-bound nodes, model
dimensions, phase timings, projection completeness, and source hashes.
`horizon_source` is `serial_heuristic` or `explicit`; `heuristic_swap_count`
records `S_h` for an automatic horizon, and `heuristic_seconds` records its
construction time. `safe_horizon` contains details of that construction.
`optimality_scope` is `global` for an automatic safe cap and `fixed_horizon`
for an explicit deadline. `proven_optimal` records a certificate within that
scope; `globally_optimal` is true only when the automatic-horizon solve proves
optimality. `symmetry` records the requested Gurobi setting.
`routing` and `movement_bound` describe the two MIPs separately; top-level
`nodes` sums their node counts. Top-level `bound` is rounded using integer
SWAP costs, while each MIP reports its raw bound. The optional `lp` diagnostic
reports the core **before** importing the computed lower bound;
`movement_bound.lp` reports the auxiliary LP. Neither diagnostic is included
in production timing.

In `schedule`, `placement[t][q]` gives the physical node of logical qubit `q`
at the start of slot `t`, `gate_times` maps gate IDs to execution slots, and
`swaps` contains `[slot, edge_index]` entries. Edge indices refer to the
reported `hardware_edges` list. A returned schedule is independently replayed
before it is marked validated. `objective` is the feasible SWAP count;
`bound` is the certified lower bound, never an auxiliary MIP incumbent.

For use from Python:

```python
from load_instance import load_instance
from our_formulation import solve

gates, dag, edges, n, _ = load_instance("benchmark/instance_01_opt2.txt")
result = solve(gates, dag, edges, n, time_limit=10, symmetry=2)
print(result["horizon"], result["objective"], result["bound"])
print(result["globally_optimal"], result["optimality_scope"])
```

The flow approach is related to
[Nannicini et al., *Optimal qubit assignment and routing via integer programming*](https://arxiv.org/abs/2106.06446).
Here gate times remain free subject to the supplied DAG, instead of fixing
input gate layers. No novelty claim is made for the formulation or the
movement-bound technique. The [hardware-topology benchmark guide](../benchmark/README.md)
describes the topology variants and their assumptions.
