# Fixed-layer SWAP-count binary integer program

This is a SWAP-count adaptation of Section 3.1 of Nannicini, Bishop,
Günlük, and Jurcevic,
[Optimal Qubit Assignment and Routing via Integer Programming](https://doi.org/10.1145/3544563).
The original paper permits a gate followed by a SWAP of the same two qubits in
one model timestep. Here, a qubit participating in a gate cannot participate in
any SWAP in that timestep. The sole objective is the number of SWAPs.

## Pregrouping and time convention

The input has `N` logical and `N` physical qubits, an undirected connected
hardware graph `H=(V,E)`, two-qubit gates indexed by `g`, and a gate DAG. All
operations have unit duration. Hardware-edge direction does not constrain gates.
At the start of each layer, collect gates whose DAG predecessors are complete.
Scan those gates in ascending index order and take a maximal subset with
pairwise disjoint logical qubits. Mark this entire subset complete only after
finishing the layer. Repeat until every gate is assigned; reject cyclic DAGs.
This is deterministic greedy grouping, not optimization over possible groups.
The gate list does not impose additional precedence beyond the supplied DAG.

Insert `D` empty slots between each pair of consecutive layers. For `L >= 1`
layers, there are `T = L + (L-1)D` total slots, indexed `0,...,T-1`. Gates of a
layer execute simultaneously at their assigned slot; they cannot move into
dummy slots. The DAG is therefore enforced by preprocessing. Each group is a
logical-qubit matching, an essential condition for the formulation below.

Let `G_t` be the gates at slot `t`, and `Q_t` their logical qubits. Placement at
`t` is the placement **before** operations in slot `t`. Gate qubits stay still;
SWAPs on other qubits determine the placement at `t+1`. The last slot has gates
but no outgoing movement variables. Initial and final placement are free.
There is no need for leading routing slots with free initial placement or
trailing routing slots with free final placement. Idle qubits may swap during
a gate slot, including the first slot. Empty input uses one empty placement slot.

## Variables and constraints

Let `A` contain both orientations of every edge, and let `A+` also include
`(i,i)` for every node. All variables below are **binary**, including `x` and `y`.
Use gate IDs, not logical-qubit pairs, so repeated gates are kept distinct.

- `w[t,q,i]`: logical qubit `q` occupies physical node `i` at slot `t`.
- `y[t,g,i,j]`: gate `g=(p,q)` is executed on oriented hardware arc `(i,j)`.
  Defined only when `g in G_t` and `(i,j) in A`.
- `x[t,q,i,j]`: qubit `q` goes from node `i` at `t` to node `j` at `t+1`.
  Defined for `0 <= t < T-1` and `(i,j) in A+`.

Bijection at every slot:

```text
sum_i w[t,q,i] = 1                     for every t,q
sum_q w[t,q,i] = 1                     for every t,i
```

Place every gate on exactly one hardware arc:

```text
sum_(i,j in A) y[t,g,i,j] = 1          for every t, g in G_t
y[t,g,i,j] >= w[t,p,i] + w[t,q,j] - 1
```

At nonfinal slots use movement-linked upper inequalities:

```text
y[t,g,i,j] <= x[t,p,i,i]
y[t,g,i,j] <= x[t,q,j,j]
```

At the **final** slot, there is no `x[t,...]`, so use ordinary McCormick bounds:

```text
y[t,g,i,j] <= w[t,p,i]
y[t,g,i,j] <= w[t,q,j]
```

The lower bound `y >= 0` is supplied by the binary variable domain. Because gate
qubits are explicitly frozen below, the strengthened upper bounds equal the
ordinary location upper bounds even in the LP relaxation; they are not an
additional LP strengthening once the freezing equations have been imposed.

Conserve every qubit through the time-expanded network:

```text
x[t,q,i,i] + sum_(j in N(i)) x[t,q,i,j] = w[t,q,i]
x[t,q,i,i] + sum_(j in N(i)) x[t,q,j,i] = w[t+1,q,i]
```

Freeze all gate qubits during that slot:

```text
x[t,q,i,i] = w[t,q,i]                 for q in Q_t, every i, t < T-1
```

Nonnegativity and outflow conservation then force all their off-diagonal flows
to zero. Adding those zero equations separately would be redundant.

Require reciprocal movement of idle qubits:

```text
sum_(q not in Q_t) x[t,q,i,j]
    = sum_(q not in Q_t) x[t,q,j,i]   for every {i,j} in E, t < T-1
```

Only one equality per undirected edge is necessary: reversing it gives the
same equation. Do not omit this constraint: flow balance alone permits
directed cycles of length three or more, which are not simultaneous SWAPs.
There are no separate `s` variables; edge usage is derived from reciprocal flow.

The objective is precisely:

```text
minimize  (1/2) sum_(t=0,...,T-2) sum_q sum_(i,j in A) x[t,q,i,j]
```

Self-loops are excluded. Each SWAP contributes two off-diagonal movements,
so the objective counts each SWAP once. Two concurrent SWAPs cost two, not one.

## Why these constraints are sufficient

At an integer solution each node has exactly one occupant. Its occupant has
exactly one outgoing flow. If it stays, the node is unused by a SWAP. If it
moves from `i` to `j`, reciprocal flow requires the occupant of `j` to move
back to `i`; unique occupancy leaves no other possible returning qubit. Thus
all nontrivial movements are disjoint transpositions on hardware edges.
Three-cycles, overlapping swaps, and teleportation are impossible. A frozen
gate qubit cannot move, and reciprocal movement prevents an idle qubit from
swapping into its node. Gate variables place each gate's two qubits on an edge.
Therefore every integer feasible solution gives a valid schedule of the fixed
layers. Conversely, every such schedule within the chosen horizon sets `w`,
`y`, and `x` satisfying these equations, with objective equal to its SWAP count.

## Sufficient dummy slots, including preservation of the minimum cost

The original instances use path hardware; the general-hardware suite adds
other topologies. On a path with `N >= 3`, use
**`D=N`**, with `D=1` for two nodes and `D=0` for one node. Number each token by
its target position along the path. Odd-even transposition sorting routes any
permutation in at most `N` parallel matching rounds; see
[Alon, Chung, and Graham, Routing Permutations on Graphs via Matchings](https://collaborate.princeton.edu/en/publications/routing-permutations-on-graphs-via-matchings-extended-abstract/).
Each performed adjacent swap removes exactly one inversion. Every route needs
at least the inversion count in swaps, so this construction also has minimum
SWAP count for those endpoints. Consequently the bound preserves a globally
minimum-SWAP solution for the fixed groups, not merely feasibility.

In particular, moving any gate-slot idle swaps into the subsequent gap and
replacing the whole gap route by this minimum route cannot increase the cost.
There are enough dummy slots for that replacement. No extra allowance for a
gate-slot swap is needed in the bound.

The source paper's Sections 4.3/4.4 use five dummy steps for its six-qubit
experiments and allow gate/SWAP merging. Five is **not** a universal six-node
path rerouting guarantee under the modified operation rules: reversing the six
tokens requires six parallel matching rounds. A fully occupied gate layer
provides no additional routing round here. The tests enumerate all 720 path
permutations and verify the six-round worst case, and also check odd-even
sorting and minimum inversion counts for every permutation through seven nodes.
This does not claim that every particular benchmark needs all six slots.

For other connected hardware graphs, the conservative default is
**`D=N(N-1)/2`**. Choose a spanning tree, bring the correct token to a leaf in at
most `N-1` swaps, remove the leaf, and repeat. The remaining tree is connected,
so this takes at most `(N-1)+(N-2)+...+1` sequential swaps. A minimum-cost route
is no longer than that constructive route, hence can be serialized in `D`
slots. This bound can be loose, but preserves the unrestricted optimum.

An explicit `--dummy-slots` override is supported. If it is below the safe
bound, results are marked `unrestricted_rerouting_guaranteed=false`; optimality
then applies only to that restricted horizon unless separately verified.

## Dummy-slot symmetry

By default, each dummy slot has a binary `z[t]` and constraints

```text
2 z[t] <= sum_(q,(i,j) in A) x[t,q,i,j] <= 2 floor(N/2) z[t]
z[t] >= z[t+1]                         within a single dummy block only
```

These are exact activity indicators and the paper's Eq. (4) compaction rule.
Moving an entirely idle slot to the end of its dummy block leaves the order of
nonempty swap matchings, the next layer placement, and cost unchanged. Thus
this removes symmetric schedules without losing an optimum. It does not
reorder nonempty matchings. Use `--no-compact` to disable it. No other cuts,
external optimal-value bounds, MIP starts, placement fixings, or secondary
objectives are enabled in the benchmark.

## Verification

`validation.py` replays each incumbent using explicit physical-edge swaps,
checking every placement, gate, dependency, resource conflict, and objective.
`exact.py` independently enumerates all placements for `N <= 7` and computes
minimum swap costs between consecutive admissible layer placements. It runs
**after** Gurobi and its time is recorded separately. It supplies no information
to the BIP. The tests also compare the BIP against a separate finite-slot
matching-enumeration dynamic program, including intentionally short horizons.

Fixed simultaneous layers restrict the original DAG-scheduling problem. The
`optimal_swaps` field in a benchmark is recorded only as a reference label;
it is never used to constrain, start, or certify this model.
