from collections import deque

import gurobipy as gp
from gurobipy import GRB

import load_instance

model = gp.Model("QMRS")

# T = # of Timesteps to consider in the model for placement and routing
T = 12

G, gate_dag, hw_edges, N, optimal_swaps = load_instance.load_instance("benchmark/qv6_depth6_01_opt7.txt")
# G, gate_dag, hw_edges, N, optimal_swaps = load_instance.load_instance("benchmark/qv6_depth6_04_opt13.txt")
# G, gate_dag, hw_edges, N, optimal_swaps = load_instance.load_instance("benchmark/instance_04_opt4.txt")
# G, gate_dag, hw_edges, N, optimal_swaps = load_instance.load_instance("benchmark/instance_01_opt2.txt")

for u, v in hw_edges:
    assert u < v


a = model.addVars(T, len(G), vtype=GRB.BINARY, name="a")
s = model.addVars(T - 1, len(hw_edges), vtype=GRB.BINARY, name="s")
w = model.addVars(T, N, N, vtype=GRB.BINARY, name="w")
# activation variables for times
z = model.addVars(T, vtype=GRB.BINARY, name="active")
"""Variable fixing"""

# Longest predecessor paths give the earliest possible execution times.
successors = [[] for _ in G]
in_degree = [0] * len(G)
for g1, g2 in gate_dag:
    successors[g1].append(g2)
    in_degree[g2] += 1

# ES is the earliest time to schedule that gate and LS is the latest schedule time.
ES = [0] * len(G)
ready = deque(g for g in range(len(G)) if in_degree[g] == 0)
topological_order = []
while ready:
    g = ready.popleft()
    topological_order.append(g)
    for successor in successors[g]:
        ES[successor] = max(ES[successor], ES[g] + 1)
        in_degree[successor] -= 1
        if in_degree[successor] == 0:
            ready.append(successor)

if len(topological_order) != len(G):
    raise ValueError("gate_dag must be acyclic to compute scheduling windows.")

# LS counts successor-path edges, so the latest execution time is T - 1 - LS[g].
LS = [0] * len(G)
for g in reversed(topological_order):
    for successor in successors[g]:
        LS[g] = max(LS[g], LS[successor] + 1)

# Fix impossible gate times before either the MIP or its relaxation is solved.
for g in range(len(G)):
    for t in range(T):
        if t < ES[g] or t > T - 1 - LS[g]:
            a[t, g].UB = 0


for e in range(len(hw_edges)):
    s[0, e].UB = 0

"""End of variable fixing"""

gates_of_q = {
    q: [g for g, (p, r) in enumerate(G) if q == p or q == r]
    for q in range(N)
}
edges_of_i = {
    i: [e for e, (u, v) in enumerate(hw_edges) if i == u or i == v]
    for i in range(N)
}

neighbors = {i: [] for i in range(N)}
for i, j in hw_edges:
    neighbors[i].append(j)
    neighbors[j].append(i)


directed_edges = []
for i, j in hw_edges:
    directed_edges.append((i, j))
    directed_edges.append((j, i))

A = len(directed_edges)


out_arcs = {
    i: [k for k, (u, v) in enumerate(directed_edges) if u == i]
    for i in range(N)
}

in_arcs = {
    i: [k for k, (u, v) in enumerate(directed_edges) if v == i]
    for i in range(N)
}

# bijection of physical and logical qubits at all times
for t in range(T):
    for p in range(N):
        model.addConstr(gp.quicksum(w[t, p, i] for i in range(N)) == 1)

for t in range(T):
    for i in range(N):
        model.addConstr(gp.quicksum(w[t, p, i] for p in range(N)) == 1)

# When a gate is scheduled, its qubits must be adjacent
# for t in range(T):
#     for g, (p, q) in enumerate(G):
#         for i in range(N):
#
#             model.addConstr(
#                 a[t, g] + w[t, p, i]
#                 <=
#                 1 + gp.quicksum(w[t, q, j] for j in neighbors[i])
#             )
#
#             model.addConstr(
#                 a[t, g] + w[t, q, i]
#                 <=
#                 1 + gp.quicksum(w[t, p, j] for j in neighbors[i])
#             )

x = model.addVars(
    T, len(G), A,
    lb=0.0, ub=1.0,
    vtype=GRB.CONTINUOUS,
    name="x"
)

# for t in range(T):
#     for g, (p, q) in enumerate(G):
#         for i in range(N):
#             for k in out_arcs[i]:
#                 u, v = directed_edges[k]
#                 model.addConstr(x[t, g, k] >= w[t, p, u] + w[t, q, v] - a[t, g])

# for t in range(T):
#     for g, (p, q) in enumerate(G):
#
#         # All executed gate mass must use a hardware edge
#         model.addConstr(
#             gp.quicksum(x[t, g, k] for k in range(A))
#             == a[t, g]
#         )

        # # Logical p is the tail of the directed physical edge
        # for i in range(N):
        #     model.addConstr(
        #         gp.quicksum(
        #             x[t, g, k] for k in out_arcs[i]
        #         )
        #         <= w[t, p, i]
        #     )
        #
        # # Logical q is the head
        # for j in range(N):
        #     model.addConstr(
        #         gp.quicksum(
        #             x[t, g, k] for k in in_arcs[j]
        #         )
        #         <= w[t, q, j]
        #     )

for t in range(T):
    # Every scheduled gate must be allocated to hardware arcs.
    for g in range(len(G)):
        model.addConstr(
            gp.quicksum(x[t, g, k] for k in range(A)) == a[t, g]
        )

    # All gates using logical qubit q at physical node i
    # must share the available placement w[t, q, i].
    for q in range(N):
        for i in range(N):
            gate_use = gp.quicksum(
                x[t, g, k]
                for g in gates_of_q[q]
                for k in (
                    out_arcs[i] if G[g][0] == q else in_arcs[i]
                )
            )

            model.addConstr(
                gate_use <= w[t, q, i],
                name=f"logical_location_capacity[{t},{q},{i}]",
            )

# Swap-location logic

# If node i was not swapped then it should hold its content
for t in range(T - 1):
    for i in range(N):

        r_ti = gp.quicksum(s[t, e] for e in edges_of_i[i])

        for q in range(N):
            model.addConstr(
                w[t + 1, q, i] - w[t, q, i] <= r_ti
            )
            model.addConstr(
                w[t, q, i] - w[t + 1, q, i] <= r_ti
            )

# If a swap occurs on a hardware edge then their content should be swapped
for t in range(T - 1):
    for e, (i, j) in enumerate(hw_edges):
        for q in range(N):

            model.addConstr(
                w[t + 1, q, i] - w[t, q, j]
                <= 1 - s[t, e]
            )
            model.addConstr(
                w[t, q, j] - w[t + 1, q, i]
                <= 1 - s[t, e]
            )

            model.addConstr(
                w[t + 1, q, j] - w[t, q, i]
                <= 1 - s[t, e]
            )
            model.addConstr(
                w[t, q, i] - w[t + 1, q, j]
                <= 1 - s[t, e]
            )

# Gates must be assigned
for g in range(len(G)):
    model.addConstr(gp.quicksum(a[t, g] for t in range(T)) == 1)

# dependency of the gates
for t in range(T):
    for g1, g2 in gate_dag:
        model.addConstr(gp.quicksum(a[tau, g2] for tau in range(t + 1)) <= gp.quicksum(a[tau, g1] for tau in range(t)))

# Each physical qubit at any given time can accommodate at most one swap or gate
# for t in range(T - 1):
#     for q in range(N):
#
#         gate_use = gp.quicksum(
#             a[t, g] for g in gates_of_q[q]
#         )
#
#         for i in range(N):
#
#             swap_use = gp.quicksum(
#                 s[t, e] for e in edges_of_i[i]
#             )
#
#             model.addConstr(
#                 swap_use + w[t, q, i] + gate_use <= 2
#             )

# A physical node can participate in at most one gate or swap
for t in range(T):
    for i in range(N):

        gate_use_at_i = gp.quicksum(
            x[t, g, k]
            for g in range(len(G))
            for k in out_arcs[i] + in_arcs[i]
        )

        if t < T - 1:
            swap_use_at_i = gp.quicksum(
                s[t, e] for e in edges_of_i[i]
            )
        else:
            swap_use_at_i = 0

        model.addConstr(
            gate_use_at_i + swap_use_at_i <= 1
        )


""" cutting planes """

# for t in range(T):
#     for q in range(N):
#         model.addConstr(
#             gp.quicksum(a[t, g] for g in gates_of_q[q]) <= 1
#         )


# for t in range(1, T - 1):
#     for i in range(N):
#         model.addConstr(
#             gp.quicksum(s[t, e] for e in edges_of_i[i]) <= 1
#         )

model.addConstr(gp.quicksum(a[0, g] for g in range(len(G))) >= 1)




"""End of Cutting planes"""

"""Symmetry Breaking"""


for t in range(T):

    operations = gp.quicksum(a[t, g] for g in range(len(G)))

    if t < T - 1:
        operations += gp.quicksum(
            s[t, e] for e in range(len(hw_edges))
        )

    model.addConstr(operations >= z[t])
    model.addConstr(operations <= (N // 2) * z[t])

for t in range(T - 1):
    model.addConstr(z[t] >= z[t + 1])



"""End Symmetry Breaking"""


# minimize the number of swaps
model.setObjective(gp.quicksum(s[t, edge_index] for t in range(T - 1) for edge_index in range(len(hw_edges))), sense=GRB.MINIMIZE)
relaxed = False
if not relaxed:
    model.optimize()
else:
    model.update()
    relaxed_model = model.relax()
    relaxed_model.optimize()

    if relaxed_model.SolCount > 0:
        def format_value(value):
            text = f"{value:.6f}".rstrip("0").rstrip(".")
            return "0" if text == "-0" else text

        def print_vector(title, labels, values):
            print(f"  {title}")
            # Wrap long gate/edge vectors while keeping each value under its label.
            for start in range(0, len(labels), 6):
                columns = labels[start:start + 6]
                entries = values[start:start + 6]
                width = max(10, max(len(label) for label in columns) + 1)
                print("    " + "".join(f"{label:>{width}}" for label in columns))
                print("    " + "".join(f"{value:>{width}}" for value in entries))

        # Read the solved relaxation's variables, not those of the original model.
        values = {var.VarName: format_value(var.X) for var in relaxed_model.getVars()}
        gate_labels = [f"g{g}({p},{q})" for g, (p, q) in enumerate(G)]
        edge_labels = [f"e{e}({i},{j})" for e, (i, j) in enumerate(hw_edges)]
        arc_labels = [f"k{k}({i}->{j})" for k, (i, j) in enumerate(directed_edges)]

        print("\nRelaxed model solution (values rounded to 6 decimal places)")
        print(f"Objective (total swaps): {format_value(relaxed_model.ObjVal)}")

        for t in range(T):
            print(f"\n{'=' * 72}\nTime t = {t}\n{'=' * 72}")
            print(f"  z[t]: time-step activity (solver variable active[{t}]) = {values[f'active[{t}]']}\n")
            print_vector(
                "a[t,g]: gates, labeled g(logical qubit pair)",
                gate_labels,
                [values[f"a[{t},{g}]"] for g in range(len(G))],
            )
            print("\n  x[t,g,k]: gate allocation to directed physical edges, labeled k(tail->head)")
            print("    For gate g(p,q), p is at the tail and q at the head; sum_k x[t,g,k] = a[t,g].")
            for g, gate_label in enumerate(gate_labels):
                print_vector(
                    f"  {gate_label} (a = {values[f'a[{t},{g}]']})",
                    arc_labels,
                    [values[f"x[{t},{g},{k}]"] for k in range(A)],
                )
            print()
            if t < T - 1:
                print_vector(
                    f"s[t,e]: swaps on physical edges, affecting placement at t = {t + 1}",
                    edge_labels,
                    [values[f"s[{t},{e}]"] for e in range(len(hw_edges))],
                )
            else:
                print("  s[t,e]: no swap variables at the final time step.")

            print("\n  w[t,q,i]: placement (rows = logical qubits, columns = physical nodes)")
            print(f"    {'q / i':>6}" + "".join(f"{f'i{i}':>10}" for i in range(N)))
            for q in range(N):
                row = [values[f"w[{t},{q},{i}]"] for i in range(N)]
                print(f"    {f'q{q}':>6}" + "".join(f"{value:>10}" for value in row))
    else:
        print(f"No solution available for the relaxed model (status: {relaxed_model.Status}).")
