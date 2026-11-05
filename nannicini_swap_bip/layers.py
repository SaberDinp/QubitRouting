"""Deterministic DAG pregrouping and provably sufficient routing horizons."""


def validate_input(gates, dag, edges, n):
    if type(n) is not int or n < 1:
        raise ValueError("N must be a positive integer")
    gates = tuple(tuple(g) for g in gates)
    dag = tuple(dict.fromkeys(tuple(e) for e in dag))
    edges = tuple(tuple(sorted(e)) for e in edges)
    for name, pairs, size in (("gates", gates, n), ("DAG", dag, len(gates)),
                              ("hardware edges", edges, n)):
        if any(len(e) != 2 or e[0] == e[1] or any(
                type(v) is not int or not 0 <= v < size for v in e) for e in pairs):
            raise ValueError(f"Invalid {name}")
    if len(set(edges)) != len(edges):
        raise ValueError("Duplicate undirected hardware edge")
    neighbors = [set() for _ in range(n)]
    for i, j in edges:
        neighbors[i].add(j)
        neighbors[j].add(i)
    seen, stack = {0}, [0]
    while stack:
        for v in neighbors[stack.pop()] - seen:
            seen.add(v)
            stack.append(v)
    if len(seen) != n:
        raise ValueError("Hardware must be connected, as assumed in the paper")
    return gates, dag, edges


def greedy_layers(gates, dag):
    """Take a maximal disjoint subset of ready gates, in gate-index order.

    Ready is evaluated at the start of each layer. A successor cannot join its
    predecessor's layer. Only the supplied DAG specifies precedence; list order
    breaks ties, including between otherwise unordered gates sharing a qubit.
    """
    pred = [set() for _ in gates]
    for u, v in dag:
        pred[v].add(u)
    done, layers = set(), []
    while len(done) != len(gates):
        used, layer = set(), []
        for g, pair in enumerate(gates):
            if g not in done and pred[g] <= done and not used.intersection(pair):
                layer.append(g)
                used.update(pair)
        if not layer:
            raise ValueError("gate_dag must be acyclic")
        layers.append(tuple(layer))
        done.update(layer)
    return tuple(layers)


def path_order(edges, n):
    """Return the hardware path in physical-node order, or None."""
    if n == 1:
        return (0,)
    neighbors = [set() for _ in range(n)]
    for i, j in edges:
        neighbors[i].add(j)
        neighbors[j].add(i)
    if len(edges) != n - 1 or max(map(len, neighbors)) > 2:
        return None
    ends = [i for i in range(n) if len(neighbors[i]) == 1]
    if len(ends) != 2:
        return None
    result = [min(ends)]
    while len(result) < n:
        choices = neighbors[result[-1]] - set(result)
        if not choices:
            return None
        result.append(min(choices))
    return tuple(result)


def sufficient_dummy_slots(edges, n):
    """Enough slots to realize a minimum-SWAP route between ANY placements.

    Paths: odd-even transposition takes at most n matching rounds and exactly
    the inversion count (the minimum number of adjacent swaps). General
    connected graphs: fixing leaves of a spanning tree gives at most n(n-1)/2
    sequential swaps; a minimum route is no longer. See FORMULATION.md.
    """
    if n <= 2:
        return max(0, n - 1)
    if path_order(edges, n) is not None:
        return n
    return n * (n - 1) // 2


def insert_dummy_slots(layers, count):
    slots, layer_times, blocks = [], [], []
    for k, layer in enumerate(layers):
        if k:
            block = tuple(range(len(slots), len(slots) + count))
            blocks.append(block)
            slots.extend([()] * count)
        layer_times.append(len(slots))
        slots.append(layer)
    return tuple(slots or [()]), tuple(layer_times), tuple(blocks)
