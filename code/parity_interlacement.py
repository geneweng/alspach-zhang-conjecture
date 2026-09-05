#!/usr/bin/env python3
"""Parity-refined circuit counts using a two-sheeted quotient and Traldi's
extended Cohn--Lempel equality (arXiv:0903.4405, Theorem 4).

The input is a perfect matching of the pentagon quotient. No existence of a
colourable support is inferred from the two nullity computations.
"""


def binary_nullity(rows):
    basis = {}
    for row in rows:
        while row:
            p = row.bit_length() - 1
            if p not in basis:
                basis[p] = row
                break
            row ^= basis[p]
    return len(rows) - len(basis)


def circuit_nullity(n, edges, transitions):
    """Euler-system interlacement matrix for a loopless 4-regular multigraph.

    transitions[v] maps each incident edge ID to its paired edge at v.
    Edges, including parallel edges, retain distinct IDs throughout.
    Returns (binary nullity, number of graph components).
    """
    incident = [[] for _ in range(n)]
    for e, (u, v) in enumerate(edges):
        assert u != v
        incident[u].append(e)
        incident[v].append(e)
    assert all(len(es) == 4 for es in incident)
    for es, pairing in zip(incident, transitions):
        assert set(pairing) == set(es)
        assert all(pairing[e] != e and pairing[pairing[e]] == e for e in es)
    remaining = [es[:] for es in incident]
    unused = set(range(len(edges)))
    nullity, connected = 0, 0
    while unused:
        start = edges[min(unused)][0]
        stack, edge_stack, reverse_vertices, reverse_edges = [start], [], [], []
        while stack:
            v = stack[-1]
            while remaining[v] and remaining[v][-1] not in unused:
                remaining[v].pop()
            if remaining[v]:
                e = remaining[v].pop()
                unused.remove(e)
                u, w = edges[e]
                stack.append(w if v == u else u)
                edge_stack.append(e)
            else:
                reverse_vertices.append(stack.pop())
                if edge_stack:
                    reverse_edges.append(edge_stack.pop())
        vertices, walk = reverse_vertices[::-1], reverse_edges[::-1]
        assert len(vertices) == len(walk) + 1 and vertices[0] == vertices[-1]
        assert all(set(edges[e]) == {vertices[i], vertices[i + 1]}
                   for i, e in enumerate(walk))
        visits = {}
        for i, v in enumerate(vertices[:-1]):
            visits.setdefault(v, []).append(i)
        assert all(len(positions) == 2 for positions in visits.values())
        retained, diagonal = [], {}
        for v, (i, j) in visits.items():
            before_i, after_i = walk[i - 1], walk[i]
            before_j, after_j = walk[j - 1], walk[j]
            pairing = transitions[v]
            if pairing[before_i] == after_i:
                assert pairing[before_j] == after_j
                continue  # follows the Euler circuit: delete its row/column
            retained.append(v)
            diagonal[v] = int(pairing[before_i] == before_j)
            assert diagonal[v] or pairing[before_i] == after_j
        rows = []
        for i, v in enumerate(retained):
            a, b = visits[v]
            row = diagonal[v] << i
            for j, w in enumerate(retained):
                c, d = visits[w]
                if (a < c < b < d) or (c < a < d < b):
                    row |= 1 << j
            rows.append(row)
        nullity += binary_nullity(rows)
        connected += 1
    return nullity, connected


def parity_profile(face, support):
    """Return the two Cohn--Lempel nullities and their oddness defect.

    R = Q - support has a forced transition pairing at each vertex:
    (e1,e4) has sign 1; (e2,e3) has sign 0, with e0 selected.
    Put half-edge weight 1 at e1 and zero at the other three positions.
    Edge voltages are the sums of endpoint half-edge weights. The same
    transitions lifted to the two-sheeted voltage cover split a circuit
    into two exactly when its prescribed sign sum is even.
    """
    assert face.is_matching(support)
    active = [e for e in range(face.m) if not support & (1 << e)]
    index = {e: i for i, e in enumerate(active)}
    edges = [face.edges[e] for e in active]
    n = face.graph.q
    weights, transitions = {}, [{} for _ in range(n)]
    for v, cyclic in enumerate(face.graph.orders):
        p = next(i for i, e in enumerate(cyclic) if support & (1 << e))
        e1, e2, e3, e4 = [index[cyclic[(p + j) % 5]] for j in range(1, 5)]
        transitions[v] = {e1: e4, e4: e1, e2: e3, e3: e2}
        for e in (e1, e2, e3, e4):
            weights[v, e] = int(e == e1)
    lifted_edges, lifted_ids = [], {}
    for e, (u, v) in enumerate(edges):
        voltage = weights[u, e] ^ weights[v, e]
        for sheet in (0, 1):
            f = len(lifted_edges)
            lifted_edges.append((2 * u + sheet, 2 * v + (sheet ^ voltage)))
            lifted_ids[u, sheet, e] = f
            lifted_ids[v, sheet ^ voltage, e] = f
    lifted_transitions = [{} for _ in range(2 * n)]
    for v, pairing in enumerate(transitions):
        for sheet in (0, 1):
            lifted_transitions[2 * v + sheet] = {
                lifted_ids[v, sheet, e]: lifted_ids[v, sheet, f]
                for e, f in pairing.items()}
    nu, c = circuit_nullity(n, edges, transitions)
    lifted_nu, lifted_c = circuit_nullity(2 * n, lifted_edges, lifted_transitions)
    odd = 2 * (nu + c) - (lifted_nu + lifted_c)
    assert odd >= 0 and odd % 2 == 0
    return {'nullity': nu, 'components': c, 'lifted_nullity': lifted_nu,
            'lifted_components': lifted_c, 'odd_circuits': odd}
