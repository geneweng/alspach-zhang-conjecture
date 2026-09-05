#!/usr/bin/env python3
"""Test whether the 2026 cycle-double-cover construction yields a Tait flow.

For a nowhere-zero F_2^3-flow, the OpenAI construction assigns a two-point
palette set P_e to every edge.  The union of these pairs is a graph on the
eight elements of F_2^3.  A proper four-colouring of that palette graph pulls
back to a nowhere-zero F_2^2-flow, hence a 3-edge-colouring of the cubic input.

This script constructs the affine compatibility system from the proof,
enumerates every solution for several distinct full-rank flows, and computes
the exact chromatic number of each resulting eight-vertex palette graph.  It
uses K_4 as a full-rank positive control, the Petersen graph as a negative
control, and Cay(S_5,{(01),(01234),(04321)}) as a nontrivial colourable Cayley
test.  Finally it starts the latter test from a rank-two (colour-aligned) flow
to check the elementary positive direction, and directly searches for the
full-rank tetrahedral palette whose existence is equivalent to that direction.
It also finds a generator-separated palette certified by a transvection
matching while ruling out all four parallel matching certificates, and
measures the induced restriction map from the weighted cut code to the
dangerous transvection equations.

Run from the repository root with:

    python3 code/cdc_palette_experiment.py
"""

import random
from itertools import combinations, product

from pysat.solvers import Cadical153

from cayley_snark_check import closure, mul, three_edge_colourable


def nz_three_bit_flow(n, edges, seed=0, forbidden=()):
    """Find a new nowhere-zero F_2^3-flow on a loopless cubic graph."""
    inc = [[] for _ in range(n)]
    for ei, (u, v) in enumerate(edges):
        inc[u].append(ei)
        inc[v].append(ei)
    assert all(len(es) == 3 for es in inc)

    clauses = []
    for e in range(len(edges)):
        clauses.append([3 * e + 1, 3 * e + 2, 3 * e + 3])
    for es in inc:
        for bit in range(3):
            a, b, c = [3 * e + bit + 1 for e in es]
            # Even parity of the three coordinate bits.
            clauses += [
                [a, b, -c],
                [a, -b, c],
                [-a, b, c],
                [-a, -b, -c],
            ]
    for old_flow in forbidden:
        # Block precisely this complete assignment.
        clauses.append([
            -(3 * e + bit + 1) if (old_flow[e] >> bit) & 1
            else 3 * e + bit + 1
            for e in range(len(edges)) for bit in range(3)
        ])

    rng = random.Random(seed)
    with Cadical153(bootstrap_with=clauses) as solver:
        phases = [
            var if rng.randrange(2) else -var
            for var in range(1, 3 * len(edges) + 1)
        ]
        solver.set_phases(phases)
        assert solver.solve()
        model = {v for v in solver.get_model() if v > 0}

    flow = [
        sum((1 << bit) for bit in range(3) if 3 * e + bit + 1 in model)
        for e in range(len(edges))
    ]
    assert all(flow)
    assert all(flow[a] ^ flow[b] ^ flow[c] == 0 for a, b, c in inc)
    return flow, inc


def solve_affine(rows, nvars):
    """Return one solution and a nullspace basis for GF(2) augmented rows."""
    coeff_mask = (1 << nvars) - 1
    rows = rows[:]
    pivots = []
    r = 0
    for col in range(nvars):
        pivot = next(
            (i for i in range(r, len(rows)) if (rows[i] >> col) & 1), None
        )
        if pivot is None:
            continue
        rows[r], rows[pivot] = rows[pivot], rows[r]
        for i in range(len(rows)):
            if i != r and ((rows[i] >> col) & 1):
                rows[i] ^= rows[r]
        pivots.append(col)
        r += 1

    for row in rows:
        if not (row & coeff_mask) and ((row >> nvars) & 1):
            return None, []

    pivot_row = {col: rows[i] for i, col in enumerate(pivots)}
    particular = sum(
        1 << col for col in pivots if (pivot_row[col] >> nvars) & 1
    )
    free = [col for col in range(nvars) if col not in pivot_row]
    basis = []
    for f in free:
        vector = 1 << f
        for col in pivots:
            if (pivot_row[col] >> f) & 1:
                vector |= 1 << col
        basis.append(vector)

    assert all(
        ((row & particular).bit_count() & 1) == ((row >> nvars) & 1)
        for row in rows
    )
    assert all(
        all(((row & vector).bit_count() & 1) == 0 for row in rows)
        for vector in basis
    )
    return particular, basis


def binary_rank(values):
    """Return the GF(2)-rank of integers regarded as binary vectors."""
    pivots = {}
    for value in values:
        while value:
            bit = value.bit_length() - 1
            if bit in pivots:
                value ^= pivots[bit]
            else:
                pivots[bit] = value
                break
    return len(pivots)


def cdc_homogeneous_basis(n, edges, flow):
    """Return a basis of the homogeneous CDC compatibility system."""
    nvars = 3 * n + len(edges)
    rows = []
    for ei, (u, v) in enumerate(edges):
        for bit in range(3):
            row = (1 << (3 * u + bit)) | (1 << (3 * v + bit))
            if (flow[ei] >> bit) & 1:
                row |= 1 << (3 * n + ei)
            rows.append(row)

    particular, basis = solve_affine(rows, nvars)
    assert particular == 0 and len(basis) >= 3
    return basis


def cdc_kernel_dimensions(n, edges, flow):
    """Return (weighted-cut-code dimension, CDC homogeneous nullity)."""
    basis = cdc_homogeneous_basis(n, edges, flow)
    return len(basis) - 3, len(basis)


def transvection_equations(n, edges, edge_kinds, flow, local_g, matching):
    """Return the affine noncollision equations for one transvection."""
    equations = []
    for ei, ((u, _), kind) in enumerate(zip(edges, edge_kinds)):
        if kind != "x":
            continue
        allowed = []
        for translation in range(8):
            p = translation ^ local_g[u][ei]
            pair = tuple(sorted((p, p ^ flow[ei])))
            if pair not in matching:
                allowed.append(translation)
        if len(allowed) == 8:
            continue
        assert len(allowed) == 4
        equation = [
            (mask, rhs)
            for mask in range(1, 8) for rhs in range(2)
            if all((((mask & t).bit_count() & 1) == rhs) == (t in allowed)
                   for t in range(8))
        ]
        assert len(equation) == 1
        mask, rhs = equation[0]
        equations.append((mask << (3 * u), rhs))
    return equations


def transvection_restriction_dimensions(
    n, edges, edge_kinds, flow, local_g, matching
):
    """Measure how the CDC kernel changes the dangerous matching equations."""
    equations = transvection_equations(
        n, edges, edge_kinds, flow, local_g, matching
    )
    extra_coefficients = [coefficient for coefficient, _ in equations]

    def restricted_images(vectors):
        return [
            sum(
                (((coefficient & vector).bit_count() & 1) << row)
                for row, coefficient in enumerate(extra_coefficients)
            )
            for vector in vectors
        ]

    kernel = cdc_homogeneous_basis(n, edges, flow)
    translations = [
        sum(1 << (3 * v + bit) for v in range(n))
        for bit in range(3)
    ]
    kernel_rank = binary_rank(restricted_images(kernel))
    translation_rank = binary_rank(restricted_images(translations))
    assert translation_rank <= kernel_rank
    return (
        len(extra_coefficients),
        kernel_rank,
        translation_rank,
        kernel_rank - translation_rank,
    )


def cdc_system(n, edges, flow, inc, rng):
    """Build and solve equation (4) of the OpenAI CDC proof."""
    m = len(edges)
    local_g = [{} for _ in range(n)]
    for v in range(n):
        a, b, c = rng.sample(inc[v], 3)
        local_g[v][a] = 0
        local_g[v][b] = flow[a]
        local_g[v][c] = 0

    nvars = 3 * n + m
    rows = []
    for ei, (u, v) in enumerate(edges):
        d = local_g[u][ei] ^ local_g[v][ei]
        for bit in range(3):
            row = (1 << (3 * u + bit)) | (1 << (3 * v + bit))
            if (flow[ei] >> bit) & 1:
                row |= 1 << (3 * n + ei)
            if (d >> bit) & 1:
                row |= 1 << nvars
            rows.append(row)

    particular, basis = solve_affine(rows, nvars)
    assert particular is not None
    return local_g, particular, basis


def palette(n, edges, flow, local_g, solution):
    """Return the simple graph on F_2^3 whose edges are the sets P_e."""
    t = [
        sum(((solution >> (3 * v + bit)) & 1) << bit for bit in range(3))
        for v in range(n)
    ]
    adj = [0] * 8
    pairs = []
    for ei, (u, v) in enumerate(edges):
        p = t[u] ^ local_g[u][ei]
        q = p ^ flow[ei]
        p2 = t[v] ^ local_g[v][ei]
        q2 = p2 ^ flow[ei]
        assert {p, q} == {p2, q2}
        adj[p] |= 1 << q
        adj[q] |= 1 << p
        pairs.append((p, q))

    inc_pairs = [[] for _ in range(n)]
    for ei, (u, v) in enumerate(edges):
        inc_pairs[u].append(pairs[ei])
        inc_pairs[v].append(pairs[ei])
    for pairs_at_v in inc_pairs:
        counts = [sum(s in pair for pair in pairs_at_v) for s in range(8)]
        assert all(count in (0, 2) for count in counts)
    return adj


def chromatic_number(adj):
    """Compute an exact colouring by subset DP (there are only 8 vertices)."""
    independent = [True] * 256
    for mask in range(1, 256):
        vbit = mask & -mask
        v = vbit.bit_length() - 1
        independent[mask] = (
            independent[mask ^ vbit] and not (adj[v] & (mask ^ vbit))
        )

    dp = [9] * 256
    choice = [0] * 256
    dp[0] = 0
    for mask in range(1, 256):
        anchor = mask & -mask
        sub = mask
        while sub:
            if (
                (sub & anchor)
                and independent[sub]
                and 1 + dp[mask ^ sub] < dp[mask]
            ):
                dp[mask] = 1 + dp[mask ^ sub]
                choice[mask] = sub
            sub = (sub - 1) & mask

    colors = [-1] * 8
    mask = 255
    color = 0
    while mask:
        sub = choice[mask]
        for v in range(8):
            if (sub >> v) & 1:
                colors[v] = color
        mask ^= sub
        color += 1
    return dp[255], colors


def complement_perfect_matching(adj):
    """A missing perfect matching is a simple sufficient 4-colour certificate."""
    def rec(mask, pairs):
        if not mask:
            return pairs
        u = (mask & -mask).bit_length() - 1
        rest = mask & ~(1 << u)
        for v in range(u + 1, 8):
            if ((rest >> v) & 1) and not ((adj[u] >> v) & 1):
                answer = rec(rest & ~(1 << v), pairs + [(u, v)])
                if answer is not None:
                    return answer
        return None

    return rec(255, [])


def gf2_rank(values):
    basis = [0, 0, 0]
    rank = 0
    for value in values:
        x = value
        while x:
            bit = x.bit_length() - 1
            if basis[bit]:
                x ^= basis[bit]
            else:
                basis[bit] = x
                rank += 1
                break
    return rank


def find_tetrahedral_palette(name, n, edges):
    """Find a full-rank palette supported on an affine tetrahedron.

    Each graph edge receives one of the six edges of a tetrahedron.  At each
    cubic vertex the three labels must form one of its four triangular faces.
    Requiring two different face types forces the induced F_2^3-flow to have
    rank three.  This is a SAT formulation of a Tait colouring, not a new way
    to find one; it is included to make the survey's diagnostic reproducible.
    """
    inc = [[] for _ in range(n)]
    for ei, (u, v) in enumerate(edges):
        inc[u].append(ei)
        inc[v].append(ei)
    assert all(len(es) == 3 for es in inc)

    tetrahedron = (0, 1, 2, 4)
    pair_indices = tuple(
        (i, j) for i in range(4) for j in range(i + 1, 4)
    )
    m = len(edges)

    def edge_var(ei, label):
        return 6 * ei + label + 1

    def face_var(v, omitted):
        return 6 * m + 4 * v + omitted + 1

    clauses = []
    for ei in range(m):
        choices = [edge_var(ei, label) for label in range(6)]
        clauses.append(choices)
        clauses.extend(
            [-choices[i], -choices[j]]
            for i in range(6) for j in range(i + 1, 6)
        )

    for v in range(n):
        faces = [face_var(v, omitted) for omitted in range(4)]
        clauses.append(faces)
        clauses.extend(
            [-faces[i], -faces[j]]
            for i in range(4) for j in range(i + 1, 4)
        )
        for label in range(6):
            for i in range(3):
                for j in range(i + 1, 3):
                    clauses.append([
                        -edge_var(inc[v][i], label),
                        -edge_var(inc[v][j], label),
                    ])
        for omitted in range(4):
            for label, pair in enumerate(pair_indices):
                if omitted in pair:
                    for ei in inc[v]:
                        clauses.append([
                            -face_var(v, omitted), -edge_var(ei, label)
                        ])

    # Break the tetrahedral symmetry and force at least two local face types.
    clauses.append([face_var(0, 0)])
    clauses.append([
        face_var(v, omitted)
        for v in range(n) for omitted in range(1, 4)
    ])

    with Cadical153(bootstrap_with=clauses) as solver:
        assert solver.solve()
        model = {literal for literal in solver.get_model() if literal > 0}

    labels = [
        next(label for label in range(6) if edge_var(ei, label) in model)
        for ei in range(m)
    ]
    faces = [
        next(omitted for omitted in range(4) if face_var(v, omitted) in model)
        for v in range(n)
    ]
    flow = [
        tetrahedron[pair_indices[label][0]]
        ^ tetrahedron[pair_indices[label][1]]
        for label in labels
    ]
    assert all(flow[es[0]] ^ flow[es[1]] ^ flow[es[2]] == 0 for es in inc)
    assert len(set(faces)) >= 2
    assert gf2_rank(flow) == 3

    label_counts = [labels.count(label) for label in range(6)]
    face_counts = [faces.count(omitted) for omitted in range(4)]
    print(
        name,
        "TETRAHEDRAL-SAT",
        "flow_rank", gf2_rank(flow),
        "clauses", len(clauses),
        "face_counts", face_counts,
        "label_counts", label_counts,
        flush=True,
    )


def find_separated_transvection_palette(name, n, edges, edge_kinds):
    """Find a separated palette certified by a nonparallel matching.

    Bit zero is the layer coordinate and the other two bits identify a point
    in F_2^2.  The chosen cross-layer matching is induced by the transvection
    (0,1,2,3) -> (0,1,3,2).  We additionally force all four lower-coordinate
    differences to occur on x-edges, which rules out every parallel matching.
    """
    inc = [[] for _ in range(n)]
    for ei, (u, v) in enumerate(edges):
        inc[u].append(ei)
        inc[v].append(ei)
    assert all(len(es) == 3 for es in inc)
    assert len(edge_kinds) == len(edges)

    palette_pairs = tuple(combinations(range(8), 2))
    pair_index = {pair: i for i, pair in enumerate(palette_pairs)}
    matching_map = (0, 1, 3, 2)
    matching = {
        tuple(sorted((2 * z, 1 + 2 * matching_map[z])))
        for z in range(4)
    }

    # A local triangle has two points in one layer and one in the other.
    # Record its within-layer edge first and its ordered cross edges second.
    states = []
    for majority_layer in range(2):
        majority = [p for p in range(8) if (p & 1) == majority_layer]
        minority = [p for p in range(8) if (p & 1) != majority_layer]
        for p, q in combinations(majority, 2):
            for r in minority:
                pairs = [
                    tuple(sorted(pair))
                    for pair in ((p, q), (p, r), (q, r))
                ]
                if any(pair in matching for pair in pairs):
                    continue
                cross = pairs[1:]
                states.append(tuple(pair_index[pair] for pair in pairs))
                states.append((
                    pair_index[pairs[0]], pair_index[cross[1]],
                    pair_index[cross[0]],
                ))
    assert len(states) == 48 and len(set(states)) == 48

    m = len(edges)

    def edge_var(ei, label):
        return 28 * ei + label + 1

    state_offset = 28 * m

    def state_var(v, state):
        return state_offset + 48 * v + state + 1

    allowed_labels = {"a": [], "x": []}
    for label, (p, q) in enumerate(palette_pairs):
        same_layer = (p & 1) == (q & 1)
        if same_layer:
            allowed_labels["a"].append(label)
        elif (p, q) not in matching:
            allowed_labels["x"].append(label)
    assert all(len(labels) == 12 for labels in allowed_labels.values())

    clauses = []
    for ei, kind in enumerate(edge_kinds):
        choices = [edge_var(ei, label) for label in allowed_labels[kind]]
        clauses.append(choices)
        clauses.extend(
            [-choices[i], -choices[j]]
            for i in range(12) for j in range(i + 1, 12)
        )

    for v in range(n):
        a_edge = next(e for e in inc[v] if edge_kinds[e] == "a")
        x_edges = [e for e in inc[v] if edge_kinds[e] == "x"]
        clauses.append([state_var(v, state) for state in range(48)])
        for state, labels in enumerate(states):
            assignment = zip((a_edge, x_edges[0], x_edges[1]), labels)
            for ei, label in assignment:
                clauses.append([-state_var(v, state), edge_var(ei, label)])

    # No parallel cross-layer matching can certify the resulting palette.
    for lower_difference in range(4):
        clauses.append([
            edge_var(ei, label)
            for ei, kind in enumerate(edge_kinds) if kind == "x"
            for label in allowed_labels["x"]
            if ((palette_pairs[label][0] ^ palette_pairs[label][1]) >> 1)
            == lower_difference
        ])

    with Cadical153(bootstrap_with=clauses) as solver:
        assert solver.solve()
        model = {literal for literal in solver.get_model() if literal > 0}

    labels = [
        next(
            label for label in allowed_labels[kind]
            if edge_var(ei, label) in model
        )
        for ei, kind in enumerate(edge_kinds)
    ]
    pairs = [palette_pairs[label] for label in labels]
    flow = [p ^ q for p, q in pairs]
    for v in range(n):
        local_pairs = [pairs[ei] for ei in inc[v]]
        local_points = {p for pair in local_pairs for p in pair}
        assert len(local_points) == 3
        assert all(
            sum(point in pair for pair in local_pairs) == 2
            for point in local_points
        )
        assert flow[inc[v][0]] ^ flow[inc[v][1]] ^ flow[inc[v][2]] == 0
    assert gf2_rank(flow) == 3
    assert all(pair not in matching for pair in pairs)
    assert all(
        (flow[ei] & 1) == (kind == "x")
        for ei, kind in enumerate(edge_kinds)
    )
    x_differences = {
        (flow[ei] >> 1) for ei, kind in enumerate(edge_kinds) if kind == "x"
    }
    assert x_differences == set(range(4))

    # Verify explicitly that these local triangles are an output of the
    # paper's affine compatibility system, not merely a parity assignment.
    local_g = [{} for _ in range(n)]
    translations = []
    for v, es in enumerate(inc):
        local_g[v][es[0]] = 0
        local_g[v][es[1]] = flow[es[0]]
        local_g[v][es[2]] = 0
        translation = next(
            t for t in range(8)
            if all(
                {t ^ local_g[v][ei], t ^ local_g[v][ei] ^ flow[ei]}
                == set(pairs[ei])
                for ei in es
            )
        )
        translations.append(translation)
    for ei, (u, v) in enumerate(edges):
        discrepancy = (
            translations[u] ^ translations[v]
            ^ local_g[u][ei] ^ local_g[v][ei]
        )
        assert discrepancy in (0, flow[ei])

    inverse_matching = [matching_map.index(z) for z in range(4)]

    def matching_color(point):
        lower = point >> 1
        return lower if not (point & 1) else inverse_matching[lower]

    tait_flow = [matching_color(p) ^ matching_color(q) for p, q in pairs]
    assert all(tait_flow)
    assert all(
        tait_flow[es[0]] ^ tait_flow[es[1]] ^ tait_flow[es[2]] == 0
        for es in inc
    )
    weighted_cut_dim, cdc_nullity = cdc_kernel_dimensions(n, edges, flow)
    (
        dangerous_x_edges,
        restriction_rank,
        translation_rank,
        cut_restriction_rank,
    ) = transvection_restriction_dimensions(
        n, edges, edge_kinds, flow, local_g, matching
    )

    print(
        name,
        "SEPARATED-TRANSVECTION-SAT",
        "flow_rank", gf2_rank(flow),
        "clauses", len(clauses),
        "weighted_cut_dim", weighted_cut_dim,
        "cdc_nullity", cdc_nullity,
        "dangerous_x_edges", dangerous_x_edges,
        "restriction_rank", restriction_rank,
        "translation_rank", translation_rank,
        "cut_restriction_rank", cut_restriction_rank,
        "palette_edges", len(set(pairs)),
        "x_lower_differences", sorted(x_differences),
        flush=True,
    )


def affine_solutions(particular, basis):
    """Enumerate an affine GF(2)-space once each, using Gray-code updates."""
    solution = particular
    previous_gray = 0
    yield solution
    for i in range(1, 1 << len(basis)):
        gray = i ^ (i >> 1)
        changed = gray ^ previous_gray
        solution ^= basis[changed.bit_length() - 1]
        yield solution
        previous_gray = gray


def examine(name, n, edges, trials, required_rank=3):
    rng = random.Random(20260904)
    best = (9, 99, None, None)
    histogram = {}
    total = 0
    seen_flows = []

    for trial in range(trials):
        while True:
            flow, inc = nz_three_bit_flow(
                n, edges, seed=rng.randrange(1 << 30), forbidden=seen_flows
            )
            seen_flows.append(tuple(flow))
            if gf2_rank(flow) == required_rank:
                break

        local_g, particular, basis = cdc_system(n, edges, flow, inc, rng)
        for solution in affine_solutions(particular, basis):
            adj = palette(n, edges, flow, local_g, solution)
            chi, coloring = chromatic_number(adj)
            nedges = sum(row.bit_count() for row in adj) // 2
            missing_matching = complement_perfect_matching(adj)
            histogram[chi] = histogram.get(chi, 0) + 1
            total += 1
            if (chi, nedges) < best[:2]:
                best = (chi, nedges, coloring, missing_matching)

        print(
            name,
            "trial", trial,
            "flow_rank", gf2_rank(flow),
            "nullity", len(basis),
            "solutions", 1 << len(basis),
            "best", best[:2],
            "hist", histogram,
            flush=True,
        )

    print(
        name, "DONE", "solutions", total, "best", best,
        "hist", histogram, flush=True
    )
    return best[0]


def check_color_aligned_flow(name, n, edges):
    ok, edge_colors = three_edge_colourable(n, edges)
    assert ok
    flow = [(1, 2, 3)[color] for color in edge_colors]
    inc = [[] for _ in range(n)]
    for ei, (u, v) in enumerate(edges):
        inc[u].append(ei)
        inc[v].append(ei)

    rng = random.Random(20260904)
    local_g, particular, basis = cdc_system(n, edges, flow, inc, rng)
    adj = palette(n, edges, flow, local_g, particular)
    chi, coloring = chromatic_number(adj)
    assert gf2_rank(flow) == 2 and chi <= 4
    print(
        name,
        "COLOR-ALIGNED",
        "flow_rank", gf2_rank(flow),
        "nullity", len(basis),
        "chi", chi,
        "palette_edges", sum(row.bit_count() for row in adj) // 2,
        "coloring", coloring,
        flush=True,
    )


def petersen():
    edges = []
    for i in range(5):
        edges.append((i, (i + 1) % 5))
        edges.append((5 + i, 5 + (i + 2) % 5))
        edges.append((i, 5 + i))
    return 10, edges


def complete_four():
    return 4, [
        (0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)
    ]


def cycle_perm(n, points):
    p = list(range(n))
    for a, b in zip(points, points[1:] + points[:1]):
        p[a] = b
    return tuple(p)


def cayley_graph(gens, return_edge_kinds=False):
    group = list(closure(gens))
    index = {g: i for i, g in enumerate(group)}
    seen = set()
    edges = []
    edge_kinds = []
    for g in group:
        u = index[g]
        for kind, generator in enumerate(gens):
            v = index[mul(g, generator)]
            edge = (min(u, v), max(u, v))
            if edge not in seen:
                seen.add(edge)
                edges.append(edge)
                edge_kinds.append("a" if kind == 0 else "x")
    assert len(edges) * 2 == 3 * len(group)
    if return_edge_kinds:
        return len(group), edges, edge_kinds
    return len(group), edges


def random_generator_separated_flow(
    n, edges, edge_kinds, rng, forbidden_quotient_flows=()
):
    """Construct Proposition 9.5's flow from a random quotient 4-flow."""
    parent = list(range(n))

    def find(v):
        while parent[v] != v:
            parent[v] = parent[parent[v]]
            v = parent[v]
        return v

    def union(u, v):
        u, v = find(u), find(v)
        if u != v:
            parent[v] = u

    for (u, v), kind in zip(edges, edge_kinds):
        if kind == "x":
            union(u, v)
    roots = sorted({find(v) for v in range(n)})
    component = {root: i for i, root in enumerate(roots)}
    vertex_component = [component[find(v)] for v in range(n)]

    a_edges = [ei for ei, kind in enumerate(edge_kinds) if kind == "a"]
    a_position = {ei: pos for pos, ei in enumerate(a_edges)}
    quotient_inc = [[] for _ in roots]
    for ei in a_edges:
        u, v = edges[ei]
        cu, cv = vertex_component[u], vertex_component[v]
        assert cu != cv
        quotient_inc[cu].append(ei)
        quotient_inc[cv].append(ei)

    def quotient_var(ei, bit):
        return 2 * a_position[ei] + bit + 1

    clauses = []
    for ei in a_edges:
        clauses.append([quotient_var(ei, 0), quotient_var(ei, 1)])
    for incident in quotient_inc:
        for bit in range(2):
            variables = [quotient_var(ei, bit) for ei in incident]
            for assignment in product(range(2), repeat=len(variables)):
                if sum(assignment) & 1:
                    clauses.append([
                        -var if value else var
                        for var, value in zip(variables, assignment)
                    ])
    for old_flow in forbidden_quotient_flows:
        clauses.append([
            -quotient_var(ei, bit)
            if (old_flow[pos] >> bit) & 1 else quotient_var(ei, bit)
            for pos, ei in enumerate(a_edges) for bit in range(2)
        ])

    with Cadical153(bootstrap_with=clauses) as solver:
        solver.set_phases([
            var if rng.randrange(2) else -var
            for var in range(1, 2 * len(a_edges) + 1)
        ])
        assert solver.solve()
        model = {literal for literal in solver.get_model() if literal > 0}
    quotient_flow = {
        ei: sum(
            (1 << bit) for bit in range(2)
            if quotient_var(ei, bit) in model
        )
        for ei in a_edges
    }
    assert all(quotient_flow[ei] for ei in a_edges)
    assert all(
        not (sum((quotient_flow[ei] >> bit) & 1 for ei in incident) & 1)
        for incident in quotient_inc for bit in range(2)
    )
    quotient_signature = tuple(quotient_flow[ei] for ei in a_edges)

    x_edges = [ei for ei, kind in enumerate(edge_kinds) if kind == "x"]
    x_position = {ei: pos for pos, ei in enumerate(x_edges)}
    inc = [[] for _ in range(n)]
    for ei, (u, v) in enumerate(edges):
        inc[u].append(ei)
        inc[v].append(ei)

    nvars = 2 * len(x_edges)
    rows = []
    for incident in inc:
        a_edge = next(ei for ei in incident if edge_kinds[ei] == "a")
        cycle_edges = [ei for ei in incident if edge_kinds[ei] == "x"]
        for bit in range(2):
            row = sum(
                1 << (2 * x_position[ei] + bit) for ei in cycle_edges
            )
            if (quotient_flow[a_edge] >> bit) & 1:
                row |= 1 << nvars
            rows.append(row)
    particular, basis = solve_affine(rows, nvars)
    assert particular is not None
    lift = particular
    for vector in basis:
        if rng.randrange(2):
            lift ^= vector

    flow = []
    for ei, kind in enumerate(edge_kinds):
        if kind == "a":
            flow.append(quotient_flow[ei] << 1)
        else:
            lower = sum(
                ((lift >> (2 * x_position[ei] + bit)) & 1) << bit
                for bit in range(2)
            )
            flow.append(1 | (lower << 1))
    assert gf2_rank(flow) == 3
    assert all(flow[es[0]] ^ flow[es[1]] ^ flow[es[2]] == 0 for es in inc)
    return flow, inc, quotient_signature


def all_transvection_matchings():
    """Return the twelve affine matchings whose linear part has order two."""
    linear_maps = []
    for image_one in range(1, 4):
        for image_two in range(1, 4):
            if image_one == image_two:
                continue
            linear = tuple(
                (image_one if z & 1 else 0) ^ (image_two if z & 2 else 0)
                for z in range(4)
            )
            if linear != tuple(range(4)) and all(
                linear[linear[z]] == z for z in range(4)
            ):
                linear_maps.append(linear)
    assert len(linear_maps) == 3
    return [
        (
            linear,
            offset,
            {
                tuple(sorted((2 * z, 1 + 2 * (linear[z] ^ offset))))
                for z in range(4)
            },
        )
        for linear in linear_maps for offset in range(4)
    ]


def solve_transvection_target(
    n, edges, edge_kinds, flow, local_g, particular, basis, matching
):
    """Solve one transvection syndrome on the CDC affine solution space."""
    equations = transvection_equations(
        n, edges, edge_kinds, flow, local_g, matching
    )
    reduced_rows = []
    for coefficient, rhs in equations:
        row = sum(
            (((coefficient & vector).bit_count() & 1) << i)
            for i, vector in enumerate(basis)
        )
        residual = rhs ^ ((coefficient & particular).bit_count() & 1)
        if residual:
            row |= 1 << len(basis)
        reduced_rows.append(row)
    choice, _ = solve_affine(reduced_rows, len(basis))
    if choice is None:
        return None, len(equations)
    solution = particular
    for i, vector in enumerate(basis):
        if (choice >> i) & 1:
            solution ^= vector
    return solution, len(equations)


def examine_random_separated_flows(name, n, edges, edge_kinds, trials=12):
    """Try the paper's method without encoding a desired Tait colouring."""
    rng = random.Random(20260905)
    matchings = all_transvection_matchings()
    total_hits = 0
    flows_with_hit = 0
    quotient_flows = []
    for trial in range(trials):
        flow, inc, quotient_signature = random_generator_separated_flow(
            n, edges, edge_kinds, rng, quotient_flows
        )
        quotient_flows.append(quotient_signature)
        local_g, particular, basis = cdc_system(n, edges, flow, inc, rng)
        hits = 0
        dangerous_counts = []
        for _, _, matching in matchings:
            solution, dangerous = solve_transvection_target(
                n, edges, edge_kinds, flow, local_g,
                particular, basis, matching
            )
            dangerous_counts.append(dangerous)
            if solution is None:
                continue
            adj = palette(n, edges, flow, local_g, solution)
            assert all(not ((adj[p] >> q) & 1) for p, q in matching)
            assert chromatic_number(adj)[0] <= 4
            hits += 1
        total_hits += hits
        flows_with_hit += bool(hits)
        print(
            name,
            "SEPARATED-RANDOM", trial,
            "weighted_cut_dim", len(basis) - 3,
            "transvection_hits", hits,
            "dangerous_range", (min(dangerous_counts), max(dangerous_counts)),
            flush=True,
        )
    print(
        name,
        "SEPARATED-RANDOM-DONE",
        "flows_with_hit", flows_with_hit,
        "trials", trials,
        "matching_hits", total_hits,
        flush=True,
    )


def main():
    n, edges = complete_four()
    assert examine("K4", n, edges, trials=5) == 4

    n, edges = petersen()
    assert examine("Petersen", n, edges, trials=3) > 4

    a = cycle_perm(5, [0, 1])
    x = cycle_perm(5, [0, 1, 2, 3, 4])
    n, edges, edge_kinds = cayley_graph([a, x], return_edge_kinds=True)
    assert n == 120
    assert examine(
        "Cay(S5; transposition, 5-cycle)", n, edges, trials=20
    ) > 4
    find_tetrahedral_palette(
        "Cay(S5; transposition, 5-cycle)", n, edges
    )
    find_separated_transvection_palette(
        "Cay(S5; transposition, 5-cycle)", n, edges, edge_kinds
    )
    examine_random_separated_flows(
        "Cay(S5; transposition, 5-cycle)", n, edges, edge_kinds
    )
    check_color_aligned_flow("Cay(S5; transposition, 5-cycle)", n, edges)


if __name__ == "__main__":
    main()
