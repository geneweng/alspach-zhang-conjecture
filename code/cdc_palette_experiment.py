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
to check the elementary positive direction.

Run from the repository root with:

    python3 code/cdc_palette_experiment.py
"""

import random

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


def cayley_graph(gens):
    group = list(closure(gens))
    index = {g: i for i, g in enumerate(group)}
    seen = set()
    edges = []
    for g in group:
        u = index[g]
        for generator in gens:
            v = index[mul(g, generator)]
            edge = (min(u, v), max(u, v))
            if edge not in seen:
                seen.add(edge)
                edges.append(edge)
    assert len(edges) * 2 == 3 * len(group)
    return len(group), edges


def main():
    n, edges = complete_four()
    assert examine("K4", n, edges, trials=5) == 4

    n, edges = petersen()
    assert examine("Petersen", n, edges, trials=3) > 4

    a = cycle_perm(5, [0, 1])
    x = cycle_perm(5, [0, 1, 2, 3, 4])
    n, edges = cayley_graph([a, x])
    assert n == 120
    assert examine(
        "Cay(S5; transposition, 5-cycle)", n, edges, trials=20
    ) > 4
    check_color_aligned_flow("Cay(S5; transposition, 5-cycle)", n, edges)


if __name__ == "__main__":
    main()
