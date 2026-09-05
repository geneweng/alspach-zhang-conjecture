#!/usr/bin/env python3
"""Check the quotient bilinear form against the original vertex equations.

Run: python3 code/layer_potential_experiment.py

The input quotient flow is sampled without a colouring constraint.  The
matrix search changes only its lifts, then tests all transvection targets.
"""

import random
from collections import Counter

from cdc_palette_experiment import (
    all_transvection_matchings, alternating_form,
    binary_rank, cayley_graph, cdc_system, chromatic_number, cycle_perm,
    layer_potential_basis, local_interlacement, palette, petersen,
    quotient_weighted_cut_dimension, random_generator_separated_flow,
    solve_affine, solve_transvection_target, xor_selected,
)


class QuotientForm:
    """An affine family of symmetric forms on the quotient subflow space."""

    def __init__(self, n, edges, kinds, flow):
        self.n, self.edges, self.kinds = n, edges, kinds
        self.a_edges = [e for e, kind in enumerate(kinds) if kind == 'a']
        self.a_position = {e: i for i, e in enumerate(self.a_edges)}
        self.values = [flow[e] >> 1 for e in self.a_edges]
        self.inc = [[] for _ in range(n)]
        x_adj = [[] for _ in range(n)]
        a_at = [None] * n
        for e, ((u, v), kind) in enumerate(zip(edges, kinds)):
            self.inc[u].append(e)
            self.inc[v].append(e)
            if kind == 'x':
                x_adj[u].append((v, e))
                x_adj[v].append((u, e))
            else:
                a_at[u] = a_at[v] = self.a_position[e]

        self.cycles, self.cycle_edges = [], []
        unseen = set(range(n))
        while unseen:
            start = min(unseen)
            order, x_edges = [], []
            previous, current = None, start
            while current not in order:
                order.append(current)
                unseen.remove(current)
                following, e = next(pair for pair in x_adj[current]
                                    if pair[0] != previous)
                x_edges.append(e)
                previous, current = current, following
            assert current == start
            self.cycles.append(order)
            self.cycle_edges.append(x_edges)

        self.q = len(self.cycles)
        self.component = [None] * n
        self.orders = [[a_at[v] for v in cycle] for cycle in self.cycles]
        self.reference = [0] * len(edges)
        constraints = []
        for C, (vertices, order, x_edges) in enumerate(zip(
                self.cycles, self.orders, self.cycle_edges)):
            y = 0
            for v, position, e in zip(vertices, order, x_edges):
                self.component[v] = C
                y ^= self.values[position]
                self.reference[e] = 1 | (y << 1)
            assert y == 0
            for bit in range(2):
                constraints.append(sum(1 << p for p in order
                                       if (self.values[p] >> bit) & 1))
        for e, value in zip(self.a_edges, self.values):
            self.reference[e] = value << 1
        _, self.basis = solve_affine(constraints, len(self.a_edges))
        self.d = len(self.basis)
        self.full = (1 << len(self.a_edges)) - 1
        self.base = [0] * self.d
        self.changes = [[0] * self.d for _ in range(2 * self.q)]
        for i, r in enumerate(self.basis):
            for j, t in enumerate(self.basis):
                self.base[i] |= self.omega(r, t) << j
                for C, order in enumerate(self.orders):
                    p = xor_selected(self.values, r & t & sum(1 << e for e in order))
                    for bit in range(2):
                        self.changes[2 * C + bit][i] |= alternating_form(p, 1 << bit) << j
        assert all(((self.base[i] >> j) & 1) == ((self.base[j] >> i) & 1)
                   for i in range(self.d) for j in range(self.d))
        assert all(not ((change[i] >> i) & 1)
                   for change in self.changes for i in range(self.d))
        assert self.omega(self.full, self.full) == 0
        for C, order in enumerate(self.orders):
            if len(order) != 5:
                continue
            values = [self.values[e] for e in order]
            repeated = next(p for p in (1, 2, 3) if values.count(p) == 3)
            first, second = self.changes[2 * C:2 * C + 2]
            assert binary_rank(first) <= 2 and binary_rank(second) <= 2
            assert all(((a if repeated & 1 else 0)
                        ^ (b if repeated & 2 else 0)) == 0
                       for a, b in zip(first, second))

    def omega(self, r, t):
        result = 0
        for order in self.orders:
            values = [self.values[e] for e in order]
            local_r = sum(((r >> e) & 1) << i for i, e in enumerate(order))
            local_t = sum(((t >> e) & 1) << i for i, e in enumerate(order))
            result = result ^ local_interlacement(values, local_r, local_t)
        return result

    def coordinates(self, r):
        rows = [sum(((b >> e) & 1) << i for i, b in enumerate(self.basis))
                | (((r >> e) & 1) << self.d)
                for e in range(len(self.a_edges))]
        answer, _ = solve_affine(rows, self.d)
        assert answer is not None
        return answer

    def matrix(self, corrections):
        rows = self.base[:]
        for bit, change in enumerate(self.changes):
            if (corrections >> bit) & 1:
                rows = [a ^ b for a, b in zip(rows, change)]
        return rows

    def flow(self, corrections):
        result = self.reference[:]
        for C, x_edges in enumerate(self.cycle_edges):
            c = (corrections >> (2 * C)) & 3
            for e in x_edges:
                result[e] ^= c << 1
        assert all(result[es[0]] ^ result[es[1]] ^ result[es[2]] == 0
                   for es in self.inc)
        return result

    def choose_lifts(self, indicators):
        """Solve the small matrix equations placing indicators in its kernel."""
        rows = []
        for r in indicators:
            vector = self.coordinates(r)
            for j in range(self.d):
                row = sum(((change[j] & vector).bit_count() & 1) << bit
                          for bit, change in enumerate(self.changes))
                row |= ((self.base[j] & vector).bit_count() & 1) << (2 * self.q)
                rows.append(row)
        return solve_affine(rows, 2 * self.q)

    def direct_lifts(self, indicators):
        """Independent check using w_u+w_v equations on the cubic graph."""
        nvars = 2 * self.q + 2 * self.n * len(indicators)
        rows = []
        for k, r in enumerate(indicators):
            vertex_r = [0] * self.n
            for position, e in enumerate(self.a_edges):
                for v in self.edges[e]:
                    vertex_r[v] = (r >> position) & 1
            offset = 2 * self.q + 2 * self.n * k
            for e, ((u, v), kind) in enumerate(zip(self.edges, self.kinds)):
                if kind == 'a':
                    p = self.reference[e] >> 1
                    row = 0
                    for bit in range(2):
                        if alternating_form(p, 1 << bit):
                            row ^= (1 << (offset + 2 * u + bit)) ^ (1 << (offset + 2 * v + bit))
                    rows.append(row)
                else:
                    delta = vertex_r[u] ^ vertex_r[v]
                    y = self.reference[e] >> 1
                    for bit in range(2):
                        row = (1 << (offset + 2 * u + bit)) ^ (1 << (offset + 2 * v + bit))
                        if delta:
                            row ^= 1 << (2 * self.component[u] + bit)
                            row ^= ((y >> bit) & 1) << nvars
                        rows.append(row)
        return solve_affine(rows, nvars)[0]


def verify_form(model, samples=40):
    rng = random.Random(173)
    statuses = Counter()
    for count in (1, 2, 3):
        for _ in range(samples if count == 1 else 6):
            indicators = [xor_selected(model.basis, rng.randrange(1 << model.d))
                          for _ in range(count)]
            corrections, _ = model.choose_lifts(indicators)
            direct = model.direct_lifts(indicators)
            assert (corrections is None) == (direct is None)
            statuses[count, corrections is not None] += 1
            if corrections is not None:
                lifted = model.flow(corrections)
                actual = layer_potential_basis(model.n, model.edges, model.kinds, lifted)
                expected = model.d - binary_rank(model.matrix(corrections)) - 1
                assert len(actual) == expected
                all_vertices = (1 << model.n) - 1
                for r in indicators:
                    mask = sum(((r >> position) & 1) << v
                               for position, e in enumerate(model.a_edges)
                               for v in model.edges[e])
                    if mask & 1:
                        mask ^= all_vertices
                    assert binary_rank(actual + [mask]) == len(actual)
    for _ in range(8):
        corrections = rng.randrange(1 << (2 * model.q))
        actual = layer_potential_basis(model.n, model.edges, model.kinds,
                                       model.flow(corrections))
        assert len(actual) == model.d - binary_rank(model.matrix(corrections)) - 1
        _, radical = solve_affine(model.matrix(corrections), model.d)
        if len(radical) >= 2:
            indicators = [xor_selected(model.basis, vector) for vector in radical[:3]]
            assert model.choose_lifts(indicators)[0] is not None
            assert model.direct_lifts(indicators) is not None
    print('FORM-CHECK', 'subflow_dim', model.d, 'direct_checks', dict(statuses), flush=True)


def improve_lifts(model, restarts=40, steps=60):
    """Minimise matrix rank; no colouring equation participates in this search."""
    rng = random.Random(541)
    active = [bit for bit, change in enumerate(model.changes) if any(change)]
    best = None
    for _ in range(restarts):
        corrections = rng.randrange(1 << (2 * model.q))
        matrix = model.matrix(corrections)
        for _ in range(steps):
            rank = binary_rank(matrix)
            if best is None or rank < best[0]:
                best = rank, corrections
            candidates = []
            for bit in active:
                candidate = [a ^ b for a, b in zip(matrix, model.changes[bit])]
                candidates.append((binary_rank(candidate), bit, candidate))
            target = min(item[0] for item in candidates)
            choices = [item for item in candidates if item[0] == target]
            if target > rank:
                break
            _, bit, matrix = rng.choice(choices)
            corrections ^= 1 << bit
    return best


def main():
    pn, pedges = petersen()
    pkinds = ['a' if (u < 5) != (v < 5) else 'x' for u, v in pedges]
    pflow, _, _ = random_generator_separated_flow(
        pn, pedges, pkinds, random.Random(132))
    verify_form(QuotientForm(pn, pedges, pkinds, pflow), samples=12)

    a = cycle_perm(5, [0, 1])
    x = cycle_perm(5, [0, 1, 2, 3, 4])
    n, edges, kinds = cayley_graph([a, x], return_edge_kinds=True)
    flow, _, _ = random_generator_separated_flow(
        n, edges, kinds, random.Random(20260904))
    model = QuotientForm(n, edges, kinds, flow)
    verify_form(model)
    rank, corrections = improve_lifts(model)
    lifted = model.flow(corrections)
    layer_dim = len(layer_potential_basis(n, edges, kinds, lifted))
    quotient_dim = quotient_weighted_cut_dimension(n, edges, kinds, lifted)
    assert layer_dim == model.d - rank - 1
    _, radical = solve_affine(model.matrix(corrections), model.d)
    indicators = [xor_selected(model.basis, vector) for vector in radical]
    assert model.direct_lifts(indicators) is not None
    local_g, particular, basis = cdc_system(n, edges, lifted, model.inc, random.Random(901))
    hits = 0
    for _, _, matching in all_transvection_matchings():
        solution, _ = solve_transvection_target(
            n, edges, kinds, lifted, local_g, particular, basis, matching)
        if solution is not None:
            assert chromatic_number(palette(n, edges, lifted, local_g, solution))[0] <= 4
            hits += 1
    print('FORM-SEARCH', 'rank', rank, 'layer_dim', layer_dim,
          'quotient_cut_dim', quotient_dim, 'transvection_hits', hits, flush=True)


if __name__ == '__main__':
    main()
