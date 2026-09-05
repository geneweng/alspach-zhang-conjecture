#!/usr/bin/env python3
"""Explore C(a)-invariant minimum-face matchings for two-cycle generators.

x has two coprime odd cycles of lengths m,n and a joins one point of each.
The C(a)-quotient has the two-subsets as vertices. We enumerate its minimum
face independently as perfect matchings of the contracted x-orbit pregraph.
"""

import argparse
import math
from collections import Counter
from itertools import combinations

from cayley_snark_check import inv, mul, order
from cdc_palette_experiment import cycle_perm


class TwoCycleQuotient:
    def __init__(self, m, n):
        assert 1 <= m <= n and m % 2 == n % 2 == 1 and math.gcd(m, n) == 1
        self.m, self.n, self.N = m, n, m + n
        self.a = cycle_perm(self.N, [0, m])
        self.x = mul(cycle_perm(self.N, list(range(m))),
                     cycle_perm(self.N, list(range(m, m + n))))
        self.xinv = inv(self.x)
        self.identity = tuple(range(self.N))
        self.vertices = list(combinations(range(self.N), 2))
        self.index = {v: i for i, v in enumerate(self.vertices)}
        self.neighbours = [[self.index[tuple(sorted(p[i] for i in pair))]
                            for p in (self.a, self.x, self.xinv)]
                           for pair in self.vertices]
        self.edges, self.kinds, self.inc = [], [], [[] for _ in self.vertices]
        edge_index = {}
        for v, neighbours in enumerate(self.neighbours):
            for kind, w in enumerate(neighbours[:2]):
                key = (kind, min(v, w), max(v, w))
                if key not in edge_index:
                    e = len(self.edges)
                    edge_index[key] = e
                    self.edges.append((v, None if v == w else w))
                    self.kinds.append(kind)
                    self.inc[v].append(e)
                    if w != v:
                        self.inc[w].append(e)
        assert all(len(es) == 3 for es in self.inc)
        self.a_at = [next(e for e in es if self.kinds[e] == 0) for es in self.inc]
        self.x_at = [edge_index[1, min(v, ns[1]), max(v, ns[1])]
                     for v, ns in enumerate(self.neighbours)]
        unseen, self.cycles = set(range(len(self.vertices))), []
        while unseen:
            start, cycle = min(unseen), []
            current = start
            while current not in cycle:
                cycle.append(current)
                unseen.remove(current)
                current = self.neighbours[current][1]
            assert current == start and len(cycle) % 2 == 1
            self.cycles.append(cycle)
        self.cycle_of = {v: c for c, cycle in enumerate(self.cycles) for v in cycle}
        assert len(self.cycles) == self.N // 2

    def supports(self):
        """Exact covers of x-cycle IDs by quotient a-edges/semi-edges."""
        rows = []
        for e, ((u, v), kind) in enumerate(zip(self.edges, self.kinds)):
            if kind:
                continue
            covered = 1 << self.cycle_of[u]
            if v is not None:
                if self.cycle_of[v] == self.cycle_of[u]:
                    continue
                covered |= 1 << self.cycle_of[v]
            rows.append((covered, e))
        incident = [[row for row in rows if row[0] & (1 << c)]
                    for c in range(len(self.cycles))]

        def visit(todo, selected):
            if not todo:
                yield selected
                return
            choices = min(([row for row in incident[c] if row[0] & todo == row[0]]
                           for c in range(len(self.cycles)) if todo & (1 << c)),
                          key=len)
            for covered, edge in choices:
                yield from visit(todo ^ covered, selected | (1 << edge))
        yield from visit((1 << len(self.cycles)) - 1, 0)

    def matching(self, support):
        matching = {e for e in range(len(self.edges)) if support & (1 << e)}
        for cycle in self.cycles:
            positions = [i for i, v in enumerate(cycle) if self.a_at[v] in matching]
            assert len(positions) == 1
            p = positions[0]
            matching.update(self.x_at[cycle[(p + d) % len(cycle)]]
                            for d in range(1, len(cycle), 2))
        assert all(len(matching.intersection(es)) == 1 for es in self.inc)
        return matching

    def canonical_support(self):
        """The circuit-free support used in the proof, for m,n > 1."""
        assert self.m >= 3
        selected_pairs = [(0, 1)]
        selected_pairs += [tuple(sorted((1, (1 + d) % self.m)))
                           for d in range(2, (self.m + 1) // 2)]
        selected_pairs += [tuple(sorted((self.m + 1,
                                         self.m + (1 + d) % self.n)))
                           for d in range(1, (self.n + 1) // 2)]
        support = 0
        for pair in selected_pairs:
            support |= 1 << self.a_at[self.index[pair]]
        # One ordinary edge covers the cross orbit and the first internal
        # orbit; every remaining selected a-dart is a semi-edge.
        assert support.bit_count() == len(self.cycles) - 1
        return support

    def complement_profile(self, matching):
        remaining = [set(es) - matching for es in self.inc]
        unseen, profile = set(range(len(self.vertices))), []
        while unseen:
            start = min(unseen)
            component, todo, component_edges = set(), [start], set()
            while todo:
                v = todo.pop()
                if v in component:
                    continue
                component.add(v)
                for e in remaining[v]:
                    component_edges.add(e)
                    todo.extend(w for w in self.edges[e] if w is not None)
            unseen -= component
            semis = [e for e in component_edges if self.edges[e][1] is None]
            if semis:
                assert len(semis) == 2
                profile.append(('path', len(component), 0))
                continue
            current, previous, word = start, None, self.identity
            while True:
                e = next(e for e in remaining[current] if e != previous)
                u, v = self.edges[e]
                target = v if u == current else u
                generator = self.a if self.kinds[e] == 0 else (
                    self.x if self.neighbours[current][1] == target else self.xinv)
                word = mul(word, generator)
                current, previous = target, e
                if current == start:
                    break
            profile.append(('cycle', len(component), order(word)))
        return tuple(sorted(profile))


def check(m, n):
    quotient = TwoCycleQuotient(m, n)
    profiles, bad, seen = Counter(), [], set()
    for support in quotient.supports():
        assert support not in seen
        seen.add(support)
        profile = quotient.complement_profile(quotient.matching(support))
        profiles[tuple(item for item in profile if item[0] == 'cycle')] += 1
        if any(length % 2 and monodromy % 2 for kind, length, monodromy in profile
               if kind == 'cycle'):
            bad.append((support, profile))
    print('TWO-CYCLE-MONODROMY', (m, n), 'x_order', order(quotient.x),
          'quotient_vertices', len(quotient.vertices), 'x_orbits', len(quotient.cycles),
          'matchings', len(seen), 'bad', len(bad), 'profiles', dict(profiles), flush=True)
    return bad


def check_canonical(limit):
    tested = 0
    for N in range(8, limit + 1, 2):
        for m in range(3, N // 2 + 1, 2):
            n = N - m
            if math.gcd(m, n) != 1:
                continue
            quotient = TwoCycleQuotient(m, n)
            matching = quotient.matching(quotient.canonical_support())
            assert not any(kind == 'cycle' for kind, _, _ in
                           quotient.complement_profile(matching))
            tested += 1
    print('TWO-CYCLE-CANONICAL', 'coprime odd pairs through degree', limit,
          'checked', tested, 'quotient_cycles', 0, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-degree', type=int, default=12)
    parser.add_argument('--canonical-limit', type=int, default=30)
    args = parser.parse_args()
    check_canonical(args.canonical_limit)
    for N in range(4, args.max_degree + 1, 2):
        for m in range(1, N // 2 + 1, 2):
            n = N - m
            if math.gcd(m, n) == 1:
                check(m, n)


if __name__ == '__main__':
    main()
