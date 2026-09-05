#!/usr/bin/env python3
"""Minimum-face matchings in the two-subset quotient of S_(s+1).

a=(0,infinity), x=(0,...,s-1), s odd. Enumerate quotient matchings
without a colouring constraint and check the monodromy of every closed
complementary component. Capped paths always lift to even circuits.
The quotient has only s(s+1)/2 vertices; S_(s+1) is not enumerated.
"""

import argparse
from collections import Counter
from itertools import combinations, product

from cayley_snark_check import inv, mul, order
from cdc_palette_experiment import cycle_perm


class PairQuotient:
    def __init__(self, s, a=None):
        assert s >= 3 and s % 2 == 1
        self.s = s
        self.a = cycle_perm(s + 1, [0, s]) if a is None else tuple(a)
        assert sorted(self.a) == list(range(s + 1)) and order(self.a) == 2
        self.x = cycle_perm(s + 1, list(range(s)))
        self.identity = tuple(range(s + 1))
        self.vertices = list(combinations(range(s + 1), 2))
        self.index = {v: i for i, v in enumerate(self.vertices)}
        self.n = len(self.vertices)
        self.neighbours = [[self.index[tuple(sorted(p[i] for i in v))]
                            for p in (self.a, self.x, inv(self.x))]
                           for v in self.vertices]
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
        self.cycles = [[self.index[tuple(sorted((s, j)))] for j in range(s)]]
        self.cycles += [[self.index[tuple(sorted((j, (j + d) % s)))]
                         for j in range(s)] for d in range(1, (s + 1) // 2)]
        assert sorted(v for c in self.cycles for v in c) == list(range(self.n))
        self.cycle_of = {v: c for c, cyc in enumerate(self.cycles) for v in cyc}

    def matchings(self):
        """Choose one a-dart per x-cycle, then its forced internal matching."""
        assert self.a == cycle_perm(self.s + 1, [0, self.s])
        for central in range(self.s):
            central_v = self.cycles[0][central]
            selected_a = self.a_at[central_v]
            opposite = self.neighbours[central_v][0]
            options = []
            for c, cycle in enumerate(self.cycles[1:], 1):
                if central and self.cycle_of[opposite] == c:
                    options.append([cycle.index(opposite)])
                else:
                    options.append([i for i, v in enumerate(cycle)
                                    if self.edges[self.a_at[v]][1] is None])
            for leaves in product(*options):
                positions = (central,) + leaves
                matching = {selected_a}
                for cycle, p in zip(self.cycles, positions):
                    matching.add(self.a_at[cycle[p]])
                    matching.update(self.x_at[cycle[(p + d) % self.s]]
                                    for d in range(1, self.s, 2))
                assert all(len(matching.intersection(es)) == 1 for es in self.inc)
                yield positions, matching

    def exact_cover_supports(self):
        """Independent enumeration: cover cycle IDs by a-edges and semi-edges."""
        rows = []
        for e, ((u, v), kind) in enumerate(zip(self.edges, self.kinds)):
            if kind == 0:
                covered = 1 << self.cycle_of[u]
                if v is not None:
                    if self.cycle_of[v] == self.cycle_of[u]:
                        continue  # two selected positions on one cycle are forbidden
                    covered |= 1 << self.cycle_of[v]
                rows.append((covered, e))

        def visit(todo, support):
            if not todo:
                yield support
                return
            first = todo & -todo
            for covered, e in rows:
                if covered & first and covered & todo == covered:
                    yield from visit(todo ^ covered, support | (1 << e))
        yield from visit((1 << len(self.cycles)) - 1, 0)

    def lift_support(self, support):
        """Force the internal matching from an exact-cover a-edge support."""
        matching = {e for e in range(len(self.edges)) if support & (1 << e)}
        assert all(self.kinds[e] == 0 for e in matching)
        for cycle in self.cycles:
            positions = [p for p, v in enumerate(cycle) if self.a_at[v] in matching]
            assert len(positions) == 1
            p = positions[0]
            matching.update(self.x_at[cycle[(p + d) % self.s]]
                            for d in range(1, self.s, 2))
        assert all(len(matching.intersection(es)) == 1 for es in self.inc)
        return matching

    def complement(self, matching):
        """Traverse quotient components and multiply each closed walk's word."""
        remaining = [set(es) - matching for es in self.inc]
        assert all(len(es) == 2 for es in remaining)
        unseen, result = set(range(self.n)), []
        while unseen:
            todo, vertices, edges = [min(unseen)], set(), set()
            while todo:
                v = todo.pop()
                if v in vertices:
                    continue
                vertices.add(v)
                for e in remaining[v]:
                    edges.add(e)
                    todo.extend(w for w in self.edges[e] if w is not None)
            unseen -= vertices
            semis = [e for e in edges if self.edges[e][1] is None]
            if semis:
                assert len(semis) == 2
                result.append(('path', len(vertices), None))
                continue
            start = current = min(vertices)
            previous, word = None, self.identity
            while True:
                e = next(e for e in sorted(remaining[current]) if e != previous)
                u, v = self.edges[e]
                target = v if u == current else u
                if self.kinds[e] == 0:
                    generator = self.a
                else:
                    generator = (self.x if self.neighbours[current][1] == target
                                 else inv(self.x))
                word = mul(word, generator)
                current, previous = target, e
                if current == start:
                    break
            assert tuple(sorted(word[i] for i in self.vertices[start])) == (
                self.vertices[start])
            result.append(('cycle', len(vertices), order(word)))
        return result


def check(s):
    quotient = PairQuotient(s)
    counts, profiles, seen, supports = Counter(), Counter(), set(), set()
    for positions, matching in quotient.matchings():
        frozen = sum(1 << e for e in matching)
        assert frozen not in seen
        seen.add(frozen)
        supports.add(sum(1 << e for e in matching if quotient.kinds[e] == 0))
        counts['matchings'] += 1
        profile = []
        for kind, length, monodromy in quotient.complement(matching):
            if kind == 'cycle':
                assert length * monodromy % 2 == 0, (s, positions, length, monodromy)
                profile.append((length, monodromy))
                counts['closed_components'] += 1
                counts['odd_closed_components'] += length % 2
                if s >= 7:
                    assert positions[0] == 0
                    assert (length, monodromy) == (5, s - 1)
        if s >= 7:
            assert len(profile) <= 1
        if positions[0] != 0:
            assert not profile
        profiles[tuple(sorted(profile))] += 1
    expected = (s - 2) ** ((s - 1) // 2) + (
        (s - 1) * (s - 2) ** ((s - 3) // 2))
    assert counts['matchings'] == expected
    independent = list(quotient.exact_cover_supports())
    assert len(independent) == len(set(independent)) == expected
    assert supports == set(independent)
    if s == 3:
        assert profiles == {(): 2, ((6, 2),): 1}
    if s == 5:
        assert profiles == {(): 16, ((5, 4),): 4, ((5, 4), (5, 4)): 1}
    if s >= 7:
        expected_closed = (s - 3) // 2 * (s - 2) ** ((s - 3) // 2)
        assert profiles == {(): expected - expected_closed,
                            ((5, s - 1),): expected_closed}
    print('TRANSPOSITION-QUOTIENT', s, dict(counts), 'cycle_profiles',
          dict(sorted(profiles.items())), 'independent_exact_cover_agreement', True,
          flush=True)


def check_words():
    for s in range(5, 102, 2):
        a = cycle_perm(s + 1, [0, s])
        x = cycle_perm(s + 1, list(range(s)))
        word = tuple(range(s + 1))
        for generator in (x, x, a, inv(x), a):
            word = mul(word, generator)
        expected = mul(cycle_perm(s + 1, list(range(s - 1))),
                       cycle_perm(s + 1, [s - 1, s]))
        assert word == expected and order(word) == s - 1
    q = PairQuotient(5)
    word = q.identity
    for generator in (q.x, q.a, q.x, q.x, q.a):
        word = mul(word, generator)
    assert word == mul(cycle_perm(6, [0, 3, 1, 4]), cycle_perm(6, [2, 5]))
    assert order(word) == 4
    print('TRANSPOSITION-WORDS', 'all odd s from 5 through 101', 'verified', flush=True)


def check_upstairs_six():
    """Independent exact-cover enumeration and direct traversal on all 720 vertices."""
    from quotient_face_experiment import QuotientFace, independent_matchings
    from translated_repair_audit import RepairGraph

    quotient = PairQuotient(5)
    graph = RepairGraph([quotient.a, quotient.x], 720)
    face = QuotientFace(graph)
    subgroup = [h for h in graph.group if mul(h, quotient.a) == mul(quotient.a, h)]
    assert len(subgroup) == 48
    pair_at = [quotient.index[tuple(sorted((g[0], g[5])))] for g in graph.group]
    item_index = {(kind, u if v is None else min(u, v),
                   u if v is None else max(u, v)): e
                  for e, ((u, v), kind) in enumerate(zip(quotient.edges, quotient.kinds))}
    edge_image = [item_index[0 if kind == 'a' else 1,
                            min(pair_at[u], pair_at[v]), max(pair_at[u], pair_at[v])]
                  for (u, v), kind in zip(graph.edges, graph.kinds)]
    lifted, counts = set(), Counter()
    for _, matching in quotient.matchings():
        upstairs = sum(1 << e for e, image in enumerate(edge_image) if image in matching)
        lifted.add(upstairs)
        graph.check(upstairs)
        lengths = graph.factor_lengths(upstairs)
        assert all(length % 2 == 0 for length in lengths)
        counts[len(lengths)] += 1
    assert lifted == set(independent_matchings(face, subgroup))
    assert counts == {16: 4, 28: 16, 40: 1}
    # H already contains a Sylow 2-subgroup of S6. The previous bound is five.
    assert min(counts) > 5
    print('TRANSPOSITION-UPSTAIRS-S6', len(lifted), dict(sorted(counts.items())),
          'independent_exact_cover_agreement', True, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-s', type=int, default=9)
    args = parser.parse_args()
    if args.max_s < 3:
        parser.error('--max-s must be at least 3')
    check_words()
    check_upstairs_six()
    for s in range(3, args.max_s + 1, 2):
        check(s)


if __name__ == '__main__':
    main()
