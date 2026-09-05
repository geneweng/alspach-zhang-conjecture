#!/usr/bin/env python3
"""Check the two-subset quotient for an odd cycle and a chord transposition.

On Z_N, x is translation by one and a=(0,k), where N is odd and gcd(N,k)=1.
The C(a)-quotient has two-subsets as vertices and one x-orbit for every cyclic
distance.  For every nonadjacent chord, the script constructs the matching
used in the proof and verifies its quotient complement.  Apart from one
harmless even quotient circuit for (N,k)=(5,2), every complement is circuit-free.
"""

import argparse
import math
from itertools import combinations

from cayley_snark_check import closure, inv, mul, order
from cdc_palette_experiment import cycle_perm


class OddTranspositionQuotient:
    def __init__(self, N, k):
        assert N >= 5 and N % 2 == 1 and 1 <= k <= N // 2
        assert math.gcd(N, k) == 1
        self.N, self.k = N, k
        self.a = cycle_perm(N, [0, k])
        self.x = cycle_perm(N, list(range(N)))
        self.xinv = inv(self.x)
        self.identity = tuple(range(N))
        self.vertices = list(combinations(range(N), 2))
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
        self.cycles = [[self.index[tuple(sorted((i, (i + d) % N)))]
                        for i in range(N)] for d in range(1, (N + 1) // 2)]
        assert sorted(v for cycle in self.cycles for v in cycle) == list(
            range(len(self.vertices)))
        self.cycle_of = {v: c for c, cycle in enumerate(self.cycles) for v in cycle}

    def matching_from_positions(self, positions):
        """Select a semi-edge at each listed position and force the x-edges."""
        assert len(positions) == len(self.cycles)
        selected_edges = set()
        for cycle, p in zip(self.cycles, positions):
            selected = self.a_at[cycle[p]]
            assert self.edges[selected][1] is None
            selected_edges.add(selected)
        return self.matching_from_a_edges(selected_edges)

    def matching_from_a_edges(self, selected_edges):
        """Force x-edges after an exact cover by a-edges and semi-edges."""
        matching = set(selected_edges)
        for cycle in self.cycles:
            positions = [p for p, v in enumerate(cycle)
                         if self.a_at[v] in selected_edges]
            assert len(positions) == 1
            p = positions[0]
            matching.update(self.x_at[cycle[(p + d) % self.N]]
                            for d in range(1, self.N, 2))
        assert all(len(matching.intersection(es)) == 1 for es in self.inc)
        return matching

    def ordinary_positions(self, cycle):
        return [p for p, v in enumerate(cycle)
                if self.edges[self.a_at[v]][1] is not None]

    def local_blockers(self, cycle):
        """Semi-edge positions leaving at most one ordinary dart per x-fragment."""
        ordinary = set(self.ordinary_positions(cycle))
        blockers = []
        for p, v in enumerate(cycle):
            if self.edges[self.a_at[v]][1] is not None:
                continue
            fragments = [{(p - 1) % self.N, p, (p + 1) % self.N}]
            fragments.extend({(p + d) % self.N, (p + d + 1) % self.N}
                             for d in range(2, self.N - 1, 2))
            if all(len(fragment & ordinary) <= 1 for fragment in fragments):
                blockers.append(p)
        return blockers

    def canonical_matching(self):
        """Return the proof's matching for every nonadjacent chord k >= 2."""
        h, k = self.N // 2, self.k
        assert k >= 2
        if (self.N, k) == (5, 2):
            return self.matching_from_positions([3, 1])
        if (self.N, k) == (7, 2):
            return self.matching_from_positions([4, 0, 1])
        if (self.N, k) == (7, 3):
            return self.matching_from_positions([1, 2, 0])
        if k == 2:
            assert self.N >= 9
            return self.matching_from_positions([4, 0] + [3] * (h - 2))
        if k == 3:
            assert self.N >= 11
            # B_(1,3)--B_(4,0) is an ordinary a-edge and covers both orbits.
            selected_edges = {self.a_at[self.cycles[0][3]]}
            selected_edges.add(self.a_at[self.cycles[1][4]])
            selected_edges.add(self.a_at[self.cycles[2][0]])
            selected_edges.update(self.a_at[self.cycles[d - 1][1]]
                                  for d in range(5, h + 1))
            return self.matching_from_a_edges(selected_edges)

        positions = []
        for d in range(1, h + 1):
            if d == k:
                p = (k + 1) % 2 if k == h else 0
            elif d == k - 1:
                p = 2
            elif k == 4 and d == 2:
                p = 5
            elif d == k - 2:
                p = 3
            elif d == 1 and k % 2 == 0:
                p = k + 1
            else:
                p = 1
            assert p in self.local_blockers(self.cycles[d - 1])
            positions.append(p)
        return self.matching_from_positions(positions)

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

    def contracted_edges(self):
        result = []
        for e, ((u, v), kind) in enumerate(zip(self.edges, self.kinds)):
            if kind == 0 and v is not None:
                result.append((self.cycle_of[u] + 1, self.cycle_of[v] + 1))
        return sorted(result)


def check(N, k, verbose=False):
    quotient = OddTranspositionQuotient(N, k)
    blockers = [quotient.local_blockers(cycle) for cycle in quotient.cycles]
    failed = [d + 1 for d, choices in enumerate(blockers) if not choices]
    if k == 1:
        if verbose:
            print('ODD-TRANSPOSITION', 'N', N, 'k', k,
                  'adjacent chord handled by the sigma-tau theorem', flush=True)
        return None
    profile = quotient.complement_profile(quotient.canonical_matching())
    circuits = tuple(item for item in profile if item[0] == 'cycle')
    if (N, k) == (5, 2):
        assert circuits == (('cycle', 2, 6),)
    else:
        assert not circuits
    if verbose:
        print('ODD-TRANSPOSITION', 'N', N, 'k', k,
              'locally_unblocked_orbits', failed,
              'complementary_quotient_circuits', circuits, flush=True)
    return circuits


def check_upstairs_small():
    """Independently lift the three exceptional matchings to S5 and S7."""
    results = []
    for N, k in ((5, 2), (7, 2), (7, 3)):
        quotient = OddTranspositionQuotient(N, k)
        matching = quotient.canonical_matching()
        group = sorted(closure([quotient.a, quotient.x]))
        assert len(group) == math.factorial(N)
        index = {g: i for i, g in enumerate(group)}
        labelled = {}
        for u, g in enumerate(group):
            for kind, generator in enumerate((quotient.a, quotient.x, quotient.xinv)):
                v = index[mul(g, generator)]
                labelled[(min(u, v), max(u, v))] = min(kind, 1)
        edges = sorted(labelled)
        kinds = [labelled[edge] for edge in edges]
        inc = [[] for _ in group]
        for e, (u, v) in enumerate(edges):
            inc[u].append(e)
            inc[v].append(e)
        assert all(len(es) == 3 for es in inc)
        pair_at = [quotient.index[tuple(sorted((g[0], g[k])))] for g in group]
        item_index = {(kind, u if v is None else min(u, v),
                       u if v is None else max(u, v)): e
                      for e, ((u, v), kind) in enumerate(
                          zip(quotient.edges, quotient.kinds))}
        edge_image = [item_index[kind,
                                 min(pair_at[u], pair_at[v]),
                                 max(pair_at[u], pair_at[v])]
                      for (u, v), kind in zip(edges, kinds)]
        upstairs = {e for e, image in enumerate(edge_image) if image in matching}
        assert all(len(upstairs.intersection(es)) == 1 for es in inc)
        remaining = [set(es) - upstairs for es in inc]
        unseen, lengths = set(range(len(group))), []
        while unseen:
            start = current = min(unseen)
            previous, length = None, 0
            while True:
                unseen.remove(current)
                e = next(e for e in remaining[current] if e != previous)
                u, v = edges[e]
                current, previous = (v if u == current else u), e
                length += 1
                if current == start:
                    break
            lengths.append(length)
        assert all(length % 2 == 0 for length in lengths)
        results.append((N, k, len(lengths)))
    print('ODD-TRANSPOSITION-UPSTAIRS', results, 'all_cycles_even', True,
          flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-N', type=int, default=51)
    parser.add_argument('--show-contracted', action='store_true')
    parser.add_argument('--verbose', action='store_true')
    parser.add_argument('--skip-upstairs', action='store_true')
    args = parser.parse_args()
    if not args.skip_upstairs:
        check_upstairs_small()
    tested = constructed = 0
    for N in range(5, args.max_N + 1, 2):
        for k in range(1, N // 2 + 1):
            if math.gcd(N, k) != 1:
                continue
            check(N, k, args.verbose)
            tested += 1
            constructed += k >= 2
            if args.show_contracted:
                quotient = OddTranspositionQuotient(N, k)
                print('CONTRACTED', N, k, quotient.contracted_edges(), flush=True)
    print('ODD-TRANSPOSITION-SUMMARY', 'through', args.max_N,
          'coprime_chord_cases', tested,
          'constructed_nonadjacent_cases', constructed,
          'construction_failures', 0, flush=True)


if __name__ == '__main__':
    main()
