#!/usr/bin/env python3
"""The one-a-edge matching face: exact exchanges and an averaging barrier.

No colouring condition is used to enumerate starts. The default run checks
local odd-cycle lifts, complete small matching spaces, and a fixed PSL(2,11)
certificate with positive mean oddness drift under two natural random moves.

python3 code/quotient_face_experiment.py
python3 code/quotient_face_experiment.py --psl11-centralizers
"""

import argparse
import json
import random
from collections import Counter
from fractions import Fraction
from types import SimpleNamespace

from cayley_snark_check import closure, inv, mul, order, psl2
from colour_support_experiment import matching_for_support
from parity_interlacement import parity_profile
from translated_block_repair import GENERATORS, independent_factor_lengths
from translated_repair_audit import RepairGraph, bits, cases, sample_matchings


DRIFT_SUPPORT = [
    1, 7, 14, 19, 30, 32, 37, 38, 45, 48, 49, 57, 60, 61, 64, 66, 72,
    73, 91, 94, 97, 108, 127, 131, 135, 138, 144, 147, 150, 151, 152,
    154, 156, 161, 163, 167, 177, 181, 188, 196, 198, 202, 205, 210,
    223, 226, 240, 258, 259, 261, 270, 276, 282, 284, 285, 288, 291,
    292, 295, 296, 306, 310, 314, 319, 320, 321,
]


def components(edges, mask):
    """Edge components by vertex-set traversal, preserving parallel edge IDs."""
    incident = {}
    for e in bits(mask):
        for v in edges[e]:
            incident.setdefault(v, set()).add(e)
    assert all(len(es) == 2 for es in incident.values())
    unseen, answer = set(incident), []
    while unseen:
        todo, vertices, chosen = [min(unseen)], set(), set()
        while todo:
            v = todo.pop()
            if v in vertices:
                continue
            vertices.add(v)
            for e in incident[v]:
                chosen.add(e)
                todo.extend(edges[e])
        unseen -= vertices
        answer.append(sum(1 << e for e in chosen))
    return sorted(answer)


class QuotientFace:
    def __init__(self, graph):
        self.graph = graph
        self.m = len(graph.a_edges)
        self.masks = [sum(1 << e for e in cyclic) for cyclic in graph.orders]
        ends = [[] for _ in graph.a_edges]
        for v, cyclic in enumerate(graph.orders):
            for e in cyclic:
                ends[e].append(v)
        assert all(len(vs) == 2 and vs[0] != vs[1] for vs in ends)
        self.edges = [tuple(vs) for vs in ends]
        positions = {e: i for i, e in enumerate(graph.a_edges)}
        self.actions = [[positions[action[e]] for e in graph.a_edges]
                        for action in graph.actions]

    def project(self, matching):
        return sum(1 << i for i, e in enumerate(self.graph.a_edges)
                   if matching & (1 << e))

    def is_matching(self, support):
        return all((support & mask).bit_count() == 1 for mask in self.masks)

    def lift(self, support):
        """Linear lift, independently checked against the forced-path routine."""
        assert self.is_matching(support)
        graph = self.graph
        matching = sum(1 << graph.a_edges[i] for i in bits(support))
        for cyclic, x_edges in zip(graph.orders, graph.cycle_edges):
            for i, e in enumerate(x_edges):
                value = sum(bool(support & (1 << cyclic[(i - d) % 5]))
                            for d in (1, 3))
                assert value in (0, 1)
                matching |= value << e
        assert set(bits(matching)) == matching_for_support(graph, set(bits(support)))
        return matching

    def check_pair(self, first, second):
        """Verify component-by-component quotient/lift correspondence."""
        p, q = self.project(first), self.project(second)
        assert self.is_matching(p) and self.is_matching(q)
        downstairs = components(self.edges, p ^ q)
        upstairs = components(self.graph.edges, first ^ second)
        assert sorted(self.project(d) for d in upstairs) == downstairs
        for d in upstairs:
            assert self.lift(p ^ self.project(d)) == first ^ d
        return len(upstairs)


def independent_matchings(face, subgroup=()):
    """Exact cover of quotient vertices by compatible subgroup edge orbits.

    This is independent of the SAT enumeration. With no subgroup each row is
    one quotient edge. A subgroup-invariant 1-factor is a union of edge orbits;
    an orbit meeting any vertex twice cannot occur and is discarded.
    """
    indices = [face.graph.index[h] for h in subgroup]
    unseen, rows = set(range(face.m)), []
    while unseen:
        e = min(unseen)
        orbit = {face.actions[i][e] for i in indices} if indices else {e}
        assert e in orbit and orbit <= unseen
        unseen -= orbit
        vertices = [v for f in orbit for v in face.edges[f]]
        if len(vertices) == len(set(vertices)):
            rows.append((sum(1 << v for v in vertices), sum(1 << f for f in orbit)))
    incident = [[row for row in rows if row[0] & (1 << v)]
                for v in range(face.graph.q)]

    def visit(todo, selected):
        if not todo:
            assert face.is_matching(selected)
            yield face.lift(selected)
            return
        best = None
        for v in bits(todo):
            choices = [row for row in incident[v] if row[0] & todo == row[0]]
            if not choices:
                return
            if best is None or len(choices) < len(best):
                best = choices
        for covered, edges in best:
            yield from visit(todo ^ covered, selected | edges)

    yield from visit((1 << face.graph.q) - 1, 0)


def check_local_lifts():
    tested = 0
    for s in (3, 5, 7, 9, 11):
        states = [{(p + d) % s for d in range(1, s, 2)} for p in range(s)]
        for p, first in enumerate(states):
            for q, second in enumerate(states):
                # Edge i means {i,i+1}. The differing internal edges form
                # the unique even path joining the two unmatched vertices.
                difference = first ^ second
                degrees = Counter(v for e in difference for v in (e, (e + 1) % s))
                assert difference == {i for i in range(s) if (
                    sum(p == (i - d) % s for d in range(1, s, 2)) !=
                    sum(q == (i - d) % s for d in range(1, s, 2)))}
                assert len(difference) % 2 == 0
                assert {v for v, degree in degrees.items() if degree == 1} == (
                    {p, q} if p != q else set())
                assert all(degree in (1, 2) for degree in degrees.values())
                tested += 1
    assert tested == 285
    print('LOCAL-ODD-CYCLE-LIFTS', tested, flush=True)


def check_descent(face, starts):
    graph, counts, paths = face.graph, Counter(), Counter()
    for matching in starts:
        assert face.is_matching(face.project(matching))
        assert face.lift(face.project(matching)) == matching
        odd = sum(n % 2 for n in graph.factor_lengths(matching))
        counts['starts'] += 1
        counts['obstructed' if odd else 'already_even'] += 1
        steps = 0
        while odd:
            witness, outcomes = graph.audit(matching)
            assert not outcomes.get('invalid', 0)
            if witness is None:
                raise AssertionError(('restricted descent stalls', hex(matching), outcomes))
            g = graph.index[tuple(witness['translation'])]
            translated = sum(1 << graph.actions[g][e] for e in bits(matching))
            face.check_pair(matching, translated)
            matching = int(witness['matching'], 16)
            assert face.is_matching(face.project(matching))
            odd = witness['new_oddness']
            steps += 1
        graph.check(matching)
        paths[steps] += 1
    return counts, paths


def check_small():
    pair_count = 0
    for name, expected in [('A5', 125), ('A5_alt', 125), ('W50', 120), ('F80', 705)]:
        generators, n = cases()[name]
        face = QuotientFace(RepairGraph(generators, n))
        independent = list(independent_matchings(face))
        status = {}
        sat = list(sample_matchings(face.graph, expected + 1, random.Random(509),
                                    True, status=status))
        assert status['exhausted'] and len(sat) == expected
        assert len(independent) == len(set(independent)) == expected
        assert set(independent) == set(sat)
        if name in ('A5', 'A5_alt'):
            for i, first in enumerate(independent):
                for second in independent[i:]:
                    face.check_pair(first, second)
                    pair_count += 1
        factor_counts = Counter()
        for matching in independent:
            lengths = face.graph.factor_lengths(matching)
            odd = sum(n % 2 for n in lengths)
            profile = parity_profile(face, face.project(matching))
            assert profile['odd_circuits'] == odd
            assert profile['nullity'] + profile['components'] == len(lengths)
            factor_counts[len(lengths), odd] += 1
        counts, paths = check_descent(face, independent)
        print('FACE-EXHAUSTIVE', name, dict(counts), 'paths', dict(sorted(paths.items())),
              'independent_set_agreement', True, flush=True)
        print('PARITY-INTERLACEMENT', name, 'checked', len(independent),
              'factor_counts', dict(sorted(factor_counts.items())), flush=True)
    assert pair_count == 15750
    print('PAIRWISE-EXCHANGE-CHECK', pair_count, flush=True)
    # The pentagon quotient of the Petersen graph has two vertices and five
    # parallel spokes, in these two cyclic orders. All five starts must fail.
    petersen = SimpleNamespace(
        m=5, edges=[(0, 1)] * 5, is_matching=lambda p: p.bit_count() == 1,
        graph=SimpleNamespace(q=2, orders=[list(range(5)), [0, 2, 4, 1, 3]]))
    for i in range(5):
        profile = parity_profile(petersen, 1 << i)
        assert profile['odd_circuits'] == 2
        assert profile['nullity'] + profile['components'] == 2
    print('PARITY-INTERLACEMENT-PETERSEN', 'five_starts', 'two_odd_circuits_each', flush=True)


def check_averaging_barrier():
    graph = RepairGraph(GENERATORS, 660)
    assert set(graph.group) == psl2(11)
    face = QuotientFace(graph)
    matching = face.lift(sum(1 << e for e in DRIFT_SUPPORT))
    graph.check(matching)
    assert independent_factor_lengths(graph, matching) == [15, 27, 618]
    assert graph.matching_stabilizer_order(matching) == 1
    profile = parity_profile(face, face.project(matching))
    assert profile['odd_circuits'] == 2
    assert profile['nullity'] + profile['components'] == 3
    counts, single, block = Counter(), Fraction(0), Fraction(0)
    max_components, subset_cases = 0, 0
    for action in graph.actions:
        translated = sum(1 << action[e] for e in bits(matching))
        face.check_pair(matching, translated)
        circuits = tuple(graph.difference_circuits(matching, action))
        # Independently constructed vertex-set components must agree.
        assert sorted(circuits) == components(graph.edges, matching ^ translated)
        r = len(circuits)
        max_components = max(max_components, r)
        assert r <= 7
        subtotal = 0
        for circuit in circuits:
            new = matching ^ circuit
            odd = sum(n % 2 for n in independent_factor_lengths(graph, new))
            counts[odd] += 1
            subtotal += odd - 2
        if r:
            single += Fraction(subtotal, r)
        current, subtotal = matching, 0
        for k in range(1 << r):
            # Gray-code enumeration visits every component subset once.
            if k:
                current ^= circuits[(k & -k).bit_length() - 1]
            assert face.is_matching(face.project(current))
            odd = sum(n % 2 for n in independent_factor_lengths(graph, current))
            assert odd == sum(n % 2 for n in graph.factor_lengths(current))
            subtotal += odd - 2
            subset_cases += 1
        block += Fraction(subtotal, 1 << r)
    assert counts == Counter({0: 215, 2: 1021, 4: 272, 6: 10})
    assert single / 660 == Fraction(197, 2475)
    assert block / 660 == Fraction(677, 21120)
    assert max_components == 7
    assert subset_cases == 4349
    witness, _ = graph.audit(matching)
    assert witness['new_oddness'] == 0
    repaired = int(witness['matching'], 16)
    face.lift(face.project(repaired))
    repaired_profile = parity_profile(face, face.project(repaired))
    assert repaired_profile['odd_circuits'] == 0
    assert repaired_profile['nullity'] + repaired_profile['components'] == 5
    print('POSITIVE-DRIFT-CERTIFICATE', json.dumps({
        'fingerprint': graph.fingerprint, 'support': DRIFT_SUPPORT,
        'initial_lengths': [15, 27, 618], 'single_flip_counts': dict(counts),
        'mean_delta_uniform_translate_then_circuit': str(single / 660),
        'mean_delta_uniform_translate_then_subset': str(block / 660),
        'subset_cases': subset_cases, 'max_components': max_components,
        'parity_profile_before': profile, 'parity_profile_after': repaired_profile,
        'colouring_witness': witness}), flush=True)


def check_psl11_centralizers():
    group = psl2(11)
    x = min(h for h in group if order(h) == 5)
    cx = [h for h in group if mul(h, x) == mul(x, h)]
    unseen, representatives = {h for h in group if order(h) == 2}, []
    while unseen:
        a = min(unseen)
        unseen -= {mul(mul(inv(h), a), h) for h in cx}
        if len(closure([a, x])) == 660:
            representatives.append(a)
    assert len(representatives) == 6
    totals, total_paths = Counter(), Counter()
    for i, a in enumerate(representatives):
        for power, y in [(1, x), (2, mul(x, x))]:
            face = QuotientFace(RepairGraph([a, y], 660))
            subgroup = [h for h in group if mul(h, a) == mul(a, h)]
            status = {}
            sat = list(sample_matchings(face.graph, 513, random.Random(509), True,
                                        subgroup, status))
            independent = list(independent_matchings(face, subgroup))
            assert status['exhausted']
            assert len(independent) == len(set(independent)) == len(sat)
            assert set(independent) == set(sat)
            counts, paths = check_descent(face, sat)
            print('CENTRALIZER-FACE', i, power, 'a', a, 'x', y, dict(counts),
                  'paths', dict(sorted(paths.items())), 'independent_set_agreement', True,
                  flush=True)
            totals.update(counts)
            total_paths.update(paths)
    assert totals == Counter(starts=4336, already_even=3906, obstructed=430)
    assert max(total_paths) == 3
    print('CENTRALIZER-FACE-SUMMARY', dict(totals), 'paths',
          dict(sorted(total_paths.items())), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--psl11-centralizers', action='store_true')
    args = parser.parse_args()
    check_local_lifts()
    check_small()
    check_averaging_barrier()
    if args.psl11_centralizers:
        check_psl11_centralizers()


if __name__ == '__main__':
    main()
