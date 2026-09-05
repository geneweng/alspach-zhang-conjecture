#!/usr/bin/env python3
"""The one-a-edge matching face, regular-map duals, and repair barriers.

No colouring condition is used to enumerate starts. The default run checks
local odd-cycle lifts, dual face colourings, defect-incidence degrees, complete
small matching spaces, and a fixed PSL(2,11) certificate with positive mean
oddness drift under two natural random moves.

python3 code/quotient_face_experiment.py
python3 code/quotient_face_experiment.py --psl11-centralizers
python3 code/quotient_face_experiment.py --psl11-unrestricted 512
"""

import argparse
import json
import random
from collections import Counter
from fractions import Fraction
from itertools import combinations
from types import SimpleNamespace

from pysat.solvers import Cadical153

from cayley_snark_check import closure, inv, mul, order, psl2
from colour_support_experiment import matching_for_support
from parity_interlacement import defect_degrees, parity_profile
from transvection_switch_experiment import partial_sums
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


def dual_face_colouring(graph, generators, colours, avoid_all_four=False):
    """Colour faces of the regular quotient map and verify the induced flow."""
    a, x = generators
    subgroup = closure([mul(a, x)])
    unseen, face_of, faces = set(graph.group), {}, []
    while unseen:
        g = min(unseen)
        face = {mul(g, h) for h in subgroup}
        identifier = len(faces)
        faces.append(face)
        for h in face:
            face_of[h] = identifier
        unseen -= face
    dual_edges = {tuple(sorted((face_of[g], face_of[mul(g, a)])))
                  for g in graph.group}
    assert all(u != v for u, v in dual_edges)
    clauses = []
    for v in range(len(faces)):
        variables = [colours * v + c + 1 for c in range(colours)]
        clauses.append(variables)
        clauses.extend([-p, -q] for p, q in combinations(variables, 2))
    for u, v in dual_edges:
        clauses.extend([-(colours * u + c + 1), -(colours * v + c + 1)]
                       for c in range(colours))
    if avoid_all_four:
        assert colours == 4
        base = colours * len(faces)
        for vertex, cycle in enumerate(graph.cycles):
            omitted = [base + 4 * vertex + c + 1 for c in range(4)]
            clauses.append(omitted)
            for c, selector in enumerate(omitted):
                clauses.extend([-selector, -(4 * face_of[graph.group[v]] + c + 1)]
                               for v in cycle)
    clauses.append([1])  # break the global colour-permutation symmetry
    with Cadical153(bootstrap_with=clauses) as solver:
        if not solver.solve():
            return None, len(faces), len(dual_edges)
        positive = set(v for v in solver.get_model() if v > 0)
    face_colours = [next(c for c in range(colours)
                         if colours * v + c + 1 in positive)
                    for v in range(len(faces))]
    if not (colours == 3 or (colours == 4 and avoid_all_four)):
        return face_colours, len(faces), len(dual_edges)
    # Embed three colours in F_2^2.  Translating every face colour does not
    # change edge differences; the same verification applies to a locally
    # constrained four-colouring.
    values = []
    for edge in graph.a_edges:
        u = graph.edges[edge][0]
        g = graph.group[u]
        value = face_colours[face_of[g]] ^ face_colours[face_of[mul(g, a)]]
        assert value in (1, 2, 3)
        values.append(value)
    colouring = [0] * len(graph.edges)
    for edge, value in zip(graph.a_edges, values):
        colouring[edge] = value
    for order_edges, x_edges in zip(graph.orders, graph.cycle_edges):
        word = [values[e] for e in order_edges]
        total = 0
        for value in word:
            total ^= value
        assert total == 0
        missing = set(range(4)) - partial_sums(word)
        assert len(missing) == 1
        current = missing.pop()
        for value, edge in zip(word, x_edges):
            current ^= value
            colouring[edge] = current
    assert all(set(colouring[e] for e in incident) == {1, 2, 3}
               for incident in graph.inc)
    return face_colours, len(faces), len(dual_edges)


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
    dual_three_expected = {'A5': True, 'A5_alt': False, 'W50': True, 'F80': False}
    for name, expected in [('A5', 125), ('A5_alt', 125), ('W50', 120), ('F80', 705)]:
        generators, n = cases()[name]
        face = QuotientFace(RepairGraph(generators, n))
        dual_three, dual_vertices, dual_edges = dual_face_colouring(
            face.graph, generators, 3)
        dual_local_four, vertices_again, edges_again = dual_face_colouring(
            face.graph, generators, 4, True)
        assert (dual_three is not None) == dual_three_expected[name]
        assert dual_local_four is not None
        assert (dual_vertices, dual_edges) == (vertices_again, edges_again)
        print('REGULAR-MAP-DUAL', name, 'vertices', dual_vertices, 'edges', dual_edges,
              'three_colourable', dual_three is not None,
              'local_four_colourable', True, flush=True)
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
            degrees = defect_degrees(face, face.project(matching))
            assert len(degrees) == len(lengths)
            assert sum(degree % 2 for degree in degrees) == odd
            factor_counts[len(lengths), odd] += 1
        minimum = min(circuits for circuits, _ in factor_counts)
        minimum_oddness = {odd: number for (circuits, odd), number in factor_counts.items()
                           if circuits == minimum}
        counts, paths = check_descent(face, independent)
        print('FACE-EXHAUSTIVE', name, dict(counts), 'paths', dict(sorted(paths.items())),
              'independent_set_agreement', True, flush=True)
        print('PARITY-INTERLACEMENT', name, 'checked', len(independent),
              'factor_counts', dict(sorted(factor_counts.items())),
              'defect_degree_formula', True, flush=True)
        print('LEXICOGRAPHIC-MINIMUM', name, 'circuits', minimum,
              'oddness_counts', dict(sorted(minimum_oddness.items())), flush=True)
        if name == 'A5':
            # A stronger shortcut would colour two disjoint quotient perfect
            # matchings differently and the other three incident edges alike.
            # It requires their two selected darts to be consecutive at every
            # quotient vertex.  Exact enumeration shows that this need not exist.
            supports = [face.project(matching) for matching in independent]
            consecutive_pairs = []
            for first in supports:
                for second in supports:
                    if first & second:
                        continue
                    positions = [[i for i, e in enumerate(order)
                                  if (first | second) & (1 << e)]
                                 for order in face.graph.orders]
                    if all(len(pair) == 2 and (pair[0] - pair[1]) % 5 in (1, 4)
                           for pair in positions):
                        consecutive_pairs.append((first, second))
            assert not consecutive_pairs
            print('CONSECUTIVE-QUOTIENT-MATCHINGS', name, 'ordered_pairs', 0,
                  'checked', len(supports) ** 2, flush=True)
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
            dual_four, dual_vertices, dual_edges = dual_face_colouring(
                face.graph, [a, y], 4)
            expected_four = order(mul(a, y)) != 11
            assert (dual_four is not None) == expected_four
            dual_five = None
            if dual_four is None:
                dual_five, vertices_again, edges_again = dual_face_colouring(
                    face.graph, [a, y], 5)
                assert dual_five is not None
                assert (dual_vertices, dual_edges) == (vertices_again, edges_again)
            subgroup = [h for h in group if mul(h, a) == mul(a, h)]
            status = {}
            sat = list(sample_matchings(face.graph, 513, random.Random(509), True,
                                        subgroup, status))
            independent = list(independent_matchings(face, subgroup))
            assert status['exhausted']
            assert len(independent) == len(set(independent)) == len(sat)
            assert set(independent) == set(sat)
            profiles = [(len(face.graph.factor_lengths(matching)),
                         sum(n % 2 for n in face.graph.factor_lengths(matching)))
                        for matching in sat]
            minimum = min(circuits for circuits, _ in profiles)
            minimum_oddness = Counter(odd for circuits, odd in profiles
                                      if circuits == minimum)
            assert set(minimum_oddness) == {0}
            counts, paths = check_descent(face, sat)
            print('CENTRALIZER-FACE', i, power, 'a', a, 'x', y, dict(counts),
                  'paths', dict(sorted(paths.items())), 'minimum_circuits', minimum,
                  'oddness_at_minimum', dict(minimum_oddness),
                  'dual_vertices', dual_vertices,
                  'dual_chromatic_barrier', 5 if dual_five is not None else 'at_most_4',
                  'independent_set_agreement', True,
                  flush=True)
            totals.update(counts)
            total_paths.update(paths)
    assert totals == Counter(starts=4336, already_even=3906, obstructed=430)
    assert max(total_paths) == 3
    print('CENTRALIZER-FACE-SUMMARY', dict(totals), 'paths',
          dict(sorted(total_paths.items())), flush=True)


def check_psl11_unrestricted(samples):
    """Sample arbitrary quotient matchings and test one-step strict descent."""
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
    totals = Counter()
    for i, a in enumerate(representatives):
        for power, y in [(1, x), (2, mul(x, x))]:
            face = QuotientFace(RepairGraph([a, y], 660))
            starts = list(sample_matchings(face.graph, samples, random.Random(1847), True))
            counts = Counter(starts=len(starts))
            for matching in starts:
                odd = sum(n % 2 for n in face.graph.factor_lengths(matching))
                if not odd:
                    counts['even'] += 1
                    continue
                counts['obstructed'] += 1
                witness, _ = face.graph.audit(matching)
                assert witness is not None and witness['new_oddness'] < odd
                counts['decreases'] += 1
            assert counts['obstructed'] == counts['decreases']
            totals.update(counts)
            print('UNRESTRICTED-FACE', i, power, dict(counts), flush=True)
    assert totals['starts'] == 12 * samples
    if samples == 512:
        assert totals['obstructed'] == totals['decreases'] == 3563
        assert totals['even'] == 2581
    print('UNRESTRICTED-FACE-SUMMARY', dict(totals), 'no_stall', True, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--psl11-centralizers', action='store_true')
    parser.add_argument('--psl11-unrestricted', type=int, metavar='SAMPLES', default=0)
    args = parser.parse_args()
    check_local_lifts()
    check_small()
    check_averaging_barrier()
    if args.psl11_centralizers:
        check_psl11_centralizers()
    if args.psl11_unrestricted:
        check_psl11_unrestricted(args.psl11_unrestricted)


if __name__ == '__main__':
    main()
