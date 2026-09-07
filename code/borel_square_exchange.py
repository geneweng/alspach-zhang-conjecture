"""Certify the residual two-chord theorem and its 13/16 density bound.

For both inverse traces one, a unique folded alternating square exists
exactly when Tr(1/Delta)=0. It colours at least one of c,c+t, although it
need not span or colour both. The construction uses quadratic fixed
points; independent graph searches and projective words audit the result.
"""

import argparse
import json
from collections import Counter
from itertools import product
from math import gcd
from pathlib import Path
from time import perf_counter

from borel_alternating_exchange import alternating_exchanges
from borel_exchange_algebra import C, T, ONE, DELTA, short_word
from borel_exchange_certificate import audit_profile_permutations, nonroot_cover_certificate
from borel_exchange_stress import LogField, canonical_gain
from borel_folded_experiment import from_point
from borel_reversal_experiment import PointDiagram


def symbolic_controls():
    u = C + T
    assert short_word((1, 1), True) == (ONE + T * u, u + T * (ONE + u**2), u, ONE + u**2)
    records = []
    for signs in product((1, -1), repeat=2):
        for reflected in (False, True):
            word = short_word(signs, reflected)
            expected = DELTA if reflected else (C + T)**2 if signs[0] == signs[1] else C**2
            assert word[0] + word[3] == expected
            assert word[2] == (C + T if signs[0] == signs[1] else C)
            records.append({'signs': signs, 'reflected': reflected,
                            'trace': repr(expected), 'lower_left': repr(word[2])})
    return records


def artin_schreier_roots(field):
    roots = {}
    for z in range(field.q):
        roots.setdefault(field.mul(z, z) ^ z, z)
    assert set(roots) == {v for v in range(field.q) if field.traces[v] == 0}
    return roots


def construct_square(field, t, c, cycle, position, labels, roots):
    """O(1) operations after field/trace/root/cycle preprocessing."""
    delta = field.mul(c, c ^ t)
    assert c not in (0, t)
    assert field.traces[field.inverse[c]] == field.traces[field.inverse[c ^ t]] == 1
    if field.traces[field.inverse[delta]]:
        return None
    # W = R (XA)^2, with u=c+t and trace Delta.
    u = c ^ t
    u2 = field.mul(u, u)
    matrix = (1 ^ field.mul(t, u), u ^ field.mul(t, 1 ^ u2), u, 1 ^ u2)
    inverse_delta = field.inverse[delta]
    y = roots[field.mul(inverse_delta, inverse_delta)]
    eigenvalue = field.mul(delta, y)
    z = field.mul(eigenvalue ^ matrix[3], field.inverse[u])
    assert field.act(matrix, z) == z
    points = [z, z ^ c]
    points.append(field.act((0, 1, 1, t), points[-1]))
    assert points[-1] != field.q
    points.append(points[-1] ^ c)
    assert field.act((0, 1, 1, t), points[-1]) == z ^ t

    def fold(point):
        i = position[point]
        return min(i, field.q + 1 - i) - 1

    vertices = [fold(p) for p in points]
    assert min(vertices) >= 0 and len(set(vertices)) == 4
    assert abs(vertices[1] - vertices[2]) == abs(vertices[3] - vertices[0]) == 1
    cuts = sorted((min(vertices[1:3]), min(vertices[3], vertices[0])))
    i, j = cuts
    assert i + 1 < j
    chords = {tuple(sorted(vertices[:2])), tuple(sorted(vertices[2:]))}

    def ratio(a, b):
        return field.mul(labels[a], field.inverse[labels[b]])

    if chords == {(i, j), (i + 1, j + 1)}:
        gain = field.mul(ratio(i, j), ratio(i + 1, j + 1))
        result = {'kind': 'crossing', 'spanning': True, 'root_vertices': field.q // 2,
                  'offroot': [], 'colours': True}
    else:
        assert chords == {(i, j + 1), (i + 1, j)}
        gain = ratio(i, j + 1)
        target = position[cycle[i + 2] ^ c]
        assert target in (j + 1, field.q - j)
        sheet = int(target != j + 1)
        length = j - i
        result = {'kind': 'nested', 'spanning': False, 'root_vertices': field.q // 2 - length,
                  'offroot': [{'folded_length': length, 'sheet_parity': sheet,
                               'point_lengths': [2 * length] if sheet else [length, length]}],
                  'colours': length % 2 == 0 or sheet == 1}
    assert gain != 1
    return {**result, 'cuts': [v + 1 for v in cuts], 'gain': gain,
            'beta': field.mul(t, gain ^ 1), 'fixed_point': z, 'short_word': matrix}


def verify_square(point, folded, square):
    graph = list(alternating_exchanges(folded, 2))
    assert len(graph) == int(square is not None)
    if square is None:
        return None
    exchange = graph[0]
    assert sorted(v + 1 for v in exchange['inserted']) == square['cuts']
    root = folded.root_path(exchange['matching'])
    assert root['gain'] == square['gain'] and len(root['vertices']) == square['root_vertices']
    cover = nonroot_cover_certificate(point, folded, exchange['matching'], root['vertices'])
    assert cover == square['offroot']
    profile = point.profile(folded.lift(point, exchange['matching']))
    assert point.good(profile) == square['colours']
    assert tuple(profile[0]['word_matrix']) == (1, square['beta'], 0, 1)
    assert sorted(p['length'] for p in profile) == sorted(
        [2 * square['root_vertices'] + 1] + [ell for c in cover for ell in c['point_lengths']])
    audit_profile_permutations(point, profile)
    return profile


def kloosterman(field, argument):
    assert argument
    return sum(1 - 2 * field.traces[u ^ field.mul(argument, field.inverse[u])]
               for u in range(1, field.q))


def full_cayley_controls():
    """Directly traverse both 32736-vertex Cayley graphs, without word orders."""
    field, t, records = LogField(32), 7, []
    for c in (26, 29):
        point = PointDiagram(field, t, c, field.cycle(t))
        folded = from_point(point)
        square = construct_square(field, t, c, point.cycle, point.position,
                                  folded.labels, artin_schreier_roots(field))
        exchange, = alternating_exchanges(folded, 2)
        selected = folded.lift(point, exchange['matching'])
        generators = (point.a, point.x, point.x_inverse)
        elements, index = [(1, 0, 0, 1)], {(1, 0, 0, 1): 0}
        for g in elements:
            for s in generators:
                other = field.multiply(s, g)
                if other not in index:
                    index[other] = len(elements)
                    elements.append(other)
        assert len(elements) == 32736
        chord_at = {v: e for e in range(point.n, len(point.edges)) for v in point.edges[e]}
        matching, complement, colours = [], [], {}
        for vertex, g in enumerate(elements):
            i = point.position[field.act(g, field.q)]
            chosen = (i == 0 or chord_at[i] in selected,
                      i in selected, (i - 1) % point.n in selected)
            assert sum(chosen) == 1
            neighbours = [index[field.multiply(s, g)] for s in generators]
            assert vertex not in neighbours and len(set(neighbours)) == 3
            mate = next(v for v, keep in zip(neighbours, chosen) if keep)
            matching.append(mate)
            complement.append([v for v, keep in zip(neighbours, chosen) if not keep])
            colours[tuple(sorted((vertex, mate)))] = 0
        assert all(matching[matching[v]] == v for v in range(len(elements)))
        assert all(v in complement[w] for v, inc in enumerate(complement) for w in inc)
        unseen, lengths = set(range(len(elements))), Counter()
        while unseen:
            start = current = min(unseen)
            previous, length = -1, 0
            while True:
                unseen.remove(current)
                other = next(v for v in complement[current] if v != previous)
                edge = tuple(sorted((current, other)))
                assert edge not in colours
                colours[edge] = 1 + length % 2
                previous, current = current, other
                length += 1
                if current == start:
                    break
            lengths[length] += 1
        good = all(length % 2 == 0 for length in lengths)
        assert good == square['colours']
        if good:
            assert all({colours[tuple(sorted((v, w)))] for w in inc + [matching[v]]} == {0, 1, 2}
                       for v, inc in enumerate(complement))
        records.append({'q': 32, 't': t, 'c': c, 'vertices': len(elements),
                        'cycle_length_counts': dict(lengths), 'colours': good})
    assert records[0]['cycle_length_counts'] == {38: 496, 434: 32}
    assert records[1]['cycle_length_counts'] == {38: 496, 217: 64}
    return records


def audit(q, canonical_audit):
    start = perf_counter()
    field, seen, totals, weighted, traces, controls = LogField(q), set(), Counter(), Counter(), [], []
    roots = artin_schreier_roots(field)
    for t in range(1, q):
        if t in seen:
            continue
        orbit = field.frobenius_orbit(t)
        seen.update(orbit)
        cycle = field.cycle(t)
        if len(cycle) != q + 1:
            continue
        position = {z: i for i, z in enumerate(cycle)}
        labels = [field.mul(z, z) ^ field.mul(t, z) ^ 1 for z in cycle[1:q // 2 + 1]]
        label_position = {k: i for i, k in enumerate(labels)}
        sampled, counts = set(), Counter()
        for c in range(1, q):
            if c == t:
                continue
            counts['generating_parameters'] += 1
            if not field.traces[field.inverse[c]] or not field.traces[field.inverse[c ^ t]]:
                counts['original_trace_successes'] += 1
                continue
            counts['residual_parameters'] += 1
            if c > (c ^ t):
                continue
            squares = [construct_square(field, t, v, cycle, position, labels, roots) for v in (c, c ^ t)]
            assert (squares[0] is None) == (squares[1] is None)
            if squares[0] is not None:
                assert squares[0]['cuts'] == squares[1]['cuts']
                assert squares[0]['gain'] == squares[1]['gain']
                assert squares[0]['colours'] or squares[1]['colours']
                counts['square_parameters'] += 2
                counts['at_least_one_per_pair'] += 1
                counts['square_successes'] += sum(s['colours'] for s in squares)
                counts['square_failures'] += sum(not s['colours'] for s in squares)
                for v, square in zip((c, c ^ t), squares):
                    signature = (square['kind'], tuple((r['folded_length'] % 2, r['sheet_parity'])
                                                       for r in square['offroot']))
                    counts['spanning_square_parameters' if square['spanning'] else 'nonspanning_square_parameters'] += 1
                    canonical = None
                    if canonical_audit and not square['colours']:
                        canonical = canonical_gain(field, labels, label_position, field.mul(v, v ^ t))
                        counts['failed_squares_with_canonical_audit'] += 1
                        if canonical == 1:
                            counts['square_and_canonical_failures'] += 1
                            print('SQUARE AND CANONICAL FAIL', q, t, v, flush=True)
                    if q <= 128 or signature not in sampled or canonical == 1:
                        point = PointDiagram(field, t, v, cycle)
                        profile = verify_square(point, from_point(point), square)
                        counts['independent_graph_and_word_audits'] += 1
                        if signature not in sampled or canonical == 1:
                            controls.append({'t': t, 'c': v, **square, 'profile': profile,
                                             'canonical_gain_if_tested': canonical})
                        sampled.add(signature)
            elif q <= 128 or 'no_square' not in sampled:
                point = PointDiagram(field, t, c, cycle)
                verify_square(point, from_point(point), None)
                counts['independent_no_square_audits'] += 1
                sampled.add('no_square')
        b = field.inverse[t]
        b2 = field.mul(b, b)
        b3, b4 = field.mul(b2, b), field.mul(b2, b2)
        k0, k1, k2 = [kloosterman(field, arg) for arg in (b2, b3 ^ b4, b2 ^ b4)]
        assert 4 * counts['original_trace_successes'] == 3 * q - 5 - k0
        assert 8 * counts['square_parameters'] == q - 3 + 2 * k0 + 2 * k1 + k2
        guaranteed = counts['original_trace_successes'] + counts['square_parameters'] // 2
        assert 16 * guaranteed == 13 * q - 23 - 2 * k0 + 2 * k1 + k2
        counts['combined_successes'] = counts['original_trace_successes'] + counts['square_successes']
        counts['guaranteed_successes'] = guaranteed
        assert counts['combined_successes'] >= guaranteed
        traces.append({'t': t, 'frobenius_weight': len(orbit), 'kloosterman': [k0, k1, k2],
                       'counts': dict(counts)})
        totals.update(counts)
        weighted.update({k: len(orbit) * v for k, v in counts.items()})
    full_traces = sum(gcd(v, q + 1) == 1 for v in range(1, q + 1)) // 2
    assert sum(r['frobenius_weight'] for r in traces) == full_traces
    assert weighted['generating_parameters'] == (q - 2) * full_traces
    result = {'q': q, 'modulus': field.modulus, 'status': 'COMPLETE',
              'scope': 'all full-projective-cycle generating parameters, modulo Frobenius on t; both c partners tested',
              'canonical_audit_on_failed_squares': canonical_audit, 'counts': dict(totals),
              'frobenius_weighted_counts': dict(weighted), 'traces': traces, 'controls': controls,
              'seconds': perf_counter() - start}
    print('SQUARE AUDIT', q, dict(totals), 'seconds', round(result['seconds'], 2), flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('q', type=int, nargs='*', default=[4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096])
    parser.add_argument('--audit-canonical-failures', action='store_true')
    parser.add_argument('--json-output', type=Path)
    args = parser.parse_args()
    output = {'symbolic_controls': symbolic_controls(), 'full_cayley_controls': full_cayley_controls(), 'fields': []}
    print('FULL CAYLEY CONTROLS', output['full_cayley_controls'], flush=True)
    for q in args.q:
        output['fields'].append(audit(q, args.audit_canonical_failures))
        if args.json_output:
            args.json_output.write_text(json.dumps(output, indent=2) + '\n')
