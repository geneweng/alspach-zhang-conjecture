"""Independent short-word certificate for the five-chord obstruction.

Enumerate fixed points of every R^eta X^eps_r A ... X^eps_1 A,
not paths in the folded graph, and compare the resulting simple cycles
with the existing graph DFS. Also audit quotient words as permutations
of projective cycle indices, independently of matrix multiplication.
"""

import argparse
import json
from collections import Counter
from itertools import product
from pathlib import Path

from borel_alternating_exchange import alternating_exchanges, inspect_exchange
from borel_exchange_stress import LogField, cut_certificate
from borel_folded_experiment import from_point
from borel_reflected_family import half_word_certificate
from borel_reversal_experiment import PointDiagram


def audit_profile_permutations(point, profile):
    """Compose x-index shifts and the explicit translation permutation."""
    a = [point.position[z ^ point.c] if z != point.field.q else 0
         for z in point.cycle]
    for component in profile:
        permutation, shift = list(range(point.n)), 0
        for symbol in component['labels']:
            if symbol == 'a':
                permutation = [a[(v + shift) % point.n] for v in permutation]
                shift = 0
            else:
                shift += 1 if symbol == 'x' else -1
        permutation = [(v + shift) % point.n for v in permutation]
        matrix = component['word_matrix']
        assert permutation == [point.position[point.field.act(matrix, z)] for z in point.cycle]
    return len(profile)


def word_cycles(point, folded, maximum):
    """All simple alternating cycles, independently via quadratic equations.

    Every folded cycle lifts from either point over any of its 2r vertices.
    Its path steps specify the signs, and the last sheet specifies eta.
    Thus each cycle must occur exactly 4r times in this enumeration.
    """
    field, occurrences = point.field, Counter()
    chord_at = {v: e for e in range(folded.h - 1, len(folded.edges))
                for v in folded.edges[e][:2]}

    def fold(z):
        i = point.position[z]
        return min(i, point.n - i) - 1

    for r in range(1, maximum + 1):
        for signs in product((1, -1), repeat=r):
            word = (1, 0, 0, 1)
            for sign in signs:
                word = field.multiply(point.a, word)
                word = field.multiply(point.x if sign == 1 else point.x_inverse, word)
            for reflected in (False, True):
                matrix = field.multiply((1, point.t, 0, 1), word) if reflected else word
                p, b, ell, d = matrix
                for start in range(field.q):
                    if field.mul(ell, field.mul(start, start)) ^ field.mul(p ^ d, start) ^ b:
                        continue
                    z, visited, cuts, chords = start, set(), set(), set()
                    for j, sign in enumerate(signs):
                        u, v = fold(z), fold(z ^ point.c)
                        if u < 0 or v < 0 or u == v or u in visited or v in visited:
                            break
                        visited.update((u, v))
                        e = chord_at[u]
                        assert set(folded.edges[e][:2]) == {u, v}
                        chords.add(e)
                        z = field.act(point.x if sign == 1 else point.x_inverse, z ^ point.c)
                        w = fold(z)
                        if w < 0 or abs(w - v) != 1:
                            break
                        cuts.add(min(v, w))
                        if w in visited and not (j == r - 1 and w == fold(start)):
                            break
                    else:
                        assert fold(z) == fold(start)
                        assert z == (start ^ point.t if reflected else start)
                        assert len(visited) == 2 * r and len(cuts) == len(chords) == r
                        occurrences[tuple(sorted(cuts))] += 1
    assert all(count == 4 * len(cuts) for cuts, count in occurrences.items())
    return occurrences


def small_controls():
    records = []
    for q, t, c in ((4, 2, 1), (8, 2, 1), (16, 8, 4), (32, 7, 1)):
        field = LogField(q)
        cycle = field.cycle(t)
        assert len(cycle) == q + 1
        point = PointDiagram(field, t, c, cycle)
        folded = from_point(point)
        words = word_cycles(point, folded, 6)
        graph = {tuple(sorted(e['inserted'])) for e in alternating_exchanges(folded, 6)}
        assert set(words) == graph
        records.append({'q': q, 't': t, 'c': c, 'cycles': len(words)})
    return records


def nonroot_cover_certificate(point, folded, selected, root_vertices):
    """Predict point-circuit lengths using the two-sheeted folded cover."""
    unseen, records = set(range(folded.h)) - set(root_vertices), []
    while unseen:
        start = current = min(unseen)
        previous, length, sheet = -1, 0, 0
        while True:
            unseen.remove(current)
            length += 1
            options = [e for e in folded.inc[current] if e not in selected and e != previous]
            e = options[0]
            u, v, kind = folded.edges[e]
            other = v if u == current else u
            if kind == 'a':
                target = point.position[point.cycle[current + 1] ^ point.c]
                assert target in (other + 1, point.n - other - 1)
                sheet ^= int(target != other + 1)
            current, previous = other, e
            if current == start:
                break
        records.append({'folded_length': length, 'sheet_parity': sheet,
                        'point_lengths': [2 * length] if sheet else [length, length]})
    return records


def obstruction():
    q, t, c = 4096, 681, 1207
    field = LogField(q)
    assert field.modulus == 4105
    assert field.traces[field.inverse[c]] == field.traces[field.inverse[c ^ t]] == 1
    cycle = field.cycle(t)
    assert len(cycle) == q + 1 and c != t
    point = PointDiagram(field, t, c, cycle)
    folded = from_point(point)
    canonical = point.profile(set(range(1, q, 2)))
    assert tuple(canonical[0]['word_matrix']) == (1, 0, 0, 1)
    assert not point.good(canonical)
    assert not point.good(json.loads(json.dumps(canonical)))
    permutation_audits = audit_profile_permutations(point, canonical)
    words = word_cycles(point, folded, 6)
    graph = {tuple(sorted(e['inserted'])): e for e in alternating_exchanges(folded, 6)}
    assert set(words) == set(graph)
    records, histogram = [], Counter()
    for cuts, exchange in sorted(graph.items(), key=lambda pair: (len(pair[0]), pair[0])):
        contracted = cut_certificate(folded, exchange)
        direct = inspect_exchange(folded, exchange)
        assert contracted['spanning'] == direct['spanning']
        assert contracted['gain'] == direct['gain']
        assert contracted['root_vertices'] == len(direct['path'])
        profile = point.profile(folded.lift(point, exchange['matching']))
        permutation_audits += audit_profile_permutations(point, profile)
        cover = nonroot_cover_certificate(point, folded, exchange['matching'], direct['path'])
        assert sorted(p['length'] for p in profile) == sorted(
            [2 * len(direct['path']) + 1] + [ell for r in cover for ell in r['point_lengths']])
        even_off_root = all(r['folded_length'] % 2 == 0 or r['sheet_parity'] for r in cover)
        if even_off_root and direct['gain'] != 1:
            assert point.good(profile)
        half = half_word_certificate(point, profile)
        assert half['beta'] == field.mul(t, contracted['gain'] ^ 1)
        if contracted['spanning']:
            assert len(profile) == 1 and profile[0]['length'] == q + 1
        kind = 'nonspanning' if not direct['spanning'] else 'unit' if direct['gain'] == 1 else 'good'
        histogram[f'{len(cuts)}_{kind}'] += 1
        records.append({**contracted, 'cuts': [i + 1 for i in cuts],
                        'fixed_point_word_occurrences': words[cuts],
                        'profile': profile, 'half_word': half, 'nonroot_cover': cover,
                        'even_off_root': even_off_root, 'colours_upstairs': point.good(profile)})
    assert dict(histogram) == {'3_nonspanning': 1, '4_nonspanning': 3,
                               '5_nonspanning': 2, '6_nonspanning': 1, '6_good': 1}
    good = next(r for r in records if r['spanning'])
    assert good['cuts'] == [358, 725, 1060, 1073, 1624, 2041]
    assert good['gain'] == 2930
    return {'q': q, 'modulus': field.modulus, 't': t, 'c': c,
            'delta': folded.delta, 'canonical_profile': canonical,
            'scope': 'all single simple alternating cycles with at most six chords',
            'histogram': dict(histogram), 'cycles': records,
            'independent_permutation_component_audits': permutation_audits}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--json-output', type=Path)
    args = parser.parse_args()
    controls = small_controls()
    print('SMALL WORD/GRAPH CONTROLS', controls, flush=True)
    case = obstruction()
    print('CERTIFIED OBSTRUCTION', case['q'], case['t'], case['c'], case['histogram'], flush=True)
    print('SIX-CHORD REPAIR', next(r['half_word'] for r in case['cycles'] if r['spanning']), flush=True)
    if args.json_output:
        args.json_output.write_text(json.dumps({'small_controls': controls, 'obstruction': case}, indent=2) + '\n')
