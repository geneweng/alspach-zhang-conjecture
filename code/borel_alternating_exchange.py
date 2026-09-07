"""Explore alternating-cycle exchanges from the all-chord folded matching.

An exchange removes r a-chords and selects r disjoint path edges. Its
complement need not span the folded graph. Spanning and non-unit gain
are checked separately, without assuming either follows from cycle length.
"""

from borel_folded_experiment import from_point
from borel_reflected_family import BinaryField
from borel_reversal_experiment import PointDiagram
from borel_spanning_exchange import absolute_trace
from borel_exchange_algebra import EXCEPTION, short_word, verify_identities


def alternating_exchanges(folded, maximum):
    """Enumerate simple alternating cycles with at most maximum chords.

    Starting at the smallest cycle vertex and traversing its chord first
    makes each undirected cycle occur once. Parallel edges stay distinct.
    """
    chords = frozenset(range(folded.h - 1, len(folded.edges)))
    chord_at = {v: e for e in chords for v in folded.edges[e][:2]}
    partner = {}
    for e in chords:
        u, v, _ = folded.edges[e]
        partner[u], partner[v] = v, u

    def walk(start, current, visited, removed, inserted):
        other = partner[current]
        if other < start or other in visited:
            return
        removed = removed | {chord_at[current]}
        visited = visited | {current, other}
        for next_vertex in (other - 1, other + 1):
            if not start <= next_vertex < folded.h:
                continue
            path_edge = min(other, next_vertex)
            added = inserted | {path_edge}
            if next_vertex == start:
                yield {'removed': removed, 'inserted': added,
                       'matching': (chords - removed) | added}
            elif next_vertex not in visited and len(removed) < maximum:
                yield from walk(start, next_vertex, visited, removed, added)

    for start in range(folded.h):
        yield from walk(start, start, frozenset(), frozenset(), frozenset())


def inspect_exchange(folded, exchange):
    root = folded.root_path(exchange['matching'])
    return {'chords': len(exchange['removed']),
            'spanning': len(root['vertices']) == folded.h,
            'gain': root['gain'], 'path': root['vertices']}


def cycle_certificate(point, folded, exchange):
    """Traverse the short alternating cycle and check its fixed-point slope."""
    cycle_edges = exchange['removed'] | exchange['inserted']
    vertices = {v for e in cycle_edges for v in folded.edges[e][:2]}
    start = min(vertices)
    current, previous, z = start, -1, point.cycle[start + 1]
    start_z = z
    word, gain, symbols, oriented = (1, 0, 0, 1), 1, [], {}

    def fold(z):
        i = point.position[z]
        return min(i, point.n - i) - 1

    while True:
        options = [e for e in folded.inc[current] if e in cycle_edges and e != previous]
        e = next(e for e in options if folded.edges[e][2] == 'a') if previous == -1 else options[0]
        u, v, kind = folded.edges[e]
        other = v if current == u else u
        if kind == 'a':
            symbol, matrix = 'a', point.a
            oriented[e] = (current, other)
            gain = point.field.times[gain][point.field.times[folded.labels[current]][
                point.field.inverse[folded.labels[other]]]]
        else:
            options = [(symbol, matrix) for symbol, matrix in (('x', point.x), ('X', point.x_inverse))
                       if fold(point.field.act(matrix, z)) == other]
            assert len(options) == 1
            symbol, matrix = options[0]
        symbols.append(symbol)
        z = point.field.act(matrix, z)
        assert fold(z) == other
        word = point.field.multiply(matrix, word)
        previous, current = e, other
        if current == start:
            break
    assert len(symbols) == 2 * len(exchange['removed'])
    assert z in (start_z, start_z ^ point.t)
    reflected = z != start_z
    if reflected:
        word = point.field.multiply((1, point.t, 0, 1), word)
    assert point.field.act(word, start_z) == start_z
    denominator = point.field.times[word[2]][start_z] ^ word[3]
    derivative = point.field.inverse[point.field.times[denominator][denominator]]
    assert derivative == gain
    assert (gain == 1) == (word[0] == word[3])
    root = folded.root_path(exchange['matching'])
    path_edges = {(u, v) for u, v in zip(root['vertices'], root['vertices'][1:])}
    signs = [1 if (u, v) in path_edges else -1 if (v, u) in path_edges else 0
             for u, v in oriented.values()]
    coherent = len(root['vertices']) == folded.h and len(set(signs)) == 1 and signs[0] != 0
    if coherent:
        assert root['gain'] == (gain if signs[0] == 1 else point.field.inverse[gain])
    return {'labels': ''.join(symbols), 'reflected': reflected,
            'matrix': word, 'trace': word[0] ^ word[3], 'gain': gain,
            'coherent': coherent, 'orientation_signs': signs}


def validate_enumerator(folded, maximum):
    """Independent exhaustive perfect-matchings check of the cycle DFS."""
    all_chords = frozenset(range(folded.h - 1, len(folded.edges)))
    expected = set()
    for matching in folded.matchings():
        difference = matching ^ all_chords
        if not difference or len(difference) > 2 * maximum:
            continue
        vertices = {v for e in difference for v in folded.edges[e][:2]}
        reached, pending = set(), [min(vertices)]
        while pending:
            v = pending.pop()
            if v in reached:
                continue
            reached.add(v)
            pending.extend(w for e in folded.inc[v] if e in difference
                           for w in folded.edges[e][:2] if w != v)
        if reached == vertices:
            expected.add(matching)
    actual = [e['matching'] for e in alternating_exchanges(folded, maximum)]
    assert len(actual) == len(set(actual))
    assert set(actual) == expected
    return len(actual)


def audit_three(q, field):
    from collections import Counter
    counts, seen, controls = Counter(), set(), []
    for t in range(1, q):
        if t in seen:
            continue
        seen.update(field.frobenius_orbit(t))
        cycle = field.cycle(t)
        if len(cycle) != q + 1:
            continue
        sampled = set()
        for c in range(1, q):
            if c == t or c > (c ^ t):
                continue
            # For larger fields concentrate on the class left by the trace test.
            if q >= 128 and (absolute_trace(field, field.inverse[c]) == 0 or
                             absolute_trace(field, field.inverse[c ^ t]) == 0):
                continue
            point = PointDiagram(field, t, c, cycle)
            folded = from_point(point)
            counts['folded_parameter_representatives'] += 1
            if q <= 16:
                counts['enumerator_control_cycles'] += validate_enumerator(folded, 5)
            exception = EXCEPTION.evaluate(field, c, t)
            pentagon = c in (1, t ^ 1)
            for exchange in alternating_exchanges(folded, 3):
                if len(exchange['removed']) != 3:
                    continue
                counts['three_chord_cycles'] += 1
                result = inspect_exchange(folded, exchange)
                if not result['spanning']:
                    continue
                counts['spanning'] += 1
                certificate = cycle_certificate(point, folded, exchange)
                signs = tuple(1 if s == 'x' else -1 for s in certificate['labels'] if s != 'a')
                polynomial_word = short_word(signs, certificate['reflected'])
                assert tuple(v.evaluate(field, c, t) for v in polynomial_word) == certificate['matrix']
                assert not certificate['coherent'] or pentagon or result['gain'] != 1
                if not pentagon and exception:
                    assert result['gain'] != 1
                    counts['quartic_certified'] += 1
                if result['gain'] == 1:
                    assert pentagon or exception == 0
                    counts['unit_gain'] += 1
                key = ('unit' if result['gain'] == 1 else 'nonunit', pentagon)
                audit = q <= 32 or key not in sampled or (result['gain'] == 1 and not pentagon)
                if audit:
                    profile = point.profile(folded.lift(point, exchange['matching']), independent=True)
                    beta = field.times[t][1 ^ result['gain']]
                    assert len(profile) == 1 and profile[0]['length'] == q + 1
                    assert tuple(profile[0]['word_matrix']) == (1, beta, 0, 1)
                    counts['independent_word_audits'] += 1
                    sampled.add(key)
                    if result['gain'] == 1 and not pentagon:
                        controls.append({'t': t, 'c': c, 'delta': folded.delta,
                                         'cuts': sorted(e + 1 for e in exchange['inserted']),
                                         'quartic': exception, 'profile': profile,
                                         'short_cycle': certificate})
    result = {'q': q, 'modulus': field.modulus, 'counts': dict(counts),
              'scope': 'all generating parameters' if q <= 64 else 'double-inverse-trace-one parameters',
              'symmetry': 'Frobenius on t and the folded c <-> c+t pairing',
              'nonpentagon_unit_controls': controls}
    print('THREE-EXCHANGE-AUDIT', q, dict(counts), flush=True)
    return result


def audit_residual(record, field, maximum):
    from collections import Counter
    results = []
    for case in record['residual_canonical_failures']:
        t, c = case['t'], case['c']
        point = PointDiagram(field, t, c, field.cycle(t))
        folded = from_point(point)
        counts, witness, controls = Counter(), None, []
        assert folded.root_path(folded.canonical)['gain'] == 1
        for exchange in alternating_exchanges(folded, maximum):
            result = inspect_exchange(folded, exchange)
            kind = 'nonspanning' if not result['spanning'] else 'unit' if result['gain'] == 1 else 'good'
            counts[f'{result["chords"]}_{kind}'] += 1
            if not result['spanning']:
                continue
            if kind == 'unit' or witness is None or result['chords'] < witness['chords']:
                profile = point.profile(folded.lift(point, exchange['matching']), independent=True)
                assert len(profile) == 1 and profile[0]['length'] == field.q + 1
                assert tuple(profile[0]['word_matrix']) == (1, field.times[t][1 ^ result['gain']], 0, 1)
                certified = {**result, 'cuts': sorted(e + 1 for e in exchange['inserted']),
                             'matching': sorted(folded.lift(point, exchange['matching'])),
                             'profile': profile, 'short_cycle': cycle_certificate(point, folded, exchange)}
                if kind == 'unit':
                    controls.append(certified)
                elif witness is None or result['chords'] < witness['chords']:
                    witness = certified
        results.append({'t': t, 'c': c, 'orbit_weight': case['orbit_weight'],
                        'counts': dict(counts), 'witness': witness, 'unit_controls': controls})
    print('RESIDUAL-EXCHANGE-AUDIT', field.q, 'cases', len(results),
          'successes', sum(r['witness'] is not None for r in results), flush=True)
    return {'q': field.q, 'maximum': maximum, 'cases': results}


if __name__ == '__main__':
    import argparse
    import json
    from pathlib import Path
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('q', type=int, nargs='*', default=[4, 8, 16, 32, 64, 128, 256, 512, 1024])
    parser.add_argument('--maximum', type=int, default=5)
    parser.add_argument('--spanning-results', type=Path, default=Path('code/borel_spanning_results.json'))
    parser.add_argument('--json-output', type=Path)
    args = parser.parse_args()
    assert args.maximum >= 3
    stored = {r['q']: r for r in json.loads(args.spanning_results.read_text())['fields']}
    output = {'symbolic_identities': verify_identities(), 'three_chord_audits': [], 'residual_audits': []}
    for q in args.q:
        field = BinaryField(q)
        output['three_chord_audits'].append(audit_three(q, field))
        if q in stored:
            output['residual_audits'].append(audit_residual(stored[q], field, args.maximum))
        if args.json_output:
            args.json_output.write_text(json.dumps(output, indent=2) + '\n')
