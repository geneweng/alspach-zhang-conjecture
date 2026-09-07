"""Audit spanning-path exchanges and their uniform trace-test family.

The starting folded matching selects ALL a-chords, unlike the older
canonical matching, which selects only path edges. A parallel-edge exchange
leaves a spanning complementary path using one chord. A crossing square
exchange leaves such a path using two chords. Either gain is non-unit.

python3 code/borel_spanning_exchange.py --json-output code/borel_spanning_results.json

All normalized generating parameters are tested, modulo Frobenius only.
Larger-field matrix audits use deterministic samples and every previously
stored canonical failure; the JSON distinguishes these from scalar checks.
"""

import argparse
import json
from collections import Counter
from math import gcd
from pathlib import Path

from borel_folded_experiment import FoldedDiagram, from_point
from borel_reflected_family import BinaryField, half_word_certificate
from borel_reversal_experiment import PointDiagram


def absolute_trace(field, value):
    result = 0
    for _ in range(field.degree):
        result ^= value
        value = field.times[value][value]
    assert result in (0, 1)
    return result


def exchange_matchings(folded):
    """Return all one-chord and crossing two-chord spanning exchanges."""
    chords = frozenset(range(folded.h - 1, len(folded.edges)))
    edge_at = {v: e for e in chords for v in folded.edges[e][:2]}
    partner = {}
    for e in chords:
        u, v, _ = folded.edges[e]
        partner[u], partner[v] = v, u
    parallel, squares = [], []
    for i in range(folded.h - 1):
        j = partner[i]
        if j == i + 1:
            parallel.append({'kind': 'parallel', 'cuts': [i],
                             'matching': (chords - {edge_at[i]}) | {i}})
        if i + 1 < j < folded.h - 1 and partner[i + 1] == j + 1:
            squares.append({'kind': 'square', 'cuts': [i, j],
                            'matching': (chords - {edge_at[i], edge_at[i + 1]}) | {i, j}})
    return parallel, squares


def spanning_certificate(folded, candidate):
    """Check topology and calculate gain without the point-quotient code."""
    selected = candidate['matching']
    assert all(sum(e in selected for e in inc) == 1 for inc in folded.inc)
    if candidate['kind'] == 'parallel':
        path = list(range(folded.h))
        chord_steps = [(candidate['cuts'][0], candidate['cuts'][0] + 1)]
    else:
        i, j = candidate['cuts']
        path = list(range(i + 1)) + list(range(j, i, -1)) + list(range(j + 1, folded.h))
        chord_steps = [(i, j), (i + 1, j + 1)]
    assert len(path) == len(set(path)) == folded.h
    steps = set(chord_steps)
    gain = 1
    for u, v in zip(path, path[1:]):
        kind = 'a' if (u, v) in steps else 'x'
        possible = [e for e in folded.inc[u] if e not in selected
                    and v in folded.edges[e][:2] and folded.edges[e][2] == kind]
        assert len(possible) == 1
        if kind == 'a':
            ratio = folded.field.times[folded.labels[u]][folded.field.inverse[folded.labels[v]]]
            gain = folded.field.times[gain][ratio]
    assert gain != 1
    return {'vertices': path, 'gain': gain, 'traversed_chords': len(chord_steps)}


def matrix_audit(point, folded, candidate, scalar, independent=False):
    assert folded.root_path(candidate['matching']) == {
        'vertices': scalar['vertices'], 'gain': scalar['gain']}
    lifted = folded.lift(point, candidate['matching'])
    profile = point.profile(lifted, independent=independent)
    assert len(profile) == 1 and profile[0]['length'] == point.n
    beta = point.field.times[point.t][1 ^ scalar['gain']]
    assert tuple(profile[0]['word_matrix']) == (1, beta, 0, 1)
    assert half_word_certificate(point, profile)['beta'] == beta
    assert point.good(profile)
    return {'matching': sorted(lifted), 'beta': beta,
            'quotient_length': point.n, 'lifted_cycle_length': 2 * point.n}


def full_cayley_audit(point, selected):
    """Construct and colour the FULL Cayley graph, not only its quotient.

    Use left generator multiplication and right Borel cosets, represented
    by g(infinity). Inversion identifies this with the usual right Cayley
    convention. No circuit-word order calculation is used in this audit.
    """
    field = point.field
    generators = (point.a, point.x, point.x_inverse)
    elements, index = [(1, 0, 0, 1)], {(1, 0, 0, 1): 0}
    for g in elements:
        for s in generators:
            other = field.multiply(s, g)
            if other not in index:
                index[other] = len(elements)
                elements.append(other)
    assert len(elements) == field.q * (field.q * field.q - 1)
    chord_at = {v: e for e in range(point.n, len(point.edges)) for v in point.edges[e]}
    matching, complement, colours = [], [], {}
    for v, g in enumerate(elements):
        i = point.position[field.act(g, field.q)]
        chosen = (i == 0 or chord_at[i] in selected,
                  i in selected, (i - 1) % point.n in selected)
        assert sum(chosen) == 1
        neighbours = [index[field.multiply(s, g)] for s in generators]
        assert v not in neighbours and len(set(neighbours)) == 3
        mate = next(w for w, keep in zip(neighbours, chosen) if keep)
        matching.append(mate)
        complement.append([w for w, keep in zip(neighbours, chosen) if not keep])
        colours[tuple(sorted((v, mate)))] = 0
    assert all(matching[matching[v]] == v for v in range(len(elements)))
    assert all(v in complement[w] for v, inc in enumerate(complement) for w in inc)
    unseen, lengths = set(range(len(elements))), []
    while unseen:
        start = current = min(unseen)
        previous, length = -1, 0
        while True:
            unseen.remove(current)
            other = next(w for w in complement[current] if w != previous)
            edge = tuple(sorted((current, other)))
            assert edge not in colours
            colours[edge] = 1 + length % 2
            previous, current = current, other
            length += 1
            if current == start:
                break
        assert length == 2 * point.n
        lengths.append(length)
    assert len(lengths) == field.q * (field.q - 1) // 2
    for v, inc in enumerate(complement):
        assert {colours[tuple(sorted((v, w)))] for w in inc + [matching[v]]} == {0, 1, 2}
    return {'vertices': len(elements), 'cycles': len(lengths), 'cycle_length': 2 * point.n}


def old_root_contained_repair(folded):
    for matching in folded.one_chord_matchings():
        root = folded.root_path(matching)
        endpoints = {v for e in matching if folded.edges[e][2] == 'a'
                     for v in folded.edges[e][:2]}
        if root['gain'] != 1 and endpoints <= set(root['vertices']):
            return matching, root
    raise AssertionError('stored canonical failure has lost its old repair')


def failure_parameters(paths):
    result = {}
    for path in paths:
        for record in json.loads(path.read_text()):
            result[record['q']] = {(case['t'], case['c']): case
                                  for case in record['canonical_failures']}
    return result


def check(q, failures):
    field, seen, totals = BinaryField(q), set(), Counter()
    traces = [absolute_trace(field, v) for v in range(q)]
    inverse_trace = {c: traces[field.inverse[c]] for c in range(1, q)}
    records, residual_failures = [], []
    for t in range(1, q):
        if t in seen:
            continue
        orbit = field.frobenius_orbit(t)
        seen.update(orbit)
        cycle = field.cycle(t)
        if len(cycle) != q + 1:
            continue
        assert inverse_trace[t] == 1
        labels = from_point(PointDiagram(field, t, 1, cycle)).labels
        # This sum omits the added 1 in Ahmadi--Granger's convention.
        parameter = field.times[field.inverse[t]][field.inverse[t]]
        kloosterman = sum(1 - 2 * traces[u ^ field.times[parameter][field.inverse[u]]]
                         for u in range(1, q))
        assert kloosterman * kloosterman <= 4 * q
        count, examples = Counter(), {}
        for c in range(1, q):
            if c == t:
                continue
            folded = FoldedDiagram(field, labels, field.times[c][c ^ t])
            parallel, squares = exchange_matchings(folded)
            trace_pattern = (inverse_trace[c], inverse_trace[c ^ t])
            assert len(parallel) == trace_pattern.count(0)
            count['parameters'] += 1
            count['trace_test'] += bool(parallel)
            count['square_only'] += bool(squares) and not parallel
            count['no_short_exchange'] += not (parallel or squares)
            previous_failure = (t, c) in failures
            if previous_failure:
                count['canonical_failures'] += 1
            kind = 'parallel' if parallel else 'square' if squares else None
            # Check EVERY candidate's explicit spanning path and scalar gain.
            for candidate in parallel + squares:
                scalar = spanning_certificate(folded, candidate)
                count['scalar_certificates'] += 1
                audit = q <= 64 or previous_failure or candidate['kind'] not in examples
                if audit:
                    point = PointDiagram(field, t, c, cycle)
                    independent = q <= 16 or previous_failure or candidate['kind'] not in examples
                    certificate = matrix_audit(point, folded, candidate, scalar, independent)
                    count['matrix_audits'] += 1
                    count['permutation_audits'] += independent
                    if q <= 8:
                        certificate['full_cayley'] = full_cayley_audit(
                            point, folded.lift(point, candidate['matching']))
                        count['full_cayley_audits'] += 1
                    examples.setdefault(candidate['kind'], {
                        'c': c, 'cuts_one_based': [i + 1 for i in candidate['cuts']],
                        'gain': scalar['gain'], **certificate})
            if previous_failure:
                count['canonical_repaired_' + (kind or 'neither')] += 1
                if kind is None:
                    # Keep a genuine residual case, with a successful OLD repair
                    # so that failure of this new subfamily is not misreported.
                    point = PointDiagram(field, t, c, cycle)
                    canonical = point.profile(set(range(1, q, 2)), independent=True)
                    assert tuple(canonical[0]['word_matrix']) == (1, 0, 0, 1)
                    repair, root = old_root_contained_repair(folded)
                    profile = point.profile(folded.lift(point, repair), independent=True)
                    assert point.good(profile)
                    residual_failures.append({
                        't': t, 'c': c, 'orbit_weight': len(orbit),
                        'inverse_trace_pattern': list(trace_pattern),
                        'canonical_length': canonical[0]['length'],
                        'old_repair_chord': [[u + 1, v + 1] for e in repair
                                             for u, v, k in [folded.edges[e]] if k == 'a'][0],
                        'old_repair_gain': root['gain'],
                        'old_repair_profile': profile})
        expected = (3 * q - 5 - kloosterman) // 4
        assert 4 * count['trace_test'] == 3 * q - 5 - kloosterman
        assert count['trace_test'] == expected >= q // 2
        assert count['parameters'] == q - 2
        for key, value in count.items():
            totals[key + '_representatives'] += value
            totals[key + '_covered'] += len(orbit) * value
        records.append({'t': t, 'orbit_weight': len(orbit), 'kloosterman': kloosterman,
                        'counts': dict(count), 'examples': examples})
    assert totals['canonical_failures_representatives'] == len(failures)
    # Full traces occur as {lambda, lambda^-1} for primitive norm-one lambda.
    trace_count = sum(gcd(j, q + 1) == 1 for j in range(1, q + 1)) // 2
    assert sum(record['orbit_weight'] for record in records) == trace_count
    assert totals['parameters_covered'] == (q - 2) * trace_count
    result = {'q': q, 'modulus': field.modulus, 'status': 'COMPLETE',
              'counts': dict(totals), 'traces': records,
              'residual_canonical_failures': residual_failures}
    print('SPANNING', q, dict(totals), flush=True)
    return result


def all_matching_control():
    field = BinaryField(16)
    folded = FoldedDiagram(field, [1, 14, 6, 10, 13, 9, 5, 2], 11)
    records = []
    for matching in folded.matchings():
        root = folded.root_path(matching)
        endpoints = {v for e in matching if folded.edges[e][2] == 'a'
                     for v in folded.edges[e][:2]}
        contained = endpoints <= set(root['vertices'])
        assert not (contained and root['gain'] != 1)
        records.append({'selected_chords': [[u + 1, v + 1] for e in matching
                                           for u, v, k in [folded.edges[e]] if k == 'a'],
                        'path': [v + 1 for v in root['vertices']],
                        'gain': root['gain'], 'root_contained': contained})
    assert len(records) == 3
    return {'scope': 'abstract reordered labels, NOT a projective-cycle instance',
            'q': 16, 'modulus': field.modulus, 'labels': folded.labels,
            'delta': folded.delta, 'all_matchings': records}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('q', type=int, nargs='*', default=[4, 8, 16, 32, 64, 128, 256, 512, 1024])
    parser.add_argument('--reflected-results', type=Path, nargs='*', default=[
        Path('code/borel_reflected_results.json'), Path('code/borel_reflected_1024.json')])
    parser.add_argument('--json-output', type=Path)
    args = parser.parse_args()
    failures = failure_parameters(args.reflected_results)
    output = {'abstract_all_matching_control': all_matching_control(), 'fields': []}
    for q in args.q:
        output['fields'].append(check(q, failures.get(q, {})))
        if args.json_output:
            args.json_output.write_text(json.dumps(output, indent=2) + '\n')
