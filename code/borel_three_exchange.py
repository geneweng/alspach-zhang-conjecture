"""Exact six-cycle counts, a short-exchange obstruction, and affine voltages.

The density result counts absence of short alternating cycles from the
all-chord matching. It does NOT count uncolourable Cayley graphs, or imply
failure of the canonical matching. Large-field trace counts are complete
modulo Frobenius; graph and word audits have the separately reported scope.
"""

import argparse
import json
from collections import Counter
from itertools import product
from pathlib import Path

from borel_alternating_exchange import alternating_exchanges
from borel_exchange_algebra import C, T, ONE, DELTA, short_word, verify_identities
from borel_exchange_certificate import audit_profile_permutations, word_cycles
from borel_exchange_stress import LogField
from borel_folded_experiment import from_point
from borel_reversal_experiment import PointDiagram


def symbolic_controls():
    assert len(verify_identities()) == 16
    # The loop--bridge--triangle correction, with k_j=t+Delta.
    assert DELTA * (DELTA + T + ONE)**2 + T**2 * (DELTA + T) == (
        DELTA**3 + DELTA + T**3)
    # Partial fractions: 1/[d^2(d+v)] = v^-2(1/d+1/(d+v))+v^-1/d^2.
    d, v = C + ONE, T + ONE
    assert d * (d + v) + d**2 + v * (d + v) == v**2
    return {'short_word_identities': 16, 'additional_polynomial_identities': 2}


def trace_data(field, t, c):
    u = c ^ t
    assert c not in (0, 1, t, t ^ 1)
    delta = field.mul(c, u)
    tau1 = field.mul(u, field.mul(c ^ 1, c ^ 1))
    tau2 = field.mul(c, field.mul(u ^ 1, u ^ 1))
    traces = [field.traces[field.inverse[v]] for v in (c, u, delta, tau1, tau2)]
    exceptional = field.mul(field.mul(delta, delta), delta) ^ delta ^ field.mul(field.mul(t, t), t)
    return {'delta': delta, 'traces': traces, 'degenerate': exceptional == 0,
            'three_cycles': 2 - traces[3] - traces[4] - int(exceptional == 0)}


def affine_certificate(point, folded, selected, audit_matrices=False):
    """Exact lift-colourability certificate for ANY folded perfect matching.

    In the frame X^i, a directed chord i->j has affine voltage
       y -> (k_j/k_i)y + w_j + (k_j/k_i)w_i,  w_i=t+z_i.
    Path-edge voltages are identity. Sheet labels are kept separately.
    """
    field, root = point.field, folded.root_path(selected)
    unseen = set(range(folded.h)) - set(root['vertices'])
    records = []
    powers = [(1, 0, 0, 1)]
    if audit_matrices:
        for _ in range(folded.h):
            powers.append(field.multiply(point.x, powers[-1]))

    def inverse(matrix):
        a, b, c, d = matrix
        return d, b, c, a

    while unseen:
        start = current = min(unseen)
        previous, alpha, nu, sheet, chords = -1, 1, 0, 0, 0
        vertices, actual_word = [], (1, 0, 0, 1)
        while True:
            unseen.remove(current)
            vertices.append(current)
            edge = next(e for e in folded.inc[current] if e not in selected and e != previous)
            u, v, kind = folded.edges[edge]
            other = v if current == u else u
            z, w = point.cycle[current + 1], point.cycle[other + 1]
            if kind == 'a':
                chords += 1
                step_alpha = field.mul(folded.labels[other], field.inverse[folded.labels[current]])
                step_nu = (point.t ^ w) ^ field.mul(step_alpha, point.t ^ z)
                alpha, nu = field.mul(step_alpha, alpha), field.mul(step_alpha, nu) ^ step_nu
                target = point.position[z ^ point.c]
                assert target in (other + 1, point.n - other - 1)
                step_sheet = int(target != other + 1)
                if audit_matrices:
                    corrected = (1, z ^ w, 0, 1)
                    bmat = field.multiply(inverse(powers[other + 1]),
                                          field.multiply(corrected, powers[current + 1]))
                    assert bmat[2] == 0
                    assert field.mul(bmat[0], bmat[0]) == step_alpha
                    assert field.mul(bmat[0], bmat[1]) == step_nu
                    actual_word = field.multiply(point.a, actual_word)
                sheet ^= step_sheet
            elif audit_matrices:
                direction = (1 if other > current else -1) * (-1 if sheet else 1)
                actual_word = field.multiply(point.x if direction == 1 else point.x_inverse, actual_word)
            previous, current = edge, other
            if current == start:
                break
        if audit_matrices:
            corrected_word = field.multiply((1, point.t, 0, 1), actual_word) if sheet else actual_word
            bmat = field.multiply(inverse(powers[start + 1]),
                                  field.multiply(corrected_word, powers[start + 1]))
            assert bmat[2] == 0
            assert (field.mul(bmat[0], bmat[0]), field.mul(bmat[0], bmat[1])) == (alpha, nu)
        length = len(vertices)
        good = length % 2 == 0 or sheet == 1 or (alpha == 1 and nu != 0)
        records.append({'vertices': [v + 1 for v in vertices], 'folded_length': length,
                        'chords': chords, 'sheet': sheet, 'alpha': alpha, 'nu': nu,
                        'point_lengths': [2 * length] if sheet else [length, length],
                        'lifts_even': good})
    return {'root_vertices': len(root['vertices']), 'root_gain': root['gain'],
            'root_beta': field.mul(point.t, root['gain'] ^ 1), 'offroot': records,
            'colours': root['gain'] != 1 and all(r['lifts_even'] for r in records)}


def audit_matching(point, folded, selected, permutations=False):
    certificate = affine_certificate(point, folded, selected, audit_matrices=True)
    profile = point.profile(folded.lift(point, selected))
    assert certificate['colours'] == point.good(profile)
    assert tuple(profile[0]['word_matrix']) == (1, certificate['root_beta'], 0, 1)
    assert sorted(p['length'] for p in profile) == sorted(
        [2 * certificate['root_vertices'] + 1] +
        [ell for r in certificate['offroot'] for ell in r['point_lengths']])
    if permutations:
        audit_profile_permutations(point, profile)
    return certificate, profile


def all_matching_controls():
    result = []
    for q in (4, 8, 16, 32):
        field, seen, counts = LogField(q), set(), Counter()
        for t in range(1, q):
            if t in seen:
                continue
            seen.update(field.frobenius_orbit(t))
            cycle = field.cycle(t)
            if len(cycle) != q + 1:
                continue
            for c in range(1, q):
                if c == t:
                    continue
                point = PointDiagram(field, t, c, cycle)
                folded = from_point(point)
                counts['parameters'] += 1
                for selected in folded.matchings():
                    cert, profile = audit_matching(point, folded, selected, permutations=True)
                    counts['matchings'] += 1
                    counts['permutation_component_audits'] += len(profile)
                    counts['colouring_matchings'] += int(cert['colours'])
                    counts['odd_point_unipotent_rescues'] += sum(
                        r['folded_length'] % 2 == 1 and r['sheet'] == 0 and
                        r['alpha'] == 1 and r['nu'] != 0 for r in cert['offroot'])
        result.append({'q': q, 'counts': dict(counts)})
        print('ALL-MATCHING AFFINE AUDIT', q, dict(counts), flush=True)
    return result


def graph_audit(field, t, c, cycle, data, words=False):
    point = PointDiagram(field, t, c, cycle)
    folded = from_point(point)
    exchanges = {tuple(sorted(e['inserted'])): e for e in alternating_exchanges(folded, 3)}
    threes = {cuts: e for cuts, e in exchanges.items() if len(cuts) == 3}
    assert len(threes) == data['three_cycles']
    assert all(len(cuts) != 1 for cuts in exchanges)
    assert sum(len(cuts) == 2 for cuts in exchanges) == 1 - data['traces'][2]
    h_partner = next(u if v == folded.h - 1 else v for u, v, kind in folded.edges
                     if kind == 'a' and folded.h - 1 in (u, v))
    triangle = any(kind == 'a' and {u, v} == {h_partner - 1, h_partner + 1}
                   for u, v, kind in folded.edges)
    assert triangle == data['degenerate']
    raw = 0
    matching_audits = permutation_components = 0
    if field.q <= 128:
        partner = PointDiagram(field, t, c ^ t, cycle)
        for exchange in threes.values():
            for diagram in (point, partner):
                cert, profile = audit_matching(diagram, folded, exchange['matching'], permutations=True)
                matching_audits += 1
                permutation_components += len(profile)
                assert all(r['chords'] <= 2 and r['alpha'] != 1 for r in cert['offroot'])
                assert cert['colours'] == (cert['root_gain'] != 1 and all(
                    r['folded_length'] % 2 == 0 or r['sheet'] == 1 for r in cert['offroot']))
    if words:
        occurrences = word_cycles(point, folded, 3)
        assert set(occurrences) == set(exchanges)
        for reflected in (False, True):
            for signs in product((1, -1), repeat=3):
                matrix = tuple(p.evaluate(field, c, t) for p in short_word(signs, reflected))
                assert matrix[0] != matrix[3]
                raw += sum(field.act(matrix, z) == z for z in range(field.q + 1))
        assert raw == 12 * (len(threes) + int(data['degenerate']))
    return {'graph_cycles': len(threes), 'word_audited': words,
            'raw_three_word_fixed_points': raw if words else None,
            'matching_audits': matching_audits, 'permutation_components': permutation_components}


def character_bound_audit(field, t, histogram):
    """Audit all 31 simple-pole character sums, without fitting an error term."""
    b = field.inverse[t]
    a = field.mul(field.inverse[t ^ 1], field.inverse[t ^ 1])
    root = field.exp[(field.log[field.inverse[t ^ 1]] * (field.q // 2)) % (field.q - 1)]
    d = a ^ root
    assert b not in (0, 1) and d != 0
    poles = (0, t, 1, t ^ 1)
    assert sum(histogram.values()) == field.q - 4
    maxima = Counter()
    for mask in range(1, 32):
        e1, e2, e3, e4, e5 = [(mask >> k) & 1 for k in range(5)]
        residues = (e1 ^ (b if e3 else 0) ^ (a if e5 else 0),
                    e2 ^ (b if e3 else 0) ^ (a if e4 else 0),
                    d if e4 else 0, d if e5 else 0)
        k = sum(r != 0 for r in residues)
        assert k >= 1
        upper = (int(bool(e1 or e3 or e5)) + int(bool(e2 or e3 or e4)) + e4 + e5)
        assert k <= upper
        maxima[upper] += 1
        value = sum(count * (-1 if (pattern & mask).bit_count() % 2 else 1)
                    for pattern, count in histogram.items())
        missing = 0
        for p, residue in zip(poles, residues):
            if residue:
                continue
            evaluation = 0
            for pole, r in zip(poles, residues):
                if r:
                    evaluation ^= field.mul(r, field.inverse[p ^ pole])
            missing += -1 if field.traces[evaluation] else 1
        assert (value + 1 + missing)**2 <= 4 * (k - 1)**2 * field.q
    assert dict(maxima) == {1: 2, 2: 9, 3: 12, 4: 8}
    m = histogram[31]
    # |32M-(q-3)| <= 114 sqrt(q)+36; no floating point in the check.
    excess = max(0, abs(32 * m - (field.q - 3)) - 36)
    assert excess**2 <= 114**2 * field.q
    return m


def audit_field(q, exhaustive_graph_limit=512, first_trace_only=False):
    field, seen, counts, records = LogField(q), set(), Counter(), []
    for t in range(1, q):
        if t in seen:
            continue
        orbit = field.frobenius_orbit(t)
        seen.update(orbit)
        cycle = field.cycle(t)
        if len(cycle) != q + 1:
            continue
        histogram, sampled, row = Counter(), set(), Counter()
        a = field.mul(field.inverse[t ^ 1], field.inverse[t ^ 1])
        d = a ^ field.exp[(field.log[field.inverse[t ^ 1]] * (q // 2)) % (q - 1)]
        for c in range(1, q):
            if c in (1, t, t ^ 1):
                continue
            data = trace_data(field, t, c)
            traces = data['traces']
            # Independent simple-pole representatives of the last two traces.
            assert traces[3] == field.traces[field.mul(a, field.inverse[c ^ t]) ^
                                             field.mul(d, field.inverse[c ^ 1])]
            assert traces[4] == field.traces[field.mul(a, field.inverse[c]) ^
                                             field.mul(d, field.inverse[c ^ t ^ 1])]
            pattern = sum(v << k for k, v in enumerate(traces))
            histogram[pattern] += 1
            row['parameters'] += 1
            if traces[:2] != [1, 1]:
                continue
            row['residual'] += 1
            assert data['three_cycles'] in (0, 1, 2)
            row['three_cycles'] += data['three_cycles']
            row['degenerate_parameters'] += int(data['degenerate'])
            if traces[2] == 1 and data['three_cycles'] == 0:
                row['no_exchange_up_to_three'] += 1
                row['degenerate_additions'] += int(pattern != 31)
            if c > (c ^ t):
                continue
            signature = tuple(traces[2:]) + (data['degenerate'],)
            if q <= exhaustive_graph_limit or signature not in sampled:
                checked = graph_audit(field, t, c, cycle, data, words=q <= 32 or data['degenerate'])
                sampled.add(signature)
                row['graph_parameter_audits'] += 1
                row['word_parameter_audits'] += int(q <= 32 or data['degenerate'])
                row['affine_matching_audits'] += checked['matching_audits']
                row['permutation_component_audits'] += checked['permutation_components']
        m = character_bound_audit(field, t, histogram)
        assert row['no_exchange_up_to_three'] == m + row['degenerate_additions']
        assert 0 <= row['degenerate_additions'] <= 6
        row['five_trace_obstructions'] = m
        for key, value in row.items():
            counts[key] += value
            if not key.endswith('_audits'):
                counts['weighted_' + key] += value * len(orbit)
        counts['trace_representatives'] += 1
        records.append({'t': t, 'frobenius_orbit': orbit, 'counts': dict(row),
                        'trace_histogram': [histogram[i] for i in range(32)]})
        if first_trace_only:
            assert q < 16384 or m > 0
            break
    print('THREE-CYCLE COUNT AUDIT', q, dict(counts), flush=True)
    return {'q': q, 'modulus': field.modulus, 'counts': dict(counts), 'traces': records,
            'scope': 'first full-cycle t only' if first_trace_only else 'all full-cycle t modulo Frobenius',
            'graph_scope': f'exhaustive for q<={exhaustive_graph_limit}; one per trace/degeneracy signature per t above that',
            'word_scope': 'all graph-audited parameters through q=32, and audited degeneracies'}


def paired_obstruction():
    q, t, c, cuts = 128, 8, 84, (18, 33, 46)
    field, results = LogField(q), []
    cycle = field.cycle(t)
    for translation in (c, c ^ t):
        point = PointDiagram(field, t, translation, cycle)
        folded = from_point(point)
        exchange = next(e for e in alternating_exchanges(folded, 3)
                        if tuple(i + 1 for i in sorted(e['inserted'])) == cuts)
        cert, profile = audit_matching(point, folded, exchange['matching'], permutations=True)
        assert cert['root_gain'] != 1 and not cert['colours']
        assert [r['folded_length'] for r in cert['offroot']] == [15, 13]
        assert all(r['chords'] == 1 and r['alpha'] != 1 for r in cert['offroot'])
        assert [r['sheet'] for r in cert['offroot']] == ([0, 1] if translation == c else [1, 0])
        results.append({'c': translation, 'certificate': cert, 'profile': profile})
    assert results[0]['certificate']['root_gain'] == results[1]['certificate']['root_gain'] == 124
    for left, right in zip(results[0]['certificate']['offroot'], results[1]['certificate']['offroot']):
        assert (left['alpha'], left['nu']) == (right['alpha'], right['nu'])
    print('PAIRED THREE-CHORD OBSTRUCTION', q, t, c, cuts, flush=True)
    return {'q': q, 'modulus': field.modulus, 't': t, 'cuts': cuts, 'partners': results,
            'scope': 'this exchange fails for both partners; neither graph is claimed uncolourable'}


def unipotent_rescue():
    """A real matching accepted by the affine test but not the parity test."""
    field = LogField(64)
    point = PointDiagram(field, 2, 6, field.cycle(2))
    folded = from_point(point)
    selected = frozenset((0, 2, 4, 6, 8, 10, 12, 14, 16, 19, 21, 24, 27, 30, 44, 46))
    cert, profile = audit_matching(point, folded, selected, permutations=True)
    assert cert['colours'] and cert['root_gain'] == 30
    assert len(cert['offroot']) == 1
    offroot = cert['offroot'][0]
    assert (offroot['folded_length'], offroot['sheet'], offroot['alpha'], offroot['nu']) == (21, 0, 1, 48)
    assert [r['length'] for r in profile] == [23, 21, 21]
    assert all(tuple(r['word_matrix']) != (1, 0, 0, 1) and
               field.multiply(r['word_matrix'], r['word_matrix']) == (1, 0, 0, 1) for r in profile)
    print('ODD-POINT UNIPOTENT RESCUE', field.q, point.t, point.c, flush=True)
    return {'q': field.q, 'modulus': field.modulus, 't': point.t, 'c': point.c,
            'selected_edge_ids_zero_based': sorted(selected),
            'selected_path_starts': [e + 1 for e in sorted(selected) if folded.edges[e][2] == 'x'],
            'selected_chords': [[u + 1, v + 1] for e in sorted(selected)
                                for u, v, kind in [folded.edges[e]] if kind == 'a'],
            'certificate': cert, 'profile': profile,
            'scope': 'this matching colours despite odd off-root point circuits; not a new parameter family'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('q', nargs='*', type=int, default=[8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096])
    parser.add_argument('--json-output', type=Path)
    args = parser.parse_args()
    result = {'symbolic_controls': symbolic_controls(),
              'affine_all_matching_controls': all_matching_controls(),
              'paired_obstruction': paired_obstruction(),
              'unipotent_rescue': unipotent_rescue(),
              'fields': [audit_field(q) for q in args.q],
              'large_positive_control': audit_field(16384, first_trace_only=True)}
    if args.json_output:
        args.json_output.write_text(json.dumps(result, indent=2) + '\n')
