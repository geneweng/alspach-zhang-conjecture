"""Signed root-path lattices for the seven-change family.

A row records the oriented folded chords of an eligible root path. Its dot
product with chord discrete logarithms is the root logarithm. An exact
integer combination equal to a power of two times one coordinate therefore
forces a non-unit root gain in any odd-order gain group with non-unit gain
on that chord. This is a conditional certificate, not a uniform existence
theorem for the required row-lattice relation.
The uniform coordinate-inclusion target for this fixed family is disproved
at q=2048 in borel_gain_rank.py; the conditional certificate remains valid.
"""

import argparse
import json
import re
from collections import Counter
from math import gcd, prod
from pathlib import Path

from borel_exchange_stress import LogField
from borel_folded_experiment import FoldedDiagram, from_point
from borel_kempe_exchange import closed_graph, eligible_certificates, normalized
from borel_reversal_experiment import PointDiagram
from borel_three_exchange import audit_matching
from borel_two_buffer import seven_switches


def root_row(folded, graph, colours, excluded):
    """Signed chord incidences, in orientation from folded 1 to h."""
    assert excluded != colours[-1]
    row = [0] * (folded.h // 2)
    current, previous, seen = 0, -1, set()
    while current != folded.h - 1:
        assert current not in seen
        seen.add(current)
        edge = next(e for e in graph.inc[current]
                    if e != previous and e != len(graph.edges) - 1 and colours[e] != excluded)
        u, v, kind = graph.edges[edge]
        other = v if current == u else u
        if kind == 'a':
            coordinate = edge - (folded.h - 1)
            assert row[coordinate] == 0
            row[coordinate] = 1 if current == u else -1
        previous, current = edge, other
    return row


def chord_logs(folded):
    field = folded.field
    return [(field.log[folded.labels[u]] - field.log[folded.labels[v]]) % (field.q - 1)
            for u, v, kind in folded.edges if kind == 'a']


def family_rows(folded, audit=False, initial=None):
    """Deduplicate nonzero rows up to sign, retaining replayable source seeds."""
    graph, natural = closed_graph(folded)
    start = natural if initial is None else initial
    assert graph.proper(start) and set(start) == {0, 1, 2}
    records, seen, counts = [], set(), Counter()
    logs = chord_logs(folded)

    def insert(colours, source):
        for excluded in sorted({0, 1, 2} - {colours[-1]}):
            raw = root_row(folded, graph, colours, excluded)
            counts['eligible_paths'] += 1
            if audit:
                selected = frozenset(e for e, c in enumerate(colours) if c == excluded)
                root = folded.root_path(selected)
                predicted = folded.field.exp[sum(a*b for a, b in zip(raw, logs)) % (folded.field.q - 1)]
                assert predicted == root['gain']
            if not any(raw):
                continue
            sign = 1 if next(v for v in raw if v) > 0 else -1
            row = tuple(sign * v for v in raw)
            if row in seen:
                continue
            seen.add(row)
            records.append({'source': source, 'excluded_colour': excluded,
                            'orientation': sign, 'row': list(row)})

    insert(start, {'kind': 'natural' if initial is None else 'initial'})
    for pivot, colour in enumerate(start):
        if colour == 0:
            continue
        result = seven_switches(graph, start, pivot)
        counts[result['status']] += 1
        for branch, candidate in enumerate(result['candidates']):
            insert(candidate['colours'], {'kind': 'seven', 'pivot_edge_id': pivot, 'branch': branch})
    return records, dict(counts)


def extended_gcd(a, b):
    old_r, r, old_s, s, old_t, t = a, b, 1, 0, 0, 1
    while r:
        q = old_r // r
        old_r, r = r, old_r - q*r
        old_s, s = s, old_s - q*s
        old_t, t = t, old_t - q*t
    if old_r < 0:
        return -old_r, -old_s, -old_t
    return old_r, old_s, old_t


def row_lattice(rows, track=False):
    """Integer row-echelon lattice basis using unimodular Euclidean steps.

    If track=True, also retain its exact expression in the input rows.
    Reduction of entries above each pivot controls coefficient growth.
    """
    width, total = len(rows[0]), len(rows)
    basis, transforms = {}, {}
    for index, source in enumerate(rows):
        row = list(source)
        coeff = [int(i == index) for i in range(total)] if track else None
        while any(row):
            column = next(i for i, v in enumerate(row) if v)
            if column not in basis:
                if row[column] < 0:
                    row = [-v for v in row]
                    if track:
                        coeff = [-v for v in coeff]
                basis[column] = row
                if track:
                    transforms[column] = coeff
                break
            old = basis[column]
            a, b = old[column], row[column]
            if b % a == 0:
                quotient = b // a
                row = [v - quotient*w for v, w in zip(row, old)]
                if track:
                    coeff = [v - quotient*w for v, w in zip(coeff, transforms[column])]
            else:
                common, s, t = extended_gcd(a, b)
                basis[column] = [s*x + t*y for x, y in zip(old, row)]
                row = [-(b//common)*x + (a//common)*y for x, y in zip(old, row)]
                if track:
                    before = transforms[column]
                    transforms[column] = [s*x + t*y for x, y in zip(before, coeff)]
                    coeff = [-(b//common)*x + (a//common)*y for x, y in zip(before, coeff)]
        # Reduce the basis to row Hermite form after every insertion.
        for column in sorted(basis):
            pivot = basis[column][column]
            for earlier in sorted(k for k in basis if k < column):
                quotient = basis[earlier][column] // pivot
                if not quotient:
                    continue
                basis[earlier] = [v - quotient*w for v, w in zip(basis[earlier], basis[column])]
                if track:
                    transforms[earlier] = [v - quotient*w for v, w in
                                           zip(transforms[earlier], transforms[column])]
    result = {'rank': len(basis), 'columns': width,
              'pivot_columns': sorted(basis), 'basis': [basis[k] for k in sorted(basis)]}
    if len(basis) == width:
        result['index'] = prod(basis[k][k] for k in range(width))
    if track:
        result['transforms'] = [transforms[k] for k in sorted(basis)]
    return result


def basis_expression(target, lattice):
    """Exact integer membership and expression in the echelon basis."""
    target, factors = list(target), []
    for pivot, basis in zip(lattice['pivot_columns'], lattice['basis']):
        if target[pivot] % basis[pivot]:
            return None
        factor = target[pivot] // basis[pivot]
        factors.append(factor)
        target = [a-factor*b for a, b in zip(target, basis)]
    return None if any(target) else factors


def doubled_coordinates(lattice):
    width = lattice['columns']
    return [j for j in range(width) if basis_expression(
        [2*int(i == j) for i in range(width)], lattice) is not None]


def audit_lattice(rows, lattice):
    """Check both inclusions when transforms are available, without row reduction."""
    assert all(basis_expression(row, lattice) is not None for row in rows)
    for basis, transform in zip(lattice['basis'], lattice.get('transforms', [])):
        assert all(sum(c*r[j] for c, r in zip(transform, rows)) == value
                   for j, value in enumerate(basis))
    return {'input_rows_in_basis': True,
            'basis_expressions_checked': len(lattice.get('transforms', []))}


def coordinate_relation(rows, lattice, coordinate=0, multiple=2):
    """Express multiple*e_coordinate in the input rows, or return None."""
    factors = basis_expression([multiple * int(i == coordinate)
                                for i in range(len(rows[0]))], lattice)
    if factors is None:
        return None
    result = [sum(factor*transform[i] for factor, transform in
                  zip(factors, lattice['transforms'])) for i in range(len(rows))]
    # Independent direct dot-product check of the complete integer certificate.
    for column in range(len(rows[0])):
        assert sum(c*r[column] for c, r in zip(result, rows)) == multiple * int(column == coordinate)
    return result


def short_coordinate_relation(rows):
    """Find a multiple 1 or 2 of one coordinate using at most two rows."""
    width = len(rows[0])
    lookup = {bytes([1]*width): (None, 0)}
    for i, row in enumerate(rows):
        for sign in (1, -1):
            lookup[bytes(1+sign*x for x in row)] = (i, sign)
    for i, row in enumerate(rows):
        encoded = bytes(1+x for x in row)
        for j, value in enumerate(encoded):
            for multiple in (1, 2):
                if value < multiple:
                    continue
                other = encoded[:j] + bytes([value-multiple]) + encoded[j+1:]
                if other not in lookup:
                    continue
                k, sign = lookup[other]
                relation = [int(a == i) for a in range(len(rows))]
                if k is not None:
                    relation[k] -= sign
                assert all(sum(a*r[column] for a, r in zip(relation, rows)) ==
                           multiple*int(column == j) for column in range(width))
                return {'coordinate': j, 'multiple': multiple, 'relation': relation}
    return None


def both_unit_control():
    """An actual projective colouring, not asserted to be a seven-change return."""
    field = LogField(128)
    point = PointDiagram(field, 8, 20, field.cycle(8))
    folded = from_point(point)
    graph, natural = closed_graph(folded)
    colours = [
        1,0,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,0,2,1,0,1,2,0,1,2,0,2,1,2,1,
        0,2,1,2,1,2,1,0,1,2,1,2,0,2,1,0,2,1,2,0,2,1,2,1,2,1,2,0,2,0,1,0,
        2,1,0,0,0,0,0,0,0,0,0,0,0,0,2,1,0,2,1,2,0,1,0,0,2,1,0,0,0,1,1,2]
    assert graph.proper(colours) and normalized(colours) != normalized(natural)
    assert folded.root_path(folded.canonical)['gain'] == 1
    cases = []
    for c in (20, 28):
        partner = PointDiagram(field, 8, c, field.cycle(8))
        for item in eligible_certificates(partner, folded, colours):
            cert, profile = audit_matching(partner, folded, frozenset(item['selected']),
                                           permutations=True)
            assert cert == item['certificate'] and cert['root_gain'] == 1
            cases.append({'c': c, **item, 'profile': profile})
    print('BOTH UNIT', [(v['c'], v['colour'], v['certificate']['root_vertices'])
                        for v in cases], flush=True)
    return {'q': 128, 'modulus': field.modulus, 't': 8, 'c': 20,
            'colours_by_edge_id': colours, 'cases': cases,
            'scope': 'non-natural proper colouring; not a claimed direct seven-change return'}


def abstract_control():
    field = LogField(16)
    labels = [1, 14, 6, 10, 13, 9, 5, 2]
    folded = FoldedDiagram(field, labels, 11)
    records, counts = family_rows(folded, audit=True)
    logs = chord_logs(folded)
    assert [r['row'] for r in records] == [[1, 1, -1, 1]]
    assert all(sum(a*b for a, b in zip(r['row'], logs)) % 15 == 0 for r in records)
    return {'q': 16, 'modulus': field.modulus, 'labels': labels, 'delta': 11,
            'rows': records, 'counts': counts, 'rank': 1, 'columns': 4,
            'all_eligible_gains_one': True, 'scope': 'reordered labels, not a Cayley instance'}


def lattice_control(q, t, c, track=True):
    field = LogField(q)
    cycle = field.cycle(t)
    assert len(cycle) == q + 1 and c not in (0, t)
    point = PointDiagram(field, t, c, cycle)
    folded = from_point(point)
    records, counts = family_rows(folded, audit=True)
    rows = [r['row'] for r in records]
    lattice = row_lattice(rows, track=track)
    checked = audit_lattice(rows, lattice)
    coordinates = doubled_coordinates(lattice)
    index = lattice.get('index')
    result = {'q': q, 'modulus': field.modulus, 't': t, 'c': c,
              'counts': counts, 'row_count': len(rows), 'rank': lattice['rank'],
              'canonical_gain': folded.root_path(folded.canonical)['gain'],
              'columns': lattice['columns'], 'index': index,
              'index_coprime_to_q_minus_one': index is not None and gcd(index, q-1) == 1,
              'doubled_coordinates': coordinates, 'audit': checked}
    if track:
        coordinate = coordinates[0] if coordinates else None
        relation = coordinate_relation(rows, lattice, coordinate) if coordinates else None
        result.update({'rows': records, 'basis': lattice['basis'], 'coordinate': coordinate,
                       'multiple': 2, 'relation': relation})
    print('LATTICE', q, t, c, 'rows', len(rows), 'rank', lattice['rank'],
          'index', index, 'doubled_coordinates', len(coordinates), flush=True)
    return result


def stored_controls(cache):
    directory, output = Path(__file__).resolve().parent, []
    for filename in ('borel_reflected_results.json', 'borel_reflected_1024.json'):
        for record in json.loads((directory / filename).read_text()):
            grouped = {}
            for case in record['canonical_failures']:
                key = (case['t'], min(case['c'], case['c'] ^ case['t']))
                group = grouped.setdefault(key, {'direct': 0, 'weighted': 0, 'c': case['c']})
                group['direct'] += 1
                group['weighted'] += case['orbit_weight']
            cases, totals = [], Counter()
            for (t, _), group in grouped.items():
                q, c = record['q'], group['c']
                result = cache.get((q, t, c))
                if result is None:
                    result = lattice_control(q, t, c, track=False)
                assert result['canonical_gain'] == 1
                result = {k: v for k, v in result.items()
                          if k not in ('rows', 'basis', 'relation', 'coordinate', 'multiple')}
                result.update(group)
                cases.append(result)
                for key in ('direct', 'weighted'):
                    totals[key] += group[key]
                totals['folded_diagrams'] += 1
                totals['has_doubled_coordinate'] += bool(result['doubled_coordinates'])
                totals['coprime_index'] += result['index_coprime_to_q_minus_one']
            output.append({'q': record['q'], 'counts': dict(totals), 'cases': cases})
    return output


def verify_saved(output):
    """Regenerate root rows and directly check serialized integer relations."""
    for case in output['cases']:
        field = LogField(case['q'])
        point = PointDiagram(field, case['t'], case['c'], field.cycle(case['t']))
        folded = from_point(point)
        records, _ = family_rows(folded, audit=True)
        assert records == case['rows']
        rows = [r['row'] for r in records]
        coordinate, multiple, relation = case['coordinate'], case['multiple'], case['relation']
        assert relation is not None
        assert all(sum(a*r[j] for a, r in zip(relation, rows)) == multiple*int(j == coordinate)
                   for j in range(folded.h//2))
        case['minimum_nonzero_root_chords'] = min(sum(x != 0 for x in row) for row in rows)
        case['canonical_gain'] = folded.root_path(folded.canonical)['gain']
        case['inverse_traces_c_partner_delta'] = [field.traces[field.inverse[v]]
            for v in (case['c'], case['c'] ^ case['t'], folded.delta)]
        short = short_coordinate_relation(rows)
        if short:
            graph, natural = closed_graph(folded)
            short['realizations'] = []
            for i, coefficient in enumerate(short['relation']):
                if not coefficient:
                    continue
                record, source = records[i], records[i]['source']
                colours = natural if source['kind'] == 'natural' else seven_switches(
                    graph, natural, source['pivot_edge_id'])['candidates'][source['branch']]['colours']
                selected = frozenset(e for e, c in enumerate(colours) if c == record['excluded_colour'])
                cert, profile = audit_matching(point, folded, selected, permutations=True)
                short['realizations'].append({'row_index': i, 'coefficient': coefficient,
                    'selected_edge_ids': sorted(selected), 'certificate': cert, 'profile': profile})
        case['short_coordinate_certificate'] = short
        print('VERIFIED SAVED', field.q, case['t'], case['c'], 'short', short is not None,
              'minimum_root_chords', case['minimum_nonzero_root_chords'], flush=True)
    return output


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--json-output', type=Path)
    parser.add_argument('--stored', action='store_true',
                        help='also audit all stored q<=1024 canonical failures (several minutes)')
    parser.add_argument('--verify-results', type=Path,
                        help='regenerate rows and check saved relations, without redoing lattices')
    parser.add_argument('--case', type=int, nargs=3, action='append', metavar=('Q', 'T', 'C'),
                        help='add a specified exact lattice control')
    args = parser.parse_args()
    if args.verify_results:
        output = json.loads(args.verify_results.read_text())
    else:
        output = {'scope': 'specified exact lattice controls; no uniform rank theorem', 'cases': []}
        for q, t, c in ((64, 7, 26), (128, 8, 20), (128, 15, 35), (256, 15, 151),
                        (256, 25, 13), (512, 47, 10), (512, 53, 66),
                        (1024, 49, 80), (1024, 49, 910), (1024, 62, 836)):
            output['cases'].append(lattice_control(q, t, c))
        output['both_unit_control'] = both_unit_control()
        output['abstract_control'] = abstract_control()
    for q, t, c in args.case or []:
        assert not any((r['q'], r['t'], r['c']) == (q, t, c) for r in output['cases'])
        output['cases'].append(lattice_control(q, t, c))
    if args.stored:
        cache = {(r['q'], r['t'], r['c']): r for r in output['cases']}
        output['stored'] = stored_controls(cache)
    verify_saved(output)
    if args.json_output:
        # Keep integer vectors on one line, avoiding multi-megabyte indentation.
        serialized = json.dumps(output, indent=2)
        serialized = re.sub(r'\[\s*((?:-?\d+\s*,\s*)*-?\d+)\s*\]',
                            lambda match: '[' + ', '.join(re.findall(r'-?\d+', match[1])) + ']',
                            serialized)
        args.json_output.write_text(serialized + '\n')
