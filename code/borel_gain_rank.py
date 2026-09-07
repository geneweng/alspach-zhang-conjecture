"""Odd-prime rank certificates using generation by all folded chord gains.

The q/2 distinct oriented chord ratios, together with 1, cannot lie in a
proper subgroup of F_q^*. Hence the chord gains generate F_q^*. Full root
row rank modulo ANY prime divisor of q-1 forces a non-unit candidate gain.
This is a sufficient certificate; no uniform rank theorem is asserted.
"""

import argparse
import json
import random
from math import gcd, prod
from pathlib import Path

from borel_exchange_stress import LogField
from borel_folded_experiment import from_point
from borel_gain_lattice import chord_logs, family_rows
from borel_kempe_exchange import audit_buffer_sequence, change, closed_graph, eligible_certificates
from borel_reversal_experiment import PointDiagram
from borel_three_exchange import audit_matching
from borel_two_buffer import seven_switches


def prime_divisors(n):
    out, divisor = [], 2
    while divisor*divisor <= n:
        if n % divisor == 0:
            out.append(divisor)
            while n % divisor == 0:
                n //= divisor
        divisor += 1
    return out + ([n] if n > 1 else [])


class PackedPrime:
    """Parallel modular arithmetic on independent fixed-width integer lanes."""

    def __init__(self, p, width):
        self.p, self.width = p, width
        self.bits = p.bit_length() + 1
        self.base = 1 << self.bits
        self.high = self.base // 2
        self.ones = sum(1 << (i*self.bits) for i in range(width))
        self.highs = self.high*self.ones
        self.offset = (self.high-p)*self.ones
        self.packed_p = p*self.ones

    def pack(self, row):
        assert len(row) == self.width
        return sum((v % self.p) << (i*self.bits) for i, v in enumerate(row))

    def unpack(self, packed):
        return [(packed >> (i*self.bits)) & (self.base-1) for i in range(self.width)]

    def reduce_sum(self, packed):
        # Each lane is in [0,2p-1]. Adding high-p marks exactly lanes >=p.
        flags = ((packed+self.offset) & self.highs) >> (self.bits-1)
        return packed-self.p*flags

    def add(self, a, b):
        return self.reduce_sum(a+b)

    def subtract(self, a, b):
        return self.reduce_sum(a+self.packed_p-b)

    def scale(self, a, multiplier):
        result = 0
        while multiplier:
            if multiplier & 1:
                result = self.add(result, a)
            multiplier >>= 1
            if multiplier:
                a = self.add(a, a)
        return result


def packed_rank(rows, p, retain_basis=False):
    """Exact echelon elimination over F_p; retain a nonsingular minor's rows."""
    width = len(rows[0])
    arithmetic = PackedPrime(p, width)
    basis, multiples, selected, pivots, factors = {}, {}, [], [], []
    for index, source in enumerate(rows):
        row = arithmetic.pack(source)
        while row:
            pivot = ((row & -row).bit_length()-1) // arithmetic.bits
            value = (row >> (pivot*arithmetic.bits)) & (arithmetic.base-1)
            if pivot not in basis:
                row = arithmetic.scale(row, pow(value, -1, p))
                basis[pivot] = row
                multiples[pivot] = {1: row}
                selected.append(index)
                pivots.append(pivot)
                factors.append(value)
                break
            table = multiples[pivot]
            if value not in table:
                table[value] = arithmetic.scale(basis[pivot], value)
            row = arithmetic.subtract(row, table[value])
        if len(basis) == width:
            break
    result = {'prime': p, 'rank': len(basis), 'columns': width,
              'selected_row_indices': selected, 'pivot_columns_in_insertion_order': pivots}
    if len(basis) == width:
        inversions = sum(a > b for i, a in enumerate(pivots) for b in pivots[i+1:])
        result['minor_determinant_mod_prime'] = ((-1)**inversions * prod(factors)) % p
        assert result['minor_determinant_mod_prime'] != 0
    if retain_basis:
        result['echelon_basis'] = {pivot: arithmetic.unpack(row) for pivot, row in basis.items()}
    return result


def ordinary_rank(rows, p):
    """Independent scalar elimination, used only on bounded audit matrices."""
    matrix = [[v % p for v in row] for row in rows]
    rank, determinant = 0, 1
    for column in range(len(matrix[0])):
        pivot = next((i for i in range(rank, len(matrix)) if matrix[i][column]), None)
        if pivot is None:
            continue
        if pivot != rank:
            matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
            determinant = -determinant
        value = matrix[rank][column]
        determinant = determinant*value % p
        inverse = pow(value, -1, p)
        matrix[rank] = [inverse*v % p for v in matrix[rank]]
        for i in range(rank+1, len(matrix)):
            value = matrix[i][column]
            if value:
                matrix[i] = [(a-value*b) % p for a, b in zip(matrix[i], matrix[rank])]
        rank += 1
    return rank, determinant % p


def arithmetic_controls():
    rng, counts = random.Random(9062026), {'packed_operations': 0, 'rank_checks': 0}
    for p in (3, 5, 7, 11, 23, 31, 89, 127):
        for width in (1, 2, 7, 19):
            arithmetic = PackedPrime(p, width)
            for _ in range(20):
                a, b = [[rng.randrange(p) for _ in range(width)] for _ in range(2)]
                multiplier = rng.randrange(p)
                pa, pb = arithmetic.pack(a), arithmetic.pack(b)
                assert arithmetic.unpack(arithmetic.add(pa, pb)) == [(x+y) % p for x, y in zip(a, b)]
                assert arithmetic.unpack(arithmetic.subtract(pa, pb)) == [(x-y) % p for x, y in zip(a, b)]
                assert arithmetic.unpack(arithmetic.scale(pa, multiplier)) == [x*multiplier % p for x in a]
                counts['packed_operations'] += 3
            for height in (width, width+3):
                rows = [[rng.randrange(-p, p) for _ in range(width)] for _ in range(height)]
                packed = packed_rank(rows, p)
                assert packed['rank'] == ordinary_rank(rows, p)[0]
                if packed['rank'] == width:
                    minor = [rows[i] for i in packed['selected_row_indices']]
                    assert packed['minor_determinant_mod_prime'] == ordinary_rank(minor, p)[1]
                counts['rank_checks'] += 1
    saved = json.loads(Path(__file__).with_name('borel_gain_lattice_results.json').read_text())
    root_checks = []
    for case in saved['cases']:
        if case['q'] > 512:
            continue
        rows = [record['row'] for record in case['rows']]
        prime = prime_divisors(case['q']-1)[0]
        packed = packed_rank(rows, prime)
        assert packed['rank'] == ordinary_rank(rows, prime)[0]
        if packed['rank'] == len(rows[0]):
            minor = [rows[i] for i in packed['selected_row_indices']]
            assert packed['minor_determinant_mod_prime'] == ordinary_rank(minor, prime)[1]
        root_checks.append({'q': case['q'], 't': case['t'], 'c': case['c'],
                            'prime': prime, 'rank': packed['rank']})
    counts['saved_root_matrix_checks'] = root_checks
    return counts


def kernel_coordinate_cover(rows, p, seed=9062026):
    """Kernel vectors excluding coordinate isolation over F_p, checked by dot products.

    If every coordinate is nonzero in some recorded kernel vector, no
    nonzero multiple of ANY coordinate lies in the row space over F_p.
    This rules out every 2^s*e_j integer relation, regardless of s.
    """
    echelon = packed_rank(rows, p, retain_basis=True)
    width, basis = len(rows[0]), echelon.pop('echelon_basis')
    free = sorted(set(range(width))-basis.keys())
    rng, vectors, covered = random.Random(seed), [], set()
    for _ in range(8):
        vector = [0]*width
        for j in free:
            vector[j] = rng.randrange(p)
        for pivot in sorted(basis, reverse=True):
            vector[pivot] = -sum(a*b for a, b in zip(basis[pivot], vector)) % p
        assert all(sum(a*b for a, b in zip(row, vector)) % p == 0 for row in rows)
        vectors.append(vector)
        covered.update(j for j, value in enumerate(vector) if value)
        if len(covered) == width:
            break
    return {'prime': p, 'rank': echelon['rank'], 'columns': width,
            'seed': seed, 'kernel_vectors': vectors,
            'covered_coordinates': sorted(covered),
            'all_coordinates_excluded': len(covered) == width,
            'no_rational_coordinate_multiple': len(covered) == width and echelon['rank'] == len(rows),
            'independent_dot_products_checked': len(rows)*len(vectors)}


def rank_control(q, t, c, audit=False):
    field = LogField(q)
    cycle = field.cycle(t)
    assert len(cycle) == q+1 and c not in (0, t)
    point = PointDiagram(field, t, c, cycle)
    folded = from_point(point)
    logs = chord_logs(folded)
    oriented = {value for exponent in logs for value in (exponent, -exponent % (q-1))}
    assert len(oriented) == q//2 and 0 not in oriented
    assert gcd(q-1, *logs) == 1
    records, counts = family_rows(folded, audit=audit)
    rows = [r['row'] for r in records]
    tests = []
    for prime in prime_divisors(q-1):
        result = packed_rank(rows, prime)
        tests.append(result)
        if result['rank'] == q//4:
            break
    passed = tests[-1]['rank'] == q//4
    # Independently compute the actual candidate gain subgroup, without rank.
    root_logs = [sum(a*b for a, b in zip(row, logs)) % (q-1) for row in rows]
    root_subgroup_index = gcd(q-1, *root_logs)
    if passed:
        prime = tests[-1]['prime']
        assert any(exponent % prime for exponent in root_logs)
    result = {'q': q, 'modulus': field.modulus, 't': t, 'c': c,
              'canonical_gain': folded.root_path(folded.canonical)['gain'],
              'inverse_traces_c_partner_delta': [field.traces[field.inverse[v]]
                  for v in (c, c ^ t, folded.delta)],
              'rows': len(rows), 'columns': q//4, 'family_counts': counts,
              'prime_rank_tests': tests, 'one_prime_certificate': passed,
              'oriented_chord_gains': len(oriented), 'chord_gain_subgroup_index': 1,
              'root_gain_subgroup_index': root_subgroup_index}
    if not passed:
        result['coordinate_obstruction'] = kernel_coordinate_cover(rows, tests[0]['prime'])
        if result['coordinate_obstruction']['all_coordinates_excluded']:
            result['coordinate_obstruction']['signed_root_rows'] = [
                ''.join({-1: '-', 0: '.', 1: '+'}[value] for value in row) for row in rows]
            result['coordinate_obstruction']['row_sources'] = [
                {k: v for k, v in record.items() if k != 'row'} for record in records]
    print('RANK', q, t, c, 'rows', len(rows), 'tests', [(r['prime'], r['rank']) for r in tests],
          'passed', passed, 'root_subgroup_index', root_subgroup_index, flush=True)
    return result


def prefixed_control():
    """One ordinary Kempe prefix restores full odd-prime rank in the obstruction.

    This is not a claim that eight changes are necessary to colour the graph:
    the original seven-change family already has successful gains.
    """
    q, t, c, prime = 2048, 343, 1766, 23
    field = LogField(q)
    point = PointDiagram(field, t, c, field.cycle(t))
    folded = from_point(point)
    graph, natural = closed_graph(folded)
    component = graph.component(natural, (0, 1), 0)
    initial = natural.copy()
    change(initial, (0, 1), component)
    records, counts = family_rows(folded, initial=initial)
    rows = [r['row'] for r in records]
    rank = packed_rank(rows, prime)
    assert rank['rank'] == q//4
    # Independently audit the determinant of the retained 512-by-512 minor.
    minor = [rows[i] for i in rank['selected_row_indices']]
    scalar_rank, determinant = ordinary_rank(minor, prime)
    assert scalar_rank == q//4 and determinant == rank['minor_determinant_mod_prime']
    logs = chord_logs(folded)
    successful = next(i for i, row in enumerate(rows) if records[i]['source']['kind'] == 'seven'
                      and sum(a*b for a, b in zip(row, logs)) % prime)
    record, source = records[successful], records[successful]['source']
    moves = [((0, 1), component)]
    if source['kind'] == 'initial':
        target = initial
    else:
        switched = seven_switches(graph, initial, source['pivot_edge_id'])['candidates'][source['branch']]
        target = switched['colours']
        moves += switched['moves']
    replay = audit_buffer_sequence(graph, natural, target, moves, closed=True)
    selected = frozenset(e for e, colour in enumerate(target) if colour == record['excluded_colour'])
    partners = []
    for translation in (c, c ^ t):
        partner = PointDiagram(field, t, translation, field.cycle(t))
        certificate, profile = audit_matching(partner, folded, selected, permutations=True)
        assert certificate['colours']
        partners.append({'c': translation, 'certificate': certificate, 'point_profile': profile})
    result = {'q': q, 'modulus': field.modulus, 't': t, 'c': c, 'partner_c': c ^ t,
              'prefix': {'colours': [0, 1], 'seed_edge_id': 0, 'component_edge_ids': sorted(component)},
              'family_counts': counts, 'rows': len(rows), 'rank_certificate': rank,
              'gains_after_prefix': [item['certificate']['root_gain']
                                     for item in eligible_certificates(point, folded, initial)],
              'independent_minor_determinant_checked': True,
              'row_sources': [{k: v for k, v in r.items() if k != 'row'} for r in records],
              'selected_successful_row': successful, 'move_replay': replay,
              'partner_lift_audits': partners,
              'scope': 'specified rank-restoring prefix, not a necessary or uniform eight-change bound'}
    print('PREFIXED', 'component', len(component), 'rows', len(rows), 'rank', scalar_rank,
          'minor_determinant', determinant, 'moves_to_selected_success', replay['moves'], flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--json-output', type=Path)
    parser.add_argument('--larger', action='store_true', help='test the 20 stored q=2048,4096 cases')
    parser.add_argument('--case', type=int, nargs=3, action='append', metavar=('Q', 'T', 'C'))
    parser.add_argument('--augment-results', type=Path,
                        help='keep an existing rank scan and add the exact obstruction/prefix controls')
    args = parser.parse_args()
    output = (json.loads(args.augment_results.read_text()) if args.augment_results else
              {'scope': 'one-prime sufficient criterion; specified finite controls only', 'cases': []})
    output['arithmetic_controls'] = arithmetic_controls()
    parameters = args.case or [(256, 15, 151), (512, 47, 10), (512, 53, 66),
                               (1024, 49, 80), (1024, 49, 910), (1024, 62, 836), (1024, 287, 683)]
    if args.larger:
        stored = json.loads(Path(__file__).with_name('borel_exchange_stress_results.json').read_text())
        parameters += [(f['q'], r['t'], r['c']) for f in stored['fields'] for r in f['cases']]
    if not args.augment_results:
        for q, t, c in parameters:
            output['cases'].append(rank_control(q, t, c, audit=q <= 1024))
    if args.augment_results or args.larger:
        # Recheck and retain the negative certificate, not just a deficient rank.
        obstruction = rank_control(2048, 343, 1766, audit=True)
        assert obstruction['coordinate_obstruction']['no_rational_coordinate_multiple']
        output['cases'] = [r for r in output['cases'] if (r['q'], r['t'], r['c']) != (2048, 343, 1766)]
        output['cases'].append(obstruction)
        output['prefixed_control'] = prefixed_control()
    if args.json_output:
        args.json_output.write_text(json.dumps(output, indent=2) + '\n')
