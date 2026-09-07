"""Exact diagnostics for root-span search and unrestricted temporary colours.

These are finite certificates, not a theorem that a recolouring closure
has full rank. The colour-state graph must not be quotiented by root span.
"""

import argparse
import json
import random
from itertools import combinations
from pathlib import Path

from borel_exchange_stress import LogField
from borel_folded_experiment import from_point
from borel_gain_lattice import chord_logs, family_rows, root_row
from borel_gain_rank import PackedPrime, ordinary_rank, packed_rank, prime_divisors
from borel_kempe_exchange import audit_buffer_sequence, change, closed_graph
from borel_reversal_experiment import PointDiagram
from borel_three_exchange import audit_matching
from borel_two_buffer import seven_switches


def kernel_signature_certificate(rows, p=23, seed=7123, count=8):
    """Certify minimum support >=3 over F_p and Q, when row ranks agree.

    Nonzero, pairwise projectively distinct kernel-column signatures
    exclude all nonzero vectors supported on at most two coordinates.
    If the input rows are independent mod p, appending a primitive
    integer vector of support <=2 raises rank both mod p and over Q.
    """
    result = packed_rank(rows, p, retain_basis=True)
    basis, width = result.pop('echelon_basis'), len(rows[0])
    assert result['rank'] == len(rows)
    assert ordinary_rank(rows, p)[0] == len(rows)
    free = sorted(set(range(width))-basis.keys())
    rng, vectors = random.Random(seed), []
    for _ in range(count):
        vector = [0]*width
        for j in free:
            vector[j] = rng.randrange(p)
        for pivot in sorted(basis, reverse=True):
            vector[pivot] = -sum(a*b for a, b in zip(basis[pivot], vector)) % p
        assert all(sum(a*b for a, b in zip(row, vector)) % p == 0 for row in rows)
        vectors.append(vector)
    signatures = [tuple(v[j] for v in vectors) for j in range(width)]
    assert all(any(s) for s in signatures)
    projective = [tuple(a*pow(next(v for v in s if v), -1, p) % p for a in s)
                  for s in signatures]
    assert len(set(projective)) == width
    minimum = min(sum(a != 0 for a in row) for row in rows)
    assert minimum == 3
    return {'prime': p, 'rank': result['rank'], 'columns': width,
            'independent_scalar_rank_check': True,
            'seed': seed, 'kernel_vectors': vectors,
            'independent_dot_products_checked': len(rows)*count,
            'nonzero_pairwise_projectively_distinct_signatures': width,
            'minimum_nonzero_support_mod_prime_and_over_Q': 3}


def span_expressions(rows, targets, p):
    """Membership certificates, independently verified by direct dot products."""
    width, height = len(rows[0]), len(rows)
    arithmetic, tracking = PackedPrime(p, width), PackedPrime(p, height)
    basis = {}
    for index, source in enumerate(rows):
        row = arithmetic.pack(source)
        coeff = tracking.pack([int(j == index) for j in range(height)])
        while row:
            pivot = ((row & -row).bit_length()-1) // arithmetic.bits
            value = (row >> (pivot*arithmetic.bits)) & (arithmetic.base-1)
            if pivot not in basis:
                inverse = pow(value, -1, p)
                basis[pivot] = (arithmetic.scale(row, inverse), tracking.scale(coeff, inverse))
                break
            b, t = basis[pivot]
            row = arithmetic.subtract(row, arithmetic.scale(b, value))
            coeff = tracking.subtract(coeff, tracking.scale(t, value))
    certificates = []
    for target in targets:
        row, coeff = arithmetic.pack(target), 0
        while row:
            pivot = ((row & -row).bit_length()-1) // arithmetic.bits
            value = (row >> (pivot*arithmetic.bits)) & (arithmetic.base-1)
            assert pivot in basis
            b, t = basis[pivot]
            row = arithmetic.subtract(row, arithmetic.scale(b, value))
            coeff = tracking.add(coeff, tracking.scale(t, value))
        expression = tracking.unpack(coeff)
        assert all((sum(a*r[j] for a, r in zip(expression, rows))-target[j]) % p == 0
                   for j in range(width))
        certificates.append(expression)
    return certificates


def all_temporary_families(folded):
    """Every natural pivot and both choices of temporary ordinary colour."""
    graph, natural = closed_graph(folded)
    records, seen, counts = [], set(), []
    for a in range(3):
        permutation = list(range(3))
        permutation[0], permutation[a] = permutation[a], permutation[0]
        initial = [permutation[colour] for colour in natural]
        family, family_counts = family_rows(folded, initial=initial, audit=True)
        for record in family:
            row = tuple(record['row'])
            if row not in seen:
                seen.add(row)
                records.append({'temporary_natural_colour': a, **record})
        rank = packed_rank([r['row'] for r in records], prime_divisors(folded.field.q-1)[0])
        counts.append({'temporary_natural_colour': a, 'family_rows': len(family),
                       'family_counts': family_counts, 'accumulated_rows': len(records),
                       'accumulated_rank_at_first_prime': rank['rank']})
    return records, counts


def unrestricted_control(q, t, c):
    field = LogField(q)
    folded = from_point(PointDiagram(field, t, c, field.cycle(t)))
    records, counts = all_temporary_families(folded)
    rows = [r['row'] for r in records]
    tests = []
    for p in prime_divisors(q-1):
        rank = packed_rank(rows, p)
        tests.append(rank)
        if rank['rank'] == q//4:
            minor = [rows[i] for i in rank['selected_row_indices']]
            independent_rank, determinant = ordinary_rank(minor, p)
            assert independent_rank == q//4
            assert determinant == rank['minor_determinant_mod_prime']
            rank['independent_scalar_minor_check'] = True
            break
    print('UNRESTRICTED', q, t, c, 'rows', len(rows),
          'ranks', [(r['prime'], r['rank']) for r in tests], flush=True)
    return {'q': q, 'modulus': field.modulus, 't': t, 'c': c,
            'temporary_colour_counts': counts, 'rows': len(rows), 'rank_tests': tests,
            'row_sources': [{k: v for k, v in r.items() if k != 'row'} for r in records]}


def obstruction_diagnostics():
    q, t, c, p = 2048, 343, 1766, 23
    field = LogField(q)
    folded = from_point(PointDiagram(field, t, c, field.cycle(t)))
    graph, natural = closed_graph(folded)
    records, counts = family_rows(folded, audit=True)
    rows = [r['row'] for r in records]
    assert len(rows) == 425
    signatures = kernel_signature_certificate(rows)
    neighbours, targets = [], []
    for pair in combinations(range(3), 2):
        for component in graph.components(natural, pair):
            target = natural.copy()
            change(target, pair, component)
            excluded = sorted({0, 1, 2}-{target[-1]})
            targets.extend(root_row(folded, graph, target, a) for a in excluded)
            neighbours.append({'pair': pair, 'seed_edge_id': min(component),
                               'length': len(component), 'excluded_colours': excluded})
    expressions = span_expressions(rows, targets, p)
    assert len(neighbours) == 9 and len(expressions) == 18
    for i, neighbour in enumerate(neighbours):
        neighbour['root_span_coefficients_mod_23'] = expressions[2*i:2*i+2]
    logs = chord_logs(folded)
    witness = next(r for r in records if sum(a != 0 for a in r['row']) == 3
                   and sum(a*b for a, b in zip(r['row'], logs)) % (q-1))
    source = witness['source']
    switched = seven_switches(graph, natural, source['pivot_edge_id'])['candidates'][source['branch']]
    replay = audit_buffer_sequence(graph, natural, switched['colours'], switched['moves'], closed=True)
    selected = frozenset(e for e, a in enumerate(switched['colours']) if a == witness['excluded_colour'])
    lift_audits = []
    for translation in (c, c ^ t):
        partner = PointDiagram(field, t, translation, field.cycle(t))
        certificate, profile = audit_matching(partner, folded, selected, permutations=True)
        assert certificate['colours']
        lift_audits.append({'c': translation, 'certificate': certificate, 'point_profile': profile})
    starts, ends = [], []
    for j, value in enumerate(witness['row']):
        if value:
            u, v, _ = folded.edges[folded.h-1+j]
            starts.append(folded.labels[u if value == 1 else v])
            ends.append(folded.labels[v if value == 1 else u])
    numerator = field.mul(field.mul(starts[0], starts[1]), starts[2])
    denominator = field.mul(field.mul(ends[0], ends[1]), ends[2])
    assert numerator != denominator
    print('SUPPORT', 'minimum', 3, 'witness gain', field.mul(numerator, field.inverse[denominator]),
          'neutral ordinary neighbours', len(neighbours), flush=True)
    return {'q': q, 'modulus': field.modulus, 't': t, 'c': c,
            'fixed_family_counts': counts, 'fixed_family_signed_rows': [
                ''.join({-1: '-', 0: '.', 1: '+'}[a] for a in row) for row in rows],
            'fixed_family_row_sources': [{k: v for k, v in r.items() if k != 'row'} for r in records],
            'minimum_support_certificate': signatures, 'ordinary_neighbours': neighbours,
            'rank_plateau_scope': 'all ordinary neighbours relative to the restricted colour-0 seven-family',
            'three_chord_witness': {'record': witness, 'starts': starts, 'ends': ends,
                'numerator': numerator, 'denominator': denominator, 'move_replay': replay,
                'partner_lift_audits': lift_audits}}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--json-output', type=Path)
    args = parser.parse_args()
    output = {'scope': 'finite support, span-membership and all-temporary-colour certificates',
              'obstruction': obstruction_diagnostics(),
              'unrestricted_controls': [unrestricted_control(q, t, c) for q, t, c in
                  [(1024, 49, 80), (1024, 287, 683), (2048, 343, 1766)]]}
    if args.json_output:
        args.json_output.write_text(json.dumps(output, indent=2)+'\n')
