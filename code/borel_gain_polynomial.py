"""Polynomial nonvanishing from the projective D_i recurrence.

Root gains equal one exactly when products of the corresponding D_i agree.
This supplies polynomial gcd/Bezout certificates independent of root-row
rank. No uniform certificate-existence theorem is asserted.
"""

import argparse
import json
import random
from functools import lru_cache
from pathlib import Path
from types import SimpleNamespace

from borel_exchange_stress import LogField
from borel_folded_experiment import from_point
from borel_gain_lattice import chord_logs, family_rows, root_row
from borel_kempe_exchange import audit_buffer_sequence, closed_graph
from borel_kempe_closure import all_three_colourings
from borel_reversal_experiment import PointDiagram
from borel_three_exchange import audit_matching
from borel_two_buffer import seven_switches


def multiply(a, b):
    result = 0
    while b:
        bit = b & -b
        result ^= a << (bit.bit_length()-1)
        b ^= bit
    return result


def coefficient_product(a, b):
    """Independent, deliberately simple coefficient-array multiplication."""
    coefficients = [0]*(a.bit_length()+b.bit_length())
    for i in range(a.bit_length()):
        if (a >> i) & 1:
            for j in range(b.bit_length()):
                if (b >> j) & 1:
                    coefficients[i+j] ^= 1
    return sum(value << i for i, value in enumerate(coefficients))


def arithmetic_controls():
    rng = random.Random(907)
    for _ in range(128):
        a, b = rng.getrandbits(24), rng.getrandbits(16) or 1
        assert multiply(a, b) == coefficient_product(a, b)
        quotient, residual = divide(a, b)
        assert coefficient_product(quotient, b) ^ residual == a
        assert residual.bit_length() < b.bit_length()
        common, u, v = bezout(a, b)
        assert coefficient_product(u, a) ^ coefficient_product(v, b) == common
    return {'independent_multiplication_division_bezout_triples': 128}


def divide(a, b):
    assert b
    quotient = 0
    while a and a.bit_length() >= b.bit_length():
        shift = a.bit_length()-b.bit_length()
        quotient ^= 1 << shift
        a ^= b << shift
    return quotient, a


def remainder(a, b):
    return divide(a, b)[1]


def gcd(a, b):
    while b:
        a, b = b, remainder(a, b)
    return a


def bezout(a, b):
    original = a, b
    s, ss, t, tt = 1, 0, 0, 1
    while b:
        quotient, residual = divide(a, b)
        a, b = b, residual
        s, ss = ss, s ^ multiply(quotient, ss)
        t, tt = tt, t ^ multiply(quotient, tt)
    assert multiply(s, original[0]) ^ multiply(t, original[1]) == a
    return a, s, t


def evaluate(poly, field, value):
    result = 0
    for i in range(poly.bit_length()-1, -1, -1):
        result = field.mul(result, value) ^ ((poly >> i) & 1)
    return result


def full_trace_polynomial(field):
    """Squarefree product over ALL full-cycle traces, with binary coefficients."""
    seen, factors, product = set(), [], 1
    for t in range(1, field.q):
        if t in seen:
            continue
        orbit = field.frobenius_orbit(t)
        seen.update(orbit)
        if len(field.cycle(t)) != field.q+1:
            continue
        coefficients = [1]
        for value in orbit:
            target = [0]*(len(coefficients)+1)
            for i, coefficient in enumerate(coefficients):
                target[i] ^= field.mul(coefficient, value)
                target[i+1] ^= coefficient
            coefficients = target
        assert all(v in (0, 1) for v in coefficients)
        factor = sum(v << i for i, v in enumerate(coefficients))
        assert gcd(product, factor) == 1
        assert all(evaluate(factor, field, value) == 0 for value in orbit)
        factors.append({'representative': t, 'orbit': sorted(orbit), 'polynomial': hex(factor)})
        product = multiply(product, factor)
    assert product == primitive_trace_polynomial(field.q+1)
    return product, factors


def recurrence_polynomials(h, modulus=None):
    values = [0, 1]
    for _ in range(1, h):
        value = (values[-1] << 1) ^ values[-2]
        values.append(value if modulus is None else remainder(value, modulus))
    return values


@lru_cache(None)
def primitive_trace_polynomial(n):
    """Independent construction: sqrt(D_n)=product_{d|n,d>1} Phi_trace_d."""
    assert n >= 3 and n % 2
    previous, current = 0, 1
    for _ in range(1, n):
        previous, current = current, (current << 1) ^ previous
    root = 0
    for i in range(current.bit_length()):
        if (current >> i) & 1:
            assert i % 2 == 0
            root |= 1 << (i//2)
    assert multiply(root, root) == current
    for d in range(3, n, 2):
        if n % d == 0:
            root, residual = divide(root, primitive_trace_polynomial(d))
            assert residual == 0
    return root


def translation_polynomial(first, second, values, modulus=None):
    """Equality of (D_i+D_j)/(D_i D_j) for two chord pairs (zero-based)."""
    a, b = (values[i+1] for i in first)
    c, d = (values[i+1] for i in second)
    left = multiply(a ^ b, multiply(c, d))
    right = multiply(c ^ d, multiply(a, b))
    return left ^ right if modulus is None else remainder(left ^ right, modulus)


def translation_control(folded, modulus, factors, values):
    """Exact realizability in the full-cycle family, verified by field recurrence."""
    pairs = [edge[:2] for edge in folded.edges if edge[2] == 'a']
    common, reducers = modulus, []
    polynomials = [translation_polynomial(pairs[0], pair, values, modulus) for pair in pairs[1:]]
    for pair, polynomial in zip(pairs[1:], polynomials):
        result, u, v = bezout(common, polynomial)
        if result != common:
            reducers.append({'compared_pair_zero_based': pair, 'polynomial': hex(polynomial),
                'previous_gcd': hex(common), 'next_gcd': hex(result),
                'bezout_previous': hex(u), 'bezout_constraint': hex(v)})
            common = result
    assert all(remainder(poly, common) == 0 for poly in polynomials)
    field, realized = folded.field, []
    for factor in factors:
        t = factor['representative']
        previous, current, labels = 0, 1, []
        for _ in range(folded.h):
            assert current
            labels.append(field.mul(field.inverse[current], field.inverse[current]))
            previous, current = current, field.mul(t, current) ^ previous
        deltas = {labels[a] ^ labels[b] for a, b in pairs}
        admissible = len(deltas) == 1
        assert (evaluate(common, field, t) == 0) == admissible
        if not admissible:
            continue
        delta = next(iter(deltas))
        translations = [c for c in range(1, field.q) if field.mul(c, c ^ t) == delta]
        assert len(translations) == 2
        for c in translations:
            reconstructed = from_point(PointDiagram(field, t, c, field.cycle(t)))
            assert reconstructed.labels == labels and reconstructed.edges == folded.edges
        realized.append({'trace_representative': t, 'trace_orbit': factor['orbit'],
                         'translations_at_representative': translations})
    degree = common.bit_length()-1
    assert sum(len(r['trace_orbit']) for r in realized) == degree
    return {'admissible_trace_polynomial': hex(common), 'degree': degree,
            'normalized_generating_pairs': 2*degree, 'reducing_constraints': reducers,
            'independently_reconstructed_orbits': realized,
            'full_trace_orbits_checked_independently': len(factors)}


def forbidden_anchored_pattern():
    """Three anchored chords incompatible with the recurrence over ANY char-2 field."""
    pairs = [(0, 2), (1, 6), (3, 5)]
    values = recurrence_polynomials(7)
    first = translation_polynomial(pairs[0], pairs[1], values)
    second = translation_polynomial(pairs[0], pairs[2], values)
    common, u, v = bezout(first, second)
    assert common == 1
    assert (first, second, u, v) == (0x397, 0x4c2, 0x43, 0x3c)
    assert coefficient_product(u, first) ^ coefficient_product(v, second) == 1
    return {'chord_pairs_one_based': [[a+1, b+1] for a, b in pairs],
            'exact_translation_polynomials': [hex(first), hex(second)],
            'bezout_coefficients': [hex(u), hex(v)], 'bezout_result': '0x1',
            'scope': 'anchored indices 1,...,7 with all D_i nonzero, any characteristic-two field'}


def perfect_pairings(vertices):
    if not vertices:
        yield []
        return
    a = vertices[0]
    for i, b in enumerate(vertices[1:], 1):
        for rest in perfect_pairings(vertices[1:i]+vertices[i+1:]):
            yield [(a, b)]+rest


def recurrence_only_controls():
    """Exhaust arbitrary pairings while KEEPING the true projective label order.

    A common chord translation is deliberately not imposed. These are
    negative controls for that relaxation, not claimed Cayley diagrams.
    """
    results = []
    for q in (8, 16):
        field, seen, tested, failures = LogField(q), set(), 0, []
        for t in range(1, q):
            if t in seen:
                continue
            seen.update(field.frobenius_orbit(t))
            cycle = field.cycle(t)
            if len(cycle) != q+1:
                continue
            point = PointDiagram(field, t, next(c for c in range(1, q) if c != t), cycle)
            labels = from_point(point).labels
            for matching in perfect_pairings(list(range(q//2))):
                diagram = SimpleNamespace(field=field, h=q//2, labels=labels,
                    edges=[(i, i+1, 'x') for i in range(q//2-1)] +
                          [(a, b, 'a') for a, b in matching])
                graph, _ = closed_graph(diagram)
                states = all_three_colourings(graph)
                rows = {tuple(root_row(diagram, graph, state, excluded))
                        for state in states for excluded in {0, 1, 2}-{state[-1]}}
                logs = chord_logs(diagram)
                tested += 1
                if any(sum(a*b for a, b in zip(row, logs)) % (q-1) for row in rows):
                    continue
                assert len(states) == 1
                deltas = [labels[a] ^ labels[b] for a, b in matching]
                assert len(set(deltas)) > 1
                modulus = primitive_trace_polynomial(q+1)
                values = recurrence_polynomials(q//2, modulus)
                nonzero_rows = sorted(row for row in rows if any(row))
                common, root_polynomials = modulus, []
                for row in nonzero_rows:
                    a, b = root_products(diagram, row, values, modulus)
                    polynomial = a ^ b
                    assert evaluate(polynomial, field, t) == 0
                    common = gcd(common, polynomial)
                    root_polynomials.append(hex(polynomial))
                translation_constraints = [translation_polynomial(matching[0], pair, values, modulus)
                                           for pair in matching[1:]]
                joint = common
                for polynomial in translation_constraints:
                    joint = gcd(joint, polynomial)
                assert joint == 1
                failures.append({'t': t, 'labels': labels, 'chord_pairs_zero_based': matching,
                    'chord_label_differences': deltas, 'all_colourings_mod_global_names': len(states),
                    'colouring': list(next(iter(states))), 'nonzero_root_rows': nonzero_rows,
                    'full_trace_polynomial': hex(modulus), 'root_polynomials_modulus': root_polynomials,
                    'common_root_polynomial': hex(common),
                    'translation_constraints_modulus': [hex(p) for p in translation_constraints],
                    'joint_gcd_with_translation_constraints': hex(joint)})
        print('RECURRENCE_ONLY', q, 'ordered_diagrams', tested, 'failures', len(failures), flush=True)
        results.append({'q': q, 'modulus': field.modulus, 'ordered_diagrams_tested': tested,
                        'failures': failures})
    return results


def root_products(folded, row, values, modulus):
    """A=product of source D_i; B=product of destination D_i; gamma=(B/A)^2."""
    numerator, denominator = 1, 1
    for coordinate, sign in enumerate(row):
        if not sign:
            continue
        u, v, _ = folded.edges[folded.h-1+coordinate]
        if sign < 0:
            u, v = v, u
        numerator = remainder(multiply(numerator, values[u+1]), modulus)
        denominator = remainder(multiply(denominator, values[v+1]), modulus)
    return numerator, denominator


def certificate_audit(folded, certificates, factors, t, c):
    """Independent coefficient identities, all trace orbits, and actual lifts."""
    field = folded.field
    if not certificates:
        return {'independent_coefficient_bezout_checks': 0,
                'independent_field_recurrence_checks': 0, 'lifts': []}
    for step in certificates:
        old = int(step['previous_gcd'], 16)
        polynomial = int(step['unit_gain_polynomial_modulus'], 16)
        u, v = int(step['bezout_previous'], 16), int(step['bezout_candidate'], 16)
        assert coefficient_product(u, old) ^ coefficient_product(v, polynomial) == int(step['next_gcd'], 16)
    full_trace_checks = 0
    for factor in factors:
        trace = factor['representative']
        previous, current, values = 0, 1, [0]
        for _ in range(folded.h):
            values.append(current)
            previous, current = current, field.mul(trace, current) ^ previous
        nonunit = False
        for step in certificates:
            source, destination = 1, 1
            for j, sign in enumerate(step['record']['row']):
                if sign:
                    u, v, _ = folded.edges[folded.h-1+j]
                    if sign < 0:
                        u, v = v, u
                    source = field.mul(source, values[u+1])
                    destination = field.mul(destination, values[v+1])
            assert source and destination
            assert source ^ destination == evaluate(int(step['unit_gain_polynomial_modulus'], 16), field, trace)
            nonunit |= source != destination
            full_trace_checks += 1
        if certificates[-1]['next_gcd'] == '0x1':
            assert nonunit
    graph, natural = closed_graph(folded)
    lifts = []
    for step in certificates:
        record = step['record']
        source = record['source']
        if source['kind'] == 'natural':
            target, moves = natural, []
        else:
            candidate = seven_switches(graph, natural, source['pivot_edge_id'])['candidates'][source['branch']]
            target, moves = candidate['colours'], candidate['moves']
        replay = audit_buffer_sequence(graph, natural, target, moves, closed=True)
        selected = frozenset(e for e, colour in enumerate(target) if colour == record['excluded_colour'])
        partners = []
        for translation in (c, c ^ t):
            point = PointDiagram(field, t, translation, field.cycle(t))
            certificate, profile = audit_matching(point, folded, selected, permutations=True)
            partners.append({'c': translation, 'certificate': certificate, 'point_profile': profile})
        lifts.append({'record_index': step['record_index'], 'move_replay': replay, 'partners': partners})
    if certificates[-1]['next_gcd'] == '0x1':
        assert any(item['partners'][0]['certificate']['colours'] for item in lifts)
    return {'independent_coefficient_bezout_checks': len(certificates),
            'independent_field_recurrence_checks': full_trace_checks, 'lifts': lifts}


def polynomial_control(q, t, c):
    field = LogField(q)
    folded = from_point(PointDiagram(field, t, c, field.cycle(t)))
    modulus, factors = full_trace_polynomial(field)
    values = recurrence_polynomials(folded.h, modulus)
    assert all(gcd(value, modulus) == 1 for value in values[1:])
    records, counts = family_rows(folded, audit=True)
    ordered = sorted(enumerate(records), key=lambda item: (sum(a != 0 for a in item[1]['row']), item[0]))
    current, certificates = modulus, []
    for index, record in ordered:
        a, b = root_products(folded, record['row'], values, modulus)
        polynomial = a ^ b
        common, u, v = bezout(current, polynomial)
        if common == current:
            continue
        certificates.append({'record_index': index, 'record': record,
            'source_product_modulus': hex(a), 'destination_product_modulus': hex(b),
            'unit_gain_polynomial_modulus': hex(polynomial),
            'previous_gcd': hex(current), 'next_gcd': hex(common),
            'bezout_previous': hex(u), 'bezout_candidate': hex(v)})
        current = common
        if current == 1:
            break
    actual_gains = []
    for step in certificates:
        a = evaluate(int(step['source_product_modulus'], 16), field, t)
        b = evaluate(int(step['destination_product_modulus'], 16), field, t)
        ratio = field.mul(b, field.inverse[a])
        actual_gains.append(field.mul(ratio, ratio))
    print('POLYNOMIAL', q, t, c, 'trace_degree', modulus.bit_length()-1,
          'rows', len(records), 'cert_rows', len(certificates),
          'gcd_degree', current.bit_length()-1, 'actual_gains', actual_gains, flush=True)
    return {'q': q, 'modulus': field.modulus, 't': t, 'c': c,
            'full_trace_polynomial': hex(modulus), 'trace_factors': factors,
            'translation_realizability': translation_control(folded, modulus, factors, values),
            'independent_certificate_audit': certificate_audit(folded, certificates, factors, t, c),
            'root_family_counts': counts, 'root_rows': len(records),
            'certificates': certificates, 'final_gcd': hex(current),
            'certified_for_every_full_trace_on_this_fixed_ordered_diagram': current == 1,
            'actual_gains_of_certificate_rows': actual_gains}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', type=int, nargs=3, action='append')
    parser.add_argument('--json-output', type=Path)
    args = parser.parse_args()
    output = {'scope': 'specified fixed ordered diagrams, not a uniform existence theorem',
              'arithmetic_controls': arithmetic_controls(),
              'forbidden_anchored_pattern': forbidden_anchored_pattern(),
              'recurrence_only_controls': recurrence_only_controls(),
              'cases': [polynomial_control(q, t, c) for q, t, c in
                  (args.case or [(64, 7, 26), (128, 8, 20), (256, 15, 151),
                                 (512, 47, 10), (512, 53, 66), (2048, 343, 1766)])]}
    if args.json_output:
        args.json_output.write_text(json.dumps(output, indent=2)+'\n')
