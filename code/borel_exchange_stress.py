"""Stress the bounded alternating-exchange existence target beyond q=1024.

Field arithmetic uses O(q) space. The candidate checker contracts the
unchanged path segments, so a size-r exchange is checked in O(r log r)
work independently of the number of quotient vertices. Both optimisations
are checked against the earlier independent implementations.
"""

import argparse
import json
from collections import Counter
from math import gcd
from pathlib import Path
from time import perf_counter

from borel_alternating_exchange import alternating_exchanges, inspect_exchange
from borel_folded_experiment import from_point
from borel_reflected_family import (BinaryField, irreducible, polynomial_product,
                                    polynomial_remainder, half_word_certificate)
from borel_reversal_experiment import PointDiagram


class ProductRow:
    def __init__(self, field, a):
        self.field, self.a = field, a

    def __getitem__(self, b):
        return self.field.mul(self.a, b)


class LogField(BinaryField):
    def __init__(self, q):
        assert q >= 4 and q & (q - 1) == 0
        self.q, self.degree = q, q.bit_length() - 1
        self.modulus = next(poly for poly in range(q + 1, 2 * q, 2) if irreducible(poly))

        def slow(a, b):
            return polynomial_remainder(polynomial_product(a, b), self.modulus)

        for generator in range(2, q):
            powers, seen, value = [], set(), 1
            while value not in seen:
                seen.add(value)
                powers.append(value)
                value = slow(value, generator)
            assert value == 1
            if len(powers) == q - 1:
                break
        assert len(powers) == q - 1
        self.generator = generator
        self.exp = powers + powers
        self.log = [-1] * q
        for i, v in enumerate(powers):
            self.log[v] = i
        self.inverse = {a: powers[(-self.log[a]) % (q - 1)] for a in range(1, q)}
        self.times = [ProductRow(self, a) for a in range(q)]
        self.traces = []
        for a in range(q):
            value, trace = a, 0
            for _ in range(self.degree):
                trace ^= value
                value = self.mul(value, value)
            assert trace in (0, 1)
            self.traces.append(trace)
        # Every element is checked against the independent polynomial product,
        # with a basis and several varying partners. Small fields are exhaustive.
        partners = list(range(q)) if q <= 32 else [1 << i for i in range(self.degree)]
        for a in range(q):
            for b in partners + [a, (a * 37 + 11) % q]:
                assert self.mul(a, b) == slow(a, b)
            if a:
                assert self.mul(a, self.inverse[a]) == 1

    def mul(self, a, b):
        return self.exp[self.log[a] + self.log[b]] if a and b else 0


def cut_certificate(folded, exchange):
    """Check the small reconnection graph on the cut path's endpoints."""
    cuts = sorted(exchange['inserted'])
    r = len(cuts)
    assert r == len(exchange['removed'])
    assert all(b > a + 1 for a, b in zip(cuts, cuts[1:]))
    slots = [0] + [v for i in cuts for v in (i, i + 1)] + [folded.h - 1]
    position = {slots[i]: i for i in range(1, 2 * r + 1)}
    assert len(position) == 2 * r
    partner = {}
    for e in exchange['removed']:
        u, v, kind = folded.edges[e]
        assert kind == 'a'
        partner[position[u]], partner[position[v]] = position[v], position[u]
    assert len(partner) == 2 * r
    current, gain, length, segments, chords = 0, 1, 0, [], []
    seen_segments = set()
    while True:
        segment = current // 2
        assert segment not in seen_segments
        seen_segments.add(segment)
        segments.append(segment)
        length += slots[2 * segment + 1] - slots[2 * segment] + 1
        current ^= 1
        if current == 2 * r + 1:
            break
        other = partner[current]
        u, v = slots[current], slots[other]
        chords.append((u, v))
        gain = folded.field.times[gain][folded.field.times[folded.labels[u]][
            folded.field.inverse[folded.labels[v]]]]
        current = other
    return {'spanning': len(segments) == r + 1, 'gain': gain,
            'root_vertices': length, 'root_segments': segments,
            'oriented_chords': chords, 'chords': r}


def canonical_gain(field, labels, position, delta):
    """Linear-time alternating traversal, independently of FoldedDiagram."""
    partner = [position[k ^ delta] for k in labels]
    current, gain, steps = 0, 1, 0
    while True:
        other = partner[current]
        gain = field.mul(gain, field.mul(labels[current], field.inverse[labels[other]]))
        steps += 1
        assert steps <= len(labels) // 2
        if other == len(labels) - 1:
            return gain
        current = other - 1 if other % 2 == 0 else other + 1


def small_controls():
    counts = Counter()
    for q in (4, 8, 16, 32):
        field, reference, seen = LogField(q), BinaryField(q), set()
        assert field.modulus == reference.modulus
        assert all(field.mul(a, b) == reference.times[a][b] for a in range(q) for b in range(q))
        for t in range(1, q):
            if t in seen:
                continue
            seen.update(field.frobenius_orbit(t))
            cycle = field.cycle(t)
            if len(cycle) != q + 1:
                continue
            for c in range(1, q):
                if c == t or c > (c ^ t):
                    continue
                folded = from_point(PointDiagram(field, t, c, cycle))
                fast = canonical_gain(field, folded.labels, {k: i for i, k in enumerate(folded.labels)}, folded.delta)
                assert fast == folded.root_path(folded.canonical)['gain']
                counts['canonical_diagrams'] += 1
                for exchange in alternating_exchanges(folded, 6):
                    small, direct = cut_certificate(folded, exchange), inspect_exchange(folded, exchange)
                    assert small['gain'] == direct['gain']
                    assert small['spanning'] == direct['spanning']
                    assert small['root_vertices'] == len(direct['path'])
                    counts['cut_diagrams'] += 1
    return dict(counts)


def check(q, maximum, trace_limit=0, stop_on_failure=False):
    start = perf_counter()
    field, seen, counts, failures = LogField(q), set(), Counter(), []
    print('FIELD', q, field.modulus, 'generator', field.generator, flush=True)
    records = []
    for t in range(1, q):
        if t in seen:
            continue
        orbit = field.frobenius_orbit(t)
        seen.update(orbit)
        cycle = field.cycle(t)
        if len(cycle) != q + 1:
            continue
        labels = from_point(PointDiagram(field, t, 1, cycle)).labels
        position = {k: i for i, k in enumerate(labels)}
        per_trace = Counter()
        for c in range(1, q):
            if c == t or c > (c ^ t):
                continue
            # Only this class is left by the proved parallel-edge trace test.
            if not field.traces[field.inverse[c]] or not field.traces[field.inverse[c ^ t]]:
                continue
            per_trace['folded_parameters'] += 1
            delta = field.mul(c, c ^ t)
            if canonical_gain(field, labels, position, delta) != 1:
                continue
            per_trace['canonical_failures'] += 1
            point = PointDiagram(field, t, c, cycle)
            folded = from_point(point)
            canonical = point.profile(set(range(1, q, 2)), independent=not failures)
            assert tuple(canonical[0]['word_matrix']) == (1, 0, 0, 1)
            assert folded.root_path(folded.canonical)['gain'] == 1
            histogram, witness = Counter(), None
            for exchange in alternating_exchanges(folded, maximum):
                candidate = cut_certificate(folded, exchange)
                kind = 'nonspanning' if not candidate['spanning'] else 'unit' if candidate['gain'] == 1 else 'good'
                histogram[f'{candidate["chords"]}_{kind}'] += 1
                if kind == 'good' and (witness is None or candidate['chords'] < witness['chords']):
                    lifted = folded.lift(point, exchange['matching'])
                    profile = point.profile(lifted)
                    assert len(profile) == 1 and profile[0]['length'] == q + 1
                    half = half_word_certificate(point, profile)
                    assert half['beta'] == field.mul(t, candidate['gain'] ^ 1) != 0
                    witness = {**candidate, 'cuts': sorted(i + 1 for i in exchange['inserted']),
                               'profile': profile, 'half_word': half}
            case = {'t': t, 'c': c, 'delta': delta, 'orbit_weight': 2 * len(orbit),
                    'canonical_profile': canonical, 'histogram': dict(histogram), 'witness': witness}
            failures.append(case)
            if witness is not None:
                per_trace['bounded_repairs'] += 1
            else:
                print('BOUNDED-EXCHANGE-FAILURE', q, t, c, dict(histogram), flush=True)
                if stop_on_failure:
                    return {'q': q, 'modulus': field.modulus, 'maximum': maximum,
                            'status': 'SPANNING_BOUND_COUNTEREXAMPLE',
                            'scope': 'partial field scan, stopped at the reported spanning-exchange obstruction',
                            'target': 'one simple alternating-cycle exchange with at most maximum chords, spanning complement and non-unit gain',
                            'cases': failures,
                            'seconds': perf_counter() - start}
        counts.update(per_trace)
        records.append({'t': t, 'orbit_weight': len(orbit), 'counts': dict(per_trace)})
        print('TRACE', q, t, dict(per_trace), 'total', dict(counts), flush=True)
        if trace_limit and len(records) >= trace_limit:
            break
    status = 'TRACE_LIMIT' if trace_limit else 'COMPLETE'
    if not trace_limit:
        trace_count = sum(gcd(j, q + 1) == 1 for j in range(1, q + 1)) // 2
        assert sum(r['orbit_weight'] for r in records) == trace_count
    return {'q': q, 'modulus': field.modulus, 'maximum': maximum, 'status': status,
            'scope': 'double-inverse-trace-one class, modulo Frobenius and c <-> c+t',
            'counts': dict(counts), 'traces': records, 'cases': failures,
            'seconds': perf_counter() - start}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('q', type=int, nargs='*', default=[2048])
    parser.add_argument('--maximum', type=int, default=5)
    parser.add_argument('--trace-limit', type=int, default=0)
    parser.add_argument('--stop-on-failure', action='store_true')
    parser.add_argument('--json-output', type=Path)
    args = parser.parse_args()
    output = {'small_controls': small_controls(), 'fields': []}
    print('SMALL-CONTROLS', output['small_controls'], flush=True)
    for q in args.q:
        output['fields'].append(check(q, args.maximum, args.trace_limit, args.stop_on_failure))
        if args.json_output:
            args.json_output.write_text(json.dumps(output, indent=2) + '\n')
