"""Stress the reflected-two-chord target on larger characteristic-two fields.

The candidate family is fixed before computing any circuit word. We test all
c for one representative of each Frobenius orbit of full-cycle traces t;
field automorphisms preserve the indexed diagram, matching, and word orders.
These are parameter tests, not generating-pair orbits. No generation claim
is inferred from a successful matching or from a failed candidate family.

python3 code/borel_reflected_family.py 64 128 256 512 --json-output code/borel_reflected_results.json
"""

import argparse
import json
from collections import Counter
from math import gcd
from pathlib import Path
from time import perf_counter

from borel_reversal_experiment import Field, IDENTITY, PointDiagram
from cayley_snark_check import GF


def polynomial_remainder(value, divisor):
    while value.bit_length() >= divisor.bit_length():
        value ^= divisor << (value.bit_length() - divisor.bit_length())
    return value


def polynomial_product(a, b):
    result = 0
    while b:
        if b & 1:
            result ^= a
        a <<= 1
        b >>= 1
    return result


def polynomial_gcd(a, b):
    while b:
        a, b = b, polynomial_remainder(a, b)
    return a


def irreducible(poly):
    degree = poly.bit_length() - 1
    power = 2
    for k in range(1, degree + 1):
        power = polynomial_remainder(polynomial_product(power, power), poly)
        if k <= degree // 2 and polynomial_gcd(power ^ 2, poly) != 1:
            return False
    return power == 2


class BinaryField(Field):
    def __init__(self, q):
        assert q >= 4 and q & (q - 1) == 0
        self.q = q
        self.degree = q.bit_length() - 1
        if q in GF.IRRED:
            self.modulus = sum(bit << i for i, bit in enumerate(GF.IRRED[q][1]))
        else:
            self.modulus = next(poly for poly in range(q + 1, 2*q, 2)
                                if irreducible(poly))
        assert irreducible(self.modulus)
        self.times = [[0]*q for _ in range(q)]
        for a in range(q):
            for b in range(q):
                result, u, v = 0, a, b
                while v:
                    if v & 1:
                        result ^= u
                    v >>= 1
                    u <<= 1
                    if u & q:
                        u ^= self.modulus
                self.times[a][b] = result
        self.inverse = {a: self.times[a].index(1) for a in range(1, q)}
        # A second multiplication algorithm checks the entire table.
        assert all(self.times[a][b] == polynomial_remainder(
                   polynomial_product(a, b), self.modulus)
                   for a in range(q) for b in range(q))
        if q <= 64:
            original = GF(q)
            assert all(self.times[a][b] == original.mulf(a, b)
                       for a in range(q) for b in range(q))

    def frobenius_orbit(self, t):
        orbit = []
        while t not in orbit:
            orbit.append(t)
            t = self.times[t][t]
        assert t == orbit[0]
        return orbit

    def cycle(self, t):
        cycle, seen, z = [], set(), self.q
        while z not in seen:
            seen.add(z)
            cycle.append(z)
            z = self.act((0, 1, 1, t), z)
        assert z == self.q
        return cycle


def generated_order(field, t, c):
    """Exact two-level Schreier calculation using affine point stabilisers.

    Powers of x are transversals for infinity. The resulting stabiliser
    generators are affine. Orbit--stabiliser at 0 then leaves a subgroup
    of the cyclic multiplicative group, whose order is an lcm of orders.
    No enumeration of the (potentially billion-element) group is needed.
    """
    q, x, a = field.q, (0, 1, 1, t), (1, c, 0, 1)
    cycle = field.cycle(t)
    assert len(cycle) == q + 1
    position = {z: i for i, z in enumerate(cycle)}
    powers, current = [], IDENTITY
    for _ in cycle:
        powers.append(current)
        current = field.multiply(x, current)
    assert current == IDENTITY
    affine = set()
    for i, z in enumerate(cycle):
        j = position[field.act(a, z)]
        r = powers[j]
        r_inverse = (r[3], r[1], r[2], r[0])
        word = field.multiply(r_inverse, field.multiply(a, powers[i]))
        assert word[2] == 0
        divisor = field.inverse[word[3]]
        affine.add((field.times[word[0]][divisor], field.times[word[1]][divisor]))
    affine = sorted(affine - {(1, 0)})
    # Store an affine representative z -> slopes[v]*z+v for every orbit point.
    slopes, orbit = {0: 1}, [0]
    for v in orbit:
        for lam, beta in affine:
            w = field.times[lam][v] ^ beta
            if w not in slopes:
                slopes[w] = field.times[lam][slopes[v]]
                orbit.append(w)
        if len(orbit) == q:
            break

    primes, remaining, divisor = [], q - 1, 2
    while divisor * divisor <= remaining:
        if remaining % divisor == 0:
            primes.append(divisor)
            while remaining % divisor == 0:
                remaining //= divisor
        divisor += 1
    if remaining > 1:
        primes.append(remaining)

    def power(base, exponent):
        result = 1
        while exponent:
            if exponent & 1:
                result = field.times[result][base]
            base = field.times[base][base]
            exponent >>= 1
        return result

    def scalar_order(value):
        n = q - 1
        for prime in primes:
            while n % prime == 0 and power(value, n // prime) == 1:
                n //= prime
        assert power(value, n) == 1
        return n

    stabiliser = 1
    for v in orbit:
        for lam, beta in affine:
            w = field.times[lam][v] ^ beta
            scalar = field.times[field.inverse[slopes[w]]][field.times[lam][slopes[v]]]
            n = scalar_order(scalar)
            stabiliser = stabiliser * n // gcd(stabiliser, n)
            if stabiliser == q - 1:
                return (q + 1) * len(orbit) * stabiliser
    return (q + 1) * len(orbit) * stabiliser


def half_word_certificate(diagram, profile):
    """Check beta=t+d^(-2), with d the lower-left half-path entry.

    For c != t, the reversal fixes no chord. Its unique fixed ordinary
    edge is the middle x-edge, which cannot belong to an invariant matching.
    This helper is intended only for the reversal-invariant matchings used
    in this experiment.
    """
    assert diagram.c != diagram.t
    field = diagram.field
    labels = profile[0]['labels']
    length = len(labels)
    assert length % 2
    halfway = length // 2
    matrices = {'a': diagram.a, 'x': diagram.x, 'X': diagram.x_inverse}
    v = IDENTITY
    for label in labels[:halfway]:
        v = field.multiply(matrices[label], v)
    center = labels[halfway]
    assert center in ('x', 'X')
    assert field.act(v, field.q) == (diagram.t ^ 1 if center == 'x' else 1)
    d = v[2]
    assert d != 0
    beta = diagram.t ^ field.times[field.inverse[d]][field.inverse[d]]
    reflection = (1, diagram.t, 0, 1)
    v_inverse = (v[3], v[1], v[2], v[0])
    predicted = field.multiply(reflection, field.multiply(v_inverse,
                field.multiply(reflection, field.multiply(matrices[center], v))))
    assert predicted == (1, beta, 0, 1) == tuple(profile[0]['word_matrix'])
    return {'half_length': halfway, 'half_word': v, 'center_label': center,
            'denominator': d, 'beta': beta}


def root_vertices(diagram, profile):
    vertices, current = set(), 0
    for symbol in profile[0]['labels']:
        vertices.add(current)
        if symbol == 'a':
            current = diagram.position[diagram.cycle[current] ^ diagram.c]
        else:
            current = (current + (1 if symbol == 'x' else -1)) % diagram.n
    assert current == 0
    return vertices


def classify(diagram, profile):
    odd = [p for p in profile if p['length'] % 2]
    assert profile[0]['length'] % 2
    w = profile[0]['word_matrix']
    assert w[0] == w[3] == 1 and w[2] == 0
    if diagram.c != diagram.t:
        half_word_certificate(diagram, profile)
    if diagram.good(profile):
        return 'good'
    if w == IDENTITY:
        return 'root_identity'
    assert len(odd) >= 3
    return 'other_odd_obstruction'


def check(q, trace_limit=0):
    started = perf_counter()
    field = BinaryField(q)
    seen, totals, weighted, failures, traces = set(), Counter(), Counter(), [], []
    for t in range(1, q):
        if t in seen:
            continue
        orbit = field.frobenius_orbit(t)
        seen.update(orbit)
        cycle = field.cycle(t)
        if len(cycle) != q + 1:
            continue
        # Explicitly check that Frobenius preserves cycle indices.
        assert field.cycle(field.times[t][t]) == [
            q if z == q else field.times[z][z] for z in cycle]
        traces.append({'t': t, 'orbit': orbit})
        for c in range(1, q):
            diagram = PointDiagram(field, t, c, cycle)
            canonical = set(range(1, q, 2))
            profile = diagram.profile(canonical, independent=(c == 1))
            assert sum(p['length'] % 2 for p in profile) == 1
            status = classify(diagram, profile)
            totals['parameters_tested'] += 1
            weighted['parameters_covered'] += len(orbit)
            if status == 'good':
                totals['canonical_good'] += 1
                weighted['canonical_good'] += len(orbit)
                continue
            assert status == 'root_identity'
            diagram.profile(canonical, independent=True)
            totals['canonical_identity'] += 1
            weighted['canonical_identity'] += len(orbit)
            histogram, odd_counts, candidates, witness = Counter(), Counter(), [], None
            strong_good = 0
            for selected in diagram.reflected_two_chord_matchings():
                candidate_profile = diagram.profile(selected)
                candidate_status = classify(diagram, candidate_profile)
                histogram[candidate_status] += 1
                odd_counts[sum(p['length'] % 2 for p in candidate_profile)] += 1
                chord_endpoints = {v for e in selected if e >= diagram.n
                                   for v in diagram.edges[e]}
                endpoints_on_root = chord_endpoints <= root_vertices(diagram, candidate_profile)
                if endpoints_on_root:
                    assert sum(p['length'] % 2 for p in candidate_profile) == 1
                    strong_good += candidate_status == 'good'
                candidates.append({'chords': [diagram.edges[e] for e in sorted(selected)
                                              if e >= diagram.n],
                                   'status': candidate_status, 'profile': candidate_profile,
                                   'all_endpoints_on_root': endpoints_on_root})
                if candidate_status == 'good' and witness is None:
                    diagram.profile(selected, independent=True)
                    witness = {'matching': sorted(selected), 'profile': candidate_profile,
                               'half_word_certificate': half_word_certificate(diagram, candidate_profile)}
            status = 'two_chord_repaired' if witness is not None else 'two_chord_failed'
            totals[status] += 1
            weighted[status] += len(orbit)
            failure = {'t': t, 'c': c, 'orbit_weight': len(orbit),
                       'canonical_profile': profile, 'candidate_histogram': dict(histogram),
                       'canonical_half_word': half_word_certificate(diagram, profile),
                       'odd_count_histogram': dict(odd_counts), 'witness': witness,
                       'good_with_all_endpoints_on_root': strong_good}
            if witness is None:
                # Retain every rejected candidate, not just the absent witness.
                failure['all_candidates'] = candidates
            failures.append(failure)
            print('REPAIR', q, t, c, dict(histogram), 'odd_counts', dict(odd_counts),
                  'strong_good', strong_good, flush=True)
        print('TRACE', q, t, 'totals', dict(totals), flush=True)
        if trace_limit and len(traces) >= trace_limit:
            break
    if not trace_limit:
        full_trace_count = sum(gcd(i, q + 1) == 1 for i in range(1, q + 1)) // 2
        assert sum(len(trace['orbit']) for trace in traces) == full_trace_count
        assert weighted['parameters_covered'] == (q - 1) * full_trace_count
    result = {'q': q, 'modulus': field.modulus,
              'status': 'TRACE_LIMIT' if trace_limit else 'EXHAUSTIVE_PARAMETERS',
              'traces': traces, 'tested': dict(totals), 'covered': dict(weighted),
              'canonical_failures': failures}
    print('SUMMARY', q, dict(totals), 'covered', dict(weighted),
          'seconds', round(perf_counter() - started, 3), flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('q', nargs='*', type=int, default=[64, 128, 256, 512])
    parser.add_argument('--trace-limit', type=int, default=0)
    parser.add_argument('--json-output', type=Path)
    args = parser.parse_args()
    results = []
    for q in args.q:
        results.append(check(q, args.trace_limit))
        if args.json_output:
            args.json_output.write_text(json.dumps(results, indent=2) + '\n')
