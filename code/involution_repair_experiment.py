#!/usr/bin/env python3
"""Check involution exchanges and a barrier to monotone symmetry growth.

python3 code/involution_repair_experiment.py --exhaustive-a5
python3 code/involution_repair_experiment.py --groups S5 A6 S6 --samples 512

The first command independently enumerates all locally valid A5 supports.
Only the explicit 60-vertex barrier and the completed enumerations are
exhaustive results. No general decreasing-move theorem is assumed.
The universal single-circuit rule, including this paired restriction, is
refuted by the PSL(2,11) certificate in translated_block_repair.py.
"""

import argparse
import json
import random
from collections import Counter

from cayley_snark_check import inv, mul, order
from cdc_palette_experiment import cycle_perm
from translated_repair_audit import (
    RepairGraph, bits, cases, independent_supports, sample_matchings,
)


def translate(mask, action):
    return sum(1 << action[e] for e in bits(mask))


def exchange_structure(graph, matching, g):
    """Check the half-turn lemma and return fixed circuits and paired circuits."""
    t, action = graph.group[g], graph.actions[g]
    assert order(t) == 2
    assert all(mul(t, h) != h for h in graph.group)
    assert all(action[action[e]] == e for e in range(len(graph.edges)))
    translated = translate(matching, action)
    circuits = set(graph.difference_circuits(matching, action))
    assert sum(circuits) == matching ^ translated  # disjoint edge masks
    fixed, pairs, seen = [], [], set()
    for circuit in sorted(circuits):
        if circuit in seen:
            continue
        other = translate(circuit, action)
        assert other in circuits
        if other == circuit:
            assert circuit.bit_count() % 4 == 2
            fixed.append(circuit)
        else:
            assert not circuit & other
            pairs.append((circuit, other))
            for d in (circuit, other):
                candidate = matching ^ d
                assert candidate ^ translate(candidate, action) == (
                    (matching ^ translated) ^ circuit ^ other)
                assert (candidate ^ translate(candidate, action)).bit_count() == (
                    (matching ^ translated).bit_count() - 2 * d.bit_count())
        seen.update((circuit, other))
    if not fixed:
        # Construct one of the 2^len(pairs) invariant mixtures. This need not
        # satisfy the additional no-all-a-pentagon condition.
        mixture = matching
        for circuit, _ in pairs:
            mixture ^= circuit
        assert translate(mixture, action) == mixture
        assert all(sum(bool(mixture & (1 << e)) for e in incident) == 1
                   for incident in graph.inc)
    return fixed, pairs


def check_mixture_count(graph, matching, g):
    """Independently enumerate every component choice in a small union."""
    action = graph.actions[g]
    circuits = list(graph.difference_circuits(matching, action))
    assert len(circuits) <= 12
    fixed, pairs = exchange_structure(graph, matching, g)
    invariant_count = 0
    for choices in range(1 << len(circuits)):
        candidate = matching
        for i in bits(choices):
            candidate ^= circuits[i]
        assert all(sum(bool(candidate & (1 << e)) for e in incident) == 1
                   for incident in graph.inc)
        invariant_count += translate(candidate, action) == candidate
    assert invariant_count == (0 if fixed else 1 << len(pairs))
    return invariant_count


def explicit_barrier():
    """A short permutation certificate; no SAT solver or matching search."""
    graph = RepairGraph(*cases()['A5'])
    # H fixes 4. The orbit coordinate g[4] transforms by the right generator.
    # On the five-vertex quotient choose a(0,1), the a-semi-edge at 4,
    # and x(2,3). Lift this matching to all 60 vertices.
    matching = 0
    for e, (u, v) in enumerate(graph.edges):
        i, j = graph.group[u][4], graph.group[v][4]
        if ((graph.kinds[e] == 'a' and i in (0, 1, 4)) or
                (graph.kinds[e] == 'x' and {i, j} == {2, 3})):
            matching |= 1 << e
    graph.check(matching)
    assert len(graph.selected(matching)) == 18
    assert graph.factor_lengths(matching) == [15] * 4
    a, x = cases()['A5'][0]
    monodromy = mul(mul(mul(mul(x, x), a), x), x)
    assert order(monodromy) == 3
    stabilizer = [g for g, action in enumerate(graph.actions)
                  if translate(matching, action) == matching]
    assert {graph.group[g] for g in stabilizer} == {
        h for h in graph.group if h[4] == 4}
    assert len(stabilizer) == 12
    t = mul(cycle_perm(5, [0, 4]), cycle_perm(5, [1, 2]))
    g = graph.index[t]
    fixed, pairs = exchange_structure(graph, matching, g)
    assert [d.bit_count() for d in fixed] == [10]
    assert sorted(d.bit_count() for d, _ in pairs) == [6, 10]
    assert check_mixture_count(graph, matching, g) == 0
    outcomes = Counter()
    for d in graph.difference_circuits(matching, graph.actions[g]):
        candidate = matching ^ d
        outcomes[d.bit_count(), graph.valid(candidate),
                 tuple(graph.factor_lengths(candidate))] += 1
        if graph.valid(candidate):
            graph.check_flip(matching, g, d, candidate)
            assert graph.matching_stabilizer_order(candidate) == 3
    assert outcomes == Counter({
        (6, True, (15, 45)): 1,
        (6, False, (5, 5, 5, 15, 15, 15)): 1,
        (10, False, (5, 15, 15, 25)): 3,
    })
    # H-conjugacy reduces every involution outside H to this representative.
    conjugates = set()
    for h in (graph.group[i] for i in stabilizer):
        conjugates.add(mul(mul(inv(h), t), h))
    outside = {h for h in graph.group if order(h) == 2 and h[4] != 4}
    assert conjugates == outside and len(outside) == 12
    # Every different translate of M is H-conjugate to tM. Thus the table
    # covers arbitrary translations, not just involutions.
    different = {translate(matching, a) for a in graph.actions} - {matching}
    from_conjugates = {translate(matching, graph.actions[graph.index[h]])
                       for h in conjugates}
    assert different == from_conjugates and len(different) == 4
    # Directly check the orbit reduction as well.
    _, all_outcomes = graph.audit(matching, exhaustive=True)
    assert all_outcomes == {'invalid': 16, '2': 4}
    # A separate positive control checks the exact count, not just a witness.
    for d in graph.difference_circuits(matching, graph.actions[g]):
        candidate = matching ^ d
        if not graph.valid(candidate):
            continue
        for j, h in enumerate(graph.group):
            if order(h) != 2:
                continue
            fixed2, pairs2 = exchange_structure(graph, candidate, j)
            if not fixed2 and pairs2:
                assert check_mixture_count(graph, candidate, j) == 1 << len(pairs2)
                break
        else:
            raise AssertionError('missing positive invariant-mixture control')
        break
    print('EXPLICIT-SYMMETRY-BARRIER', json.dumps({
        'fingerprint': graph.fingerprint, 'matching': hex(matching),
        'translation': t, 'stabilizer_order': 12,
        'factor_lengths': [15] * 4, 'fixed_difference_lengths': [10],
        'paired_difference_lengths': [6, 10],
        'valid_flip_factor_lengths': [15, 45],
        'valid_flip_stabilizer_order': 3,
        'all_distinct_flip_outcomes': all_outcomes,
    }), flush=True)


def audit_group(name, count, exhaustive):
    graph = RepairGraph(*cases()[name])
    involutions = [g for g, h in enumerate(graph.group) if order(h) == 2]
    status, counts, barriers = {}, Counter(), []
    if exhaustive:
        assert name in ('A5', 'A5_alt')
        starts = independent_supports(graph)
    else:
        starts = sample_matchings(graph, count, random.Random(509), False, status=status)
    for matching in starts:
        counts['starts'] += 1
        odd = sum(n % 2 for n in graph.factor_lengths(matching))
        if not odd:
            counts['already_even'] += 1
            continue
        counts['obstructed'] += 1
        witness, can_impose_new_involution = None, False
        for g in involutions:
            fixed, pairs = exchange_structure(graph, matching, g)
            if not fixed and pairs:
                can_impose_new_involution = True
            if witness is not None:
                continue  # still verify the lemma for every involution
            for pair in pairs:
                for d in pair:
                    candidate = matching ^ d
                    if not graph.valid(candidate):
                        continue
                    new_odd = sum(n % 2 for n in graph.factor_lengths(candidate))
                    if new_odd < odd:
                        graph.check_flip(matching, g, d, candidate)
                        witness = (g, d)
                        break
                if witness is not None:
                    break
        if witness is None:
            counts['paired_rule_stalled'] += 1
            # This would refute only the stronger involution-paired rule.
            unrestricted, outcomes = graph.audit(matching)
            print('PAIRED-RULE-STALL', name, hex(matching), odd,
                  'unrestricted_witness', unrestricted, 'outcomes', outcomes, flush=True)
        else:
            counts['paired_decrease'] += 1
        if not can_impose_new_involution:
            counts['cannot_impose_new_involution'] += 1
            barriers.append(matching)
    completed = exhaustive or status.get('exhausted', False)
    if exhaustive:
        expected = {'A5': (5875, 2400, 10), 'A5_alt': (4625, 1690, 0)}[name]
        assert (counts['starts'], counts['obstructed'], len(barriers)) == expected
        assert counts['paired_decrease'] == counts['obstructed']
        for matching in barriers:
            assert graph.factor_lengths(matching) == [15] * 4
            assert graph.matching_stabilizer_order(matching) == 12
    print('INVOLUTION-AUDIT', name, graph.fingerprint, dict(counts),
          'exhausted', completed, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--groups', nargs='+', choices=list(cases()), default=['A5', 'A5_alt'])
    parser.add_argument('--samples', type=int, default=512)
    parser.add_argument('--exhaustive-a5', action='store_true',
                        help='use the non-SAT complete local-pattern enumerator on A5 cases')
    args = parser.parse_args()
    if args.samples < 1:
        parser.error('--samples must be positive')
    if args.exhaustive_a5 and any(g not in ('A5', 'A5_alt') for g in args.groups):
        parser.error('--exhaustive-a5 requires only A5 or A5_alt')
    explicit_barrier()
    for name in args.groups:
        audit_group(name, args.samples, args.exhaustive_a5)


if __name__ == '__main__':
    main()
