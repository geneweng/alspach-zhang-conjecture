#!/usr/bin/env python3
"""Exhaust the remaining seven-transposition boundary degrees 19--27.

Write x=(0,1,...,N-1), where N is one of 19,21,23,25,27.  Configurations
in which two fixed points of the involution are consecutive are covered by
Proposition 4.6 of the survey.  For every remaining independent fixed set,
this script reduces by the dihedral normalizer of <x>, then enumerates every
perfect matching of the other fourteen points modulo the residual stabilizer.

The fast point-quotient test enumerates proper colourings of the Hamilton
cycle whose missing colours on the fourteen chord endpoints use at most two
values.  It stores and directly checks a cycle-colouring witness for every
projected pattern.  A chord matching has only 2^6 endpoint assignments after
global colour complementation.  Any matching not covered by this restricted
family is checked by CaDiCaL on the full three-colour point quotient; genuine
point failures are sent to the two-set quotient.  Every SAT witness is checked
by the shared direct pregraph verifier.
"""

import argparse
import multiprocessing as mp
from itertools import combinations
from math import comb
from time import perf_counter

from sparse_involution_cycle_experiment import (
    colour_pregraph,
    dihedral_image,
    perfect_matchings,
    point_pregraph,
    subset_pregraph,
)


SUPPORT_SIZE = 14
DEGREES = (19, 21, 23, 25, 27)
INDEPENDENT_COUNTS = {19: 2717, 21: 5148, 23: 3289, 25: 650, 27: 27}
FIXED_ORBIT_COUNTS = {19: 79, 21: 133, 23: 79, 25: 16, 27: 1}
STABILIZER_COUNTS = {
    19: {1: 64, 2: 15},
    21: {1: 113, 2: 19, 14: 1},
    23: {1: 64, 2: 15},
    25: {1: 10, 2: 6},
    27: {2: 1},
}
TYPE_COUNTS = {
    19: 9671925,
    21: 16576723,
    23: 9671925,
    25: 1760664,
    27: 68219,
}
EXPECTED_TOTALS = {
    19: (9671473, 412, 40),
    21: (16576721, 2, 0),
    23: (9671925, 0, 0),
    25: (1760664, 0, 0),
    27: (68219, 0, 0),
}


def set_image(n, points, sign, shift):
    return tuple(sorted((sign * point + shift) % n for point in points))


def canonical_fixed_set(n, points):
    return min(set_image(n, points, sign, shift)
               for sign in (1, -1) for shift in range(n))


def fixed_set_representatives(n):
    fixed_count = n - SUPPORT_SIZE
    independent = []
    for points in combinations(range(n), fixed_count):
        point_set = set(points)
        if all((point + 1) % n not in point_set for point in point_set):
            independent.append(points)
    assert len(independent) == INDEPENDENT_COUNTS[n]
    representatives = sorted(set(canonical_fixed_set(n, points)
                                 for points in independent))
    assert len(representatives) == FIXED_ORBIT_COUNTS[n]
    return representatives


def stabilizer(n, fixed):
    return [(sign, shift) for sign in (1, -1) for shift in range(n)
            if set_image(n, fixed, sign, shift) == fixed]


def decode(code):
    values = [0] * SUPPORT_SIZE
    for index in range(SUPPORT_SIZE - 1, -1, -1):
        values[index] = code % 3
        code //= 3
    assert code == 0
    return tuple(values)


def binary_cycle_colour_patterns(n, fixed):
    """Return two-valued missing-colour patterns with checked witnesses.

    Cycle edge i joins i to i+1.  Colour symmetry fixes edges n-1 and 0 to
    colours 0 and 1.  A state consists of the projected base-three word, its
    set of used missing colours, and the outgoing cycle-edge colour.  States
    with the same data have identical continuations, so one packed witness
    per state is enough.
    """
    support = set(range(n)) - set(fixed)
    initial_code = 2 if 0 in support else 0
    initial_seen = 1 << 2 if 0 in support else 0
    states = {(initial_code, initial_seen, 1): 1}
    for vertex in range(1, n - 1):
        new_states = {}
        for (code, seen, incoming), witness in states.items():
            for outgoing in range(3):
                if outgoing == incoming:
                    continue
                missing = 3 - incoming - outgoing
                if vertex in support:
                    new_seen = seen | (1 << missing)
                    if new_seen.bit_count() > 2:
                        continue
                    new_code = 3 * code + missing
                else:
                    new_seen, new_code = seen, code
                key = new_code, new_seen, outgoing
                new_states.setdefault(
                    key, witness | (outgoing << (2 * vertex)))
        states = new_states

    witnesses = {}
    for (code, seen, incoming), witness in states.items():
        if incoming == 0:
            continue
        missing = 3 - incoming
        if n - 1 in support:
            final_seen = seen | (1 << missing)
            if final_seen.bit_count() > 2:
                continue
            final_code = 3 * code + missing
        else:
            final_seen, final_code = seen, code
        witnesses.setdefault(final_code, witness)

    patterns = []
    for code in sorted(witnesses):
        pattern = decode(code)
        assert len(set(pattern)) <= 2
        edges = tuple((witnesses[code] >> (2 * index)) & 3
                      for index in range(n))
        assert edges[n - 1] == 0 and edges[0] == 1
        assert all(edges[vertex - 1] != edges[vertex]
                   for vertex in range(n))
        observed = tuple(3 - edges[vertex - 1] - edges[vertex]
                         for vertex in sorted(support))
        assert pattern == observed
        patterns.append(pattern)
    return patterns


def binary_partitions(patterns):
    """Map unlabelled endpoint bipartitions to witnessed patterns."""
    full = (1 << SUPPORT_SIZE) - 1
    result = {}
    for pattern in patterns:
        first_value = pattern[0]
        mask = sum((value != first_value) << index
                   for index, value in enumerate(pattern))
        result.setdefault(min(mask, full ^ mask), pattern)
    return result


def binary_witness(matching, position, partitions):
    """Find a witnessed bipartition constant on every chord."""
    full = (1 << SUPPORT_SIZE) - 1
    edge_masks = [(1 << position[first]) | (1 << position[second])
                  for first, second in matching]
    # Put the first matched pair in class zero; complementation loses nothing.
    candidates = [0]
    for edge_mask in edge_masks[1:]:
        candidates += [candidate | edge_mask for candidate in candidates]
    for candidate in candidates:
        key = min(candidate, full ^ candidate)
        if key not in partitions:
            continue
        pattern = partitions[key]
        assert all(pattern[position[first]] == pattern[position[second]]
                   for first, second in matching)
        return pattern
    return None


def check_block(task):
    n, fixed = task
    started = perf_counter()
    symmetries = stabilizer(n, fixed)
    support = tuple(point for point in range(n) if point not in fixed)
    position = {point: index for index, point in enumerate(support)}
    patterns = binary_cycle_colour_patterns(n, fixed)
    partitions = binary_partitions(patterns)

    total = binary_colourable = three_colour_rescue = 0
    point_failed = twoset_failed = 0
    failed_examples = []
    for matching in perfect_matchings(support):
        if len(symmetries) > 1 and matching != min(
                dihedral_image(n, matching, sign, shift)
                for sign, shift in symmetries):
            continue
        total += 1
        if binary_witness(matching, position, partitions) is not None:
            binary_colourable += 1
            continue

        # This complete three-colour test also directly verifies its model.
        if colour_pregraph(n, point_pregraph(n, matching)) is not None:
            three_colour_rescue += 1
            continue

        point_failed += 1
        vertices, items = subset_pregraph(n, matching, 2)
        semiedges = 7 + comb(n - SUPPORT_SIZE, 2)
        assert vertices == comb(n, 2)
        assert len(items) == (3 * vertices + semiedges) // 2
        assert sum(second is None for _, second in items) == semiedges
        if colour_pregraph(vertices, items) is None:
            twoset_failed += 1
            if len(failed_examples) < 10:
                failed_examples.append(matching)

    return {
        'degree': n,
        'fixed': fixed,
        'stabilizer': len(symmetries),
        'binary_patterns': len(patterns),
        'types': total,
        'binary_colourable': binary_colourable,
        'three_colour_rescue': three_colour_rescue,
        'point_failed': point_failed,
        'twoset_failed': twoset_failed,
        'failed_examples': failed_examples,
        'seconds': round(perf_counter() - started, 2),
    }


def check(degrees, workers):
    tasks = [(n, fixed) for n in degrees
             for fixed in fixed_set_representatives(n)]
    started = perf_counter()
    results = []
    if workers == 1:
        iterator = map(check_block, tasks)
        pool = None
    else:
        pool = mp.Pool(workers)
        iterator = pool.imap_unordered(check_block, tasks)
    try:
        for result in iterator:
            results.append(result)
            print('BLOCK', result, flush=True)
    finally:
        if pool is not None:
            pool.close()
            pool.join()

    for n in degrees:
        degree_results = [result for result in results
                          if result['degree'] == n]
        observed_stabilizers = {}
        for result in degree_results:
            size = result['stabilizer']
            observed_stabilizers[size] = observed_stabilizers.get(size, 0) + 1
        assert observed_stabilizers == STABILIZER_COUNTS[n]
        totals = {
            key: sum(result[key] for result in degree_results)
            for key in ('types', 'binary_colourable', 'three_colour_rescue',
                        'point_failed', 'twoset_failed')
        }
        assert totals['types'] == TYPE_COUNTS[n]
        assert (totals['binary_colourable'] + totals['three_colour_rescue']
                + totals['point_failed'] == totals['types'])
        assert (totals['binary_colourable'], totals['three_colour_rescue'],
                totals['point_failed']) == EXPECTED_TOTALS[n]
        assert totals['twoset_failed'] == 0, [
            example for result in degree_results
            for example in result['failed_examples']]
        print('SEPTUPLE-TRANSPOSITION N', n,
              'independent-fixed-set-orbits', FIXED_ORBIT_COUNTS[n],
              'dihedral-types', totals['types'],
              'binary-pattern-colourable', totals['binary_colourable'],
              'three-colour-point-rescue', totals['three_colour_rescue'],
              'rescued-by-two-set', totals['point_failed'],
              'two-set-failed', totals['twoset_failed'], flush=True)
    print('TOTAL-SECONDS', round(perf_counter() - started, 2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--degrees', type=int, nargs='+', default=list(DEGREES))
    arguments = parser.parse_args()
    if arguments.workers < 1:
        parser.error('--workers must be positive')
    if any(degree not in DEGREES for degree in arguments.degrees):
        parser.error('degrees must be chosen from 19, 21, 23, 25, 27')
    check(tuple(dict.fromkeys(arguments.degrees)), arguments.workers)
