#!/usr/bin/env python3
"""Exhaust the seven-transposition, degree-17 boundary case.

Write x=(0,1,...,16).  If the three fixed points of the involution contain
two consecutive points, Proposition 4.6 of the survey applies.  The remaining
442 independent triples form 16 orbits under the dihedral normalizer of <x>.
For one representative fixed triple from each orbit, this script enumerates
the perfect matchings of the other 14 points modulo the triple's stabilizer.

It first tries the 17-vertex point-stabilizer quotient.  Only failures are
sent to the 136-vertex two-set-stabilizer quotient.  Every positive SAT model
is checked directly by ``colour_pregraph`` in the shared certificate module.

The 16 independent blocks may be checked in parallel.  For a deliberately
serial run, use ``--workers 1``.
"""

import argparse
import multiprocessing as mp
from itertools import combinations
from time import perf_counter

from sparse_involution_cycle_experiment import (
    colour_pregraph,
    dihedral_image,
    perfect_matchings,
    point_pregraph,
    subset_pregraph,
)


N = 17


def set_image(points, sign, shift):
    return tuple(sorted((sign * point + shift) % N for point in points))


def canonical_fixed_set(points):
    return min(set_image(points, sign, shift)
               for sign in (1, -1) for shift in range(N))


def fixed_set_representatives():
    independent = []
    for points in combinations(range(N), 3):
        point_set = set(points)
        if all((point + 1) % N not in point_set for point in point_set):
            independent.append(points)
    assert len(independent) == 442
    representatives = sorted(set(map(canonical_fixed_set, independent)))
    assert len(representatives) == 16
    return representatives


def stabilizer(fixed):
    return [(sign, shift) for sign in (1, -1) for shift in range(N)
            if set_image(fixed, sign, shift) == fixed]


def check_block(fixed):
    """Return counts for one fixed-set orbit."""
    started = perf_counter()
    symmetries = stabilizer(fixed)
    nonidentity = [element for element in symmetries if element != (1, 0)]
    assert len(symmetries) in (1, 2)
    endpoints = tuple(point for point in range(N) if point not in fixed)

    total = point_failed = twoset_failed = 0
    failed_examples = []
    for matching in perfect_matchings(endpoints):
        # Each symmetric fixed triple has a single residual reflection.
        if nonidentity:
            sign, shift = nonidentity[0]
            if matching > dihedral_image(N, matching, sign, shift):
                continue
        total += 1
        if colour_pregraph(N, point_pregraph(N, matching)) is not None:
            continue
        point_failed += 1

        vertices, items = subset_pregraph(N, matching, 2)
        assert vertices == 136 and len(items) == 209
        assert sum(second is None for _, second in items) == 10
        if colour_pregraph(vertices, items) is None:
            twoset_failed += 1
            if len(failed_examples) < 10:
                failed_examples.append(matching)

    return {
        'fixed': fixed,
        'stabilizer': len(symmetries),
        'types': total,
        'point_failed': point_failed,
        'twoset_failed': twoset_failed,
        'failed_examples': failed_examples,
        'seconds': round(perf_counter() - started, 2),
    }


def check(workers):
    representatives = fixed_set_representatives()
    expected = {
        (0, 2, 4): (2, 68219, 0),
        (0, 2, 5): (1, 135135, 0),
        (0, 2, 6): (1, 135135, 0),
        (0, 2, 7): (1, 135135, 0),
        (0, 2, 8): (1, 135135, 0),
        (0, 2, 9): (1, 135135, 0),
        (0, 3, 6): (2, 68219, 1146),
        (0, 3, 7): (1, 135135, 427),
        (0, 3, 8): (1, 135135, 1381),
        (0, 3, 9): (1, 135135, 442),
        (0, 3, 10): (2, 68219, 634),
        (0, 4, 8): (2, 68219, 165),
        (0, 4, 9): (1, 135135, 330),
        (0, 4, 10): (1, 135135, 307),
        (0, 5, 10): (2, 68219, 544),
        (0, 5, 11): (2, 68219, 188),
    }
    assert set(representatives) == set(expected)
    started = perf_counter()
    results = []
    if workers == 1:
        iterator = map(check_block, representatives)
        pool = None
    else:
        pool = mp.Pool(workers)
        iterator = pool.imap_unordered(check_block, representatives)
    try:
        for result in iterator:
            results.append(result)
            print('BLOCK', result, flush=True)
    finally:
        if pool is not None:
            pool.close()
            pool.join()

    results.sort(key=lambda result: result['fixed'])
    assert all((result['stabilizer'], result['types'],
                result['point_failed']) == expected[result['fixed']]
               for result in results)
    totals = {
        key: sum(result[key] for result in results)
        for key in ('types', 'point_failed', 'twoset_failed')
    }
    assert totals == {
        'types': 1760664, 'point_failed': 5564, 'twoset_failed': 0,
    }, [
        example for result in results for example in result['failed_examples']]
    print('SEPTUPLE-TRANSPOSITION N 17',
          'independent-fixed-set-orbits', len(representatives),
          'dihedral-types', totals['types'],
          'point-colourable', totals['types'] - totals['point_failed'],
          'rescued-by-two-set', totals['point_failed'],
          'two-set-failed', totals['twoset_failed'],
          'seconds', round(perf_counter() - started, 2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--workers', type=int, default=4)
    arguments = parser.parse_args()
    if arguments.workers < 1:
        parser.error('--workers must be positive')
    check(arguments.workers)
