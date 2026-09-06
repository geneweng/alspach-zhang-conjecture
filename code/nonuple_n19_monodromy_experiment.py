#!/usr/bin/env python3
"""Exhaust the nine-transposition, degree-19 boundary by monodromy.

Write x=(0,1,...,18), and let a fix 0 and pair the other eighteen
points.  In the point-stabiliser quotient the a-dart at 0 is the unique
semi-edge.  Every quotient perfect matching must select it, and its
complement is therefore a 2-factor on nineteen vertices.

For each of the 17!! labelled choices of a, this script enumerates quotient
perfect matchings until it finds one for which every odd complementary
circuit has even-order monodromy.  The matching-monodromy lemma in the
survey then lifts it to a Tait colouring of the Cayley graph.  Before
computing a word, the checker uses the stronger elementary certificate that
an odd circuit containing an odd number of a-edges has odd sign and hence
even order.  Word orders are computed only for the remaining circuits.

The seventeen blocks prescribe the mate of point 1 and can be checked in
parallel.  A single block contains 15!! = 2,027,025 chord matchings.  Use
``--first-mate`` for a quick reproducible partial run; omitting it performs
the theorem's full exhaustive calculation.
"""

import argparse
import multiprocessing as mp
from collections import Counter
from math import prod
from time import perf_counter

from sparse_involution_cycle_experiment import perfect_matchings


N = 19
CHORD_OFFSET = N
FIRST_MATES = tuple(range(2, N))
MATCHINGS_PER_BLOCK = prod(range(1, 17, 2))
TOTAL_MATCHINGS = prod(range(1, N, 2))


def has_even_order(permutation):
    """Whether a permutation has even order, from its cycle lengths."""
    seen = 0
    for start in range(N):
        if seen & (1 << start):
            continue
        current = start
        length = 0
        while not seen & (1 << current):
            seen |= 1 << current
            current = permutation[current]
            length += 1
        if length % 2 == 0:
            return True
    return False


class PointQuotient:
    """The 19-cycle, its nine chords, and the selected semi-edge at 0."""

    def __init__(self, matching):
        self.partner = list(range(N))
        self.chord_at = [-1] * N
        self.ends = [(edge, (edge + 1) % N) for edge in range(N)]
        for index, (first, second) in enumerate(matching):
            edge = CHORD_OFFSET + index
            self.partner[first] = second
            self.partner[second] = first
            self.chord_at[first] = self.chord_at[second] = edge
            self.ends.append((first, second))
        assert self.chord_at[0] == -1
        assert all(self.chord_at[vertex] >= CHORD_OFFSET
                   for vertex in range(1, N))

    def incident(self, vertex):
        """Ordinary quotient items at a vertex; the semi-edge is selected."""
        cycle_items = ((vertex - 1) % N, vertex)
        if vertex == 0:
            return cycle_items
        return cycle_items + (self.chord_at[vertex],)

    def target(self, edge, vertex):
        first, second = self.ends[edge]
        return second if first == vertex else first

    def word_has_even_order(self, start, selected):
        """Compute the monodromy around one specified complement circuit."""
        word = list(range(N))
        current, previous = start, -1
        while True:
            edge = next(item for item in self.incident(current)
                        if item != previous and not selected & (1 << item))
            target = self.target(edge, current)
            if edge >= CHORD_OFFSET:
                word = [self.partner[value] for value in word]
            elif target == (current + 1) % N:
                word = [(value + 1) % N for value in word]
            else:
                word = [(value - 1) % N for value in word]
            current, previous = target, edge
            if current == start:
                break
        return has_even_order(word)

    def certificate_kind(self, selected):
        """Return ``parity``, ``order``, or None for one quotient matching."""
        unseen = (1 << N) - 1
        used_order = False
        while unseen:
            start = (unseen & -unseen).bit_length() - 1
            current, previous = start, -1
            length = chord_count = 0
            while True:
                unseen &= ~(1 << current)
                edge = next(item for item in self.incident(current)
                            if item != previous and not selected & (1 << item))
                chord_count += edge >= CHORD_OFFSET
                length += 1
                current, previous = self.target(edge, current), edge
                if current == start:
                    break
            if length % 2 == 0 or chord_count % 2 == 1:
                continue
            used_order = True
            if not self.word_has_even_order(start, selected):
                return None
        return 'order' if used_order else 'parity'

    def find_certificate(self):
        """Depth-first enumeration of all quotient perfect matchings.

        The unique semi-edge at 0 is implicit and forced.  The first leaf is
        the canonical matching of cycle edges (1,2),(3,4),...,(17,18).
        Parallel cycle and chord items remain distinct choices.
        """
        attempts = 0

        def search(uncovered, selected):
            nonlocal attempts
            if not uncovered:
                attempts += 1
                kind = self.certificate_kind(selected)
                return (kind, attempts) if kind is not None else None

            first = (uncovered & -uncovered).bit_length() - 1
            for edge in self.incident(first):
                second = self.target(edge, first)
                if not uncovered & (1 << second):
                    continue
                result = search(
                    uncovered ^ (1 << first) ^ (1 << second),
                    selected | (1 << edge))
                if result is not None:
                    return result
            return None

        result = search(((1 << N) - 1) ^ 1, 0)
        return result if result is not None else ('failed', attempts)


def check_examples():
    """Fixed controls requiring both adaptation and actual word orders."""
    adaptive = (
        (1, 7), (2, 8), (3, 17), (4, 5), (6, 9),
        (10, 12), (11, 18), (13, 15), (14, 16),
    )
    canonical_order = (
        (1, 3), (2, 4), (5, 14), (6, 8), (7, 17),
        (9, 16), (10, 12), (11, 13), (15, 18),
    )
    assert PointQuotient(adaptive).find_certificate() == ('order', 2)
    assert PointQuotient(canonical_order).find_certificate() == ('order', 1)


def check_block(first_mate):
    started = perf_counter()
    remaining = tuple(vertex for vertex in range(2, N)
                      if vertex != first_mate)
    counts = Counter()
    attempt_sum = max_attempts = 0
    failures = []
    for rest in perfect_matchings(remaining):
        matching = ((1, first_mate),) + rest
        kind, attempts = PointQuotient(matching).find_certificate()
        position = 'canonical' if attempts == 1 else 'adaptive'
        counts[f'{position}_{kind}'] += 1
        attempt_sum += attempts
        max_attempts = max(max_attempts, attempts)
        if kind == 'failed' and len(failures) < 5:
            failures.append(matching)

    assert sum(counts.values()) == MATCHINGS_PER_BLOCK
    return {
        'first_mate': first_mate,
        'types': sum(counts.values()),
        'counts': dict(sorted(counts.items())),
        'attempt_sum': attempt_sum,
        'max_attempts': max_attempts,
        'failures': failures,
        'seconds': round(perf_counter() - started, 2),
    }


def check(workers, first_mates):
    check_examples()
    started = perf_counter()
    results = []
    if workers == 1:
        iterator = map(check_block, first_mates)
        pool = None
    else:
        pool = mp.Pool(workers)
        iterator = pool.imap_unordered(check_block, first_mates)
    try:
        for result in iterator:
            results.append(result)
            print('BLOCK', result, flush=True)
    finally:
        if pool is not None:
            pool.close()
            pool.join()

    totals = Counter()
    for result in results:
        totals.update(result['counts'])
    checked = sum(totals.values())
    failures = [matching for result in results
                for matching in result['failures']]
    assert totals.get('canonical_failed', 0) == 0
    assert totals.get('adaptive_failed', 0) == 0, failures
    if tuple(sorted(first_mates)) == FIRST_MATES:
        assert checked == TOTAL_MATCHINGS

    print('NONUPLE-TRANSPOSITION N 19',
          'labelled-types', checked,
          'certificates', dict(sorted(totals.items())),
          'quotient-matchings-tried', sum(r['attempt_sum'] for r in results),
          'max-attempts', max(r['max_attempts'] for r in results),
          'failed', 0,
          'seconds', round(perf_counter() - started, 2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--first-mate', type=int, action='append',
                        choices=FIRST_MATES)
    arguments = parser.parse_args()
    if arguments.workers < 1:
        parser.error('--workers must be positive')
    selected = tuple(arguments.first_mate) if arguments.first_mate else FIRST_MATES
    if len(set(selected)) != len(selected):
        parser.error('--first-mate values must be distinct')
    check(arguments.workers, selected)
