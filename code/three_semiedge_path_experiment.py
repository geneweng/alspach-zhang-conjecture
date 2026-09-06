#!/usr/bin/env python3
"""Test the three-semi-edge alternating-path question from the survey.

Let C be an odd cycle and let A consist of three semi-edges and a matching
of all other vertices by chords.  Starting with the perfect pregraph
matching A, flip one simple A-alternating path between a prescribed pair of
semi-edges.  The resulting matching is a sign certificate if every odd
circuit of its complement contains an odd number of chords.

The sign condition says that the corresponding circuit has non-trivial
voltage in the two-sheeted cover.  It therefore lifts to an even circuit;
the capped complementary path always lifts to an even circuit as well.

By default the script exhausts every labelled configuration through order
13, and checks all three pairs of semi-edges separately.  The optional
random mode is intended for larger stress tests.
"""

import argparse
import random
from collections import Counter
from itertools import combinations
from math import comb, prod
from time import perf_counter


def perfect_matchings(vertices):
    """Generate all perfect matchings of a sorted tuple."""
    if not vertices:
        yield ()
        return
    first = vertices[0]
    for index in range(1, len(vertices)):
        second = vertices[index]
        remaining = vertices[1:index] + vertices[index + 1:]
        for rest in perfect_matchings(remaining):
            yield ((first, second),) + rest


class ThreeSemiedgeDiagram:
    """An odd cycle, three semi-edges, and chords at all other vertices."""

    def __init__(self, n, semiedges, matching):
        self.n = n
        self.semiedges = tuple(sorted(semiedges))
        self.semi_set = set(self.semiedges)
        self.matching = tuple(matching)
        self.partner = [-1] * n
        self.chord_at = [-1] * n
        for chord, (first, second) in enumerate(self.matching):
            self.partner[first] = second
            self.partner[second] = first
            self.chord_at[first] = self.chord_at[second] = chord
        assert len(self.semiedges) == 3
        assert set(vertex for edge in self.matching for vertex in edge) == (
            set(range(n)) - self.semi_set
        )

    def is_certificate(self, cycle_mask, chord_mask, endpoints):
        """Check one flipped path by an independent complement traversal."""
        adjacency = [[] for _ in range(self.n)]
        for edge in range(self.n):
            if cycle_mask & (1 << edge):
                continue
            first, second = edge, (edge + 1) % self.n
            adjacency[first].append((second, False))
            adjacency[second].append((first, False))
        for chord, (first, second) in enumerate(self.matching):
            if not chord_mask & (1 << chord):
                continue
            adjacency[first].append((second, True))
            adjacency[second].append((first, True))

        # Verify directly that the flipped set is a perfect matching.
        matching_degree = [0] * self.n
        for edge in range(self.n):
            if cycle_mask & (1 << edge):
                matching_degree[edge] += 1
                matching_degree[(edge + 1) % self.n] += 1
        for chord, (first, second) in enumerate(self.matching):
            if not chord_mask & (1 << chord):
                matching_degree[first] += 1
                matching_degree[second] += 1
        for vertex in self.semiedges:
            if vertex not in endpoints:
                matching_degree[vertex] += 1
        assert matching_degree == [1] * self.n

        endpoint_mask = sum(1 << vertex for vertex in endpoints)
        unseen = (1 << self.n) - 1
        capped_components = 0
        while unseen:
            root_bit = unseen & -unseen
            stack = [root_bit.bit_length() - 1]
            component = degree_sum = chord_twice = 0
            while stack:
                vertex = stack.pop()
                if component & (1 << vertex):
                    continue
                component |= 1 << vertex
                for neighbour, is_chord in adjacency[vertex]:
                    degree_sum += 1
                    chord_twice += is_chord
                    stack.append(neighbour)
            unseen &= ~component
            if component & endpoint_mask:
                assert component & endpoint_mask == endpoint_mask
                capped_components += 1
                continue
            length = degree_sum // 2
            chords = chord_twice // 2
            if length % 2 == 1 and chords % 2 == 0:
                return False
        assert capped_components == 1
        return True

    def first_certificate(self, first, second):
        """Return leaves tried before a one-path certificate, or None."""
        assert first in self.semi_set and second in self.semi_set
        leaves = 0

        def search(vertex, used, cycle_mask, chord_mask):
            nonlocal leaves
            choices = (
                ((vertex - 1) % self.n, (vertex - 1) % self.n),
                (vertex, (vertex + 1) % self.n),
            )
            for edge, neighbour in choices:
                if used & (1 << neighbour):
                    continue
                next_used = used | (1 << neighbour)
                next_cycle = cycle_mask | (1 << edge)
                if neighbour in self.semi_set:
                    if neighbour != second:
                        continue
                    leaves += 1
                    if self.is_certificate(
                        next_cycle, chord_mask, {first, second}
                    ):
                        return leaves
                    continue

                partner = self.partner[neighbour]
                if next_used & (1 << partner):
                    continue
                result = search(
                    partner,
                    next_used | (1 << partner),
                    next_cycle,
                    chord_mask | (1 << self.chord_at[neighbour]),
                )
                if result is not None:
                    return result
            return None

        return search(first, 1 << first, 0, 0)

    def check_all_pairs(self):
        attempts = []
        for first, second in combinations(self.semiedges, 2):
            leaves = self.first_certificate(first, second)
            if leaves is None:
                return None
            attempts.append(leaves)
        return tuple(attempts)


def labelled_configuration_count(n):
    return comb(n, 3) * prod(range(1, n - 3, 2))


def exhaustive(n):
    started = perf_counter()
    configurations = pairs = 0
    attempt_distribution = Counter()
    maximum = 0
    for semiedges in combinations(range(n), 3):
        remaining = tuple(vertex for vertex in range(n)
                          if vertex not in semiedges)
        for matching in perfect_matchings(remaining):
            configurations += 1
            attempts = ThreeSemiedgeDiagram(
                n, semiedges, matching
            ).check_all_pairs()
            assert attempts is not None, (n, semiedges, matching)
            attempt_distribution.update(attempts)
            pairs += len(attempts)
            maximum = max(maximum, *attempts)
    assert configurations == labelled_configuration_count(n)
    result = {
        "n": n,
        "labelled_configurations": configurations,
        "terminal_pairs": pairs,
        "max_leaves": maximum,
        "attempt_distribution": dict(sorted(attempt_distribution.items())),
        "seconds": round(perf_counter() - started, 3),
    }
    print("EXACT", result, flush=True)


def random_matching(vertices, rng):
    vertices = list(vertices)
    rng.shuffle(vertices)
    return tuple(zip(vertices[::2], vertices[1::2]))


def random_checks(n, samples, seed, independent):
    started = perf_counter()
    rng = random.Random(seed)
    maximum = 0
    attempt_sum = pairs = 0
    for sample in range(samples):
        while True:
            semiedges = tuple(sorted(rng.sample(range(n), 3)))
            if not independent or all(
                (vertex + 1) % n not in semiedges for vertex in semiedges
            ):
                break
        remaining = tuple(vertex for vertex in range(n)
                          if vertex not in semiedges)
        matching = random_matching(remaining, rng)
        attempts = ThreeSemiedgeDiagram(
            n, semiedges, matching
        ).check_all_pairs()
        assert attempts is not None, (n, semiedges, matching, sample)
        maximum = max(maximum, *attempts)
        attempt_sum += sum(attempts)
        pairs += len(attempts)
    result = {
        "n": n,
        "samples": samples,
        "terminal_pairs": pairs,
        "independent_semiedges": independent,
        "mean_leaves": round(attempt_sum / pairs, 6),
        "max_leaves": maximum,
        "seed": seed,
        "seconds": round(perf_counter() - started, 3),
    }
    print("RANDOM", result, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--max-exact", type=int, default=13,
        help="exhaust odd orders 5 through this value (default: 13)",
    )
    parser.add_argument("--random-n", type=int)
    parser.add_argument("--samples", type=int, default=20_000)
    parser.add_argument("--seed", type=int, default=20260905)
    parser.add_argument(
        "--allow-adjacent", action="store_true",
        help="do not require the three random semi-edges to be independent",
    )
    args = parser.parse_args()
    assert args.max_exact < 5 or args.max_exact % 2 == 1
    for n in range(5, args.max_exact + 1, 2):
        exhaustive(n)
    if args.random_n is not None:
        assert args.random_n >= 5 and args.random_n % 2 == 1
        random_checks(
            args.random_n, args.samples, args.seed, not args.allow_adjacent
        )


if __name__ == "__main__":
    main()
