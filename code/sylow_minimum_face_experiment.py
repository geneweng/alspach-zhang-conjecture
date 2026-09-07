"""Exact and uniformly sampled minimum-face tests for Sylow quotients.

Each quotient x-cycle must select exactly one involution dart. Contracting
these cycles turns this condition into a pregraph perfect-matching problem;
ordinary loops cannot be selected. Completion counts give an exact uniform
sampler, with no colouring or monodromy condition in its distribution.

Run: python3 code/sylow_minimum_face_experiment.py 5 7 8 9
"""

import argparse
import json
import random
from collections import Counter
from fractions import Fraction
from pathlib import Path
from time import perf_counter

from cayley_snark_check import closure, inv, mul, order, psl2
from sylow_monodromy_review import Quotient, sylow


class CountLimit(Exception):
    pass


class MinimumFace:
    def __init__(self, quotient, max_states=500000):
        self.quotient = quotient
        self.cycles, self.cycle_of = [], {}
        unseen = set(range(quotient.n))
        while unseen:
            start, cycle = min(unseen), []
            current = start
            while current not in cycle:
                cycle.append(current)
                self.cycle_of[current] = len(self.cycles)
                current = quotient.edges[quotient.x_at[current]][1]
            assert current == start
            self.cycles.append(cycle)
            unseen -= set(cycle)
        s = order(quotient.x)
        assert all(len(cycle) == s for cycle in self.cycles)
        self.n = len(self.cycles)
        self.options = [[] for _ in self.cycles]
        for item in sorted(set(quotient.a_at.values())):
            u, v = quotient.edges[item]
            first = self.cycle_of[u]
            if v is None:
                self.options[first].append((1 << first, item))
            else:
                second = self.cycle_of[v]
                if first != second:
                    mask = (1 << first) | (1 << second)
                    self.options[first].append((mask, item))
                    self.options[second].append((mask, item))
        self.cache = {0: 1}
        self.max_states = max_states
        self.full = (1 << self.n) - 1

    def choices(self, mask):
        best = None
        todo = mask
        while todo:
            bit = todo & -todo
            vertex = bit.bit_length() - 1
            todo ^= bit
            choices = [(covered, item) for covered, item in self.options[vertex]
                       if covered & mask == covered]
            if best is None or len(choices) < len(best):
                best = choices
                if len(best) <= 1:
                    break
        return best

    def count(self, mask=None):
        if mask is None:
            mask = self.full
        if mask in self.cache:
            return self.cache[mask]
        if len(self.cache) >= self.max_states:
            raise CountLimit
        value = sum(self.count(mask ^ covered)
                    for covered, _ in self.choices(mask))
        self.cache[mask] = value
        return value

    def lift(self, chosen):
        selected = set(chosen)
        q = self.quotient
        for cycle in self.cycles:
            positions = [i for i, u in enumerate(cycle) if q.a_at[u] in chosen]
            assert len(positions) == 1
            start = positions[0]
            for offset in range(1, len(cycle), 2):
                selected.add(q.x_at[cycle[(start + offset) % len(cycle)]])
        assert all(sum(item in selected for item in inc) == 1 for inc in q.inc)
        return selected

    def enumerate(self, mask=None, chosen=frozenset()):
        if mask is None:
            mask = self.full
        if not mask:
            yield self.lift(chosen)
            return
        for covered, item in self.choices(mask):
            if self.count(mask ^ covered):
                yield from self.enumerate(mask ^ covered, chosen | {item})

    def sample(self, rng):
        mask, chosen = self.full, set()
        while mask:
            ticket = rng.randrange(self.count(mask))
            for covered, item in self.choices(mask):
                count = self.count(mask ^ covered)
                if ticket < count:
                    chosen.add(item)
                    mask ^= covered
                    break
                ticket -= count
            else:
                raise AssertionError('sampling weights disagree')
        return self.lift(chosen)

    def statistics(self, exact_limit, samples, seed):
        total = self.count()
        if total == 0:
            return {'status': 'EMPTY', 'quotient_cycles': self.n}
        exact = total <= exact_limit
        rng = random.Random(seed)
        starts = self.enumerate() if exact else (self.sample(rng) for _ in range(samples))
        histogram = Counter()
        witnessed_kinds = set()
        for selected in starts:
            bad = sum(length % 2 and word == self.quotient.identity
                      for length, word, _ in self.quotient.components(selected))
            histogram[bad] += 1
            # Check all small lifts and at least one positive and negative
            # example in each larger test against a separate full traversal.
            if self.quotient.n <= 63 or (bad == 0) not in witnessed_kinds:
                lengths = self.quotient.verify_lift(selected, require_even=False)
                assert sum(length % 2 for length in lengths) == len(self.quotient.subgroup)*bad
                witnessed_kinds.add(bad == 0)
        tested = sum(histogram.values())
        if exact:
            assert tested == total
        mean = Fraction(sum(bad * count for bad, count in histogram.items()), tested)
        return {'status': 'EXACT' if exact else 'SAMPLED', 'quotient_cycles': self.n,
                'face_size': total, 'states': len(self.cache), 'tested': tested,
                'bad_histogram': dict(sorted(histogram.items())),
                'mean_bad': str(mean), 'mean_bad_decimal': float(mean)}


def check(q, args):
    started = perf_counter()
    group = sorted(psl2(q))
    a = next(g for g in group if order(g) == 2)
    assert {mul(mul(inv(g), a), g) for g in group} == {
        g for g in group if order(g) == 2}
    subgroup = sylow(group, a)
    centralizer = [g for g in group if mul(a, g) == mul(g, a)]
    seen, totals, records = set(), Counter(), []
    for x in group:
        s = order(x)
        if x in seen or s <= 1 or s % 2 == 0:
            continue
        seen |= {mul(mul(inv(g), y), g) for g in centralizer for y in (x, inv(x))}
        if args.x_order and s != args.x_order:
            continue
        if len(closure([a, x])) != len(group):
            continue
        quotient = Quotient(group, subgroup, a, x)
        face = MinimumFace(quotient, args.max_states)
        try:
            result = face.statistics(args.exact_limit, args.samples, args.seed)
        except CountLimit:
            result = {'status': 'COUNT_LIMIT', 'quotient_cycles': face.n,
                      'states': len(face.cache)}
        totals[result['status']] += 1
        records.append({'q': q, 'group_order': len(group), 'sylow_order': len(subgroup),
                        'a': a, 'x': x, 'x_order': s, 'ax_order': order(mul(a, x)),
                        'seed': args.seed, **result})
        print('FACE', q, 'x_order', s, 'ax_order', order(mul(a, x)), result, flush=True)
        if result['status'] == 'EXACT' and Fraction(result['mean_bad']) >= 1:
            print('MOMENT-OBSTRUCTION', q, 'a', a, 'x', x, flush=True)
        if args.pairs and sum(totals.values()) >= args.pairs:
            break
    print('SUMMARY', q, dict(totals), 'seconds', round(perf_counter()-started, 3), flush=True)
    return records


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('q', nargs='*', type=int, default=[5, 7, 8, 9])
    parser.add_argument('--x-order', type=int)
    parser.add_argument('--pairs', type=int, default=0)
    parser.add_argument('--exact-limit', type=int, default=100000)
    parser.add_argument('--samples', type=int, default=1000)
    parser.add_argument('--max-states', type=int, default=500000)
    parser.add_argument('--seed', type=int, default=20260905)
    parser.add_argument('--json-output', type=Path)
    args = parser.parse_args()
    records = []
    for q in args.q:
        records.extend(check(q, args))
    if args.json_output:
        args.json_output.write_text(json.dumps(records, indent=2) + '\n')
