"""Review experiment for the Sylow-2 monodromy proposal.

Run from the repository root:
    python3 code/sylow_monodromy_review.py

This checks one representative per generating-pair orbit under inner
conjugation and x-inversion for the listed PSL(2,q) groups. It chooses a
Sylow subgroup containing a at each representative. A positive result
therefore proves existence for some Sylow subgroup at every pair in that
orbit, not for every Sylow subgroup.

The search uses perfect-matching clauses, then excludes any complementary
odd circuit whose word is the identity. Positive witnesses are independently
lifted and traversed on the full Cayley graph. A second search forbids all
odd quotient circuits to distinguish quotient colouring from monodromy.
No existing census result files are read or changed.
"""
import argparse
from collections import Counter
from fractions import Fraction
from itertools import combinations
from time import perf_counter

from cayley_snark_check import closure, inv, mul, order, psl2
from cayley_snark_check2 import colour_pregraph
from pysat.solvers import Cadical153


def sylow(group, a):
    target = len(group) & -len(group)
    generators = [a]
    subgroup = closure(generators)
    candidates = [g for g in group if order(g) & (order(g) - 1) == 0]
    while len(subgroup) < target:
        for g in candidates:
            if g in subgroup:
                continue
            gi = inv(g)
            if not all(mul(mul(gi, h), g) in subgroup for h in subgroup):
                continue
            extension = closure(generators + [g])
            if len(extension) & (len(extension) - 1) == 0:
                generators.append(g)
                subgroup = extension
                break
        else:
            raise AssertionError('Sylow construction stalled')
    assert len(subgroup) == target
    return subgroup


class Quotient:
    def __init__(self, group, subgroup, a, x):
        self.group, self.subgroup, self.a, self.x = group, subgroup, a, x
        self.identity = tuple(range(len(a)))
        self.coset, self.reps = {}, []
        for g in group:
            if g not in self.coset:
                cid = len(self.reps)
                self.reps.append(g)
                for h in subgroup:
                    self.coset[mul(h, g)] = cid
        self.n = len(self.reps)
        self.edges, self.labels = [], []
        self.a_at, self.x_at = {}, {}
        for i, g in enumerate(self.reps):
            if i in self.a_at:
                continue
            j = self.coset[mul(g, a)]
            self.a_at[i] = self.a_at[j] = len(self.edges)
            self.edges.append((i, None if i == j else j))
            self.labels.append(a)
        for i, g in enumerate(self.reps):
            j = self.coset[mul(g, x)]
            assert i != j
            self.x_at[i] = len(self.edges)
            self.edges.append((i, j))
            self.labels.append(x)
        self.inc = [[] for _ in self.reps]
        for e, (u, v) in enumerate(self.edges):
            self.inc[u].append(e)
            if v is not None:
                self.inc[v].append(e)
        assert all(len(es) == 3 for es in self.inc)

    def components(self, selected):
        unseen, answer = set(range(self.n)), []
        while unseen:
            start = min(unseen)
            stack, vertices, edges = [start], set(), set()
            while stack:
                u = stack.pop()
                if u in vertices:
                    continue
                vertices.add(u)
                for e in self.inc[u]:
                    if e in selected:
                        continue
                    edges.add(e)
                    stack.extend(v for v in self.edges[e] if v is not None)
            unseen -= vertices
            semis = sum(self.edges[e][1] is None for e in edges)
            if semis:
                assert semis == 2
                answer.append((len(edges), None, edges))
                continue
            word, current, previous = self.identity, start, -1
            for _ in range(len(edges)):
                e = next(e for e in self.inc[current]
                         if e not in selected and e != previous)
                u, v = self.edges[e]
                label = self.labels[e] if current == u else inv(self.labels[e])
                word = mul(word, label)
                current, previous = (v if current == u else u), e
            assert current == start
            r = self.reps[start]
            assert mul(mul(r, word), inv(r)) in self.subgroup
            answer.append((len(edges), word, edges))
        return answer

    def verify_lift(self, selected, require_even=True):
        adjacency = {g: [] for g in self.group}
        degree = Counter()
        for g in self.group:
            h = mul(g, self.a)
            if g < h:
                if self.a_at[self.coset[g]] in selected:
                    degree[g] += 1
                    degree[h] += 1
                else:
                    adjacency[g].append(h)
                    adjacency[h].append(g)
            h = mul(g, self.x)
            if self.x_at[self.coset[g]] in selected:
                degree[g] += 1
                degree[h] += 1
            else:
                adjacency[g].append(h)
                adjacency[h].append(g)
        assert all(degree[g] == 1 and len(adjacency[g]) == 2 for g in self.group)
        unseen, lengths = set(self.group), []
        while unseen:
            stack, component = [min(unseen)], set()
            while stack:
                g = stack.pop()
                if g not in component:
                    component.add(g)
                    stack.extend(adjacency[g])
            unseen -= component
            lengths.append(len(component))
        if require_even:
            assert all(length % 2 == 0 for length in lengths)
        return sorted(lengths)

    def uniform_statistics(self):
        """Exact first moment, with independent full lifts of every matching."""
        histogram = Counter()
        costs = {}
        a_items = set(self.a_at.values())

        def visit(uncovered, selected):
            if not uncovered:
                bad = sum(length % 2 and word == self.identity
                          for length, word, _ in self.components(selected))
                lengths = self.verify_lift(selected, require_even=False)
                assert sum(length % 2 for length in lengths) == (
                    len(self.subgroup) * bad)
                histogram[bad] += 1
                cost = sum(1 if self.edges[e][1] is None else 2
                           for e in selected & a_items)
                costs.setdefault(cost, Counter())[bad] += 1
                return
            u = min(uncovered)
            for e in self.inc[u]:
                covered = set(self.edges[e]) - {None}
                if covered <= uncovered:
                    visit(uncovered - covered, selected | {e})

        visit(set(range(self.n)), set())
        total = sum(histogram.values())
        assert total > 0
        expectation = Fraction(sum(bad * count for bad, count in histogram.items()), total)
        minimum = min(costs)
        minimal = costs[minimum]
        minimal_expectation = Fraction(sum(bad * count for bad, count in minimal.items()),
                                       sum(minimal.values()))
        return {'matchings': total, 'bad_circuit_histogram': dict(sorted(histogram.items())),
                'expected_bad': str(expectation), 'first_moment_certifies': expectation < 1,
                'minimum_a_darts': minimum, 'minimal_matchings': sum(minimal.values()),
                'minimal_expected_bad': str(minimal_expectation)}

    def solve(self, max_rounds=20000, require_even=False):
        clauses = []
        for es in self.inc:
            clauses.append([e + 1 for e in es])
            clauses.extend([-e - 1, -f - 1] for e, f in combinations(es, 2))
        with Cadical153(bootstrap_with=clauses) as solver:
            for attempts in range(1, max_rounds + 1):
                if not solver.solve():
                    return {'status': 'UNSAT', 'attempts': attempts - 1}
                selected = {v - 1 for v in solver.get_model() if v > 0}
                comps = self.components(selected)
                bad = [es for length, word, es in comps
                       if length % 2 and word is not None
                       and (require_even or word == self.identity)]
                if not bad:
                    lifted = self.verify_lift(selected)
                    return {'status': 'SAT', 'attempts': attempts,
                            'odd_nonidentity': sum(length % 2 for length, word, _
                               in comps if word is not None),
                            'lift_circuits': len(lifted)}
                for es in bad:
                    solver.add_clause([e + 1 for e in es])
        return {'status': 'LIMIT', 'attempts': max_rounds}


def check(q):
    started = perf_counter()
    group = sorted(psl2(q))
    a = next(g for g in group if order(g) == 2)
    # PSL(2,q) has one involution class in these test groups. Verify it here.
    assert {mul(mul(inv(g), a), g) for g in group} == {
        g for g in group if order(g) == 2}
    subgroup = sylow(group, a)
    centralizer = [g for g in group if mul(a, g) == mul(g, a)]
    seen, counts = set(), Counter()
    max_attempts = 0
    for x in group:
        if x in seen or order(x) <= 1 or order(x) % 2 == 0:
            continue
        orbit = {mul(mul(inv(g), y), g)
                 for g in centralizer for y in (x, inv(x))}
        seen |= orbit
        if len(closure([a, x])) != len(group):
            continue
        quotient = Quotient(group, subgroup, a, x)
        if quotient.n <= 21:
            print('EXACT-MOMENT', q, 'x_order', order(x),
                  quotient.uniform_statistics(), flush=True)
        result = quotient.solve()
        counts[result['status']] += 1
        max_attempts = max(max_attempts, result['attempts'])
        if result['status'] != 'SAT':
            print('FAILURE', q, a, x, result, flush=True)
        else:
            counts['found_using_odd_quotient_circuit'] += bool(result['odd_nonidentity'])
            ordinary = quotient.solve(require_even=True)
            counts['ordinary_' + ordinary['status']] += 1
            if ordinary['status'] != 'SAT':
                if ordinary['status'] == 'UNSAT':
                    # Different encoding: solve for quotient colours directly.
                    assert colour_pregraph(quotient.n, quotient.edges) == (False, None)
                print('MONODROMY-ONLY', q, 'x_order', order(x), 'a', a,
                      'x', x, 'ordinary', ordinary, 'monodromy', result,
                      flush=True)
    print('PSL2', q, 'order', len(group), 'sylow_order', len(subgroup),
          'quotient_order', len(group) // len(subgroup), 'pair_orbits', dict(counts),
          'max_attempts', max_attempts, 'seconds', round(perf_counter()-started, 3),
          flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('q', nargs='*', type=int,
                        default=[5, 7, 8, 9, 11, 13, 17, 19, 23])
    args = parser.parse_args()
    for q in args.q:
        check(q)
