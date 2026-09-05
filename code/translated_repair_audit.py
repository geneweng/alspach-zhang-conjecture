#!/usr/bin/env python3
"""Stress-test strict odd-circuit descent by translated alternating-cycle flips.

Run: python3 code/translated_repair_audit.py --samples 128

Complete two-enumerator A5 check:
    python3 code/translated_repair_audit.py --samples 100000 --groups A5 A5_alt --independent

Sampling imposes only a perfect matching and no all-a pentagon. A reported
stall exhausts every translate and every individual difference circuit;
finding a decreasing move is only a positive finite test, not a theorem.
"""

import argparse
import hashlib
import json
import random
from collections import Counter
from itertools import combinations

from pysat.solvers import Cadical153

from cayley_snark_check import closure, inv, mul, order
from cdc_palette_experiment import cycle_perm, petersen
from colour_support_experiment import (
    complementary_factor, matching_for_support, support_system, verify_colouring,
)


def bits(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


class RepairGraph:
    """Canonical Cayley geometry without any auxiliary quotient flow."""

    def __init__(self, generators, expected_order):
        self.group = sorted(closure(generators))
        self.n = len(self.group)
        assert self.n == expected_order
        self.index = {g: i for i, g in enumerate(self.group)}
        labelled = {}
        for i, g in enumerate(self.group):
            for j, generator in enumerate(generators):
                v = self.index[mul(g, generator)]
                labelled[tuple(sorted((i, v)))] = 'a' if j == 0 else 'x'
        self.edges = sorted(labelled)
        self.kinds = [labelled[e] for e in self.edges]
        self.edge_index = {edge: e for e, edge in enumerate(self.edges)}
        self.inc = [[] for _ in range(self.n)]
        for e, (u, v) in enumerate(self.edges):
            self.inc[u].append(e)
            self.inc[v].append(e)
        assert all(len(incident) == 3 for incident in self.inc)
        self.a_edges = [e for e, kind in enumerate(self.kinds) if kind == 'a']
        a_position = {e: i for i, e in enumerate(self.a_edges)}
        a_at = [next(a_position[e] for e in incident if self.kinds[e] == 'a')
                for incident in self.inc]
        self.cycles, self.cycle_edges, self.orders = [], [], []
        unseen = set(range(self.n))
        while unseen:
            start = min(unseen)
            current, previous = start, None
            vertices, x_edges = [], []
            while True:
                unseen.remove(current)
                vertices.append(current)
                e = next(e for e in self.inc[current]
                         if self.kinds[e] == 'x' and e != previous)
                x_edges.append(e)
                u, v = self.edges[e]
                current, previous = (v if u == current else u), e
                if current == start:
                    break
            assert len(vertices) == 5
            self.cycles.append(vertices)
            self.cycle_edges.append(x_edges)
            self.orders.append([a_at[v] for v in vertices])
        self.q = len(self.cycles)
        self.pentagon_masks = [sum(1 << self.a_edges[i] for i in order)
                              for order in self.orders]
        self.actions = []
        for h in self.group:
            vertices = [self.index[mul(h, g)] for g in self.group]
            action = [self.edge_index[tuple(sorted((vertices[u], vertices[v])))]
                      for u, v in self.edges]
            assert all(self.kinds[e] == self.kinds[f] for e, f in enumerate(action))
            self.actions.append(action)
        self.fingerprint = hashlib.sha256(json.dumps(
            [self.group, self.edges, self.kinds], separators=(',', ':')).encode()).hexdigest()

    def selected(self, matching):
        return {i for i, e in enumerate(self.a_edges) if matching & (1 << e)}

    def matching_stabilizer_order(self, matching):
        edges = list(bits(matching))
        return sum(sum(1 << action[e] for e in edges) == matching for action in self.actions)

    def valid(self, matching):
        return all(matching & mask != mask for mask in self.pentagon_masks)

    def factor_lengths(self, matching):
        unseen = set(range(self.n))
        lengths = []
        while unseen:
            start = min(unseen)
            current, previous, length = start, None, 0
            while True:
                unseen.remove(current)
                length += 1
                e = next(e for e in self.inc[current]
                         if not matching & (1 << e) and e != previous)
                u, v = self.edges[e]
                current, previous = (v if u == current else u), e
                if current == start:
                    break
            lengths.append(length)
        return sorted(lengths)

    def difference_circuits(self, matching, action):
        translated = sum(1 << action[e] for e in bits(matching))
        unseen = matching ^ translated
        while unseen:
            first = (unseen & -unseen).bit_length() - 1
            start = self.edges[first][0]
            current, previous, circuit = start, None, 0
            while True:
                e = next(e for e in self.inc[current]
                         if (matching ^ translated) & (1 << e) and e != previous)
                unseen ^= 1 << e
                circuit |= 1 << e
                u, v = self.edges[e]
                current, previous = (v if u == current else u), e
                if current == start:
                    break
            assert circuit.bit_count() % 2 == 0
            yield circuit

    def check_flip(self, matching, g, circuit, candidate):
        """Check the move certificate using sets and an independent traversal."""
        translated = {self.actions[g][e] for e in bits(matching)}
        difference = set(bits(matching)) ^ translated
        changed = set(bits(circuit))
        assert changed and changed <= difference
        assert set(bits(candidate)) == set(bits(matching)) ^ changed
        touched = {v for e in changed for v in self.edges[e]}
        for v in touched:
            assert set(self.inc[v]) & difference == set(self.inc[v]) & changed
            assert len(set(self.inc[v]) & changed) == 2
        reached, todo = set(), [min(touched)]
        while todo:
            v = todo.pop()
            if v in reached:
                continue
            reached.add(v)
            for e in self.inc[v]:
                if e in changed:
                    todo.extend(self.edges[e])
        assert reached == touched
        assert self.valid(candidate)
        assert sum(n % 2 for n in self.factor_lengths(candidate)) < sum(
            n % 2 for n in self.factor_lengths(matching))
        self.check(candidate)

    def audit(self, matching, exhaustive=False):
        initial = sum(length % 2 for length in self.factor_lengths(matching))
        outcomes, distinct = Counter(), set()
        witness = None
        for g, action in enumerate(self.actions):
            for circuit in self.difference_circuits(matching, action):
                candidate = matching ^ circuit
                if candidate in distinct:
                    continue
                distinct.add(candidate)
                if not self.valid(candidate):
                    outcomes['invalid'] += 1
                    continue
                odd = sum(length % 2 for length in self.factor_lengths(candidate))
                outcomes[str(odd)] += 1
                if odd < initial and witness is None:
                    witness = {'translation': self.group[g], 'circuit': hex(circuit),
                               'new_oddness': odd, 'matching': hex(candidate)}
                    self.check_flip(matching, g, circuit, candidate)
                    if not exhaustive:
                        return witness, dict(outcomes)
        return witness, dict(outcomes)

    def check(self, matching):
        assert all(sum(bool(matching & (1 << e)) for e in incident) == 1
                   for incident in self.inc)
        selected = self.selected(matching)
        assert matching_for_support(self, selected) == set(bits(matching))
        solution, _, cycles = support_system(self, selected)
        assert sorted(cycles) == complementary_factor(self, selected)
        assert sum(sign for _, sign in cycles) == sum(
            length % 2 for length in self.factor_lengths(matching))
        if solution is not None:
            verify_colouring(self, selected, solution)


def sample_matchings(graph, count, rng, quotient_only, invariant=(), status=None):
    if status is not None:
        status['exhausted'] = False
    clauses = []
    for incident in graph.inc:
        clauses.append([e + 1 for e in incident])
        clauses.extend([-(e + 1), -(f + 1)] for e, f in combinations(incident, 2))
    for order in graph.orders:
        edges = [graph.a_edges[i] for i in order]
        if quotient_only:
            clauses.extend([-(e + 1), -(f + 1)] for e, f in combinations(edges, 2))
        else:
            clauses.append([-(e + 1) for e in edges])
    for h in invariant:
        action = graph.actions[graph.index[h]]
        for e, f in enumerate(action):
            clauses.extend([[-(e + 1), f + 1], [e + 1, -(f + 1)]])
    with Cadical153(bootstrap_with=clauses) as solver:
        for _ in range(count):
            solver.set_phases([v if rng.randrange(2) else -v
                               for v in range(1, len(graph.edges) + 1)])
            if not solver.solve():
                if status is not None:
                    status['exhausted'] = True
                return
            matching = sum(1 << (v - 1) for v in solver.get_model() if v > 0)
            assert graph.valid(matching)
            yield matching
            solver.add_clause([-(e + 1) for e in bits(matching)])


def independent_supports(graph):
    """Backtrack over the ten local patterns; no SAT solver is used."""
    states = []
    for cyclic in graph.orders:
        edges = [graph.a_edges[i] for i in cyclic]
        states.append([1 << e for e in edges] + [sum(
            1 << edges[(i + j) % 5] for j in range(3)) for i in range(5)])

    def visit(todo, known, selected):
        if not todo:
            support = graph.selected(selected)
            matching = sum(1 << e for e in matching_for_support(graph, support))
            graph.check(matching)
            yield matching
            return
        best = None
        for C in bits(todo):
            mask = graph.pentagon_masks[C]
            compatible = [state for state in states[C] if not (state ^ selected) & known & mask]
            if not compatible:
                return
            if best is None or len(compatible) < len(best[1]):
                best = C, compatible
        C, compatible = best
        for state in compatible:
            yield from visit(todo ^ (1 << C), known | graph.pentagon_masks[C], selected | state)

    yield from visit((1 << graph.q) - 1, 0, 0)


def classify_a5_pairs():
    """Check that the two A5 examples cover all involution/5-cycle pairs."""
    examples = cases()
    first = examples['A5'][0]
    group, x = sorted(closure(first)), first[1]
    involutions = [a for a in group if order(a) == 2]
    generating = {a for a in involutions if len(closure([a, x])) == 60}
    centralizer = [h for h in group if mul(h, x) == mul(x, h)]
    orbits = []
    for name in ['A5', 'A5_alt']:
        a = examples[name][0][0]
        orbit = {mul(mul(inv(h), a), h) for h in centralizer}
        orbits.append(orbit)
    assert len(involutions) == 15 and len(generating) == 10
    assert len(orbits[0]) == len(orbits[1]) == 5
    assert not orbits[0] & orbits[1] and orbits[0] | orbits[1] == generating
    print('A5-PAIR-CLASSIFICATION', 'involutions', 15, 'generating', 10,
          'centralizer_orbits', [5, 5], flush=True)


def check_petersen_stalls():
    """A transitive but nonregular action where every allowed start stalls."""
    graph = object.__new__(RepairGraph)
    graph.n, original_edges = petersen()
    graph.edges = sorted(tuple(sorted(e)) for e in original_edges)
    graph.edge_index = {e: i for i, e in enumerate(graph.edges)}
    graph.kinds = ['a' if (u < 5) != (v < 5) else 'x' for u, v in graph.edges]
    graph.inc = [[] for _ in range(graph.n)]
    for e, (u, v) in enumerate(graph.edges):
        graph.inc[u].append(e)
        graph.inc[v].append(e)
    graph.a_edges = [e for e, kind in enumerate(graph.kinds) if kind == 'a']
    a_position = {e: i for i, e in enumerate(graph.a_edges)}
    graph.cycles = [list(range(5)), [5, 7, 9, 6, 8]]
    graph.q = 2
    graph.cycle_edges = [[graph.edge_index[tuple(sorted((v, cycle[(i + 1) % 5])))]
                         for i, v in enumerate(cycle)] for cycle in graph.cycles]
    graph.orders = [[next(a_position[e] for e in graph.inc[v] if graph.kinds[e] == 'a')
                     for v in cycle] for cycle in graph.cycles]
    graph.pentagon_masks = [sum(1 << graph.a_edges[i] for i in order) for order in graph.orders]
    rotation = tuple((v // 5) * 5 + (v + 1) % 5 for v in range(10))
    exchange = tuple((1 - v // 5) * 5 + (2 * v) % 5 for v in range(10))
    graph.group = sorted(closure([rotation, exchange]))
    assert len(graph.group) == 20
    assert {g[0] for g in graph.group} == set(range(10))
    graph.actions = [[graph.edge_index[tuple(sorted((g[u], g[v])))]
                      for u, v in graph.edges] for g in graph.group]
    assert all(graph.kinds[e] == graph.kinds[f] for action in graph.actions
               for e, f in enumerate(action))
    for i in range(5):
        matching = sum(1 << e for e in matching_for_support(graph, {i}))
        graph.check(matching)
        assert graph.matching_stabilizer_order(matching) == 4
        assert graph.factor_lengths(matching) == [5, 5]
        witness, outcomes = graph.audit(matching, exhaustive=True)
        assert witness is None and outcomes == {'2': 4}
        assert {circuit.bit_count() for action in graph.actions
                for circuit in graph.difference_circuits(matching, action)} == {8}
    print('PETERSEN-STALL-CONTROL', 'transitive_group_order', 20,
          'locally_valid_starts', 5, 'distinct_moves_each', 4,
          'matching_stabilizer_order', 4, 'factor_lengths_always', [5, 5], flush=True)


def cases():
    x5 = cycle_perm(5, list(range(5)))
    x6 = cycle_perm(6, list(range(5)))
    return {
        'A5': ([mul(cycle_perm(5, [0, 1]), cycle_perm(5, [2, 3])), x5], 60),
        'A5_alt': ([mul(cycle_perm(5, [0, 2]), cycle_perm(5, [1, 3])), x5], 60),
        'S5': ([cycle_perm(5, [0, 1]), x5], 120),
        'S5_alt': ([cycle_perm(5, [0, 2]), x5], 120),
        'A6': ([mul(cycle_perm(6, [0, 5]), cycle_perm(6, [1, 2])), x6], 360),
        'S6': ([cycle_perm(6, [0, 5]), x6], 720),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--samples', type=int, default=128)
    parser.add_argument('--groups', nargs='+', choices=list(cases()), default=['A5', 'S5', 'A6', 'S6'])
    parser.add_argument('--mode', choices=['local', 'quotient'], default='local')
    parser.add_argument('--invariant', choices=[
        'none', 'point', 'pair', 'a', 'three', 'blocks', 'centralizer'], default='none')
    parser.add_argument('--independent', action='store_true',
                        help='compare complete A5 enumeration with non-SAT backtracking')
    parser.add_argument('--descend', action='store_true',
                        help='follow decreasing moves all the way to an even 2-factor or a stall')
    args = parser.parse_args()
    if args.samples < 1:
        parser.error('--samples must be positive')
    classify_a5_pairs()
    check_petersen_stalls()
    for name in args.groups:
        generators, size = cases()[name]
        graph = RepairGraph(generators, size)
        print('GRAPH', name, graph.n, graph.fingerprint, flush=True)
        invariant = []
        if args.invariant == 'point':
            invariant = [h for h in graph.group if h[-1] == len(h) - 1]
        elif args.invariant == 'pair':
            invariant = [h for h in graph.group if h[-2:] == tuple(range(len(h) - 2, len(h)))]
        elif args.invariant == 'a':
            invariant = [generators[0]]
        elif args.invariant == 'three':
            invariant = [cycle_perm(len(generators[0]), [0, 1, 2])]
        elif args.invariant == 'centralizer':
            invariant = [h for h in graph.group if mul(h, generators[0]) == mul(generators[0], h)]
        elif args.invariant == 'blocks':
            assert len(generators[0]) == 6
            invariant = [cycle_perm(6, [0, 1, 2]), cycle_perm(6, [3, 4, 5]),
                         mul(cycle_perm(6, [0, 3, 1, 4]), cycle_perm(6, [2, 5]))]
            if name == 'S6':
                invariant.append(cycle_perm(6, [0, 1]))
        subgroup_size = len(closure(invariant)) if invariant else 1
        symmetry_divisor = subgroup_size & -subgroup_size
        print('INVARIANCE', args.invariant, 'subgroup_order', subgroup_size,
              'odd_cycle_divisor', symmetry_divisor, flush=True)
        counts, odd_counts, path_lengths = Counter(), Counter(), Counter()
        status = {}
        enumerated = set() if args.independent else None
        if args.independent:
            assert name in ('A5', 'A5_alt') and args.mode == 'local' and args.invariant == 'none'
        for trial, matching in enumerate(sample_matchings(
                graph, args.samples, random.Random(509), args.mode == 'quotient', invariant, status)):
            lengths = graph.factor_lengths(matching)
            odd = sum(length % 2 for length in lengths)
            assert odd % symmetry_divisor == 0
            odd_counts[odd] += 1
            counts['sampled'] += 1
            if enumerated is not None:
                assert matching not in enumerated
                enumerated.add(matching)
            if not odd:
                counts['already_even'] += 1
                path_lengths[0] += 1
                continue
            current, steps = matching, 0
            while True:
                witness, outcomes = graph.audit(current)
                if witness is None:
                    graph.check(current)
                    counts['stalled'] += 1
                    print('STALL', json.dumps({'group': name, 'trial': trial,
                          'fingerprint': graph.fingerprint, 'matching': hex(current),
                          'starting_matching': hex(matching), 'steps': steps,
                          'factor_lengths': graph.factor_lengths(current),
                          'outcomes': outcomes}), flush=True)
                    break
                current = int(witness['matching'], 16)
                steps += 1
                if not args.descend or witness['new_oddness'] == 0:
                    counts['decreased'] += 1
                    path_lengths[steps] += 1
                    break
            if counts['stalled']:
                break
            if trial % 512 == 511:
                print('PROGRESS', name, trial + 1, dict(counts), flush=True)
        print('DONE', name, dict(counts), 'odd_circuit_counts', dict(sorted(odd_counts.items())),
              'exhausted', status['exhausted'], 'path_lengths', dict(sorted(path_lengths.items())), flush=True)
        if args.independent:
            assert status['exhausted'] and not counts['stalled']
            independent = list(independent_supports(graph))
            assert len(independent) == len(set(independent))
            assert set(independent) == enumerated
            expected = {'A5': {0: 3475, 2: 2330, 4: 70},
                        'A5_alt': {0: 2935, 2: 1590, 4: 100}}
            assert dict(odd_counts) == expected[name]
            symmetry_counts = Counter()
            for matching in independent:
                odd = sum(length % 2 for length in graph.factor_lengths(matching))
                if odd:
                    stabilizer = graph.matching_stabilizer_order(matching)
                    assert odd % (stabilizer & -stabilizer) == 0
                    symmetry_counts[odd, stabilizer] += 1
            print('INDEPENDENT-EXHAUSTIVE', name, 'supports', len(independent),
                  'exact_set_agreement', True, flush=True)
            print('STABILIZER-DIVISIBILITY', name, dict(sorted(symmetry_counts.items())), flush=True)


if __name__ == '__main__':
    main()
