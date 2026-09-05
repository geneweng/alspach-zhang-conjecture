#!/usr/bin/env python3
"""Exact implication constraints for locally valid translated block flips.

The PSL(2,11) certificate refutes single-component local descent, not
Tait colourability. A simultaneous two-component flip repairs it.
The second certificate shows that minimal valid blocks need not improve,
although a union of two individually neutral valid flips can colour.

python3 code/translated_block_repair.py
python3 code/translated_block_repair.py --stress-psl11 512 --descend
"""

import argparse
import json
import random
from collections import Counter

from cayley_snark_check import closure, inv, mul, order, psl2
from colour_support_experiment import matching_for_support
from translated_repair_audit import RepairGraph, bits, sample_matchings


GENERATORS = [
    (2, 10, 0, 9, 5, 4, 11, 8, 7, 3, 1, 6),
    (0, 1, 3, 9, 11, 7, 2, 4, 5, 10, 6, 8),
]
MATCHING = int(
    '1c94060808010f02070554403a104d17f9a58054c2c24c4c0a2cb054d283154941'
    '155110504130c42a480a8012a2314306184a2a224992188951128c86130a4a225424'
    'c45248625422a6a9049504569a929a800a5a02116804202018594962811c6e925291'
    '00430a84c54d4492c33a81968a44528932488a5454902c', 16)
PATTERNS = [1, 1, 19, 19, 2, 28, 28, 2, 25, 25, 25]
MINIMAL_BLOCK_PATTERNS = [2, 2, 16, 14, 8, 7, 1, 14, 25, 14, 25]


def mixture(matching, circuits, choice):
    for i in bits(choice):
        matching ^= circuits[i]
    return matching


def implication_system(graph, matching, action):
    """Return difference circuits and the necessary/sufficient implications.

Variable z_i says to flip circuit i. A dangerous pentagon has one M-x-edge
and one translated-M-x-edge, disjoint in vertices. Removing the former
requires inserting the latter, hence z_i -> z_j.
"""
    translated = sum(1 << action[e] for e in bits(matching))
    assert graph.valid(matching) and graph.valid(translated)
    circuits = tuple(sorted(graph.difference_circuits(matching, action)))
    owner = {e: i for i, circuit in enumerate(circuits) for e in bits(circuit)}
    arcs, reasons = set(), {}
    for c, edges in enumerate(graph.cycle_edges):
        first = [e for e in edges if matching & (1 << e)]
        second = [e for e in edges if translated & (1 << e)]
        assert len(first) in (1, 2) and len(second) in (1, 2)
        if len(first) != 1 or len(second) != 1:
            continue
        e, f = first[0], second[0]
        if set(graph.edges[e]) & set(graph.edges[f]):
            continue
        i, j = owner[e], owner[f]
        if i != j:
            arcs.add((i, j))
            reasons.setdefault((i, j), []).append(c)
    return circuits, arcs, reasons


def sink_blocks(size, arcs):
    """Inclusion-minimal nonempty implication-closed choices (sink SCCs)."""
    reach = [{i} for i in range(size)]
    for i, j in arcs:
        reach[i].add(j)
    for k in range(size):
        for i in range(size):
            if k in reach[i]:
                reach[i].update(reach[k])
    unseen, blocks = set(range(size)), []
    while unseen:
        i = min(unseen)
        component = {j for j in reach[i] if i in reach[j]}
        unseen -= component
        if all(j in component for u, j in arcs if u in component):
            blocks.append(sum(1 << j for j in component))
    return blocks


def compact_support(graph, matching):
    """Encode the H-invariant support on eleven representative pentagons."""
    a = GENERATORS[0]
    subgroup = [i for i, h in enumerate(graph.group) if mul(h, a) == mul(a, h)]
    assert len(subgroup) == 12
    assert all(sum(1 << graph.actions[i][e] for e in bits(matching)) == matching
               for i in subgroup)
    cycle_by_edges = {frozenset(edges): i for i, edges in enumerate(graph.cycle_edges)}
    unseen, patterns, selected = set(range(graph.q)), [], set()
    while unseen:
        c = min(unseen)
        support_mask = sum(1 << j for j, pos in enumerate(graph.orders[c])
                           if matching & (1 << graph.a_edges[pos]))
        patterns.append((c, support_mask))
        for i in subgroup:
            action = graph.actions[i]
            d = cycle_by_edges[frozenset(action[e] for e in graph.cycle_edges[c])]
            unseen.discard(d)
            for j in bits(support_mask):
                selected.add(action[graph.a_edges[graph.orders[c][j]]])
    assert len(patterns) == 11
    support = {i for i, e in enumerate(graph.a_edges) if e in selected}
    assert sum(1 << e for e in matching_for_support(graph, support)) == matching
    return patterns


def matching_from_patterns(graph, patterns=PATTERNS):
    """Construct the certificate from eleven local patterns, without SAT."""
    a = GENERATORS[0]
    actions = [action for h, action in zip(graph.group, graph.actions)
               if mul(h, a) == mul(a, h)]
    selected = set()
    for c, pattern in enumerate(patterns):
        for j in bits(pattern):
            e = graph.a_edges[graph.orders[c][j]]
            selected.update(action[e] for action in actions)
    support = {i for i, e in enumerate(graph.a_edges) if e in selected}
    return sum(1 << e for e in matching_for_support(graph, support))


def check_local_implications():
    """Exhaust all 100 ordered pairs of local pentagon matching patterns."""
    edges = [{i, (i + 1) % 5} for i in range(5)]
    states = [{i} for i in range(5)] + [{i, (i + 2) % 5} for i in range(5)]
    tested = 0
    for first in states:
        for second in states:
            unseen, components = set(first ^ second), []
            while unseen:
                reached, todo = set(), [min(unseen)]
                while todo:
                    e = todo.pop()
                    if e in reached:
                        continue
                    reached.add(e)
                    todo.extend(f for f in unseen - reached if edges[e] & edges[f])
                unseen -= reached
                components.append(reached)
            owner = {e: i for i, part in enumerate(components) for e in part}
            arcs = []
            if len(first) == len(second) == 1:
                e, f = next(iter(first)), next(iter(second))
                if not edges[e] & edges[f]:
                    arcs.append((owner[e], owner[f]))
            for choice in range(1 << len(components)):
                selected = set(first)
                for i in bits(choice):
                    selected ^= components[i]
                assert len({v for e in selected for v in edges[e]}) == 2 * len(selected)
                predicted = all(not choice & (1 << i) or choice & (1 << j) for i, j in arcs)
                assert bool(selected) == bool(predicted)
                tested += 1
    print('LOCAL-IMPLICATION-CHECK', 'ordered_pattern_pairs', 100,
          'component_choices', tested, flush=True)


def independent_factor_lengths(graph, matching):
    """Check degree one, then find complementary components by set-based DFS."""
    selected = {graph.edges[e] for e in bits(matching)}
    degrees, adjacency = [0] * graph.n, [set() for _ in range(graph.n)]
    for u, v in selected:
        degrees[u] += 1
        degrees[v] += 1
    assert degrees == [1] * graph.n
    for u, v in set(graph.edges) - selected:
        adjacency[u].add(v)
        adjacency[v].add(u)
    assert all(len(ns) == 2 for ns in adjacency)
    unseen, lengths = set(range(graph.n)), []
    while unseen:
        reached, todo = set(), [min(unseen)]
        while todo:
            v = todo.pop()
            if v not in reached:
                reached.add(v)
                todo.extend(adjacency[v] - reached)
        unseen -= reached
        lengths.append(len(reached))
    return sorted(lengths)


def independent_single_flip_check(graph, matching):
    """Rebuild the graph and all translates with vertex-pair edge sets and BFS.

    No RepairGraph cycle traversal, stored edge action, validity predicate,
    or factor-length routine is used to establish the stall.
    """
    group, index = graph.group, {h: i for i, h in enumerate(graph.group)}
    a, x = GENERATORS
    edges = {tuple(sorted((i, index[mul(h, s)])))
             for i, h in enumerate(group) for s in (a, x, inv(x))}
    assert edges == set(graph.edges)
    chosen = {graph.edges[e] for e in bits(matching)}
    x_cycles, unseen = [], set(range(len(group)))
    while unseen:
        start, vertices, current = min(unseen), [], min(unseen)
        while current not in vertices:
            vertices.append(current)
            unseen.remove(current)
            current = index[mul(group[current], x)]
        assert current == start and len(vertices) == 5
        x_cycles.append({tuple(sorted((v, vertices[(j + 1) % 5])))
                         for j, v in enumerate(vertices)})

    def lengths(selected):
        incident = [set() for _ in group]
        degrees = [0] * len(group)
        for u, v in selected:
            degrees[u] += 1
            degrees[v] += 1
        assert degrees == [1] * len(group)
        for u, v in edges - selected:
            incident[u].add(v)
            incident[v].add(u)
        assert all(len(ns) == 2 for ns in incident)
        todo, sizes = set(range(len(group))), []
        while todo:
            stack, component = [min(todo)], set()
            while stack:
                v = stack.pop()
                if v in component:
                    continue
                component.add(v)
                stack.extend(incident[v] - component)
            todo -= component
            sizes.append(len(component))
        return sorted(sizes)

    initial = lengths(chosen)
    assert initial == [27] * 4 + [38] * 6 + [120, 204]
    seen, outcomes = set(), Counter()
    for h in group:
        vertices = [index[mul(h, g)] for g in group]
        translated = {tuple(sorted((vertices[u], vertices[v]))) for u, v in chosen}
        difference = chosen ^ translated
        adjacency = {}
        for u, v in difference:
            adjacency.setdefault(u, set()).add(v)
            adjacency.setdefault(v, set()).add(u)
        assert all(len(ns) == 2 for ns in adjacency.values())
        todo = set(adjacency)
        while todo:
            stack, reached = [min(todo)], set()
            while stack:
                v = stack.pop()
                if v in reached:
                    continue
                reached.add(v)
                stack.extend(adjacency[v] - reached)
            todo -= reached
            component = {e for e in difference if e[0] in reached}
            candidate = frozenset(chosen ^ component)
            if candidate in seen:
                continue
            seen.add(candidate)
            bad_pentagons = sum(not candidate & c for c in x_cycles)
            assert bad_pentagons > 0
            outcomes[sum(n % 2 for n in lengths(candidate)), bad_pentagons] += 1
    assert len(seen) == 168
    print('INDEPENDENT-BFS-STALL', len(seen), 'all_single_flips_invalid', True, flush=True)
    return outcomes


def audit_blocks(graph, matching, all_closed=False):
    """Search sink blocks, or all closed choices; return a witness or None.

    all_closed explicitly rejects large systems instead of silently claiming
    exhaustion. A sink-block stall does not refute all-closed-set descent.
    """
    odd = sum(n % 2 for n in graph.factor_lengths(matching))
    seen = {matching}
    for g, action in enumerate(graph.actions):
        translated = sum(1 << action[e] for e in bits(matching))
        if translated in seen:
            continue
        seen.add(translated)
        circuits, arcs, _ = implication_system(graph, matching, action)
        if all_closed:
            if len(circuits) > 16:
                raise ValueError('all-closed search exceeds 16 circuits; not exhaustive')
            choices = (s for s in range(1, (1 << len(circuits)) - 1)
                       if all(not s & (1 << i) or s & (1 << j) for i, j in arcs))
        else:
            choices = sink_blocks(len(circuits), arcs)
        for choice in choices:
            if choice == (1 << len(circuits)) - 1:
                continue  # this gives only the translate itself
            candidate = mixture(matching, circuits, choice)
            assert graph.valid(candidate)
            new_odd = sum(n % 2 for n in graph.factor_lengths(candidate))
            if new_odd < odd:
                graph.check(candidate)
                return {'translation': graph.group[g], 'choice': choice,
                        'circuit_lengths': [circuits[i].bit_count() for i in bits(choice)],
                        'new_oddness': new_odd, 'matching': hex(candidate)}
    return None


def check_certificate():
    def projective_map(a, b, c, d):
        assert (a * d - b * c) % 11 == 1
        values = []
        for z in range(12):
            numerator, denominator = (a, c) if z == 11 else ((a * z + b) % 11, (c * z + d) % 11)
            values.append(11 if denominator == 0 else numerator * pow(denominator, -1, 11) % 11)
        return tuple(values)

    assert GENERATORS == [projective_map(2, 7, 4, 9), projective_map(8, 0, 1, 7)]
    graph = RepairGraph(GENERATORS, 660)
    assert set(graph.group) == psl2(11)
    assert order(GENERATORS[0]) == 2 and order(GENERATORS[1]) == 5
    assert order(mul(*GENERATORS)) == 11
    assert matching_from_patterns(graph) == MATCHING
    graph.check(MATCHING)
    assert compact_support(graph, MATCHING) == list(enumerate(PATTERNS))
    assert graph.matching_stabilizer_order(MATCHING) == 12
    print('GRAPH', graph.fingerprint, 'patterns', PATTERNS, flush=True)
    seen, counts, systems, best = set(), Counter(), Counter(), None
    block_outcomes = Counter()
    for g, action in enumerate(graph.actions):
        circuits, arcs, reasons = implication_system(graph, MATCHING, action)
        blocks = sink_blocks(len(circuits), arcs)
        systems[len(circuits), tuple(sorted(c.bit_count() for c in blocks))] += 1
        assert len(circuits) <= 12
        for choice in range(1 << len(circuits)):
            candidate = mixture(MATCHING, circuits, choice)
            predicted = all(not choice & (1 << i) or choice & (1 << j) for i, j in arcs)
            assert bool(predicted) == graph.valid(candidate)
        for i, circuit in enumerate(circuits):
            candidate = MATCHING ^ circuit
            if candidate not in seen:
                seen.add(candidate)
                odd = sum(n % 2 for n in graph.factor_lengths(candidate))
                invalid = sum(candidate & mask == mask for mask in graph.pentagon_masks)
                counts[odd, invalid] += 1
            assert not graph.valid(candidate)
            assert any(u == i for u, _ in arcs)
        for choice in blocks:
            candidate = mixture(MATCHING, circuits, choice)
            assert graph.valid(candidate)
            odd = sum(n % 2 for n in graph.factor_lengths(candidate))
            block_outcomes[choice.bit_count(), odd] += 1
            if best is None or odd < best[0]:
                graph.check(candidate)
                best = (odd, graph.group[g], order(graph.group[g]),
                        [circuits[i].bit_count() for i in bits(choice)],
                        hex(candidate), graph.factor_lengths(candidate),
                        sorted(arcs), [i for i in bits(choice)])
    assert dict(systems) == {(0, ()): 12, (2, (2,)): 216, (3, (3,)): 144,
                             (4, (2,)): 144, (4, (4,)): 144}
    assert independent_single_flip_check(graph, MATCHING) == counts
    witness, original_outcomes = graph.audit(MATCHING, exhaustive=True)
    assert witness is None and original_outcomes == {'invalid': 168}
    assert best[0] == 0 and best[3] == [22, 214] and best[5] == [182, 478]
    print('SINGLE-FLIP-OUTCOMES', len(seen), dict(sorted(counts.items())), flush=True)
    print('IMPLICATION-SYSTEMS', dict(sorted(systems.items())), flush=True)
    print('SINK-BLOCK-OUTCOMES', dict(sorted(block_outcomes.items())), flush=True)
    print('BEST-SINK-BLOCK', best, flush=True)


def check_minimal_block_barrier():
    """Every minimal closed block fails; a union of two singleton blocks colours."""
    graph = RepairGraph(GENERATORS, 660)
    matching = matching_from_patterns(graph, MINIMAL_BLOCK_PATTERNS)
    assert compact_support(graph, matching) == list(enumerate(MINIMAL_BLOCK_PATTERNS))
    graph.check(matching)
    assert graph.matching_stabilizer_order(matching) == 12
    assert independent_factor_lengths(graph, matching) == [15] * 4 + [44] * 3 + [78] * 6
    seen, counts = {matching}, Counter()
    for action in graph.actions:
        translated = sum(1 << action[e] for e in bits(matching))
        if translated in seen:
            continue
        seen.add(translated)
        circuits, arcs, _ = implication_system(graph, matching, action)
        # Independently enumerate ALL choices, then determine the minimal
        # nonempty ones by inclusion, without the SCC routine.
        assert len(circuits) <= 8
        valid = {choice for choice in range(1, 1 << len(circuits))
                 if graph.valid(mixture(matching, circuits, choice))}
        minimal = {s for s in valid if not any(t != s and t & s == t for t in valid)}
        assert minimal == set(sink_blocks(len(circuits), arcs))
        for choice in minimal:
            candidate = mixture(matching, circuits, choice)
            odd = sum(n % 2 for n in independent_factor_lengths(graph, candidate))
            assert odd >= 4
            counts[choice.bit_count(), odd] += 1
    assert counts == Counter({(2, 4): 24, (1, 4): 48, (4, 4): 12, (1, 6): 12})
    g = (0, 5, 8, 10, 2, 7, 6, 3, 11, 4, 1, 9)
    circuits, arcs, _ = implication_system(graph, matching, graph.actions[graph.index[g]])
    assert order(g) == 5
    assert circuits[5].bit_count() == circuits[6].bit_count() == 100
    assert all(1 << i in sink_blocks(len(circuits), arcs) for i in (5, 6))
    for i in (5, 6):
        candidate = matching ^ circuits[i]
        assert sum(n % 2 for n in independent_factor_lengths(graph, candidate)) == 4
    coloured = matching ^ circuits[5] ^ circuits[6]
    assert independent_factor_lengths(graph, coloured) == [28, 28, 36, 184, 192, 192]
    graph.check(coloured)
    print('MINIMAL-BLOCK-BARRIER', json.dumps({
        'patterns': MINIMAL_BLOCK_PATTERNS, 'matching': hex(matching),
        'distinct_translates': len(seen), 'minimal_block_count': sum(counts.values()),
        'translation': g, 'circuit_lengths': [100, 100],
        'individual_odd_counts': [4, 4], 'joint_odd_count': 0,
        'joint_factor_lengths': [28, 28, 36, 184, 192, 192]}), flush=True)


def stress_psl11(samples, descend=False):
    """Bounded centralizer-invariant samples, with no colouring constraints."""
    group = psl2(11)
    x = min(h for h in group if order(h) == 5)
    centralizer = [h for h in group if mul(h, x) == mul(x, h)]
    unseen, representatives = {h for h in group if order(h) == 2}, []
    while unseen:
        a = min(unseen)
        unseen -= {mul(mul(inv(h), a), h) for h in centralizer}
        if len(closure([a, x])) == 660:
            representatives.append(a)
    assert len(representatives) == 6
    totals, total_paths = Counter(), Counter()
    for i, a in enumerate(representatives):
        for power, generator in [(1, x), (2, mul(x, x))]:
            graph = RepairGraph([a, generator], 660)
            subgroup = [h for h in group if mul(h, a) == mul(a, h)]
            status, counts, paths = {}, Counter(), Counter()
            for matching in sample_matchings(
                    graph, samples, random.Random(509), False, subgroup, status):
                odd = sum(n % 2 for n in graph.factor_lengths(matching))
                counts['starts'] += 1
                if not odd:
                    counts['already_even'] += 1
                    paths[0] += 1
                    continue
                current, steps = matching, 0
                while True:
                    witness = audit_blocks(graph, current)
                    if witness is None:
                        counts['sink_block_stall'] += 1
                        witness = audit_blocks(graph, current, all_closed=True)
                        print('SINK-BLOCK-STALL', json.dumps({
                            'pair': i, 'power': power, 'a': a, 'x': generator,
                            'fingerprint': graph.fingerprint, 'matching': hex(current),
                            'factor_lengths': graph.factor_lengths(current),
                            'all_closed_witness': witness}), flush=True)
                        if witness is None:
                            counts['all_closed_stall'] += 1
                            break
                    current = int(witness['matching'], 16)
                    steps += 1
                    if not descend or witness['new_oddness'] == 0:
                        counts['decreased'] += 1
                        paths[steps] += 1
                        break
                if counts['all_closed_stall']:
                    break
            print('BLOCK-STRESS', i, power, dict(counts), 'exhausted',
                  status.get('exhausted', False), 'path_lengths', dict(paths), flush=True)
            totals.update(counts)
            total_paths.update(paths)
    print('BLOCK-STRESS-SUMMARY', dict(totals), 'path_lengths',
          dict(sorted(total_paths.items())), 'descend', descend, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stress-psl11', type=int, default=0, metavar='SAMPLES')
    parser.add_argument('--descend', action='store_true',
                        help='follow strict decreases to a colouring or a closed-set stall')
    args = parser.parse_args()
    if args.stress_psl11 < 0:
        parser.error('--stress-psl11 must be nonnegative')
    check_local_implications()
    check_certificate()
    check_minimal_block_barrier()
    if args.stress_psl11:
        stress_psl11(args.stress_psl11, args.descend)


if __name__ == '__main__':
    main()
