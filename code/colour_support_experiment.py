#!/usr/bin/env python3
"""An exact signed-cycle test for a prescribed quotient colour support.

Run: python3 code/colour_support_experiment.py

The odd generator has order five. The support specifies which a-edges
receive colour 1; the other two colours are determined by signed cycles.
No CDC compatibility or quotient-subflow constraint is imposed.
"""

import random
from itertools import combinations, product

from pysat.solvers import Cadical153

from cdc_palette_experiment import (
    cayley_graph, cycle_perm, petersen, random_generator_separated_flow,
    solve_affine,
)
from cayley_snark_check import mul
from layer_potential_experiment import QuotientForm
from transvection_switch_experiment import (
    LINEAR_MAPS, partial_sums, sample_filtered_flow, switch_system,
)


def local_relations(order, selected):
    """Pairs of nonselected edges whose other-colour bits agree or differ."""
    positions = [i for i, e in enumerate(order) if e in selected]
    if len(positions) == 1:
        i = positions[0]
        return [(order[(i + 2) % 5], order[(i + 3) % 5], 0),
                (order[(i + 1) % 5], order[(i + 4) % 5], 1)]
    if len(positions) == 3 and any(
            set(positions) == {i, (i + 1) % 5, (i + 2) % 5} for i in range(5)):
        other = [e for e in order if e not in selected]
        return [(other[0], other[1], 1)]
    return None


def support_system(model, selected):
    """Return colour bits, their dimension, and (length, sign) of each cycle."""
    variables = [e for e in range(len(model.a_edges)) if e not in selected]
    index = {e: i for i, e in enumerate(variables)}
    relations = []
    for order in model.orders:
        local = local_relations(order, selected)
        if local is None:
            return None, None, None
        relations.extend(local)
    adjacency = {e: [] for e in variables}
    equations = []
    for identifier, (u, v, sign) in enumerate(relations):
        adjacency[u].append((v, identifier, sign))
        adjacency[v].append((u, identifier, sign))
        equations.append((1 << index[u]) ^ (1 << index[v]) ^ (sign << len(variables)))
    assert all(len(incident) == 2 for incident in adjacency.values())
    unseen = set(variables)
    cycles = []
    while unseen:
        start = min(unseen)
        current, previous, length, parity = start, None, 0, 0
        while True:
            unseen.remove(current)
            following, identifier, sign = next(
                entry for entry in adjacency[current] if entry[1] != previous)
            length += 1
            parity ^= sign
            current, previous = following, identifier
            if current == start:
                break
        cycles.append((length, parity))
    assert sum(sign for _, sign in cycles) % 2 == model.q % 2
    solution, kernel = solve_affine(equations, len(variables))
    assert (solution is not None) == all(sign == 0 for _, sign in cycles)
    if solution is None:
        return None, None, cycles
    assert len(kernel) == len(cycles)
    bits = {e: (solution >> index[e]) & 1 for e in variables}
    return bits, len(kernel), cycles


def verify_colouring(model, selected, bits):
    values = [1 if e in selected else 2 ^ bits[e] for e in range(len(model.a_edges))]
    colouring = [0] * len(model.edges)
    for e, value in zip(model.a_edges, values):
        colouring[e] = value
    for order, x_edges in zip(model.orders, model.cycle_edges):
        word = [values[e] for e in order]
        missing = set(range(4)) - partial_sums(word)
        assert len(missing) == 1
        current = missing.pop()
        for p, e in zip(word, x_edges):
            current ^= p
            colouring[e] = current
    assert all(set(colouring[e] for e in incident) == {1, 2, 3}
               for incident in model.inc)
    return colouring


def matching_for_support(model, selected):
    """Build the matching forced by a locally valid support."""
    matching = {model.a_edges[e] for e in selected}
    for order, x_edges in zip(model.orders, model.cycle_edges):
        start = next(i for i, e in enumerate(order) if e in selected)
        offset = 1
        while offset < 5:
            i = (start + offset) % 5
            if order[i] in selected:
                offset += 1
                continue
            assert order[(i + 1) % 5] not in selected
            matching.add(x_edges[i])
            offset += 2
    assert all(sum(e in matching for e in incident) == 1 for incident in model.inc)
    return matching


def complementary_factor(model, selected):
    """Independently build the forced matching and inspect its 2-factor."""
    matching = matching_for_support(model, selected)
    unseen = set(range(model.n))
    cycles = []
    while unseen:
        start = min(unseen)
        current, previous, length, a_edges = start, None, 0, 0
        while True:
            unseen.remove(current)
            e = next(e for e in model.inc[current] if e not in matching and e != previous)
            u, v = model.edges[e]
            following = v if current == u else u
            length += 1
            a_edges += model.kinds[e] == 'a'
            current, previous = following, e
            if current == start:
                break
        cycles.append((a_edges, length % 2))
    return sorted(cycles)


def translated_repair(model, selected, group):
    """Try one alternating-cycle flip against a left translate of the matching.

    This is a finite search, not a guarantee that such a repair exists.
    """
    matching = matching_for_support(model, selected)
    index = {g: i for i, g in enumerate(group)}
    edge_index = {tuple(sorted(edge)): e for e, edge in enumerate(model.edges)}
    tested = 0
    for h in group:
        vertices = [index[mul(h, g)] for g in group]
        translated = {
            edge_index[tuple(sorted((vertices[u], vertices[v])))]
            for e in matching for u, v in [model.edges[e]]
        }
        difference = matching ^ translated
        adjacency = [[] for _ in range(model.n)]
        for e in difference:
            u, v = model.edges[e]
            adjacency[u].append(e)
            adjacency[v].append(e)
        assert all(len(incident) in (0, 2) for incident in adjacency)
        unseen = set(difference)
        while unseen:
            first = min(unseen)
            start = model.edges[first][0]
            current, previous = start, None
            circuit = set()
            while True:
                e = next(e for e in adjacency[current] if e != previous)
                unseen.remove(e)
                circuit.add(e)
                u, v = model.edges[e]
                current, previous = (v if current == u else u), e
                if current == start:
                    break
            assert len(circuit) % 2 == 0
            candidate = matching ^ circuit
            changed = {i for i, e in enumerate(model.a_edges) if e in candidate}
            tested += 1
            bits, _, cycles = support_system(model, changed)
            if bits is not None:
                assert matching_for_support(model, changed) == candidate
                assert complementary_factor(model, changed) == sorted(cycles)
                verify_colouring(model, changed, bits)
                return {'translation': h, 'flip_length': len(circuit),
                        'candidates_tested': tested, 'new_cycles': cycles}
    return None


def sample_quotient_matching(model, rng, forbidden):
    """Sample a 1-factor of Q, with no colouring or circuit-parity constraints."""
    m = len(model.a_edges)
    clauses = []
    for order in model.orders:
        clauses.append([e + 1 for e in order])
        clauses.extend([-(e + 1), -(f + 1)] for e, f in combinations(order, 2))
    for selected in forbidden:
        clauses.append([-(e + 1) if e in selected else e + 1 for e in range(m)])
    with Cadical153(bootstrap_with=clauses) as solver:
        solver.set_phases([v if rng.randrange(2) else -v for v in range(1, m + 1)])
        assert solver.solve()
        assignment = set(solver.get_model())
    selected = {e for e in range(m) if e + 1 in assignment}
    assert all(sum(e in selected for e in order) == 1 for order in model.orders)
    return selected


def check_local_table():
    cases = 0
    for mask in range(32):
        selected = {e for e in range(5) if (mask >> e) & 1}
        relations = local_relations(list(range(5)), selected)
        variables = [e for e in range(5) if e not in selected]
        for assignment in product((0, 1), repeat=len(variables)):
            bits = dict(zip(variables, assignment))
            word = [1 if e in selected else 2 ^ bits[e] for e in range(5)]
            total = 0
            for p in word:
                total ^= p
            admissible = total == 0 and len(partial_sums(word)) < 4
            predicted = relations is not None and all(
                bits[u] ^ bits[v] == sign for u, v, sign in relations)
            assert predicted == admissible
            cases += 1
    assert cases == 243
    print('SUPPORT-LOCAL-TABLE', 'words', cases, flush=True)


def main():
    check_local_table()
    n, edges, kinds, group = cayley_graph(
        [cycle_perm(5, [0, 1]), cycle_perm(5, list(range(5)))],
        return_edge_kinds=True, return_group=True)
    flow, _, _ = random_generator_separated_flow(n, edges, kinds, random.Random(20260904))
    model = QuotientForm(n, edges, kinds, flow)
    old, rng = [], random.Random(982)
    support_outcomes = {}
    support_hits = switch_hits = 0
    for trial in range(12):
        candidate, signature = sample_filtered_flow(model, 1, rng, old)
        old.append(signature)
        current = QuotientForm(n, edges, kinds, candidate)
        selected = {e for e, p in enumerate(current.values) if p == 1}
        bits, dimension, cycles = support_system(current, selected)
        assert complementary_factor(current, selected) == sorted(cycles)
        key = tuple(sorted(selected))
        if key in support_outcomes:
            assert support_outcomes[key] == (bits is not None)
        support_outcomes[key] = bits is not None
        r, _ = switch_system(current, LINEAR_MAPS[0])
        switch_hits += r is not None
        if bits is not None:
            verify_colouring(current, selected, bits)
            support_hits += 1
        assert r is None or bits is not None
        print('SUPPORT-S5', trial, 'consistent', bits is not None,
              'dimension', dimension, 'cycles', cycles, flush=True)
    print('SUPPORT-S5-DONE', 'support_hits', support_hits,
          'subflow_switch_hits', switch_hits, 'trials', len(old),
          'distinct_supports', len(support_outcomes),
          'distinct_support_hits', sum(support_outcomes.values()), flush=True)
    for selected, successful in support_outcomes.items():
        if not successful:
            repair = translated_repair(model, set(selected), group)
            print('TRANSLATED-REPAIR', repair, flush=True)

    matchings, rng = [], random.Random(517)
    direct_hits = repairs = 0
    for trial in range(12):
        selected = sample_quotient_matching(model, rng, matchings)
        matchings.append(selected)
        bits, _, cycles = support_system(model, selected)
        assert complementary_factor(model, selected) == sorted(cycles)
        repair = None
        if bits is not None:
            verify_colouring(model, selected, bits)
            direct_hits += 1
        else:
            repair = translated_repair(model, selected, group)
            repairs += repair is not None
        print('MATCHING-START', trial, 'direct', bits is not None,
              'odd_circuits', sum(sign for _, sign in cycles), 'repair', repair, flush=True)
    print('MATCHING-START-DONE', 'trials', len(matchings),
          'direct_hits', direct_hits, 'one_flip_repairs', repairs, flush=True)

    pn, pedges = petersen()
    pkinds = ['a' if (u < 5) != (v < 5) else 'x' for u, v in pedges]
    pflow, _, _ = random_generator_separated_flow(pn, pedges, pkinds, random.Random(132))
    negative = QuotientForm(pn, pedges, pkinds, pflow)
    valid = 0
    for mask in range(1 << len(negative.a_edges)):
        selected = {e for e in range(len(negative.a_edges)) if (mask >> e) & 1}
        bits, _, cycles = support_system(negative, selected)
        if cycles is not None:
            valid += 1
            assert complementary_factor(negative, selected) == sorted(cycles)
            assert bits is None
            assert sorted(cycles) == [(2, 1), (2, 1)]
    assert valid == 5
    print('SUPPORT-PETERSEN', 'all_supports', 32, 'locally_valid', valid,
          'odd_signed_cycles_each', 2, flush=True)


if __name__ == '__main__':
    main()
