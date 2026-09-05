#!/usr/bin/env python3
"""Check a power-of-two congruence between unsigned and Tait loop counts.

python3 code/loop_congruence_experiment.py
python3 code/loop_congruence_experiment.py --s5

The unsigned count uses only numbers of complementary components. A separate
vertex-set expansion enumerates negative states and their G x C_2 orbits,
where C_2 complements the chosen set of circuits. No descent is used.
"""

import argparse
from collections import Counter
from types import SimpleNamespace

from cayley_snark_check import closure, mul
from cdc_palette_experiment import cayley_graph, cycle_perm, petersen
from colour_support_experiment import support_system
from quotient_face_experiment import QuotientFace, independent_matchings
from translated_repair_audit import RepairGraph, bits, cases


def translate(mask, permutation):
    return sum(1 << permutation[i] for i in bits(mask))


def factor_vertex_sets(graph, matching):
    """Connected vertex sets of X-M, by a set-based traversal."""
    adjacent = [set() for _ in range(graph.n)]
    degrees = [0] * graph.n
    for e, (u, v) in enumerate(graph.edges):
        if not matching & (1 << e):
            adjacent[u].add(v)
            adjacent[v].add(u)
            degrees[u] += 1
            degrees[v] += 1
    assert set(degrees) == {2}
    unseen, components = set(range(graph.n)), []
    while unseen:
        todo, vertices = [min(unseen)], set()
        while todo:
            v = todo.pop()
            if v not in vertices:
                vertices.add(v)
                todo.extend(adjacent[v] - vertices)
        unseen -= vertices
        components.append(sum(1 << v for v in vertices))
    return components


def cayley_vertex_actions(graph):
    return [[graph.index[mul(g, h)] for h in graph.group]
            for g in graph.group]


def matching_orbit_certificates(graph, matchings):
    """Check the equivalent circuit-count/stabiliser bound on every orbit."""
    unseen = set(matchings)
    group_order = len(graph.actions)
    h = group_order & -group_order
    modulus = (2 if h == 2 else 4) * h
    certifying = Counter()
    while unseen:
        matching = min(unseen)
        orbit = {translate(matching, p) for p in graph.actions}
        assert orbit <= unseen
        unseen -= orbit
        assert group_order % len(orbit) == 0
        stabiliser = group_order // len(orbit)
        r = (stabiliser & -stabiliser).bit_length() - 1
        components = factor_vertex_sets(graph, matching)
        c = len(components)
        certified = c <= r + (h >= 4)
        assert certified == (len(orbit) * (1 << c) % modulus != 0)
        if certified:
            assert all(component.bit_count() % 2 == 0
                       for component in components)
            certifying[c] += 1
    return dict(sorted(certifying.items()))


def counts_and_negative_orbits(graph, matchings, vertex_actions, require_free=True):
    """Audit the congruence and every orbit of negative (M, vertex subset) states."""
    matchings = set(matchings)
    group_order = len(vertex_actions)
    h = group_order & -group_order
    assert len(graph.actions) == group_order
    if require_free:
        assert group_order % 2 == 0
        assert all(all(p[v] != v for v in range(graph.n))
                   for p in vertex_actions if p != list(range(graph.n)))
    modulus = (4 if graph.n % 4 == 0 else 2) * h
    unsigned = signed = affine_count = hamiltonian = 0
    negative = set()
    for matching in sorted(matchings):
        components = factor_vertex_sets(graph, matching)
        unsigned += 1 << len(components)
        hamiltonian += len(components) == 1
        # Compare the upstairs components with the signed quotient equations.
        if hasattr(graph, 'orders'):
            selected = {i for i, e in enumerate(graph.a_edges)
                        if matching & (1 << e)}
            solution, dimension, _ = support_system(graph, selected)
            affine_count += 0 if solution is None else 1 << dimension
        chosen = 0
        for i in range(1 << len(components)):
            if i:
                chosen ^= components[(i & -i).bit_length() - 1]
            sign = -1 if chosen.bit_count() % 2 else 1
            signed += sign
            if sign == -1:
                negative.add((matching, chosen))
    if hasattr(graph, 'orders'):
        assert signed == affine_count
    assert unsigned - signed == 2 * len(negative)
    if require_free:
        assert (unsigned - signed) % modulus == 0

    unseen, orbit_sizes, stabilizers = set(negative), Counter(), Counter()
    all_vertices = (1 << graph.n) - 1
    while unseen:
        matching, chosen = min(unseen)
        orbit = set()
        ordinary_orbit = set()
        for edge_action, vertex_action in zip(graph.actions, vertex_actions):
            moved_matching = translate(matching, edge_action)
            moved_chosen = translate(chosen, vertex_action)
            assert moved_matching in matchings
            ordinary_orbit.add((moved_matching, moved_chosen))
            orbit.add((moved_matching, moved_chosen))
            orbit.add((moved_matching, all_vertices ^ moved_chosen))
        assert orbit <= negative
        assert group_order % len(ordinary_orbit) == 0
        assert 2 * group_order % len(orbit) == 0
        if require_free:
            assert (group_order // len(ordinary_orbit)) % 2 == 1
            if graph.n % 4 == 0:
                assert (2 * group_order // len(orbit)) % 2 == 1
            assert len(orbit) % (modulus // 2) == 0
        unseen -= orbit
        orbit_sizes[len(orbit)] += 1
        stabilizers[2 * group_order // len(orbit)] += 1
    return {
        'matchings': len(matchings), 'unsigned': unsigned, 'tait': signed,
        'negative_states': len(negative), 'modulus': modulus,
        'residue': unsigned % modulus,
        'hamiltonian_matchings': hamiltonian,
        'nonhamiltonian_residue': (unsigned - 2 * hamiltonian) % modulus,
        'negative_orbit_sizes': dict(sorted(orbit_sizes.items())),
        'negative_stabilizer_orders': dict(sorted(stabilizers.items())),
    }


def check_faces(names):
    expected = {
        'A5': (125, 860, 540, 16), 'A5_alt': (125, 860, 540, 16),
        'W50': (120, 880, 440, 4), 'F80': (705, 8320, 8320, 64),
        'S5': (26305, 191480, 115960, 32),
        'S5_alt': (26305, 191480, 115960, 32),
    }
    expected_certificates = {
        'A5': {1: 6}, 'A5_alt': {1: 6}, 'W50': {},
        'F80': {2: 6, 3: 8, 4: 6},
        'S5': {1: 92, 2: 55, 3: 8}, 'S5_alt': {1: 92, 2: 55, 3: 8},
    }
    for name in names:
        face = QuotientFace(RepairGraph(*cases()[name]))
        matchings = list(independent_matchings(face))
        result = counts_and_negative_orbits(
            face.graph, matchings, cayley_vertex_actions(face.graph))
        assert tuple(result[key] for key in (
            'matchings', 'unsigned', 'tait', 'modulus')) == expected[name]
        assert result['nonhamiltonian_residue'] == 0
        result['certifying_matching_orbits_by_circuit_count'] = (
            matching_orbit_certificates(face.graph, matchings))
        assert result['certifying_matching_orbits_by_circuit_count'] == (
            expected_certificates[name])
        print('LOOP-CONGRUENCE', name, result, flush=True)
        if name == 'F80':
            # In the full family the residue is zero. A single invariant
            # orbit can still give a certificate without a Hamiltonian state.
            unseen = {m for m in matchings
                      if len(factor_vertex_sets(face.graph, m)) == 2}
            orbit_count = 0
            while unseen:
                matching = min(unseen)
                orbit = {translate(matching, p) for p in face.graph.actions}
                assert orbit <= unseen
                unseen -= orbit
                restricted = counts_and_negative_orbits(
                    face.graph, orbit, cayley_vertex_actions(face.graph))
                assert restricted['matchings'] == 40
                assert restricted['unsigned'] == restricted['tait'] == 160
                assert restricted['residue'] == 32
                assert restricted['hamiltonian_matchings'] == 0
                orbit_count += 1
            assert orbit_count == 6
            print('LOOP-CONGRUENCE-NONHAMILTONIAN-ORBIT', name,
                  'six_orbits_with_same_counts', restricted, flush=True)


def generic_cayley(generators):
    n, edges, kinds, group = cayley_graph(
        generators, return_edge_kinds=True, return_group=True)
    graph = SimpleNamespace(n=n, edges=edges, group=group,
                            index={g: i for i, g in enumerate(group)})
    vertex_actions = cayley_vertex_actions(graph)
    edge_index = {tuple(sorted(edge)): i for i, edge in enumerate(edges)}
    graph.actions = [[edge_index[tuple(sorted((p[u], p[v])))] for u, v in edges]
                     for p in vertex_actions]
    all_a = sum(1 << e for e, kind in enumerate(kinds) if kind == 'a')
    return graph, vertex_actions, all_a


def check_sharpness():
    # The singleton family consisting of all a-edges has respectively two
    # and four odd triangles. Thus the universal moduli 4 and 16 are sharp.
    examples = [
        ('S3', [cycle_perm(3, [0, 1]), cycle_perm(3, [0, 1, 2])], 6, 4),
        ('A4', [mul(cycle_perm(4, [0, 1]), cycle_perm(4, [2, 3])),
                cycle_perm(4, [0, 1, 2])], 12, 16),
    ]
    for name, generators, n, expected in examples:
        graph, actions, matching = generic_cayley(generators)
        assert graph.n == n
        result = counts_and_negative_orbits(graph, [matching], actions)
        assert result['unsigned'] == result['modulus'] == expected
        assert result['tait'] == 0
        print('LOOP-CONGRUENCE-SHARP', name, result, flush=True)


def check_nonfree_control():
    n, edges = petersen()
    edges = sorted(tuple(sorted(edge)) for edge in edges)
    rotation = tuple(5 * (i // 5) + (i + 1) % 5 for i in range(10))
    swap = tuple(5 * (1 - i // 5) + 2 * i % 5 for i in range(10))
    vertex_actions = sorted(closure([rotation, swap]))
    assert len(vertex_actions) == 20
    assert any(p != tuple(range(10)) and any(p[i] == i for i in range(10))
               for p in vertex_actions)
    edge_index = {edge: i for i, edge in enumerate(edges)}
    actions = [[edge_index[tuple(sorted((p[u], p[v])))] for u, v in edges]
               for p in vertex_actions]
    graph = SimpleNamespace(n=n, edges=edges, actions=actions)
    matchings = []
    for chosen in range(5):
        matching = 1 << edge_index[chosen, chosen + 5]
        for base, step in ((0, 1), (5, 2)):
            for offset in (1, 3):
                edge = tuple(sorted((base + (chosen + offset * step) % 5,
                                     base + (chosen + (offset + 1) * step) % 5)))
                matching |= 1 << edge_index[edge]
        matchings.append(matching)
    result = counts_and_negative_orbits(
        graph, matchings, vertex_actions, require_free=False)
    assert (result['unsigned'], result['tait']) == (20, 0)
    assert (result['unsigned'] - result['tait']) % (2 * 4) != 0
    print('LOOP-CONGRUENCE-NONFREE-PETERSEN', result, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--s5', action='store_true')
    args = parser.parse_args()
    check_sharpness()
    check_nonfree_control()
    names = ['A5', 'A5_alt', 'W50', 'F80']
    if args.s5:
        names += ['S5', 'S5_alt']
    check_faces(names)


if __name__ == '__main__':
    main()
