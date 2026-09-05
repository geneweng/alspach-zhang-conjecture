#!/usr/bin/env python3
"""Test the affine quotient-switch criterion for odd generators of order five.

Run: python3 code/transvection_switch_experiment.py

This constructs Tait flows from a quotient flow and binary switches, checks
the affine layer equation against the original CDC system, and uses the
Petersen graph as a negative control.
"""

import random
from collections import Counter
from itertools import product

from pysat.solvers import Cadical153

from cdc_palette_experiment import (
    affine_solutions, all_transvection_matchings, binary_basis, binary_rank,
    cayley_graph, cycle_perm, find_separated_transvection_palette, palette,
    petersen, random_generator_separated_flow, solve_affine, xor_selected,
)
from layer_potential_experiment import QuotientForm, improve_lifts


LINEAR_MAPS = [L for L, offset, _ in all_transvection_matchings() if offset == 0]


def switch_values(values, layer, linear):
    return [linear[p] if (layer >> i) & 1 else p for i, p in enumerate(values)]


def partial_sums(values):
    current, visited = 0, set()
    for p in values:
        current ^= p
        visited.add(current)
    assert current == 0
    return visited


def local_switches(values, linear):
    return [r for r in range(1 << len(values))
            if xor_selected(values, r) == 0
            and len(partial_sums(switch_values(values, r, linear))) < 4]


def affine_equations(points, dimension):
    """Recover and verify an affine space from its complete point set."""
    if not points:
        return [1 << dimension]
    basis = binary_basis([p ^ points[0] for p in points])
    assert len(points) == 1 << len(basis)
    equations = [normal | (((normal & points[0]).bit_count() & 1) << dimension)
                 for normal in range(1, 1 << dimension)
                 if all((normal & b).bit_count() % 2 == 0 for b in basis)]
    start, kernel = solve_affine(equations, dimension)
    assert set(affine_solutions(start, kernel)) == set(points)
    return equations


def switch_system(model, linear):
    """An exact linear system for r in S_phi with admissible switched words."""
    m = len(model.a_edges)
    equations = []
    for order in model.orders:
        assert len(order) == 5
        values = [model.values[e] for e in order]
        allowed = local_switches(values, linear)
        for row in affine_equations(allowed, 5):
            coefficient = sum(((row >> i) & 1) << e for i, e in enumerate(order))
            equations.append(coefficient | (((row >> 5) & 1) << m))
    return solve_affine(equations, m)


def tait_from_switch(model, linear, r):
    """Construct and independently verify the resulting nowhere-zero 4-flow."""
    result = [0] * len(model.edges)
    changed = switch_values(model.values, r, linear)
    for e, value in zip(model.a_edges, changed):
        result[e] = value
    for order, x_edges in zip(model.orders, model.cycle_edges):
        values = [changed[e] for e in order]
        missing = set(range(4)) - partial_sums(values)
        assert len(missing) == 1
        current = missing.pop()
        for p, e in zip(values, x_edges):
            current ^= p
            result[e] = current
    assert all(result)
    assert all(set(result[e] for e in incident) == {1, 2, 3}
               for incident in model.inc)
    return result


def characteristic_layers(model, corrections):
    matrix = model.matrix(corrections)
    rows = [row | (((row >> i) & 1) << model.d)
            for i, row in enumerate(matrix)]
    start, kernel = solve_affine(rows, model.d)
    assert start is not None
    return start, kernel


def check_characteristic_layers(model, corrections):
    """Compare the affine first-coordinate image with the full CDC equations."""
    n = model.n
    flow = model.flow(corrections)
    local_g = [{} for _ in range(n)]
    for vertices, order, x_edges in zip(model.cycles, model.orders, model.cycle_edges):
        for i, v in enumerate(vertices):
            local_g[v] = {model.a_edges[order[i]]: 0,
                          x_edges[i - 1]: 2 * model.values[order[i]], x_edges[i]: 0}
    nvars = 3 * n + len(model.edges)
    rows = []
    for e, (u, v) in enumerate(model.edges):
        rhs = local_g[u][e] ^ local_g[v][e]
        for bit in range(3):
            row = (1 << (3 * u + bit)) ^ (1 << (3 * v + bit))
            if (flow[e] >> bit) & 1:
                row ^= 1 << (3 * n + e)
            row |= ((rhs >> bit) & 1) << nvars
            rows.append(row)
    solution, homogeneous = solve_affine(rows, nvars)
    assert solution is not None
    palette(n, model.edges, flow, local_g, solution)

    def project(vector):
        r = sum(((vector >> (3 * model.edges[e][0])) & 1) << pos
                for pos, e in enumerate(model.a_edges))
        return model.coordinates(r)

    start, radical = characteristic_layers(model, corrections)
    projected = binary_basis([project(v) for v in homogeneous])
    assert len(projected) == len(radical)
    assert binary_rank(projected + radical) == len(radical)
    assert binary_rank(radical + [project(solution) ^ start]) == len(radical)
    return start, radical


def check_local_switch_table():
    table = {
        (1, 1, 1, 2, 3): {1: None, 2: (3, 0), 3: (6, 0)},
        (1, 1, 2, 1, 3): {1: False, 2: (10, 1), 3: (9, 1)},
    }
    for values, columns in table.items():
        subflows = [r for r in range(32) if xor_selected(values, r) == 0]
        for linear in LINEAR_MAPS:
            fixed = next(p for p in (1, 2, 3) if linear[p] == p)
            condition = columns[fixed]
            if condition is False:
                expected = []
            elif condition is None:
                expected = subflows
            else:
                mask, rhs = condition
                expected = [r for r in subflows if (r & mask).bit_count() % 2 == rhs]
            assert local_switches(values, linear) == expected
    counts = Counter()
    for values in product((1, 2, 3), repeat=5):
        if xor_selected(values, 31):
            continue
        for linear in LINEAR_MAPS:
            allowed = local_switches(values, linear)
            affine_equations(allowed, 5)
            counts[len(allowed)] += 1
    assert counts == {0: 30, 4: 120, 8: 30}
    print('LOCAL-SWITCH-TABLE', 'word_map_cases', sum(counts.values()),
          'allowed_counts', dict(sorted(counts.items())), flush=True)


def repeated_obstructions(model):
    bad = Counter()
    for order in model.orders:
        values = [model.values[e] for e in order]
        repeated = next(p for p in (1, 2, 3) if values.count(p) == 3)
        positions = {i for i, p in enumerate(values) if p == repeated}
        consecutive = any(positions == {i, (i + 1) % 5, (i + 2) % 5}
                          for i in range(5))
        if not consecutive:
            bad[repeated] += 1
    return bad


def sample_filtered_flow(model, fixed_value, rng, forbidden=()):
    """Choose phi using only flow and local repeated-value restrictions."""
    m = len(model.a_edges)
    clauses = [[2 * e + 1, 2 * e + 2] for e in range(m)]
    words = []
    for values in product((1, 2, 3), repeat=5):
        if xor_selected(values, 31):
            continue
        positions = {i for i, p in enumerate(values) if p == fixed_value}
        if len(positions) == 3 and not any(
                positions == {i, (i + 1) % 5, (i + 2) % 5} for i in range(5)):
            continue
        words.append(values)
    assert len(words) == 50
    next_var = 2 * m + 1
    for order in model.orders:
        states = []
        for values in words:
            state = next_var
            next_var += 1
            states.append(state)
            for e, p in zip(order, values):
                for bit in range(2):
                    literal = 2 * e + bit + 1
                    clauses.append([-state, literal if (p >> bit) & 1 else -literal])
        clauses.append(states)
    for old in forbidden:
        clauses.append([-(2 * e + bit + 1) if (p >> bit) & 1 else 2 * e + bit + 1
                        for e, p in enumerate(old) for bit in range(2)])
    with Cadical153(bootstrap_with=clauses) as solver:
        solver.set_phases([v if rng.randrange(2) else -v for v in range(1, 2 * m + 1)])
        assert solver.solve()
        assignment = set(solver.get_model())
    values = [sum((1 << bit) for bit in range(2) if 2 * e + bit + 1 in assignment)
              for e in range(m)]
    flow = [0] * len(model.edges)
    for e, value in zip(model.a_edges, values):
        flow[e] = value << 1
    for order, x_edges in zip(model.orders, model.cycle_edges):
        y = 0
        for pos, e in zip(order, x_edges):
            y ^= values[pos]
            flow[e] = 1 | (y << 1)
        assert y == 0
    return flow, tuple(values)


def main():
    check_local_switch_table()
    a = cycle_perm(5, [0, 1])
    x = cycle_perm(5, list(range(5)))
    n, edges, kinds = cayley_graph([a, x], return_edge_kinds=True)
    flow, _, _ = random_generator_separated_flow(n, edges, kinds, random.Random(20260904))
    model = QuotientForm(n, edges, kinds, flow)
    rank, corrections = improve_lifts(model)
    start, radical = check_characteristic_layers(model, corrections)
    bad = repeated_obstructions(model)
    assert set(bad) == {1, 2, 3}
    assert all(switch_system(model, L)[0] is None for L in LINEAR_MAPS)
    print('ALL-LIFTS-OBSTRUCTED', 'form_rank', rank, 'layer_assignments', 1 << len(radical),
          'bad_vertices_by_fixed_value', dict(sorted(bad.items())),
          'excluded_lifts', 1 << (2 * model.q), flush=True)

    successful, _, _ = find_separated_transvection_palette('SWITCH-CONTROL-S5', n, edges, kinds)
    witness = QuotientForm(n, edges, kinds, successful)
    witness_corrections = sum((successful[x_edges[-1]] >> 1) << (2 * C)
                              for C, x_edges in enumerate(witness.cycle_edges))
    assert witness.flow(witness_corrections) == successful
    check_characteristic_layers(witness, witness_corrections)
    r, _ = switch_system(witness, LINEAR_MAPS[0])
    assert r is not None
    tait_from_switch(witness, LINEAR_MAPS[0], r)
    print('SWITCH-POSITIVE-CONTROL', 'verified', True, flush=True)

    rng = random.Random(982)
    old = []
    successes = 0
    for trial in range(12):
        candidate, signature = sample_filtered_flow(model, 1, rng, old)
        old.append(signature)
        current = QuotientForm(n, edges, kinds, candidate)
        assert repeated_obstructions(current)[1] == 0
        r, kernel = switch_system(current, LINEAR_MAPS[0])
        if r is not None:
            tait_from_switch(current, LINEAR_MAPS[0], r)
            successes += 1
        print('FILTERED-SWITCH', trial, 'switch_consistent', r is not None,
              'solution_dim', len(kernel) if r is not None else None, flush=True)
    print('FILTERED-SWITCH-DONE', 'successes', successes, 'trials', len(old), flush=True)

    pn, pedges = petersen()
    pkinds = ['a' if (u < 5) != (v < 5) else 'x' for u, v in pedges]
    pflow, _, _ = random_generator_separated_flow(pn, pedges, pkinds, random.Random(132))
    negative = QuotientForm(pn, pedges, pkinds, pflow)
    check_characteristic_layers(negative, 0)
    old = []
    for _ in range(12):
        pflow, signature = sample_filtered_flow(negative, 1, rng, old)
        old.append(signature)
        current = QuotientForm(pn, pedges, pkinds, pflow)
        assert repeated_obstructions(current)[1] == 0
        assert all(switch_system(current, L)[0] is None for L in LINEAR_MAPS)
    print('SWITCH-NEGATIVE-CONTROL', 'filtered_flows', len(old), 'hits', 0, flush=True)


if __name__ == '__main__':
    main()
