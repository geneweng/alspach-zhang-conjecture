#!/usr/bin/env python3
"""Certify the odd-cycle/sparse-involution proposition in the survey.

Let x=(0,1,...,N-1) and let a be a product of three disjoint
or five disjoint transpositions.  Once N >= 4r+1, a fixes enough points
that two are consecutive on the x-cycle, and Proposition 4.6 of the survey
applies.  This script checks the remaining odd degrees, up to the dihedral
normalizer of <x> that is relevant here (rotations and reversal).

For N=9 and 11 it colours the point-stabilizer quotient, an N-cycle with
three chords and N-6 semi-edges.  The only failed orbit for N=9 is the
Petersen-minus-a-vertex quotient already displayed in the survey; exact
Schreier recursion gives |<a,x>|=162.  For N=7 the point quotient has only
one semi-edge and cannot be coloured.  Instead we colour the 21-vertex
two-set quotient.  Its sole failed orbit is the reflection of the 7-cycle,
and |<a,x>|=14.

For five transpositions, only N=11,13,15,17,19 remain after the consecutive-
fixed-point argument.  The point quotient colours every configuration with
nonconsecutive fixed points except all 513 types at N=11 and 46 types at
N=13.  Two-set quotients colour all of these except the N=11 reflection,
which generates the dihedral group of order 22.

The script also checks the first seven-transposition boundary N=15.  After
normalizing its unique fixed point, 68,219 types remain modulo reflection.
The two-set quotient colours 68,144; exact Schreier recursion shows that all
75 failures generate proper subgroups, so none is relevant to S_15.
The next boundary N=17 is checked by ``septuple_n17_experiment.py``; the
remaining boundaries N=19,21,23,25,27 are checked by
``septuple_later_experiment.py``.  Both import the quotient and direct-
verification routines from this module.

Every positive SAT witness is checked directly as a pregraph edge-colouring.
The exceptional group orders are computed by an exact implementation of
Schreier's lemma, not by enumerating the symmetric groups.
"""

from itertools import combinations
from math import factorial

from pysat.solvers import Cadical153


def mul(p, q):
    """Apply p and then q."""
    return tuple(q[i] for i in p)


def inv(p):
    result = [0] * len(p)
    for i, j in enumerate(p):
        result[j] = i
    return tuple(result)


def cycle(n):
    return tuple((i + 1) % n for i in range(n))


def matching_permutation(n, matching):
    result = list(range(n))
    for u, v in matching:
        result[u], result[v] = v, u
    return tuple(result)


def group_order(generators):
    """Compute a permutation-group order by exact Schreier recursion."""
    n = len(generators[0])
    identity = tuple(range(n))

    def stabilizer_order(gens, fixed):
        gens = list(dict.fromkeys(g for g in gens if g != identity))
        if not gens:
            return 1
        base = next((i for i in range(n)
                     if i not in fixed and any(g[i] != i for g in gens)), None)
        if base is None:
            return 1

        # reps[alpha] maps base to alpha.  Schreier's generators
        # reps[alpha] s reps[alpha^s]^-1 generate the point stabilizer.
        reps = {base: identity}
        todo = [base]
        while todo:
            alpha = todo.pop()
            for generator in gens:
                beta = generator[alpha]
                if beta not in reps:
                    reps[beta] = mul(reps[alpha], generator)
                    todo.append(beta)

        schreier, seen = [], set()
        for alpha, representative in reps.items():
            for generator in gens:
                beta = generator[alpha]
                stabilizer_generator = mul(
                    mul(representative, generator), inv(reps[beta]))
                assert stabilizer_generator[base] == base
                if stabilizer_generator != identity and stabilizer_generator not in seen:
                    seen.add(stabilizer_generator)
                    schreier.append(stabilizer_generator)
        return len(reps) * stabilizer_order(schreier, fixed | {base})

    return stabilizer_order(generators, set())


def perfect_matchings(vertices):
    """All perfect matchings on a sorted tuple of vertices."""
    if not vertices:
        yield ()
        return
    first = vertices[0]
    for j in range(1, len(vertices)):
        second = vertices[j]
        remaining = vertices[1:j] + vertices[j + 1:]
        for rest in perfect_matchings(remaining):
            yield ((first, second),) + rest


def matchings(n, r):
    """All matchings of size r on n labelled points, without repetition."""
    count = 0
    for support in combinations(range(n), 2 * r):
        for matching in perfect_matchings(support):
            count += 1
            yield matching
    expected = factorial(n) // (2 ** r * factorial(r) * factorial(n - 2 * r))
    assert count == expected


def dihedral_image(n, matching, sign, shift):
    return tuple(sorted(tuple(sorted(((sign * u + shift) % n,
                                      (sign * v + shift) % n)))
                        for u, v in matching))


def canonical_dihedral(n, matching):
    return min(dihedral_image(n, matching, sign, shift)
               for sign in (1, -1) for shift in range(n))


def fixed_points_are_independent(n, matching):
    endpoints = {point for edge in matching for point in edge}
    fixed = set(range(n)) - endpoints
    return all((point + 1) % n not in fixed for point in fixed)


def orbit_representatives(n, r, independent_fixed_only=False):
    if not independent_fixed_only:
        source = matchings(n, r)
    else:
        def independent_source():
            for fixed_tuple in combinations(range(n), n - 2 * r):
                fixed = set(fixed_tuple)
                if any((point + 1) % n in fixed for point in fixed):
                    continue
                endpoints = tuple(point for point in range(n) if point not in fixed)
                yield from perfect_matchings(endpoints)
        source = independent_source()
    representatives = sorted(set(canonical_dihedral(n, matching)
                                 for matching in source))
    assert all(not independent_fixed_only or
               fixed_points_are_independent(n, matching)
               for matching in representatives)
    return representatives


def one_fixed_point_representatives(n):
    """Near-perfect matchings modulo the dihedral group, fixing point 0.

    Rotation first takes the unique fixed point to 0.  The residual dihedral
    stabilizer of 0 consists of the identity and i -> -i.
    """
    representatives = set()
    for matching in perfect_matchings(tuple(range(1, n))):
        reflected = tuple(sorted(tuple(sorted(((-u) % n, (-v) % n)))
                                 for u, v in matching))
        representatives.add(min(matching, reflected))
    return sorted(representatives)


def subset_pregraph(n, matching, k):
    """Quotient by a k-set stabilizer, on the k-subsets of {0,...,n-1}."""
    vertices = list(combinations(range(n), k))
    index = {vertex: i for i, vertex in enumerate(vertices)}
    permutations = (matching_permutation(n, matching), cycle(n))
    items, seen = [], set()
    for i, vertex in enumerate(vertices):
        for kind, permutation in enumerate(permutations):
            image_vertex = tuple(sorted(permutation[point] for point in vertex))
            image = index[image_vertex]
            key = kind, min(i, image), max(i, image)
            if key in seen:
                continue
            seen.add(key)
            items.append((i, None if i == image else image))
    return len(vertices), items


def point_pregraph(n, matching):
    vertices, items = subset_pregraph(n, matching, 1)
    assert vertices == n
    return items


def verify_colouring(n, items, colours):
    at = [set() for _ in range(n)]
    for item, ((u, v), colour) in enumerate(zip(items, colours)):
        assert colour in (0, 1, 2), (item, colour)
        assert colour not in at[u], (item, u, colour)
        at[u].add(colour)
        if v is not None:
            assert colour not in at[v], (item, v, colour)
            at[v].add(colour)
    assert all(colours_at_vertex == {0, 1, 2}
               for colours_at_vertex in at)


def colour_pregraph(n, items):
    """Find and independently verify a Tait colouring of a cubic pregraph."""
    incident = [[] for _ in range(n)]
    for item, (u, v) in enumerate(items):
        incident[u].append(item)
        if v is not None:
            incident[v].append(item)
    assert all(len(row) == 3 for row in incident)

    # m(e) says that e has colour 2.  Exactly one incident item has colour 2;
    # the Boolean c(e) distinguishes colours 0 and 1 on the other two.
    matching_var = lambda item: item + 1
    colour_var = lambda item: len(items) + item + 1
    solver = Cadical153()
    for row in incident:
        solver.add_clause([matching_var(item) for item in row])
        for first, second in combinations(row, 2):
            solver.add_clause([-matching_var(first), -matching_var(second)])
            solver.add_clause([matching_var(first), matching_var(second),
                               colour_var(first), colour_var(second)])
            solver.add_clause([matching_var(first), matching_var(second),
                               -colour_var(first), -colour_var(second)])
    satisfiable = solver.solve()
    if not satisfiable:
        solver.delete()
        return None
    model = set(literal for literal in solver.get_model() if literal > 0)
    solver.delete()
    colours = [2 if matching_var(item) in model else
               (1 if colour_var(item) in model else 0)
               for item in range(len(items))]
    verify_colouring(n, items, colours)
    return colours


def check():
    triple_orbits = {7: 11, 9: 85, 11: 350}
    triple_exception = {
        7: ((0, 1), (2, 6), (3, 5)),
        9: ((0, 4), (1, 6), (3, 7)),
    }
    triple_exception_order = {7: 14, 9: 162}

    for n in (7, 9, 11):
        representatives = orbit_representatives(n, 3)
        assert len(representatives) == triple_orbits[n]
        point_colours = {}
        for matching in representatives:
            items = point_pregraph(n, matching)
            point_colours[matching] = colour_pregraph(n, items)

        if n == 7:
            # One semi-edge on an odd number of vertices violates the parity
            # condition, and the SAT calculation independently sees this.
            assert all(colouring is None for colouring in point_colours.values())
            twoset_colours = {}
            for matching in representatives:
                vertices, items = subset_pregraph(n, matching, 2)
                assert vertices == 21 and len(items) == 33
                assert sum(v is None for _, v in items) == 3
                twoset_colours[matching] = colour_pregraph(vertices, items)
            failed = [matching for matching, colouring in twoset_colours.items()
                      if colouring is None]
            assert failed == [triple_exception[n]]
            coloured = len(twoset_colours) - len(failed)
            quotient = 'two-set'
        else:
            failed = [matching for matching, colouring in point_colours.items()
                      if colouring is None]
            assert failed == ([] if n == 11 else [triple_exception[n]])
            coloured = len(point_colours) - len(failed)
            quotient = 'point-stabilizer'

        failed_orders = {}
        for matching in failed:
            order = group_order([matching_permutation(n, matching), cycle(n)])
            failed_orders[order] = failed_orders.get(order, 0) + 1
            assert order == triple_exception_order[n]
        print('TRIPLE-TRANSPOSITION', 'N', n,
              'dihedral_orbits', len(representatives),
              quotient + '_colourable', coloured,
              'failed_group_orders', failed_orders)

    quintuple_orbits = {11: 513, 13: 5832, 15: 12152, 17: 5832, 19: 513}
    quintuple_reflection = ((0, 1), (2, 10), (3, 9), (4, 8), (5, 7))
    for n in (11, 13, 15, 17, 19):
        # Configurations with consecutive fixed points are already covered by
        # Proposition 4.6, so enumerate only the remaining configurations.
        representatives = orbit_representatives(n, 5, independent_fixed_only=True)
        assert len(representatives) == quintuple_orbits[n]
        point_failed = []
        for matching in representatives:
            if colour_pregraph(n, point_pregraph(n, matching)) is None:
                point_failed.append(matching)
        expected_point_failures = {11: 513, 13: 46, 15: 0, 17: 0, 19: 0}
        assert len(point_failed) == expected_point_failures[n]

        twoset_failed = []
        for matching in point_failed:
            vertices, items = subset_pregraph(n, matching, 2)
            if n == 11:
                assert vertices == 55 and len(items) == 85
                assert sum(v is None for _, v in items) == 5
            if n == 13:
                assert vertices == 78 and len(items) == 121
                assert sum(v is None for _, v in items) == 8
            if colour_pregraph(vertices, items) is None:
                twoset_failed.append(matching)

        if n == 11:
            assert twoset_failed == [quintuple_reflection]
            failed_order = group_order([
                matching_permutation(n, quintuple_reflection), cycle(n)])
            assert failed_order == 22
            failed_orders = {22: 1}
        else:
            assert not twoset_failed
            failed_orders = {}
        print('QUINTUPLE-TRANSPOSITION', 'N', n,
              'nonconsecutive-fixed-point_orbits', len(representatives),
              'point_colourable', len(representatives) - len(point_failed),
              'rescued_by_two-set', len(point_failed) - len(twoset_failed),
              'failed_group_orders', failed_orders)

    representatives = one_fixed_point_representatives(15)
    assert len(representatives) == 68219
    twoset_failed = []
    for matching in representatives:
        vertices, items = subset_pregraph(15, matching, 2)
        assert vertices == 105 and len(items) == 161
        assert sum(v is None for _, v in items) == 7
        if colour_pregraph(vertices, items) is None:
            twoset_failed.append(matching)
    assert len(twoset_failed) == 75
    failed_orders = {}
    for matching in twoset_failed:
        order = group_order([matching_permutation(15, matching), cycle(15)])
        failed_orders[order] = failed_orders.get(order, 0) + 1
    assert failed_orders == {
        30: 1, 360: 2, 750: 2, 2430: 4, 3000: 5,
        29160: 10, 38880: 21, 466560: 30,
    }
    print('SEPTUPLE-TRANSPOSITION', 'N', 15,
          'dihedral_orbits', len(representatives),
          'two-set_colourable', len(representatives) - len(twoset_failed),
          'failed_group_orders', failed_orders)

    # Sanity controls for the exact group-order routine.
    for n in range(3, 9):
        assert group_order([cycle(n), matching_permutation(n, ((0, 1),))]) == factorial(n)
    print('SPARSE-INVOLUTION', 'r=3,N>=13 or r=5,N>=21',
          'consecutive_fixed_points_by_pigeonhole', True)


if __name__ == '__main__':
    check()
