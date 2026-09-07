"""Target-free seven-change switches with at most two temporary-colour edges.

One temporary edge cannot change the three-colour Kempe class of a cubic
graph: completing it with the common missing colour retracts every allowed
move to a three-colour Kempe move (or no change). Two temporary edges do
suffice in the rigid q=256 control, by a seven-change, target-free rule.

The local rule has an exact endpoint-pairing return test. Return to three
colours does NOT imply non-unit root gain; both conditions are checked.
No uniform success or two-buffer connectivity theorem is asserted.
"""

import argparse
import json
from collections import Counter
from itertools import combinations, permutations
from pathlib import Path

from borel_exchange_stress import LogField, canonical_gain
from borel_folded_experiment import FoldedDiagram, from_point
from borel_kempe_exchange import (audit_buffer_sequence, change, closed_graph,
                                 eligible_certificates, normalized)
from borel_reversal_experiment import PointDiagram
from borel_three_exchange import audit_matching


def missing(graph, colours, vertex):
    return set(range(3)) - {colours[e] for e in graph.inc[vertex]}


def endpoints(graph, edges):
    counts = Counter(v for e in edges for v in graph.edges[e][:2])
    return {v for v, degree in counts.items() if degree == 1}


def seven_switches(graph, initial, pivot, a=0):
    """All two path choices of the seven-change move at a specified pivot.

    The pivot uv has colour b != a. Its incident a-edges are uw and vz.
    Temporarily recolour uw,vz, then recolour uv to a; swap the (a,c)-path
    from w to z. There are two (b,c)-paths on the four deficient endpoints.
    A swap followed by restoring uw,vz succeeds iff u's path ends != w.
    """
    assert graph.proper(initial) and set(initial) == {0, 1, 2}
    u, v, _ = graph.edges[pivot]
    b = initial[pivot]
    assert a != b
    c = next(iter({0, 1, 2} - {a, b}))
    e1 = next(e for e in graph.inc[u] if initial[e] == a)
    e2 = next(e for e in graph.inc[v] if initial[e] == a)
    if e1 == e2:
        return {'status': 'PARALLEL_SEED', 'candidates': []}
    w = next(x for x in graph.edges[e1][:2] if x != u)
    z = next(x for x in graph.edges[e2][:2] if x != v)
    assert len({u, v, w, z}) == 4
    colours, prefix = initial.copy(), []

    def perform(state, moves, pair, edges):
        assert graph.component(state, pair, min(edges)) == edges
        change(state, pair, edges)
        assert graph.proper(state) and state.count(3) <= 2
        moves.append((pair, frozenset(edges)))

    perform(colours, prefix, (a, 3), {e1})
    perform(colours, prefix, (a, 3), {e2})
    perform(colours, prefix, (a, b), {pivot})
    seed = next(e for e in graph.inc[w] if colours[e] in (a, c))
    ac_path = graph.component(colours, (a, c), seed)
    assert endpoints(graph, ac_path) == {w, z}
    perform(colours, prefix, (a, c), ac_path)
    assert missing(graph, colours, u) == missing(graph, colours, v) == {b}
    assert missing(graph, colours, w) == missing(graph, colours, z) == {c}
    paths = [comp for comp in graph.components(colours, (b, c)) if endpoints(graph, comp)]
    assert len(paths) == 2
    pairs = [endpoints(graph, path) for path in paths]
    assert set.union(*pairs) == {u, v, w, z}
    blocked = {u, w} in pairs
    candidates = []
    for path in paths:
        trial, moves = colours.copy(), prefix.copy()
        perform(trial, moves, (b, c), path)
        replacements = [missing(graph, trial, left) & missing(graph, trial, right)
                        for left, right in ((u, w), (v, z))]
        assert all(len(s) == 1 for s in replacements) == (not blocked)
        if blocked:
            continue
        for edge, replacement in zip((e1, e2), replacements):
            perform(trial, moves, (3, next(iter(replacement))), {edge})
        assert 3 not in trial and len(moves) == 7
        candidates.append({'colours': trial, 'moves': moves,
                           'changed_unlabelled_colouring': normalized(trial) != normalized(initial)})
    return {'status': 'BLOCKED_ENDPOINT_PAIRING' if blocked else 'RETURNS',
            'pivot_edge_id': pivot, 'pivot_vertices': [u + 1, v + 1],
            'temporary_edge_ids': [e1, e2], 'four_vertices_uvwz': [x + 1 for x in (u, v, w, z)],
            'ordinary_colours_abc': [a, b, c],
            'ac_path_length': len(ac_path),
            'bc_path_lengths': [len(path) for path in paths],
            'bc_endpoint_pairs': [sorted(x + 1 for x in pair) for pair in pairs],
            'candidates': candidates}


def discover(point, folded, exhaustive=False):
    """Search from the natural colouring, with no supplied target or random seed.

    Only the original chord colour is used as a, and every path edge (including
    e*) is tried as pivot. Failure exhausts this specified seven-change family,
    not all two-buffer walks or all three-colourings.
    """
    graph, initial = closed_graph(folded)
    counts, successes, first = Counter(), [], None
    for pivot, colour in enumerate(initial):
        if colour == 0:
            continue
        counts['pivots_tested'] += 1
        result = seven_switches(graph, initial, pivot)
        counts[result['status']] += 1
        for candidate in result['candidates']:
            counts['three_colour_returns'] += 1
            counts['changed_unlabelled_returns'] += int(candidate['changed_unlabelled_colouring'])
            eligible = eligible_certificates(point, folded, candidate['colours'])
            good = [item for item in eligible if item['certificate']['colours']]
            counts['returns_with_a_good_lift'] += int(bool(good))
            counts['eligible_good_matchings'] += len(good)
            counts['nontrivial_returns_with_both_gains_unit'] += int(
                candidate['changed_unlabelled_colouring'] and not good)
            if not good:
                continue
            if first is None:
                chosen = good[0]
                checked, profile = audit_matching(point, folded, frozenset(chosen['selected']),
                                                  permutations=True)
                assert checked == chosen['certificate']
                replay = audit_buffer_sequence(graph, initial, candidate['colours'],
                                               candidate['moves'], closed=True)
                assert replay['moves'] == 7 and replay['max_edges_of_buffer_colour'] == 2
                first = {k: v for k, v in result.items() if k != 'candidates'}
                first.update({'eligible_root_gains': [item['certificate']['root_gain'] for item in eligible],
                              'selected': chosen['selected'], 'certificate': checked,
                              'point_profile': profile, 'seven_change_certificate': replay})
            successes.append(pivot)
            if not exhaustive:
                return {'status': 'REPAIRED', 'counts': dict(counts), 'first': first}
    return {'status': 'REPAIRED' if first else 'FAMILY_EXHAUSTED',
            'counts': dict(counts), 'successful_pivots': sorted(set(successes)), 'first': first}


def complete_one_buffer(graph, colours):
    temporary = [e for e, c in enumerate(colours) if c == 3]
    assert len(temporary) <= 1 and graph.proper(colours)
    completed = list(colours)
    if temporary:
        edge = temporary[0]
        u, v, _ = graph.edges[edge]
        assert missing(graph, colours, u) == missing(graph, colours, v)
        completed[edge] = next(iter(missing(graph, colours, u)))
    assert graph.proper(completed) and set(completed) == {0, 1, 2}
    return completed


def audit_one_buffer_retraction(graph, initial, closed_class=False):
    """Every allowable move from every one-edge deletion of the supplied start.

    Independently tests that unique completions differ by one ordinary Kempe
    change or agree. With closed_class=True, also exhausts the 6(|E|+1) states
    over all global permutations and certifies closure under every <=1 move.
    """
    parents = {tuple(p[c] for c in initial) for p in permutations(range(3))} if closed_class else {tuple(initial)}
    states = set(parents)
    for parent in parents:
        for e in graph.all:
            state = list(parent)
            state[e] = 3
            states.add(tuple(state))
    counts = Counter(states=len(states))
    for state in states:
        parent = complete_one_buffer(graph, state)
        for pair in combinations(range(4), 2):
            for comp in graph.components(state, pair):
                nxt = list(state)
                change(nxt, pair, comp)
                if nxt.count(3) > 1:
                    continue
                after = complete_one_buffer(graph, nxt)
                difference = {e for e in graph.all if parent[e] != after[e]}
                if difference:
                    used = {parent[e] for e in difference}
                    assert len(used) == 2
                    assert graph.component(parent, used, min(difference)) == difference
                    check = parent.copy()
                    change(check, sorted(used), difference)
                    assert check == after
                if closed_class:
                    assert tuple(nxt) in states
                counts['moves_checked'] += 1
    return dict(counts)


def small_controls():
    counts = Counter()
    for q in (4, 8, 16, 32):
        field, seen = LogField(q), set()
        for t in range(1, q):
            if t in seen:
                continue
            seen.update(field.frobenius_orbit(t))
            cycle = field.cycle(t)
            if len(cycle) != q + 1:
                continue
            for c in range(1, q):
                if c == t:
                    continue
                point = PointDiagram(field, t, c, cycle)
                folded = from_point(point)
                graph, initial = closed_graph(folded)
                counts['parameters'] += 1
                audit = audit_one_buffer_retraction(graph, initial)
                counts['one_buffer_states'] += audit['states']
                counts['one_buffer_moves'] += audit['moves_checked']
                # Test both choices of temporary ordinary colour at every edge,
                # not only the natural-chord seeds used by discover().
                for pivot in graph.all:
                    for a in sorted({0, 1, 2} - {initial[pivot]}):
                        switched = seven_switches(graph, initial, pivot, a)
                        counts['seven_change_pivots'] += 1
                        counts[switched['status']] += 1
                        for candidate in switched['candidates']:
                            audit_buffer_sequence(graph, initial, candidate['colours'],
                                                  candidate['moves'], closed=True)
                            for item in eligible_certificates(point, folded, candidate['colours']):
                                cert, _ = audit_matching(point, folded, frozenset(item['selected']),
                                                         permutations=True)
                                assert cert == item['certificate']
                                counts['matching_audits'] += 1
    return dict(counts)


def abstract_return_control():
    """A genuine three-colouring change with two unit gains, after REORDERING.

    This is not a projective-cycle instance and has no claimed Cayley lift.
    It excludes automatic nonvanishing based only on the affine label set.
    """
    field = LogField(32)
    actual = from_point(PointDiagram(field, 6, 3, field.cycle(6)))
    labels = [1, 25, 27, 9, 22, 30, 4, 3, 11, 12, 14, 20, 28, 19, 17, 6]
    folded = FoldedDiagram(field, labels, 15)
    assert set(labels) == set(actual.labels) and labels != actual.labels
    assert folded.delta == actual.delta and folded.root_path(folded.canonical)['gain'] == 1
    graph, initial = closed_graph(folded)
    switched = seven_switches(graph, initial, 0)
    candidate = next(c for c in switched['candidates'] if c['changed_unlabelled_colouring'])
    roots = []
    for colour in sorted({0, 1, 2} - {candidate['colours'][-1]}):
        selected = frozenset(e for e, c in enumerate(candidate['colours']) if c == colour)
        root = folded.root_path(selected)
        assert root['gain'] == 1
        roots.append({'colour': colour, 'selected': sorted(selected), 'gain': root['gain'],
                      'vertices': [v + 1 for v in root['vertices']]})
    return {'scope': 'abstract reordered labels, NOT a projective-cycle or Cayley counterexample',
            'q': 32, 'modulus': field.modulus, 't': 6, 'c': 3, 'delta': 15,
            'actual_labels': actual.labels, 'permuted_labels': labels,
            'pivot_vertices': [1, 2], 'eligible_roots': roots,
            'seven_change_certificate': audit_buffer_sequence(
                graph, initial, candidate['colours'], candidate['moves'], closed=True)}


def single_trace_control():
    """A complete canonical-gain scan for ONE full-cycle trace at q=8192.

    Folded translation pairing halves the scan; seven-change repairs have even
    off-path folded circuits, so their successes hold for BOTH translations.
    This is not a complete q=8192 all-traces colouring audit.
    """
    field, t = LogField(8192), 13
    cycle = field.cycle(t)
    assert len(cycle) == field.q + 1
    labels = [field.mul(z, z) ^ field.mul(t, z) ^ 1 for z in cycle[1:field.q // 2 + 1]]
    position = {k: i for i, k in enumerate(labels)}
    tested, failures = 0, []
    for c in range(1, field.q):
        if c == t or c > (c ^ t):
            continue
        tested += 1
        delta = field.mul(c, c ^ t)
        if canonical_gain(field, labels, position, delta) != 1:
            continue
        point = PointDiagram(field, t, c, cycle)
        folded = from_point(point)
        assert folded.root_path(folded.canonical)['gain'] == 1
        result = discover(point, folded)
        assert result['status'] == 'REPAIRED'
        # Independently audit the partner too, rather than relying just on parity.
        partner = PointDiagram(field, t, c ^ t, cycle)
        cert, _ = audit_matching(partner, folded, frozenset(result['first']['selected']), permutations=True)
        assert cert['colours'] and cert['root_gain'] == result['first']['certificate']['root_gain']
        failures.append({'c': c, 'partner_c': c ^ t, **result})
    assert tested == 4095
    return {'scope': 'complete generating translations for the single trace t=13 only',
            'q': field.q, 'modulus': field.modulus, 't': t,
            'folded_translation_pairs': tested, 'normalized_generating_parameters': 2 * tested,
            'canonical_failure_pairs': len(failures), 'cases': failures}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--json-output', type=Path)
    parser.add_argument('--exhaustive-stored', action='store_true',
                        help='audit every seven-change pivot on stored q<=1024 failures (several minutes)')
    args = parser.parse_args()
    output = {'scope': 'exact two-buffer switch; small controls; stored failures and specified larger cases',
              'edge_ids': 'zero-based; displayed vertex labels are one-based',
              'stored_search': 'all pivots' if args.exhaustive_stored else 'first success or exhaustive failure',
              'small_controls': small_controls(), 'fields': []}
    print('SMALL CONTROLS', output['small_controls'], flush=True)
    directory = Path(__file__).resolve().parent
    for filename in ('borel_reflected_results.json', 'borel_reflected_1024.json'):
        for record in json.loads((directory / filename).read_text()):
            field, cases, counts = LogField(record['q']), [], Counter()
            for case in record['canonical_failures']:
                point = PointDiagram(field, case['t'], case['c'], field.cycle(case['t']))
                result = discover(point, from_point(point), exhaustive=args.exhaustive_stored)
                counts['direct_failures'] += 1
                counts['weighted_failures'] += case['orbit_weight']
                counts['direct_repaired'] += int(result['status'] == 'REPAIRED')
                counts['weighted_repaired'] += case['orbit_weight'] * int(result['status'] == 'REPAIRED')
                cases.append({'t': case['t'], 'c': case['c'], 'orbit_weight': case['orbit_weight'], **result})
            print('STORED FAILURES', field.q, dict(counts), flush=True)
            output['fields'].append({'q': field.q, 'counts': dict(counts), 'cases': cases})
    field = LogField(256)
    point = PointDiagram(field, 15, 151, field.cycle(15))
    folded = from_point(point)
    graph, initial = closed_graph(folded)
    closure = audit_one_buffer_retraction(graph, initial, closed_class=True)
    assert closure['states'] == 1158
    rigid = discover(point, folded, exhaustive=True)
    assert rigid['status'] == 'REPAIRED'
    output['rigid_control'] = {'q': 256, 'modulus': field.modulus, 't': 15, 'c': 151,
                               'one_buffer_closed_class': closure, **rigid}
    print('RIGID CONTROL', closure, rigid['counts'], rigid['first']['eligible_root_gains'], flush=True)
    output['large_controls'] = []
    for record in json.loads((directory / 'borel_exchange_stress_results.json').read_text())['fields']:
        field = LogField(record['q'])
        for case in record['cases']:
            t, c = case['t'], case['c']
            point = PointDiagram(field, t, c, field.cycle(t))
            result = discover(point, from_point(point))
            print('LARGE CONTROL', field.q, t, c, result['status'], result['counts']['pivots_tested'], flush=True)
            output['large_controls'].append({'q': field.q, 't': t, 'c': c, **result})
    output['abstract_return_control'] = abstract_return_control()
    output['single_trace_control'] = single_trace_control()
    print('SINGLE TRACE', {k: v for k, v in output['single_trace_control'].items() if k != 'cases'}, flush=True)
    if args.json_output:
        args.json_output.write_text(json.dumps(output, indent=2) + '\n')
