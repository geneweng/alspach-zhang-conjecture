"""Global folded recolouring, a rigid three-colour class, and a buffer escape.

Closing the folded path gives an explicitly three-edge-coloured cubic graph H.
Every colour class avoiding the closing edge is a folded perfect matching with
even off-root circuits. Its Cayley lift colours exactly when the root gain is
not one. This is a useful restriction, NOT an existence theorem.

The four-colour algorithm below implements the vertex-deletion proof of
McDonald--Mohar--Scheide (2012), Lemma 7, on F=H-e*. Its target is supplied
in advance. Intermediate states need not extend over e*, but a separately
checked lifting routine DOES extend the rigid control's sequence to H.
"""

import argparse
import json
from collections import Counter, deque
from itertools import combinations, permutations
from pathlib import Path

from borel_exchange_stress import LogField
from borel_folded_experiment import from_point
from borel_reversal_experiment import PointDiagram
from borel_rotation_experiment import exchange_components
from borel_three_exchange import affine_certificate, audit_matching


class EdgeGraph:
    def __init__(self, vertices, edges):
        self.vertices, self.edges = vertices, list(edges)
        self.inc = [set() for _ in range(vertices)]
        for e, (u, v, _) in enumerate(self.edges):
            assert u != v
            self.inc[u].add(e)
            self.inc[v].add(e)
        self.all = frozenset(range(len(edges)))

    def proper(self, colours, active=None):
        active = self.all if active is None else active
        return all(len({colours[e] for e in local & active}) == len(local & active)
                   for local in self.inc)

    def component(self, colours, pair, seed, active=None):
        active = self.all if active is None else active
        assert seed in active and colours[seed] in pair
        found, pending = {seed}, [seed]
        while pending:
            e = pending.pop()
            for v in self.edges[e][:2]:
                for f in self.inc[v] & active:
                    if f not in found and colours[f] in pair:
                        found.add(f)
                        pending.append(f)
        return frozenset(found)

    def components(self, colours, pair, active=None):
        active = self.all if active is None else active
        unseen = {e for e in active if colours[e] in pair}
        result = []
        while unseen:
            found = self.component(colours, pair, min(unseen), active)
            result.append(found)
            unseen -= found
        return result


def change(colours, pair, edges):
    a, b = pair
    for e in edges:
        assert colours[e] in pair
        colours[e] = b if colours[e] == a else a


def normalized(colours):
    names = {}
    return tuple(names.setdefault(c, len(names)) for c in colours)


def closed_graph(folded):
    graph = EdgeGraph(folded.h, folded.edges + [(0, folded.h - 1, 'close')])
    colours = [0 if kind == 'a' else 2 if kind == 'close' else 1 + u % 2
               for u, _, kind in graph.edges]
    assert graph.proper(colours)
    assert all(len(local) == 3 for local in graph.inc)
    return graph, colours


def colour_matching(graph, selected):
    """Complete a perfect matching to a three-colouring, if its complement is even."""
    assert all(len(local & selected) == 1 for local in graph.inc)
    colours = [0 if e in selected else -1 for e in range(len(graph.edges))]
    unseen = set(graph.all - selected)
    while unseen:
        edge, start = min(unseen), graph.edges[min(unseen)][0]
        current, colour = start, 1
        while True:
            assert edge in unseen
            unseen.remove(edge)
            colours[edge] = colour
            u, v, _ = graph.edges[edge]
            current = v if current == u else u
            if current == start:
                assert colour == 2
                break
            edge = next(iter(graph.inc[current] - selected - {edge}))
            colour = 3 - colour
    assert graph.proper(colours)
    return colours


def eligible_certificates(point, folded, colours):
    graph, _ = closed_graph(folded)
    assert graph.proper(colours) and set(colours) == {0, 1, 2}
    result = []
    for colour in sorted({0, 1, 2} - {colours[-1]}):
        selected = frozenset(e for e, c in enumerate(colours) if c == colour)
        cert = affine_certificate(point, folded, selected)
        assert all(r['folded_length'] % 2 == 0 for r in cert['offroot'])
        assert cert['root_vertices'] % 2 == 0
        assert cert['colours'] == (cert['root_gain'] != 1)
        result.append({'colour': colour, 'selected': sorted(selected), 'certificate': cert})
    return result


def parity_cubes(point, folded):
    """Search both monochromatic-path-parity cubes, stopping at first success.

    A cube is all subsets of natural (0,1) or (0,2) circuits not containing
    e*. No cycle-length cap. A failure exhausts BOTH cubes; success can be early.
    """
    graph, natural = closed_graph(folded)
    all_chords = frozenset(e for e, c in enumerate(natural) if c == 0)
    sizes, tested = [], 0
    for colour in (1, 2):
        circuits = [c for c in graph.components(natural, (0, colour))
                    if len(graph.edges) - 1 not in c]
        sizes.append([len(c) // 2 for c in circuits])
        for mask in range(1, 1 << len(circuits)):
            selected = set(all_chords)
            colours = natural.copy()
            for i, circuit in enumerate(circuits):
                if mask >> i & 1:
                    selected.symmetric_difference_update(circuit)
                    change(colours, (0, colour), circuit)
            assert graph.proper(colours)
            cert = affine_certificate(point, folded, selected)
            assert all(r['folded_length'] % 2 == 0 for r in cert['offroot'])
            tested += 1
            if cert['colours']:
                checked, profile = audit_matching(point, folded, selected, permutations=True)
                assert checked == cert
                return {'status': 'REPAIRED', 'tested_nonempty_subsets': tested,
                        'cube_chord_counts_examined': sizes,
                        'path_start_parity': 'odd' if colour == 1 else 'even',
                        'selected': sorted(selected), 'certificate': cert,
                        'point_lengths': [p['length'] for p in profile]}
    return {'status': 'BOTH_CUBES_EXHAUSTED', 'tested_nonempty_subsets': tested,
            'cube_chord_counts_examined': sizes}


def star_completion(graph, active, colours, target, star):
    """Small BFS version of Lemma 6 for the at-most-two-edge disagreement path."""
    keys = sorted(star)
    initial, goal = tuple(colours[e] for e in keys), tuple(target[e] for e in keys)
    pending, parents = deque([initial]), {initial: None}
    while pending:
        state = pending.popleft()
        if state == goal:
            moves = []
            while parents[state] is not None:
                state, move = parents[state]
                moves.append(move)
            return list(reversed(moves))
        trial = colours.copy()
        for e, c in zip(keys, state):
            trial[e] = c
        for pair in combinations(range(4), 2):
            for edges in graph.components(trial, pair, active):
                if not edges <= star:
                    continue
                nxt = trial.copy()
                change(nxt, pair, edges)
                assert graph.proper(nxt, active)
                key = tuple(nxt[e] for e in keys)
                if key not in parents:
                    parents[key] = state, (pair, edges)
                    pending.append(key)
    raise AssertionError('Lemma 6 star completion failed')


def buffer_sequence(graph, initial, target):
    """Construct a four-colour Kempe sequence on a 2-degenerate subcubic graph.

    Direct vertex-deletion/lifting implementation of MMS Lemma 7. Each recorded
    move swaps a MAXIMAL bichromatic component of the graph at its level.
    """
    assert graph.proper(initial) and graph.proper(target)
    assert max(initial + target) <= 3

    def recurse(active, start):
        if not active:
            return []
        v = next(v for v, inc in enumerate(graph.inc) if 0 < len(inc & active) <= 2)
        star = graph.inc[v] & active
        lower = active - star
        lower_moves = recurse(lower, start.copy())
        colours, moves = start.copy(), []

        def perform(pair, edges):
            assert graph.component(colours, pair, min(edges), active) == edges
            change(colours, pair, edges)
            assert graph.proper(colours, active)
            moves.append((pair, edges))

        for pair, low_edges in lower_moves:
            full = graph.component(colours, pair, min(low_edges), active)
            if full & lower != low_edges:
                # Two lower components were joined through the deleted vertex.
                assert len(star) == 2 and {colours[e] for e in star} == set(pair)
                endpoints = {u for e in low_edges for u in graph.edges[e][:2]}
                bridge = next(e for e in sorted(star)
                              if (set(graph.edges[e][:2]) - {v}) & endpoints)
                x = next(u for u in graph.edges[bridge][:2] if u != v)
                missing = set(range(4)) - {colours[e] for e in (graph.inc[v] | graph.inc[x]) & active}
                assert missing
                perform((colours[bridge], min(missing)), frozenset({bridge}))
                full = graph.component(colours, pair, min(low_edges), active)
            assert full & lower == low_edges
            perform(pair, full)
        assert all(colours[e] == target[e] for e in lower)
        for pair, edges in star_completion(graph, active, colours, target, star):
            perform(pair, edges)
        assert all(colours[e] == target[e] for e in active)
        return moves

    return recurse(graph.all, initial.copy())


def audit_buffer_sequence(graph, initial, target, moves, closed=False):
    """Independent edge-boundary and incidence replay (not component traversal)."""
    colours, max_buffer, unclosable, three_states = initial.copy(), 0, 0, 0
    records = []
    for pair, edges in moves:
        assert edges and all(colours[e] in pair for e in edges)
        vertices = {v for e in edges for v in graph.edges[e][:2]}
        # Maximality plus connectedness, checked independently of component().
        assert all(colours[e] not in pair for v in vertices for e in graph.inc[v] - edges)
        reached, pending = {next(iter(vertices))}, []
        pending.extend(reached)
        while pending:
            v = pending.pop()
            for e in graph.inc[v] & edges:
                for w in graph.edges[e][:2]:
                    if w not in reached:
                        reached.add(w)
                        pending.append(w)
        assert reached == vertices
        change(colours, pair, edges)
        assert all(len({colours[e] for e in inc}) == len(inc) for inc in graph.inc)
        max_buffer = max(max_buffer, colours.count(3))
        if not closed:
            missing = [set(range(4)) - {colours[e] for e in graph.inc[v]}
                       for v in (0, graph.vertices - 1)]
            unclosable += int(not missing[0] & missing[1])
        three_states += int(3 not in colours)
        records.append({'pair': list(pair), 'edge_ids': sorted(edges)})
    assert colours == target
    return {'graph': 'H, including closing edge' if closed else 'F, closing edge absent',
            'moves': len(moves),
            'max_edges_of_buffer_colour': max_buffer,
            'states_without_buffer_colour': three_states,
            'states_not_extendible_over_closing_edge_with_four_colours': None if closed else unclosable,
            'initial_colours': initial, 'target_colours': target, 'sequence': records}


def lift_buffer_sequence(graph, initial, target, moves):
    """Try to lift a specified sequence from F to H; failure is permitted.

    Before each change, recolour only e* if needed, then require its full
    bichromatic component in H to restrict to exactly the prescribed F move.
    All valid colours of e* are linked by singleton Kempe changes. Thus an
    unsuccessful step has exhausted this particular step-by-step lifting rule.
    This routine is NOT an implementation of the general cubic theorem.
    """
    colours, lifted, closing = initial.copy(), [], len(graph.edges) - 1
    for pair, edges in moves:
        candidates = []
        for c in range(4):
            trial = colours[:-1] + [c]
            if not graph.proper(trial):
                continue
            full = graph.component(trial, pair, min(edges))
            if full - {closing} == edges:
                candidates.append((c, full))
        if not candidates:
            return None
        c, full = next((item for item in candidates if item[0] == colours[-1]), candidates[0])
        if c != colours[-1]:
            recolour = (colours[-1], c)
            assert graph.component(colours, recolour, closing) == {closing}
            lifted.append((recolour, frozenset({closing})))
            change(colours, recolour, {closing})
        lifted.append((pair, full))
        change(colours, pair, full)
        assert graph.proper(colours)
    assert colours[:-1] == target[:-1]
    if colours[-1] != target[-1]:
        recolour = (colours[-1], target[-1])
        assert graph.component(colours, recolour, closing) == {closing}
        lifted.append((recolour, frozenset({closing})))
        change(colours, recolour, {closing})
    assert colours == target
    return lifted


def rigid_control():
    field = LogField(256)
    point = PointDiagram(field, 15, 151, field.cycle(15))
    folded = from_point(point)
    graph, natural = closed_graph(folded)
    assert len({frozenset((u, v)) for u, v, _ in graph.edges}) == len(graph.edges)
    profiles = {str(pair): [len(c) for c in graph.components(natural, pair)]
                for pair in combinations(range(3), 2)}
    assert all(lengths == [128] for lengths in profiles.values())
    # These six states are CLOSED under every possible three-colour Kempe move.
    orbit = {tuple(p[c] for c in natural) for p in permutations(range(3))}
    assert len(orbit) == 6
    bad_matchings = set()
    for state in orbit:
        assert graph.proper(state)
        for pair in combinations(range(3), 2):
            for edges in graph.components(state, pair):
                nxt = list(state)
                change(nxt, pair, edges)
                assert tuple(nxt) in orbit
        for item in eligible_certificates(point, folded, state):
            assert item['certificate']['root_gain'] == 1
            bad_matchings.add(tuple(item['selected']))
    assert len(bad_matchings) == 2
    for selected in bad_matchings:
        cert, profile = audit_matching(point, folded, frozenset(selected), permutations=True)
        assert len(profile) == 1 and profile[0]['length'] == 257
        assert tuple(profile[0]['word_matrix']) == (1, 0, 0, 1)
    chord = next(e for e, (u, v, kind) in enumerate(folded.edges)
                 if kind == 'a' and {u + 1, v + 1} == {3, 122})
    starts = [1] + list(range(4, 121, 2)) + [123, 125, 127]
    selected = frozenset({chord} | {s - 1 for s in starts})
    cert, profile = audit_matching(point, folded, selected, permutations=True)
    assert cert['root_gain'] == 28 and cert['root_beta'] == 187 and not cert['offroot']
    assert [p['length'] for p in profile] == [257]
    repaired = colour_matching(graph, selected)
    assert normalized(repaired) != normalized(natural)
    difference = selected ^ folded.canonical
    assert len(difference) == 120
    assert Counter(natural[e] for e in difference) == {0: 1, 1: 60, 2: 59}
    # The all-three-colour difference is one circuit, not a Kempe component.
    assert len(graph.component(natural, (0, 1, 2), min(difference), difference)) == 120
    partner = PointDiagram(field, 15, 152, field.cycle(15))
    partner_cert, _ = audit_matching(partner, folded, selected, permutations=True)
    assert partner_cert['root_gain'] == 28 and partner_cert['root_beta'] == 187
    open_graph = EdgeGraph(folded.h, folded.edges)
    moves = buffer_sequence(open_graph, natural[:-1], repaired[:-1])
    buffer = audit_buffer_sequence(open_graph, natural[:-1], repaired[:-1], moves)
    assert buffer['max_edges_of_buffer_colour'] > 0
    lifted = lift_buffer_sequence(graph, natural, repaired, moves)
    assert lifted is not None
    closed_buffer = audit_buffer_sequence(graph, natural, repaired, lifted, closed=True)
    return {'q': 256, 'modulus': field.modulus, 't': 15, 'c': 151, 'partner_c': 152,
            'vertices_H': graph.vertices, 'edges_H': len(graph.edges),
            'natural_colours_H': natural, 'natural_bicolour_cycle_lengths': profiles,
            'three_colour_labelled_orbit_size': 6, 'eligible_matchings_in_orbit': 2,
            'all_orbit_root_gains': [1, 1], 'orbit_closed_under_all_three_colour_moves': True,
            'repair': {'selected_chord': [3, 122], 'path_edge_starts': starts,
                       'certificate': cert, 'point_profile': profile,
                       'colours_H': repaired,
                       'all_chord_difference_chord_counts': exchange_components(folded, selected),
                       'canonical_difference_edges': 120,
                       'canonical_difference_natural_colour_counts': {str(k): v for k, v in
                           sorted(Counter(natural[e] for e in difference).items())}},
            'buffer_escape_on_F': buffer, 'buffer_escape_on_H': closed_buffer}


def small_controls():
    """Every even-complement matching through q=32, all generating representatives."""
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
                graph, natural = closed_graph(folded)
                counts['parameter_representatives'] += 1
                for selected in folded.matchings():
                    cert = affine_certificate(point, folded, selected)
                    if any(r['folded_length'] % 2 for r in cert['offroot']):
                        continue
                    target = colour_matching(graph, selected)
                    assert any(frozenset(r['selected']) == selected for r in
                               eligible_certificates(point, folded, target))
                    checked, _ = audit_matching(point, folded, selected, permutations=True)
                    assert checked == cert
                    # The source theorem is stated for simple graphs: only run
                    # its constructive audit when F has no parallel edges.
                    if len({frozenset((u, v)) for u, v, _ in folded.edges}) == len(folded.edges):
                        open_graph = EdgeGraph(folded.h, folded.edges)
                        moves = buffer_sequence(open_graph, natural[:-1], target[:-1])
                        audit_buffer_sequence(open_graph, natural[:-1], target[:-1], moves)
                        counts['buffer_sequences'] += 1
                        counts['buffer_moves'] += len(moves)
                    counts['even_complement_matchings'] += 1
    return dict(counts)


def stored_failures():
    directory = Path(__file__).resolve().parent
    records = []
    for name in ('borel_reflected_results.json', 'borel_reflected_1024.json'):
        records.extend(json.loads((directory / name).read_text()))
    result = []
    for record in records:
        field, counts, cases = LogField(record['q']), Counter(), []
        for case in record['canonical_failures']:
            t, c, weight = case['t'], case['c'], case['orbit_weight']
            point = PointDiagram(field, t, c, field.cycle(t))
            folded = from_point(point)
            assert folded.root_path(folded.canonical)['gain'] == 1
            tested = parity_cubes(point, folded)
            counts['direct_failures'] += 1
            counts['weighted_failures'] += weight
            counts['direct_repaired'] += int(tested['status'] == 'REPAIRED')
            counts['weighted_repaired'] += weight * int(tested['status'] == 'REPAIRED')
            counts['tested_nonempty_subsets'] += tested['tested_nonempty_subsets']
            cases.append({'t': t, 'c': c, 'orbit_weight': weight, **tested})
        print('PARITY CUBES', field.q, dict(counts), flush=True)
        result.append({'q': field.q, 'counts': dict(counts), 'cases': cases})
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--json-output', type=Path)
    args = parser.parse_args()
    output = {'scope': 'global recolouring lemmas; exhaustive small controls; stored failures only',
              'edge_ids': 'zero-based; vertices and path starts in repair descriptions are one-based',
              'source': 'https://arxiv.org/abs/1005.2248', 'small_controls': small_controls()}
    print('SMALL CONTROLS', output['small_controls'], flush=True)
    output['stored_failures'] = stored_failures()
    output['rigid_control'] = rigid_control()
    buffer = output['rigid_control']['buffer_escape_on_H']
    print('RIGID CONTROL', {k: v for k, v in buffer.items()
                            if k not in ('initial_colours', 'target_colours', 'sequence')}, flush=True)
    output['large_controls'] = []
    for q, t, c in ((2048, 278, 1567), (4096, 681, 1207)):
        field = LogField(q)
        point = PointDiagram(field, t, c, field.cycle(t))
        result = parity_cubes(point, from_point(point))
        print('LARGE CONTROL', q, t, c, result['status'], flush=True)
        output['large_controls'].append({'q': q, 't': t, 'c': c, **result})
    if args.json_output:
        args.json_output.write_text(json.dumps(output, indent=2) + '\n')
