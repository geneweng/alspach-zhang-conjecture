"""Exact small recolouring closures, independent of arithmetic rank success.

The finite-state search never discards a colouring because it is rank-neutral.
It distinguishes exhaustive reachability from a merely increasing row span.
No general connectivity or full-rank theorem is inferred from the controls.
"""

import argparse
import json
from collections import Counter, deque
from itertools import combinations
from pathlib import Path

from borel_exchange_stress import LogField
from borel_folded_experiment import FoldedDiagram, from_point
from borel_gain_lattice import chord_logs, root_row
from borel_gain_rank import packed_rank, prime_divisors
from borel_kempe_exchange import change, closed_graph, normalized
from borel_reversal_experiment import PointDiagram
from borel_three_exchange import audit_matching
from borel_two_buffer import seven_switches


def all_three_colourings(graph):
    """Independent edge-domain backtracking; quotient only global colour names."""
    colours = [-1]*len(graph.edges)
    used = [0]*graph.vertices
    for colour, edge in enumerate(sorted(graph.inc[0])):
        colours[edge] = colour
        for v in graph.edges[edge][:2]:
            used[v] |= 1 << colour

    def visit():
        chosen, allowed = None, None
        for edge, colour in enumerate(colours):
            if colour >= 0:
                continue
            u, v, _ = graph.edges[edge]
            mask = 7 & ~(used[u] | used[v])
            if not mask:
                return
            if allowed is None or mask.bit_count() < allowed.bit_count():
                chosen, allowed = edge, mask
                if mask.bit_count() == 1:
                    break
        if chosen is None:
            assert graph.proper(colours)
            yield normalized(colours)
            return
        u, v, _ = graph.edges[chosen]
        for colour in range(3):
            bit = 1 << colour
            if not (allowed & bit):
                continue
            colours[chosen] = colour
            used[u] |= bit
            used[v] |= bit
            yield from visit()
            used[u] ^= bit
            used[v] ^= bit
        colours[chosen] = -1

    result = set(visit())
    assert all(graph.proper(state) for state in result)
    return result


def closure(graph, initial, seven=True):
    pending = deque([normalized(initial)])
    seen = set(pending)
    counts = Counter()
    while pending:
        state = list(pending.popleft())
        candidates = []
        for pair in combinations(range(3), 2):
            for component in graph.components(state, pair):
                target = state.copy()
                change(target, pair, component)
                candidates.append(normalized(target))
                counts['ordinary_moves'] += 1
        if seven:
            for pivot, colour in enumerate(state):
                for temporary in sorted({0, 1, 2}-{colour}):
                    result = seven_switches(graph, state, pivot, temporary)
                    counts['seven_seeds'] += 1
                    counts[result['status']] += 1
                    candidates.extend(normalized(c['colours']) for c in result['candidates'])
        for target in candidates:
            if target not in seen:
                seen.add(target)
                pending.append(target)
    return seen, dict(counts)


def audit_diagram(folded, point=None):
    graph, natural = closed_graph(folded)
    all_states = all_three_colourings(graph)
    ordinary, ordinary_counts = closure(graph, natural, seven=False)
    reached, counts = closure(graph, natural)
    assert ordinary <= reached <= all_states
    rows = set()
    logs = chord_logs(folded)
    both_unit, lift_checks = 0, 0
    for state in all_states:
        gains = []
        for excluded in {0, 1, 2}-{state[-1]}:
            raw = root_row(folded, graph, state, excluded)
            exponent = sum(a*b for a, b in zip(raw, logs)) % (folded.field.q-1)
            gains.append(exponent)
            selected = frozenset(e for e, colour in enumerate(state) if colour == excluded)
            gain = folded.field.exp[exponent]
            assert folded.root_path(selected)['gain'] == gain
            if point is not None:
                certificate, _ = audit_matching(point, folded, selected, permutations=True)
                assert certificate['root_gain'] == gain
                assert certificate['colours'] == (gain != 1)
                lift_checks += 1
            if any(raw):
                sign = next(x for x in raw if x)
                rows.add(tuple(sign*x for x in raw))
        both_unit += all(v == 0 for v in gains)
    rows = sorted(rows)
    ranks = {p: packed_rank(rows, p)['rank'] for p in prime_divisors(folded.field.q-1)}
    return {'all_colourings_mod_global_names': len(all_states),
            'ordinary_reachable': len(ordinary), 'combined_reachable': len(reached),
            'combined_exhausts_all_colourings': reached == all_states,
            'ordinary_counts': ordinary_counts, 'combined_counts': counts,
            'all_colouring_root_rows': len(rows), 'root_ranks': ranks,
            'both_unit_colourings': both_unit,
            'independent_matrix_and_permutation_lift_checks': lift_checks}


def small_controls():
    output = []
    for q in (4, 8, 16, 32):
        field, seen, cases = LogField(q), set(), []
        for t in range(1, q):
            if t in seen:
                continue
            seen.update(field.frobenius_orbit(t))
            cycle = field.cycle(t)
            if len(cycle) != q+1:
                continue
            for c in range(1, q):
                if c == t:
                    continue
                point = PointDiagram(field, t, c, cycle)
                folded = from_point(point)
                result = audit_diagram(folded, point)
                cases.append({'t': t, 'c': c, **result})
        output.append({'q': q, 'modulus': field.modulus, 'cases': cases})
        print('CLOSURE', q, 'parameters', len(cases), 'colourings',
              sum(r['all_colourings_mod_global_names'] for r in cases), 'unreached',
              sum(not r['combined_exhausts_all_colourings'] for r in cases), flush=True)
    return output


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--json-output', type=Path)
    args = parser.parse_args()
    output = {'scope': 'complete small closures, not a uniform rank-growth theorem',
              'small': small_controls()}
    field = LogField(16)
    abstract = FoldedDiagram(field, [1,14,6,10,13,9,5,2], 11)
    output['abstract'] = {'scope': 'reordered labels, no claimed Cayley lift',
                          **audit_diagram(abstract)}
    if args.json_output:
        args.json_output.write_text(json.dumps(output, indent=2)+'\n')
