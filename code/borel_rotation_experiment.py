"""Audit endpoint rotations as an unbounded spanning-path construction.

In the simple folded graph the rotation graph of Hamiltonian paths from
one fixed endpoint has degrees one and two. Starting at a path to the
other endpoint therefore reaches a distinct such path. This is the
classical Thomason parity mechanism, not a non-unit-gain theorem.
"""

import argparse
import json
from collections import Counter
from pathlib import Path

from borel_exchange_certificate import audit_profile_permutations
from borel_exchange_stress import LogField
from borel_folded_experiment import from_point
from borel_reflected_family import half_word_certificate
from borel_reversal_experiment import PointDiagram


def adjacency(folded):
    adj = [set() for _ in range(folded.h)]
    for u, v, _ in folded.edges:
        assert v not in adj[u], 'rotation audit is for the simple residual class'
        adj[u].add(v)
        adj[v].add(u)
    assert [len(a) for a in adj] == [2] + [3] * (folded.h - 2) + [2]
    return adj


def edge_gain(folded, u, v):
    if abs(u - v) == 1:
        return 1
    assert folded.labels[u] ^ folded.labels[v] == folded.delta
    return folded.field.mul(folded.labels[u], folded.field.inverse[folded.labels[v]])


def path_gain(folded, path):
    gain = 1
    for u, v in zip(path, path[1:]):
        gain = folded.field.mul(gain, edge_gain(folded, u, v))
    return gain


def mate(folded, initial, limit=50000, audit_steps=False):
    adj, path, previous_pivot, pivots = adjacency(folded), tuple(initial), None, []
    assert path[0] in (0, folded.h - 1) and path[-1] == folded.h - 1 - path[0]
    for _ in range(limit):
        options = adj[path[-1]] - {path[-2], previous_pivot}
        assert len(options) == 1
        pivot = next(iter(options))
        i = path.index(pivot)
        assert i < len(path) - 2
        rotated = path[:i + 1] + (path[-1],) + path[-2:i:-1]
        if audit_steps:
            # Exact gain update, checked independently by recomputing both paths.
            tail_gain = path_gain(folded, path[i + 1:])
            denominator = folded.field.mul(edge_gain(folded, pivot, path[i + 1]),
                                            folded.field.mul(tail_gain, tail_gain))
            ratio = folded.field.mul(edge_gain(folded, pivot, path[-1]),
                                      folded.field.inverse[denominator])
            assert path_gain(folded, rotated) == folded.field.mul(path_gain(folded, path), ratio)
            assert len(set(rotated)) == folded.h
            assert all(v in adj[u] for u, v in zip(rotated, rotated[1:]))
        path, previous_pivot = rotated, pivot
        pivots.append(pivot)
        if path[-1] == initial[-1]:
            assert path != tuple(initial)
            return {'status': 'COMPLETE', 'path': path, 'pivots': pivots, 'steps': len(pivots)}
    return {'status': 'STEP_LIMIT', 'steps': limit}


def replay(folded, initial, pivots):
    """Separate list-edit replay of the recorded rotations."""
    adj, path = adjacency(folded), list(initial)
    for pivot in pivots:
        i = path.index(pivot)
        assert path[-1] in adj[pivot] and i < len(path) - 2
        tail = path[i + 1:-1]
        endpoint = path.pop()
        del path[i + 1:]
        path.append(endpoint)
        path.extend(reversed(tail))
    assert set(path) == set(range(folded.h))
    assert all(v in adj[u] for u, v in zip(path, path[1:]))
    return tuple(path)


def path_matching(folded, path):
    used = {frozenset((u, v)) for u, v in zip(path, path[1:])}
    selected = frozenset(e for e, (u, v, _) in enumerate(folded.edges)
                         if frozenset((u, v)) not in used)
    root = folded.root_path(selected)
    assert root['vertices'] == list(path)
    assert root['gain'] == path_gain(folded, path)
    return selected


def exchange_components(folded, selected):
    difference = set(selected ^ frozenset(range(folded.h - 1, len(folded.edges))))
    sizes = []
    while difference:
        start = folded.edges[min(difference)][0]
        pending, vertices, edges = [start], set(), set()
        while pending:
            v = pending.pop()
            if v in vertices:
                continue
            vertices.add(v)
            local = set(folded.inc[v]) & difference
            assert len(local) == 2
            edges.update(local)
            pending.extend(w for e in local for w in folded.edges[e][:2] if w != v)
        assert len(edges) == len(vertices) and len(edges) % 2 == 0
        sizes.append(len(edges) // 2)
        difference -= edges
    return sorted(sizes, reverse=True)


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
                if c == t or c > (c ^ t):
                    continue
                if not field.traces[field.inverse[c]] or not field.traces[field.inverse[c ^ t]]:
                    continue
                point = PointDiagram(field, t, c, cycle)
                folded = from_point(point)
                adjacency(folded)
                paths = [tuple(root['vertices']) for selected in folded.matchings()
                         if len((root := folded.root_path(selected))['vertices']) == folded.h]
                assert len(paths) >= 2 and len(paths) % 2 == 0
                for path in paths:
                    result = mate(folded, path, audit_steps=True)
                    assert result['status'] == 'COMPLETE' and result['path'] in paths
                    assert replay(folded, path, result['pivots']) == result['path']
                    back = mate(folded, result['path'])
                    assert back['status'] == 'COMPLETE' and back['path'] == path
                    counts['rotation_steps_with_gain_audits'] += result['steps']
                    counts['paired_spanning_paths'] += 1
                counts['folded_diagrams'] += 1
    return dict(counts)


def audit(q, t, c, limit):
    field = LogField(q)
    cycle = field.cycle(t)
    assert len(cycle) == q + 1 and c not in (0, t)
    point, records = PointDiagram(field, t, c, cycle), []
    folded = from_point(point)
    adjacency(folded)
    for reverse in (False, True):
        initial = tuple(reversed(range(folded.h))) if reverse else tuple(range(folded.h))
        result = mate(folded, initial, limit)
        if result['status'] != 'COMPLETE':
            records.append({**result, 'anchor': folded.h if reverse else 1})
            continue
        assert replay(folded, initial, result['pivots']) == result['path']
        back = mate(folded, result['path'], limit)
        assert back['status'] == 'COMPLETE' and back['path'] == initial
        path = tuple(reversed(result['path'])) if reverse else result['path']
        selected = path_matching(folded, path)
        profile = point.profile(folded.lift(point, selected))
        assert len(profile) == 1 and profile[0]['length'] == q + 1
        audit_profile_permutations(point, profile)
        half, gain = half_word_certificate(point, profile), path_gain(folded, path)
        assert half['beta'] == field.mul(t, gain ^ 1)
        records.append({'status': 'COMPLETE', 'anchor': folded.h if reverse else 1,
                        'steps': result['steps'], 'path': [v + 1 for v in path],
                        'pivots': [v + 1 for v in result['pivots']], 'gain': gain,
                        'exchange_component_chord_counts': exchange_components(folded, selected),
                        'profile': profile, 'half_word': half})
    return {'q': q, 'modulus': field.modulus, 't': t, 'c': c,
            'canonical_gain': folded.root_path(folded.canonical)['gain'], 'mates': records}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--step-limit', type=int, default=50000)
    parser.add_argument('--stress-results', type=Path)
    parser.add_argument('--json-output', type=Path)
    args = parser.parse_args()
    cases = [(16, 8, 4), (32, 7, 1), (128, 15, 2), (128, 15, 35), (4096, 681, 1207)]
    if args.stress_results:
        for record in json.loads(args.stress_results.read_text())['fields']:
            cases.extend((record['q'], c['t'], c['c']) for c in record['cases'])
    output = {'scope': 'small exhaustive parity controls and specified rotation walks',
              'small_controls': small_controls(), 'cases': []}
    print('SMALL ROTATION CONTROLS', output['small_controls'], flush=True)
    for q, t, c in dict.fromkeys(cases):
        record = audit(q, t, c, args.step_limit)
        print('ROTATION MATES', q, t, c, 'canonical', record['canonical_gain'],
              [(r['anchor'], r['status'], r['steps'], r.get('gain')) for r in record['mates']], flush=True)
        output['cases'].append(record)
        if args.json_output:
            args.json_output.write_text(json.dumps(output, indent=2) + '\n')
