"""Exchanges along full orbits of short Moebius words in the point quotient.

For a sign pattern eps = (eps_1, ..., eps_p), the alternating walk
z -> a z -> x^{eps_1} a z -> ... returns after ell periods, where ell is the
orbit length of the Moebius word W = x^{eps_p} a ... x^{eps_1} a at z.  When
all 2 p ell visited points are distinct, this is an alternating cycle of the
point quotient with respect to the all-chord matching M0.  Unlike the earlier
short exchanges (fixed points of W, one period), the cycle length here is
unbounded: it is p times the order of W.  This script exchanges M0 along one
such cycle, or along a cycle and its reversal image, and tests the lifted
colouring by the exact circuit-word criterion.  It is an experiment, not a
proof; no uniform statement is inferred.
"""

import argparse
import json
from collections import Counter
from pathlib import Path

from borel_exchange_stress import LogField
from borel_reversal_experiment import PointDiagram


def patterns(max_period):
    """Sign patterns up to cyclic rotation."""
    seen, result = set(), []
    for p in range(1, max_period+1):
        for bits in range(1 << p):
            eps = tuple(1 if bits >> i & 1 else -1 for i in range(p))
            rotations = {eps[i:]+eps[:i] for i in range(p)}
            if rotations & seen:
                continue
            seen |= rotations
            result.append(eps)
    return result


class Orbits:
    def __init__(self, point):
        self.point = point
        n = point.n
        self.n = n
        self.chord_at = [None]*n
        self.a_image = [0]*n
        for e in range(n, len(point.edges)):
            u, v = point.edges[e]
            self.chord_at[u] = e
            self.chord_at[v] = e
            self.a_image[u], self.a_image[v] = v, u
        self.all_chords = frozenset(range(n, len(point.edges)))

    def walk(self, start, eps):
        """Alternating walk from start; returns (vertices, edges) or None if not simple."""
        n, p = self.n, len(eps)
        vertices, edges, seen = [], [], set()
        z, k = start, 0
        while True:
            if z == 0 or z in seen:
                return None
            seen.add(z)
            vertices.append(z)
            edges.append(self.chord_at[z])
            z = self.a_image[z]
            if z in seen:
                return None
            seen.add(z)
            vertices.append(z)
            step = eps[k]
            edges.append(z if step == 1 else (z-1) % n)
            z = (z+step) % n
            k = (k+1) % p
            if z == start and k == 0:
                return vertices, edges

    def cycles(self, eps):
        found, done = [], set()
        for start in range(1, self.n):
            if start in done:
                continue
            result = self.walk(start, eps)
            if result is None:
                continue
            vertices, edges = result
            done.update(vertices)
            found.append((vertices, edges))
        return found

    def reverse(self, vertices):
        return [(-v) % self.n for v in vertices]

    def exchange(self, edge_sets):
        selected = set(self.all_chords)
        for edges in edge_sets:
            selected ^= set(edges)
        return selected

    def test(self, selected):
        profile = self.point.profile(selected)
        return self.point.good(profile), profile


def analyse(q, t, c, max_period, field=None):
    field = field or LogField(q)
    cycle = field.cycle(t)
    point = PointDiagram(field, t, c, cycle)
    orbits = Orbits(point)
    n = point.n
    out = {'q': q, 't': t, 'c': c, 'patterns': []}
    total = Counter()
    for eps in patterns(max_period):
        cycles = orbits.cycles(eps)
        record = {'pattern': ''.join('+' if e == 1 else '-' for e in eps),
                  'simple_cycles': len(cycles), 'period_orbits': None,
                  'single': Counter(), 'symmetrized': Counter(), 'witness': None}
        if cycles:
            record['period_orbits'] = sorted(Counter(len(v)//(2*len(eps)) for v, _ in cycles).items())
        index = {}
        for i, (vertices, edges) in enumerate(cycles):
            for v in vertices:
                index[v] = i
        for i, (vertices, edges) in enumerate(cycles):
            good, profile = orbits.test(orbits.exchange([edges]))
            record['single']['colours' if good else 'fails'] += 1
            odd = [p for p in profile if p['length'] % 2]
            total['single_tested'] += 1
            if good:
                total['single_colours'] += 1
                if record['witness'] is None:
                    record['witness'] = {'kind': 'single', 'start': vertices[0],
                                         'cycle_length': len(vertices),
                                         'point_profile': sorted(p['length'] for p in profile)}
            # reversal image: the same pattern negated; symmetrize when disjoint
            mirror = orbits.reverse(vertices)
            mirror_set = set(mirror)
            if mirror_set == set(vertices):
                record['symmetrized']['self_symmetric'] += 1
            elif mirror_set & set(vertices):
                record['symmetrized']['overlapping'] += 1
            else:
                mirror_edges = []
                for k in range(0, len(mirror), 2):
                    z = mirror[k]
                    mirror_edges.append(orbits.chord_at[z])
                    w = mirror[k+1]
                    nxt = mirror[(k+2) % len(mirror)]
                    mirror_edges.append(w if (w+1) % n == nxt else (w-1) % n)
                sel = orbits.exchange([edges, mirror_edges])
                good2, profile2 = orbits.test(sel)
                record['symmetrized']['colours' if good2 else 'fails'] += 1
                total['sym_tested'] += 1
                if good2:
                    total['sym_colours'] += 1
                    if record['witness'] is None:
                        record['witness'] = {'kind': 'symmetrized', 'start': vertices[0],
                                             'cycle_length': len(vertices),
                                             'point_profile': sorted(p['length'] for p in profile2)}
        record['single'] = dict(record['single'])
        record['symmetrized'] = dict(record['symmetrized'])
        out['patterns'].append(record)
    out['totals'] = dict(total)
    out['repaired'] = any(r['witness'] for r in out['patterns'])
    out['repairing_patterns'] = [r['pattern'] for r in out['patterns'] if r['witness']]
    return out


def stored_failures(paths):
    jobs = []
    for path in paths:
        for block in json.loads(Path(path).read_text()):
            seen = set()
            for f in block['canonical_failures']:
                key = (f['t'], min(f['c'], f['c'] ^ f['t']))
                if key in seen:
                    continue
                seen.add(key)
                jobs.append((block['q'], f['t'], f['c']))
    return jobs


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reflected', nargs='+',
                        default=['code/borel_reflected_results.json', 'code/borel_reflected_1024.json'])
    parser.add_argument('--max-period', type=int, default=4)
    parser.add_argument('--max-q', type=int, default=1024)
    parser.add_argument('--json-output', type=Path)
    args = parser.parse_args()
    jobs = [j for j in stored_failures(args.reflected) if j[0] <= args.max_q]
    fields, results = {}, []
    for q, t, c in jobs:
        fields.setdefault(q, LogField(q))
        r = analyse(q, t, c, args.max_period, fields[q])
        results.append(r)
        print('ORBIT', q, t, c, 'repaired', r['repaired'], 'by', r['repairing_patterns'],
              'totals', r['totals'], flush=True)
    summary = {'diagrams': len(results), 'repaired': sum(r['repaired'] for r in results),
               'single_tested': sum(r['totals'].get('single_tested', 0) for r in results),
               'single_colours': sum(r['totals'].get('single_colours', 0) for r in results),
               'sym_tested': sum(r['totals'].get('sym_tested', 0) for r in results),
               'sym_colours': sum(r['totals'].get('sym_colours', 0) for r in results)}
    print('SUMMARY', summary)
    if args.json_output:
        args.json_output.write_text(json.dumps({'scope': 'orbit exchanges on stored canonical failures; experiment only',
                                                'summary': summary, 'results': results}, indent=2)+'\n')
