"""Independent folded-graph and scalar-gain audit for the Borel construction.

python3 code/borel_folded_experiment.py --json-output code/borel_folded_results.json

The deterministic abstract control deliberately permutes the field labels;
it is NOT an instance of the projective-cycle problem. It rules out a proof
of the stronger root-contained repair claim using arbitrary path orderings.
"""

import argparse
import json
import random
from collections import Counter
from pathlib import Path

from borel_reflected_family import BinaryField, generated_order, half_word_certificate
from borel_reversal_experiment import PointDiagram
from cayley_snark_check import closure


class FoldedDiagram:
    def __init__(self, field, labels, delta):
        self.field, self.labels, self.delta = field, labels, delta
        self.h = len(labels)
        assert self.h % 2 == 0 and delta != 0
        position = {k: i for i, k in enumerate(labels)}
        assert len(position) == self.h and 0 not in position
        partner = [position[k ^ delta] for k in labels]
        self.edges = [(i, i + 1, 'x') for i in range(self.h - 1)]
        self.edges.extend((i, j, 'a') for i, j in enumerate(partner) if i < j)
        self.inc = [[] for _ in labels]
        for e, (u, v, _) in enumerate(self.edges):
            self.inc[u].append(e)
            self.inc[v].append(e)
        self.canonical = frozenset(range(0, self.h - 1, 2))

    def matchings(self, uncovered=None, selected=frozenset()):
        if uncovered is None:
            uncovered = frozenset(range(self.h))
        if not uncovered:
            yield selected
            return
        v = min(uncovered)
        for e in self.inc[v]:
            endpoints = frozenset(self.edges[e][:2])
            if endpoints <= uncovered:
                yield from self.matchings(uncovered - endpoints, selected | {e})

    def one_chord_matchings(self):
        for e in range(self.h - 1, len(self.edges)):
            u, v, _ = self.edges[e]
            if u % 2 or v % 2 == 0:
                continue
            chosen = {e}
            uncovered = set(range(self.h)) - {u, v}
            while uncovered:
                w = min(uncovered)
                assert w + 1 in uncovered
                chosen.add(w)
                uncovered -= {w, w + 1}
            yield frozenset(chosen)

    def root_path(self, selected):
        assert all(sum(e in selected for e in inc) == 1 for inc in self.inc)
        current, previous, gain, vertices = 0, -1, 1, []
        while True:
            assert current not in vertices
            vertices.append(current)
            if current == self.h - 1:
                return {'vertices': vertices, 'gain': gain}
            e = next(e for e in self.inc[current] if e != previous and e not in selected)
            u, v, kind = self.edges[e]
            other = v if current == u else u
            if kind == 'a':
                ratio = self.field.times[self.labels[current]][
                    self.field.inverse[self.labels[other]]]
                gain = self.field.times[gain][ratio]
            previous, current = e, other

    def lift(self, diagram, selected):
        def folded(v):
            return min(v, diagram.n - v) - 1
        original_chords = {}
        for e in range(diagram.n, len(diagram.edges)):
            key = frozenset(folded(v) for v in diagram.edges[e])
            original_chords.setdefault(key, set()).add(e)
        lifted = set()
        for e in selected:
            u, v, kind = self.edges[e]
            if kind == 'a':
                items = original_chords[frozenset((u, v))]
                assert len(items) == 2
                lifted.update(items)
            else:
                original = u + 1
                lifted.update((original, (-original - 1) % diagram.n))
        return lifted


def from_point(diagram):
    field, t, c = diagram.field, diagram.t, diagram.c
    h = field.q // 2
    labels = [field.times[z][z] ^ field.times[t][z] ^ 1
              for z in diagram.cycle[1:h + 1]]
    assert labels[0] == 1 and labels[-1] == t
    return FoldedDiagram(field, labels, field.times[c][c ^ t])


def check_small(q):
    field = BinaryField(q)
    totals, seen = Counter(), set()
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
            partner = from_point(PointDiagram(field, t, c ^ t, cycle))
            assert (folded.labels, folded.delta, folded.edges) == (
                partner.labels, partner.delta, partner.edges)
            one_chord = {tuple(sorted(folded.lift(point, m)))
                         for m in folded.one_chord_matchings()}
            assert one_chord == {tuple(sorted(m))
                                 for m in point.reflected_two_chord_matchings()}
            totals['parameter_representatives'] += 1
            for selected in folded.matchings():
                totals['matchings'] += 1
                path = folded.root_path(selected)
                lifted = folded.lift(point, selected)
                profile = point.profile(lifted, independent=q <= 16)
                certificate = half_word_certificate(point, profile)
                predicted_beta = field.times[t][1 ^ path['gain']]
                assert certificate['beta'] == predicted_beta
                if path['gain'] == 1:
                    totals['root_identity'] += 1
                chord_endpoints = {v for e in selected if folded.edges[e][2] == 'a'
                                   for v in folded.edges[e][:2]}
                if chord_endpoints <= set(path['vertices']):
                    assert all(p['length'] % 2 == 0 for p in profile[1:])
                    assert point.good(profile) == (path['gain'] != 1)
                    totals['root_contained_matchings'] += 1
    print('FOLDED-AUDIT', q, dict(totals), flush=True)
    return {'q': q, 'counts': dict(totals)}


def check_pentagon_family(q):
    field, seen, tested = BinaryField(q), set(), 0
    for t in range(1, q):
        if t in seen:
            continue
        seen.update(field.frobenius_orbit(t))
        cycle = field.cycle(t)
        if len(cycle) != q + 1:
            continue
        for c in (1, t ^ 1):
            assert c not in (0, t)
            point = PointDiagram(field, t, c, cycle)
            profile = point.profile(set(range(1, q, 2)), independent=tested == 0)
            assert profile[0]['length'] == 5
            assert tuple(profile[0]['word_matrix']) == (1, t ^ 1, 0, 1)
            assert point.good(profile)
            tested += 1
    print('PENTAGON-FAMILY', q, 'representatives', tested, flush=True)
    return {'q': q, 'representatives': tested}


def check_schreier():
    records = []
    for q in (4, 8, 16, 32):
        field = BinaryField(q)
        t = next(t for t in range(1, q) if len(field.cycle(t)) == q + 1)
        for c in (t, 1):
            a = tuple([z ^ c for z in range(q)] + [q])
            x = tuple(field.act((0, 1, 1, t), z) for z in range(q + 1))
            direct, calculated = len(closure([a, x])), generated_order(field, t, c)
            assert direct == calculated
            records.append({'q': q, 't': t, 'c': c, 'order': calculated})
    print('SCHREIER-CONTROLS', len(records), 'agree with direct closure', flush=True)
    return records


def abstract_candidates(folded):
    out = []
    for selected in folded.one_chord_matchings():
        path = folded.root_path(selected)
        e = next(e for e in selected if folded.edges[e][2] == 'a')
        u, v, _ = folded.edges[e]
        out.append({'chord': [u + 1, v + 1], 'path': [v + 1 for v in path['vertices']],
                    'gain': path['gain'],
                    'both_endpoints_on_path': u in path['vertices'] and v in path['vertices']})
    return out


def abstract_control():
    field = BinaryField(16)
    labels = [1, 14, 6, 10, 13, 9, 5, 2]
    folded = FoldedDiagram(field, labels, 11)
    actual = from_point(PointDiagram(field, 2, 6, field.cycle(2)))
    assert set(labels) == set(actual.labels) and labels != actual.labels
    assert folded.delta == actual.delta
    canonical = folded.root_path(folded.canonical)
    candidates = abstract_candidates(folded)
    assert canonical['gain'] == 1 and len(candidates) == 1
    assert not any(c['gain'] != 1 and c['both_endpoints_on_path'] for c in candidates)
    result = {'q': 16, 'modulus': field.modulus, 't': 2, 'c': 6, 'delta': 11,
              'permuted_labels': labels, 'actual_labels': actual.labels,
              'canonical': canonical, 'candidates': candidates,
              'scope': 'abstract reordered labels, not a projective-cycle instance'}
    print('ABSTRACT-CONTROL', result, flush=True)
    return result


def actual_controls(paths):
    """Audit all root-contained candidates at stored canonical failures.

    A candidate can keep both endpoints on the path and STILL have gain 1.
    Such controls are rechecked by full projective permutations; generation
    of the first control in each field is checked by exact Schreier steps.
    """
    results = []
    for path in paths:
        for record in json.loads(path.read_text()):
            q = record['q']
            field = BinaryField(q)
            total, controls = 0, []
            for case in record['canonical_failures']:
                t, c = case['t'], case['c']
                point = PointDiagram(field, t, c, field.cycle(t))
                folded = from_point(point)
                assert folded.root_path(folded.canonical)['gain'] == 1
                for selected in folded.one_chord_matchings():
                    root = folded.root_path(selected)
                    e = next(e for e in selected if folded.edges[e][2] == 'a')
                    u, v, _ = folded.edges[e]
                    if not {u, v} <= set(root['vertices']):
                        continue
                    total += 1
                    lifted = folded.lift(point, selected)
                    profile = point.profile(lifted)
                    assert all(p['length'] % 2 == 0 for p in profile[1:])
                    assert half_word_certificate(point, profile)['beta'] == field.times[t][1 ^ root['gain']]
                    if root['gain'] != 1:
                        continue
                    point.profile(lifted, independent=True)
                    control = {'t': t, 'c': c, 'folded_chord': [u + 1, v + 1],
                               'matching': sorted(lifted), 'profile': profile,
                               'root_path': root, 'half_word': half_word_certificate(point, profile)}
                    if not controls:
                        control['generated_order'] = generated_order(field, t, c)
                        assert control['generated_order'] == q*(q*q - 1)
                    controls.append(control)
                    print('ACTUAL-IDENTITY-CONTROL', q, t, c, [u + 1, v + 1],
                          'root_length', profile[0]['length'], flush=True)
            result = {'q': q, 'modulus': field.modulus,
                      'root_contained_candidates': total, 'identity_controls': controls}
            print('ACTUAL-FOLDED-AUDIT', q, 'candidates', total, 'identity', len(controls), flush=True)
            results.append(result)
    return results


def random_abstract(q, samples, seed):
    field, rng = BinaryField(q), random.Random(seed)
    t = next(t for t in range(1, q) if len(field.cycle(t)) == q + 1)
    point = PointDiagram(field, t, next(c for c in range(1, q) if c != t), field.cycle(t))
    labels = from_point(point).labels
    interior = [k for k in labels if k not in (1, t)]
    totals = Counter()
    for _ in range(samples):
        rng.shuffle(interior)
        c = rng.choice([c for c in range(1, q) if c != t])
        folded = FoldedDiagram(field, [1] + interior + [t], field.times[c][c ^ t])
        totals['diagrams'] += 1
        if folded.root_path(folded.canonical)['gain'] != 1:
            continue
        totals['canonical_identity'] += 1
        for candidate in abstract_candidates(folded):
            if candidate['both_endpoints_on_path']:
                totals['root_contained_candidates'] += 1
                if candidate['gain'] == 1:
                    result = {'q': q, 'seed': seed, 'counts': dict(totals),
                              'counterexample': {'t': t, 'c': c, 'labels': folded.labels,
                                                 'candidate': candidate}}
                    print('ABSTRACT-STRONG-IDENTITY', result, flush=True)
                    return result
    result = {'q': q, 'seed': seed, 'counts': dict(totals), 'counterexample': None}
    print('ABSTRACT-RANDOM', result, flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('q', nargs='*', type=int, default=[4, 8, 16, 32])
    parser.add_argument('--abstract-samples', type=int, default=20000)
    parser.add_argument('--seed', type=int, default=20260906)
    parser.add_argument('--json-output', type=Path)
    parser.add_argument('--reflected-results', nargs='*', type=Path, default=[])
    args = parser.parse_args()
    result = {'small_audits': [check_small(q) for q in args.q],
              'pentagon_families': [check_pentagon_family(q)
                                    for q in (4, 8, 16, 32, 64, 128, 256, 512, 1024)],
              'schreier_controls': check_schreier(),
              'abstract_control': abstract_control(),
              'actual_controls': actual_controls(args.reflected_results),
              'abstract_random': [random_abstract(q, args.abstract_samples, args.seed)
                                  for q in (16, 32, 64)]}
    if args.json_output:
        args.json_output.write_text(json.dumps(result, indent=2) + '\n')
