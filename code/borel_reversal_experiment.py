"""Test the reversal-symmetric Borel construction in PSL(2,2^m).

Use x(z)=1/(z+t), a(z)=z+c, with x a full projective-line cycle.
The canonical matching takes the semi-edge at infinity and alternate cycle
edges on the other q points. Its unique odd complementary circuit has
unipotent or identity monodromy, as proved in the survey. The latter occurs
at q=64. Adaptive quotient matchings repair all twelve failures there.

python3 code/borel_reversal_experiment.py --json-output code/borel_reversal_results.json
python3 code/borel_reversal_experiment.py 64 --verify-generation
"""

import argparse
import json
from collections import Counter
from pathlib import Path

from cayley_snark_check import GF, closure

GF.IRRED[4] = (2, [1, 1, 1])
IDENTITY = (1, 0, 0, 1)


class Field:
    def __init__(self, q):
        self.q = q
        field = GF(q)
        assert field.p == 2
        self.times = [[field.mulf(i, j) for j in range(q)] for i in range(q)]
        self.inverse = field.inv_table

    def multiply(self, m, n):
        a, b, c, d = m
        e, f, g, h = n
        t = self.times
        return (t[a][e] ^ t[b][g], t[a][f] ^ t[b][h],
                t[c][e] ^ t[d][g], t[c][f] ^ t[d][h])

    def act(self, matrix, z):
        a, b, c, d = matrix
        numerator = a if z == self.q else self.times[a][z] ^ b
        denominator = c if z == self.q else self.times[c][z] ^ d
        return self.times[numerator][self.inverse[denominator]] if denominator else self.q


class PointDiagram:
    def __init__(self, field, t, c, cycle):
        self.field, self.t, self.c = field, t, c
        self.n = field.q + 1
        self.cycle = cycle
        self.position = {z: i for i, z in enumerate(cycle)}
        self.x = (0, 1, 1, t)
        self.x_inverse = (t, 1, 1, 0)
        self.a = (1, c, 0, 1)
        self.edges = [(i, (i + 1) % self.n) for i in range(self.n)]
        for i in range(1, self.n):
            j = self.position[cycle[i] ^ c]
            if i < j:
                self.edges.append((i, j))
        self.inc = [[] for _ in cycle]
        for e, (u, v) in enumerate(self.edges):
            self.inc[u].append(e)
            self.inc[v].append(e)
        assert len(self.inc[0]) == 2
        assert all(len(self.inc[v]) == 3 for v in range(1, self.n))

    def profile(self, selected, independent=False):
        assert all(sum(e in selected for e in self.inc[v]) == int(v != 0)
                   for v in range(self.n))
        unseen, result = set(range(self.n)), []
        while unseen:
            start = current = min(unseen)
            previous, word, labels = -1, IDENTITY, []
            while True:
                unseen.remove(current)
                e = next(e for e in self.inc[current] if e != previous and e not in selected)
                u, v = self.edges[e]
                if e >= self.n:
                    label = self.a
                    symbol = 'a'
                else:
                    label = self.x if current == u else self.x_inverse
                    symbol = 'x' if current == u else 'X'
                word = self.field.multiply(label, word)
                labels.append(symbol)
                current, previous = (v if current == u else u), e
                if current == start:
                    break
            if independent:
                # Compose the induced permutations point by point separately.
                permutation = list(range(self.n))
                for symbol in labels:
                    matrix = {'a': self.a, 'x': self.x, 'X': self.x_inverse}[symbol]
                    permutation = [self.field.act(matrix, z) for z in permutation]
                assert permutation == [self.field.act(word, z) for z in range(self.n)]
            assert self.field.act(word, self.cycle[start]) == self.cycle[start]
            result.append({'length': len(labels), 'word_matrix': word,
                           'labels': ''.join(labels)})
        return result

    @staticmethod
    def good(profile):
        return all(p['length'] % 2 == 0 or
                   (tuple(p['word_matrix']) != IDENTITY and
                    p['word_matrix'][0] == p['word_matrix'][3]) for p in profile)

    def adaptive(self, max_attempts=100000):
        attempts = 0

        def search(uncovered, selected):
            nonlocal attempts
            if not uncovered:
                attempts += 1
                if attempts > max_attempts:
                    raise RuntimeError('adaptive attempt limit reached')
                profile = self.profile(selected)
                if self.good(profile):
                    self.profile(selected, independent=True)
                    return {'attempts': attempts, 'matching': sorted(selected),
                            'profile': profile}
                return None
            vertex = min(uncovered)
            for e in self.inc[vertex]:
                covered = set(self.edges[e])
                if covered <= uncovered:
                    result = search(uncovered - covered, selected | {e})
                    if result is not None:
                        return result
            return None

        return search(set(range(1, self.n)), set())

    def reflected_two_chord_matchings(self):
        """Generate every reflected chord pair completed by cycle edges.

        This is a structurally specified linear-size candidate family, not
        a search over all quotient perfect matchings. Every candidate keeps
        the semi-edge at infinity and is invariant under i -> -i.
        """
        chord_at = {frozenset(edge): e for e, edge in enumerate(self.edges)
                    if e >= self.n}
        for e in range(self.n, len(self.edges)):
            reflected = frozenset((-v) % self.n for v in self.edges[e])
            other = chord_at[reflected]
            if e >= other:
                continue
            endpoints = set(self.edges[e]) | set(self.edges[other])
            if len(endpoints) != 4:
                continue
            selected = {e, other}
            uncovered = set(range(1, self.n)) - endpoints
            while uncovered:
                v = min(uncovered)
                if v + 1 not in uncovered:
                    break
                selected.add(v)
                uncovered -= {v, v + 1}
            if uncovered:
                continue
            reflected_matching = {(-item - 1) % self.n if item < self.n
                                  else chord_at[frozenset((-v) % self.n
                                                         for v in self.edges[item])]
                                  for item in selected}
            assert reflected_matching == selected
            yield selected

    def two_chord_repairs(self, independent=True):
        tested, good, witness = 0, 0, None
        for selected in self.reflected_two_chord_matchings():
            tested += 1
            profile = self.profile(selected, independent=independent)
            # The component through the unique fixed point is reversed by r.
            # Other odd components, if any, still need their own word checks.
            first = profile[0]
            assert first['length'] % 2 == 1
            assert first['word_matrix'][0] == first['word_matrix'][3] == 1
            assert first['word_matrix'][2] == 0
            if self.good(profile):
                good += 1
                if witness is None:
                    witness = {'matching': sorted(selected), 'profile': profile}
        return {'tested': tested, 'good': good, 'witness': witness}


def check(q, verify_generation):
    field = Field(q)
    counts, failures = Counter(), []
    for t in range(1, q):
        x = (0, 1, 1, t)
        cycle, current = [], q
        while current not in cycle:
            cycle.append(current)
            current = field.act(x, current)
        if len(cycle) != q + 1:
            continue
        for c in range(1, q):
            diagram = PointDiagram(field, t, c, cycle)
            profile = diagram.profile(set(range(1, q, 2)), independent=True)
            odd = [p for p in profile if p['length'] % 2]
            assert len(odd) == 1
            word = odd[0]['word_matrix']
            assert word[0] == word[3] == 1 and word[2] == 0
            if diagram.good(profile):
                counts['canonical_good'] += 1
            else:
                assert word == IDENTITY
                counts['canonical_identity_failure'] += 1
                repaired = diagram.adaptive()
                assert repaired is not None
                counts['adaptive_repaired'] += 1
                failure = {'q': q, 't': t, 'c': c, 'canonical_odd': odd[0],
                           'repaired': repaired,
                           'two_chord_repairs': diagram.two_chord_repairs()}
                if verify_generation and not failures:
                    a_perm = tuple([z ^ c for z in range(q)] + [q])
                    x_perm = tuple(field.act(x, z) for z in range(q + 1))
                    generated_order = len(closure([a_perm, x_perm]))
                    assert generated_order == q * (q*q - 1)
                    failure['independently_checked_generated_order'] = generated_order
                failures.append(failure)
                print('IDENTITY-FAILURE', q, t, c, 'length', odd[0]['length'],
                      'repair_attempts', repaired['attempts'], flush=True)
    print('BOREL-REVERSAL', q, dict(counts), flush=True)
    return {'q': q, 'counts': dict(counts), 'failures': failures}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('q', nargs='*', type=int, default=[4, 8, 16, 32, 64])
    parser.add_argument('--verify-generation', action='store_true')
    parser.add_argument('--json-output', type=Path)
    args = parser.parse_args()
    records = [check(q, args.verify_generation) for q in args.q]
    if args.json_output:
        args.json_output.write_text(json.dumps(records, indent=2) + '\n')
