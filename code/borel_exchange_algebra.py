"""Exact characteristic-two polynomial identities for short cycle words.

Sparse polynomials are in F_2[c,t]. This file uses no symbolic-algebra
dependency; division is exact and all products are checked explicitly.
"""

from itertools import product


class Polynomial:
    def __init__(self, terms=()):
        self.terms = frozenset(terms)

    def __add__(self, other):
        return Polynomial(self.terms ^ other.terms)

    def __mul__(self, other):
        result = set()
        for i, j in self.terms:
            for k, l in other.terms:
                term = (i + k, j + l)
                if term in result:
                    result.remove(term)
                else:
                    result.add(term)
        return Polynomial(result)

    def __pow__(self, exponent):
        result = ONE
        for _ in range(exponent):
            result = result * self
        return result

    def __eq__(self, other):
        return self.terms == other.terms

    def __repr__(self):
        return ' + '.join('*'.join(([f'c^{i}'] if i else []) + ([f't^{j}'] if j else [])) or '1'
                          for i, j in sorted(self.terms, reverse=True)) or '0'

    def divide(self, divisor):
        assert divisor.terms
        quotient, remainder, lead = ZERO, self, max(divisor.terms)
        while remainder.terms:
            a, b = max(remainder.terms)
            i, j = lead
            if a < i or b < j:
                return None
            term = Polynomial([(a - i, b - j)])
            quotient = quotient + term
            remainder = remainder + term * divisor
        assert quotient * divisor == self
        return quotient

    def evaluate(self, field, c, t):
        cp, tp = [1], [1]
        for _ in range(max((i for i, _ in self.terms), default=0)):
            cp.append(field.times[cp[-1]][c])
        for _ in range(max((j for _, j in self.terms), default=0)):
            tp.append(field.times[tp[-1]][t])
        result = 0
        for i, j in self.terms:
            result ^= field.times[cp[i]][tp[j]]
        return result


ZERO, ONE = Polynomial(), Polynomial([(0, 0)])
C, T = Polynomial([(1, 0)]), Polynomial([(0, 1)])
IDENTITY = (ONE, ZERO, ZERO, ONE)
A = (ONE, C, ZERO, ONE)
X = (ZERO, ONE, ONE, T)
XI = (T, ONE, ONE, ZERO)
R = (ONE, T, ZERO, ONE)
DELTA = C * (C + T)
EXCEPTION = DELTA**4 + T**3 * DELTA**2 + T**4 * DELTA + ONE + T**3 + T**4 + T**5


def multiply(left, right):
    a, b, c, d = left
    e, f, g, h = right
    return a * e + b * g, a * f + b * h, c * e + d * g, c * f + d * h


def short_word(signs, reflected):
    word = IDENTITY
    for sign in signs:
        word = multiply(A, word)
        word = multiply(X if sign == 1 else XI, word)
    return multiply(R, word) if reflected else word


def reversed_chord_resultant(word):
    """Eliminate z from W(z)=z and (denominator(W,z)+1) Q(z)=Delta.

    W=[a b;l d], Q=z^2+tz+1. For F=l*z^2+(a+d)*z+b and
    G=(l*z+d+1)*Q+Delta, the exact identity
      l*G+(l*z+(l*t+a+1))*F = alpha*z+beta
    yields the displayed necessary resultant. Zero l is not divided out.
    """
    a, b, l, d = word
    trace = a + d
    b2 = l * T + a + ONE
    b1 = l + T * (d + ONE) + b
    b0 = d + ONE + DELTA
    alpha, beta = l * b1 + trace * b2, l * b0 + b * b2
    return l * beta**2 + trace * alpha * beta + b * alpha**2


def trial_factors(poly):
    factors = [C, T, C + T, C + ONE, T + ONE, C + T + ONE,
               C**2 + C + ONE, T**2 + T + ONE,
               DELTA + ONE, DELTA + T, DELTA + T**2,
               DELTA + T + ONE, DELTA + T**2 + ONE]
    result = []
    for factor in factors:
        count = 0
        while poly != ZERO:
            quotient = poly.divide(factor)
            if quotient is None:
                break
            poly = quotient
            count += 1
        if count:
            result.append((factor, count))
    return result, poly


def verify_identities():
    """Prove the six-row trace/resultant table by exact polynomial equality."""
    u = C + T
    h0 = (DELTA + ONE)**2 + T**2 * (u + ONE)
    h1 = (DELTA + ONE)**2 + T**2 * (C + ONE)
    assert h0 * h1 == EXCEPTION
    records = []
    for reflected in (False, True):
        for signs in product((1, -1), repeat=3):
            word = short_word(signs, reflected)
            assert word[0] * word[3] + word[1] * word[2] == ONE
            p, b, l, d = word
            tau, j = p + d, l * T + p + ONE
            f = [b, tau, l, ZERO]
            g = [d + ONE + DELTA, l + T * (d + ONE), l * T + d + ONE, l]
            reduced = [l * g[k] + j * f[k] + (l * f[k - 1] if k else ZERO)
                       for k in range(4)]
            assert reduced == [l * (d + ONE + DELTA) + b * j,
                               l * (l + T * (d + ONE) + b) + tau * j, ZERO, ZERO]
            if signs[0] != signs[2]:
                kind = 'different_ends'
                expected_trace = C * (u + ONE)**2 if reflected else u * (C + ONE)**2
                expected_resultant = (C * u**2 * (DELTA + ONE) * h1 if reflected else
                                      C**2 * u * (DELTA + ONE) * h0)
            elif signs[0] != signs[1]:
                kind = 'alternating'
                expected_trace = C * (C + ONE)**2 if reflected else u * (C + ONE)**2
                expected_resultant = (C * u**2 * (C + ONE)**6 if reflected else
                                      C**2 * u * (C + ONE)**4 * (u + ONE)**2)
            else:
                kind = 'constant'
                expected_trace = C * (u + ONE)**2 if reflected else u * (u + ONE)**2
                expected_resultant = (C * u**2 * (C + ONE)**2 * (u + ONE)**4 if reflected else
                                      C**2 * u * (u + ONE)**6)
            assert word[0] + word[3] == expected_trace
            assert reversed_chord_resultant(word) == expected_resultant
            records.append({'signs': signs, 'reflected': reflected, 'kind': kind,
                            'trace': repr(expected_trace), 'resultant': repr(expected_resultant)})
    return records


if __name__ == '__main__':
    assert len(verify_identities()) == 16
    for reflected in (False, True):
        for signs in product((1, -1), repeat=3):
            word = short_word(signs, reflected)
            trace = word[0] + word[3]
            print('WORD', signs, reflected, 'trace', trial_factors(trace),
                  'resultant', trial_factors(reversed_chord_resultant(word)))
