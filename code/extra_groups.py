#!/usr/bin/env python3
"""
The three groups of order < 12 200 admissible for a smallest Cayley snark
(by the Nedela-Skoviera theorem combined with the census) that are not
covered by cayley_snark_check.py: PSU(3,3), G2(2) = PSU(3,3).2 and
PSL(3,3).2.  Constructed as permutation groups and checked with the same
routine.  Writes results_extra.json.
"""
import itertools, json, random, sys
import cayley_snark_check as C

# ------------------------------------------------------------ GF(9) helpers
F9 = C.GF(9)
def conj9(a):            # Frobenius x -> x^3 on GF(9)
    return F9.mulf(F9.mulf(a, a), a)
def dot_h(u, v):         # Hermitian form  sum u_i * conj(v_i)
    r = 0
    for a, b in zip(u, v):
        r = F9.add(r, F9.mulf(a, conj9(b)))
    return r
def vec_add(u, v): return tuple(F9.add(a, b) for a, b in zip(u, v))
def vec_scale(c, u): return tuple(F9.mulf(c, a) for a in u)
def normalise(v):
    for c in v:
        if c:
            ci = F9.inv_table[c]
            return tuple(F9.mulf(ci, a) for a in v)
    raise ValueError

# isotropic points of PG(2,9): 28 of them
iso = []
seen = set()
for v in itertools.product(range(9), repeat=3):
    if v == (0, 0, 0): continue
    if dot_h(v, v) != 0: continue
    w = normalise(v)
    if w not in seen:
        seen.add(w); iso.append(w)
assert len(iso) == 28, len(iso)
iso_idx = {v: i for i, v in enumerate(iso)}

def mat_apply(M, v):
    return tuple(
        (lambda s: s)(
            __import__('functools').reduce(F9.add, (F9.mulf(M[i][j], v[j]) for j in range(3)), 0))
        for i in range(3))

def perm_of_matrix(M):
    return tuple(iso_idx[normalise(mat_apply(M, v))] for v in iso)

def det(M):
    a, b, c = M[0]; d, e, f = M[1]; g, h, i = M[2]
    m = F9.mulf; ad = F9.add; ng = F9.neg
    return ad(ad(m(a, ad(m(e, i), ng(m(f, h)))), ng(m(b, ad(m(d, i), ng(m(f, g)))))),
              m(c, ad(m(d, h), ng(m(e, g)))))

def random_unitary():
    """Random matrix preserving the Hermitian form (columns orthonormal)."""
    cols = []
    while len(cols) < 3:
        v = tuple(random.randrange(9) for _ in range(3))
        # orthogonalise against previous columns (they are orthonormal)
        for c in cols:
            lam = dot_h(v, c)            # coefficient so that v - lam*c is orthogonal to c
            v = vec_add(v, vec_scale(F9.neg(lam), c))
        n = dot_h(v, v)                  # lies in GF(3) = {0,1,2}
        if n == 0: continue
        if n != 1:
            # need c with c*conj(c) = n^{-1}; n in {1,2}: find c with N(c) = 2
            for c in range(1, 9):
                if F9.mulf(c, conj9(c)) == n:   # then (c^{-1} v) has norm 1
                    v = vec_scale(F9.inv_table[c], v); break
            else:
                continue
        assert dot_h(v, v) == 1
        cols.append(v)
    M = [[cols[j][i] for j in range(3)] for i in range(3)]
    return M

def psu33():
    random.seed(1)
    gens = []
    while len(gens) < 6:
        M = random_unitary()
        if det(M) == 1:
            gens.append(perm_of_matrix(M))
    G = C.closure(gens)
    assert len(G) == 6048, len(G)
    return G

def g22():
    G0 = psu33()
    frob = tuple(iso_idx[normalise(tuple(conj9(a) for a in v))] for v in iso)
    G = C.closure(list(G0)[:50] + [frob])   # a few elements of PSU(3,3) generate it
    assert len(G) == 12096, len(G)
    return G

# ------------------------------------------------------------ PSL(3,3).2
def psl33_2():
    p = 3
    pts = []
    seen = set()
    for v in itertools.product(range(p), repeat=3):
        if v == (0, 0, 0): continue
        for c in v:
            if c: break
        w = tuple((x * (1 if c == 1 else 2)) % p for x in v)
        if w not in seen:
            seen.add(w); pts.append(w)
    assert len(pts) == 13
    idx = {v: i for i, v in enumerate(pts)}
    def norm(w):
        for c in w:
            if c: break
        return tuple((x * (1 if c == 1 else 2)) % p for x in w)
    def matvec(M, v):
        return norm(tuple(sum(M[i][j] * v[j] for j in range(3)) % p for i in range(3)))
    def transpose(M): return [[M[j][i] for j in range(3)] for i in range(3)]
    def inverse(M):
        # brute force over GL(3,3) is too big; use adjugate: det in {1,2}
        a, b, c = M[0]; d, e, f = M[1]; g, h, i = M[2]
        D = (a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g)) % p
        Di = 1 if D == 1 else 2
        adj = [[(e*i - f*h), -(b*i - c*h), (b*f - c*e)],
               [-(d*i - f*g), (a*i - c*g), -(a*f - c*d)],
               [(d*h - e*g), -(a*h - b*g), (a*e - b*d)]]
        return [[(Di * adj[r][s]) % p for s in range(3)] for r in range(3)]
    # points 0..12, lines 13..25 (line with coordinates l = set of points v with l.v = 0)
    def act(M):
        Minv_t = transpose(inverse(M))
        out = [idx[matvec(M, v)] for v in pts]            # points: v -> M v
        out += [13 + idx[matvec(Minv_t, l)] for l in pts]  # lines:  l -> M^{-T} l
        return tuple(out)
    polarity = tuple(list(range(13, 26)) + list(range(13)))  # v <-> line with the same coordinates
    gens = [act([[1, 1, 0], [0, 1, 0], [0, 0, 1]]),
            act([[0, 0, 1], [1, 0, 0], [0, 1, 0]]),
            act([[1, 0, 0], [0, 1, 1], [0, 0, 1]]),
            polarity]
    G = C.closure(gens)
    assert len(G) == 11232, len(G)
    return G

if __name__ == "__main__":
    C.petersen_control()
    import multiprocessing as mp
    mp.set_start_method("fork")
    summary = []
    for name, ctor in [("PSU(3,3)", psu33), ("PSL(3,3).2", psl33_2), ("G2(2)", g22)]:
        G = ctor()
        N, ncls, tested, nbad, nund, results = C.check_group(name, G, pool=True)
        summary.append(dict(group=name, order=N, involution_classes=ncls,
                            pairs_tested=tested, non_colourable=nbad, undecided=nund,
                            methods=sorted(set(r[3].split(" (")[0] for r in results)),
                            orders_x=sorted(set(r[0] for r in results)),
                            pairs=sorted(set((r[0], r[1]) for r in results))))
        with open("results_extra.json", "w") as f:
            json.dump(summary, f, indent=1)
