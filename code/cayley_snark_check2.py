#!/usr/bin/env python3
r"""
Second-generation 3-edge-colourability check for cubic Cayley graphs

    X(G;a,x) = Cay(G, {a, x, x^-1}),   a^2 = 1,  ord(x) odd,  <a,x> = G.

Difference from cayley_snark_check.py.  There the Cayley graph X was reduced
by the quotient X -> H\X for an ABELIAN subgroup H of ODD order containing no
conjugate of x.  Here H may be ANY subgroup containing no conjugate of x:
the quotient H\X (left cosets Hg, edges Hg -- Hga and Hg -- Hgx) is then a
cubic *pregraph*: it is loopless (no conjugate of x in H) but it may have
semi-edges, namely at the cosets Hg with g a g^-1 in H.  A 3-edge-colouring
of the pregraph in which the semi-edges are coloured as well lifts to a
3-edge-colouring of X (Lemma "generalised quotient lemma" in the survey).  So
H can be a large subgroup -- a point stabiliser, a set stabiliser, a
centraliser or normaliser -- and the SAT instance is tiny.  The lifted
colouring is verified independently on the full Cayley graph.

For each group we enumerate all pairs (a, x) up to conjugation in G and
inversion of x, with a an involution, x of odd order and <a, x> = G.  The
generation test uses the socle T of G (T = G for simple G, |G:T| = 2
otherwise): <a,x> = G iff a is outside T and <x, x^a> = T (or <a,x> = T when
T = G), and a subgroup of T is T as soon as it is larger than |T|/m, where m
is the minimal index of a proper subgroup of T (the minimal degree, from the
ATLAS).

Groups are permutation groups, given by generators (constructed here or read
from the ATLAS of Finite Group Representations in MeatAxe text format).  The
group order is checked against the known value in every case.

Run:  python3 cayley_snark_check2.py [max_order] [group names ...]
"""
import sys, os, time, json, math, random, itertools, pickle, select
from pysat.solvers import Cadical153
import cayley_snark_check as C

mul, inv, order, closure = C.mul, C.inv, C.order, C.closure

# ------------------------------------------------------------------ groups
ATLAS_DIR = os.environ.get("ATLAS_DIR", os.path.join(os.path.dirname(__file__), "atlas"))

def load_meataxe(path):
    """Permutation in MeatAxe text format (header '12 1 n 1', then images, 1-based)."""
    toks = open(path).read().split()
    assert toks[0] == "12", path
    n = int(toks[2])
    imgs = [int(t) - 1 for t in toks[4:4 + n]]
    assert sorted(imgs) == list(range(n)), path
    return tuple(imgs)

def atlas_group(name, expected):
    gens = [load_meataxe(os.path.join(ATLAS_DIR, name + ".m1")),
            load_meataxe(os.path.join(ATLAS_DIR, name + ".m2"))]
    G = closure(gens)
    assert len(G) == expected, (name, len(G), expected)
    return G

def wreath_Z2(S, expected):
    """S wr Z_2 for a permutation group S on n points, acting on 2n points."""
    S = list(S)
    n = len(S[0])
    gens = []
    # a few generators of S (random elements generate S with high probability;
    # we verify the order at the end)
    rnd = random.Random(1)
    for _ in range(4):
        s = rnd.choice(S)
        gens.append(tuple(list(s) + list(range(n, 2 * n))))
    gens.append(tuple(list(range(n, 2 * n)) + list(range(n))))
    G = closure(gens)
    assert len(G) == expected, (len(G), expected)
    return G

def direct_with_frob(q, e, expected):
    return C.psl2(q, frob=e, expected=expected)

# GF(4) is needed for PSL(3,4)
C.GF.IRRED[4] = (2, [1, 1, 1])          # x^2 + x + 1
C.GF.IRRED[81] = (3, [2, 1, 0, 0, 1])   # x^4 + x + 2 (irreducible over GF(3), checked below)

def psl3(q, expected, field_aut=None, graph_aut=False):
    """PSL(3,q) (q = 4 here) acting on the points (and, if graph_aut, also
    the lines) of PG(2,q).  field_aut = e adds x -> x^e on coordinates;
    graph_aut adds the inverse-transpose (point-line duality)."""
    F = C.GF(q)
    add, mulf, neg = F.add, F.mulf, F.neg
    def normalise(v):
        for c in v:
            if c:
                ci = F.inv_table[c]
                return tuple(mulf(ci, a) for a in v)
        raise ValueError
    pts, seen = [], set()
    for v in itertools.product(range(q), repeat=3):
        if v == (0, 0, 0):
            continue
        w = normalise(v)
        if w not in seen:
            seen.add(w); pts.append(w)
    idx = {v: i for i, v in enumerate(pts)}
    npts = len(pts)
    def dot(u, v):
        r = 0
        for a, b in zip(u, v):
            r = add(r, mulf(a, b))
        return r
    # lines = points of the dual plane; line w = {v : w.v = 0}
    def apply(M, v):
        return tuple(functools_reduce(add, (mulf(M[i][j], v[j]) for j in range(3)), 0) for i in range(3))
    import functools
    def functools_reduce(f, it, init):
        return functools.reduce(f, it, init)
    def transpose(M):
        return [[M[j][i] for j in range(3)] for i in range(3)]
    def det(M):
        a, b, c = M[0]; d, e, f = M[1]; g, h, i = M[2]
        return add(add(mulf(a, add(mulf(e, i), neg(mulf(f, h)))),
                       neg(mulf(b, add(mulf(d, i), neg(mulf(f, g)))))),
                   mulf(c, add(mulf(d, h), neg(mulf(e, g)))))
    def adjugate(M):
        # cofactor matrix transposed
        def minor(i, j):
            rows = [r for r in range(3) if r != i]; cols = [c for c in range(3) if c != j]
            return add(mulf(M[rows[0]][cols[0]], M[rows[1]][cols[1]]),
                       neg(mulf(M[rows[0]][cols[1]], M[rows[1]][cols[0]])))
        A = [[None] * 3 for _ in range(3)]
        for i in range(3):
            for j in range(3):
                m = minor(j, i)
                A[i][j] = m if (i + j) % 2 == 0 else neg(m)
        return A
    def perm_of_matrix(M):
        # points: v -> Mv ; lines: w -> (M^-T) w  (so that incidence is preserved)
        P = [idx[normalise(apply(M, v))] for v in pts]
        if graph_aut:
            Mt = transpose(adjugate(M))       # proportional to M^-T
            P += [npts + idx[normalise(apply(Mt, w))] for w in pts]
        return tuple(P)
    gens = []
    z = F.zeta
    E = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    T1 = [[1, 1, 0], [0, 1, 0], [0, 0, 1]]
    T2 = [[1, 0, 0], [0, 1, 1], [0, 0, 1]]
    T3 = [[1, 0, 0], [0, 1, 0], [1, 0, 1]]
    Cy = [[0, 0, 1], [1, 0, 0], [0, 1, 0]]
    D = [[z, 0, 0], [0, F.inv_table[z], 0], [0, 0, 1]]
    for M in (T1, T2, T3, Cy, D):
        gens.append(perm_of_matrix(M))
    frob = dual = None
    if field_aut is not None:
        def fp(x):
            r = 1
            for _ in range(field_aut):
                r = mulf(r, x)
            return r
        P = [idx[normalise(tuple(fp(c) for c in v))] for v in pts]
        if graph_aut:
            P += [npts + idx[normalise(tuple(fp(c) for c in w))] for w in pts]
        frob = tuple(P)
    if graph_aut:
        # duality: point v -> line v, line w -> point w
        dual = tuple(list(range(npts, 2 * npts)) + list(range(npts)))
    if frob is not None and dual is not None:
        gens.append(mul(frob, dual))        # the graph-field automorphism only
    elif frob is not None:
        gens.append(frob)
    elif dual is not None:
        gens.append(dual)
    G = closure(gens)
    assert len(G) == expected, (len(G), expected)
    return G

def psp4_3(similitude=False):
    """PSp(4,3) = PSU(4,2) acting on the 40 points of PG(3,3); with
    similitude, PGSp(4,3) = PSU(4,2).2 = W(E_6)."""
    p = 3
    def normalise(v):
        for c in v:
            if c:
                ci = 1 if c == 1 else 2
                return tuple((ci * a) % p for a in v)
        raise ValueError
    pts, seen = [], set()
    for v in itertools.product(range(p), repeat=4):
        if v == (0, 0, 0, 0):
            continue
        w = normalise(v)
        if w not in seen:
            seen.add(w); pts.append(w)
    idx = {v: i for i, v in enumerate(pts)}
    def perm(M):
        return tuple(idx[normalise(tuple(sum(M[i][j] * v[j] for j in range(4)) % p for i in range(4)))] for v in pts)
    # symplectic form J = [[0,I],[-I,0]]; generators of Sp(4,3): transvections
    # x -> x + c<x,u>u for various u (these generate Sp), we take a few and check the order
    def transvection(u, c):
        # J u
        Ju = ((u[2]) % p, (u[3]) % p, (-u[0]) % p, (-u[1]) % p)
        # <x,u> = x^T J u ; T x = x + c <x,u> u
        M = [[(int(i == j) + c * u[i] * Ju[j]) % p for j in range(4)] for i in range(4)]
        return M
    gens = [perm(transvection(u, 1)) for u in [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1), (1, 1, 0, 0), (1, 0, 1, 0), (0, 1, 0, 1), (1, 1, 1, 1)]]
    if similitude:
        gens.append(perm([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 2, 0], [0, 0, 0, 2]]))   # multiplier -1
    G = closure(gens)
    assert len(G) == (51840 if similitude else 25920), len(G)
    return G

def psu3(q, expected, field_aut=None):
    """PSU(3,q) acting on the isotropic points of PG(2,q^2) (q^3+1 of them);
    generated by random unitary matrices (order checked)."""
    Q = q * q
    F = C.GF(Q)
    add, mulf, neg = F.add, F.mulf, F.neg
    def conj(a):
        r = 1
        for _ in range(q):
            r = mulf(r, a)
        return r
    def dot_h(u, v):
        r = 0
        for a, b in zip(u, v):
            r = add(r, mulf(a, conj(b)))
        return r
    def normalise(v):
        for c in v:
            if c:
                ci = F.inv_table[c]
                return tuple(mulf(ci, a) for a in v)
        raise ValueError
    iso, seen = [], set()
    for v in itertools.product(range(Q), repeat=3):
        if v == (0, 0, 0) or dot_h(v, v) != 0:
            continue
        w = normalise(v)
        if w not in seen:
            seen.add(w); iso.append(w)
    assert len(iso) == q ** 3 + 1, len(iso)
    idx = {v: i for i, v in enumerate(iso)}
    def apply(M, v):
        out = []
        for i in range(3):
            s = 0
            for j in range(3):
                s = add(s, mulf(M[i][j], v[j]))
            out.append(s)
        return tuple(out)
    def perm(M):
        return tuple(idx[normalise(apply(M, v))] for v in iso)
    rnd = random.Random(7)
    # norm-1 scalars for rescaling
    def random_unitary():
        cols = []
        while len(cols) < 3:
            v = tuple(rnd.randrange(Q) for _ in range(3))
            for c in cols:
                lam = dot_h(v, c)
                v = tuple(add(a, neg(mulf(lam, b))) for a, b in zip(v, c))
            nrm = dot_h(v, v)
            if nrm == 0:
                continue
            # find c with c*conj(c) = nrm, then v/c has norm 1
            for c in range(1, Q):
                if mulf(c, conj(c)) == nrm:
                    v = tuple(mulf(F.inv_table[c], a) for a in v); break
            else:
                continue
            assert dot_h(v, v) == 1
            cols.append(v)
        return [[cols[j][i] for j in range(3)] for i in range(3)]
    gens = [perm(random_unitary()) for _ in range(3)]
    psu_order = expected if field_aut is None else expected // 2
    G = closure(gens)
    if len(G) != psu_order:
        # the unitary matrices generate PGU(3,q), of order gcd(3,q+1)*|PSU(3,q)|; PSU(3,q) is
        # its derived subgroup, generated by commutators of a few random elements
        assert len(G) == 3 * psu_order, (len(G), psu_order)
        Gl = list(G)
        gens = []
        for _ in range(6):
            g, h = rnd.choice(Gl), rnd.choice(Gl)
            gens.append(mul(mul(mul(inv(g), inv(h)), g), h))
        G = closure(gens)
        assert len(G) == psu_order, (len(G), psu_order)
    if field_aut is not None:
        def fp(x):
            r = 1
            for _ in range(field_aut):
                r = mulf(r, x)
            return r
        gens = gens + [tuple(idx[normalise(tuple(fp(c) for c in v))] for v in iso)]
        G = closure(gens)
    assert len(G) == expected, (len(G), expected)
    return G

def check_psl2_min_index(q):
    return {5: 5, 7: 7, 9: 6, 11: 11}.get(q, q + 1)

# ------------------------------------------------------------------ SAT on pregraphs
def colour_pregraph(n, items):
    """3-edge-colour a cubic pregraph.  items[i] = (u, v) for an edge, (u, None)
    for a semi-edge.  Returns (True, colours) or (False, None).
    Encoding: a perfect matching M (colour 2, semi-edges allowed) such that the
    remaining items can be properly 2-coloured (colours 0/1) -- i.e. any two
    non-M items at a vertex get different colours."""
    inc = [[] for _ in range(n)]
    for i, (u, v) in enumerate(items):
        inc[u].append(i)
        if v is not None:
            inc[v].append(i)
    assert all(len(L) == 3 for L in inc)
    m = lambda i: i + 1
    c = lambda i: len(items) + i + 1
    s = Cadical153()
    for v in range(n):
        L = inc[v]
        s.add_clause([m(i) for i in L])
        for i, j in itertools.combinations(L, 2):
            s.add_clause([-m(i), -m(j)])
            s.add_clause([m(i), m(j), c(i), c(j)])
            s.add_clause([m(i), m(j), -c(i), -c(j)])
    s.add_clause([m(inc[0][0])])
    s.add_clause([-c(inc[0][1])])
    r = s.solve()
    if not r:
        s.delete()
        return False, None
    model = set(l for l in s.get_model() if l > 0)
    s.delete()
    col = [2 if m(i) in model else (1 if c(i) in model else 0) for i in range(len(items))]
    check_pregraph_colouring(n, items, col)
    return True, col

def check_pregraph_colouring(n, items, col):
    at = [set() for _ in range(n)]
    for i, (u, v) in enumerate(items):
        assert col[i] in (0, 1, 2)
        assert col[i] not in at[u]; at[u].add(col[i])
        if v is not None:
            assert col[i] not in at[v]; at[v].add(col[i])
    assert all(len(s) == 3 for s in at)

def run_with_timeout(fn, args, time_limit):
    """CaDiCaL cannot be interrupted from Python: run fn(*args) in a forked
    child and kill it after time_limit seconds.  Returns fn's result or None."""
    r, w = os.pipe()
    pid = os.fork()
    if pid == 0:
        os.close(r)
        try:
            data = pickle.dumps(fn(*args))
        except Exception as ex:
            data = pickle.dumps(("error", repr(ex)))
        while data:
            k = os.write(w, data); data = data[k:]
        os.close(w)
        os._exit(0)
    os.close(w)
    chunks, deadline = [], time.time() + time_limit
    while True:
        remaining = deadline - time.time()
        if remaining <= 0:
            os.kill(pid, 9); os.waitpid(pid, 0); os.close(r)
            return None
        rl, _, _ = select.select([r], [], [], remaining)
        if rl:
            chunk = os.read(r, 1 << 16)
            if not chunk:
                break
            chunks.append(chunk)
    os.close(r); os.waitpid(pid, 0)
    out = pickle.loads(b"".join(chunks))
    if isinstance(out, tuple) and out and out[0] == "error":
        raise RuntimeError(out[1])
    return out

# ------------------------------------------------------------------ context (shared with forked workers)
class Ctx:
    pass
K = Ctx()

def setup_group(name, G, min_index, socle_index):
    t0 = time.time()
    elems = list(G)
    index = {g: i for i, g in enumerate(elems)}
    N = len(elems)
    n = len(elems[0])
    e = tuple(range(n))
    ords = {g: order(g) for g in elems}
    K.name, K.elems, K.index, K.N, K.n, K.e, K.ords = name, elems, index, N, n, e, ords
    K.min_index = min_index
    # odd-order conjugacy classes
    odd = [g for g in elems if ords[g] % 2 == 1 and ords[g] > 1]
    class_of = {}
    classes = []
    for y in odd:
        if y in class_of:
            continue
        cls = set()
        for g in elems:
            cls.add(mul(mul(inv(g), y), g))
        cid = len(classes)
        for z in cls:
            class_of[z] = cid
        classes.append((y, sorted(cls, key=index.get)))
    K.odd_classes, K.class_of = classes, class_of
    # involution classes
    invols = [g for g in elems if ords[g] == 2]
    inv_reps, unassigned = [], set(invols)
    while unassigned:
        b = next(iter(unassigned))
        cls = {mul(mul(inv(g), b), g) for g in elems}
        unassigned -= cls
        inv_reps.append((b, len(cls)))
    K.inv_reps = inv_reps
    # socle
    if socle_index == 1:
        K.T, K.Tsize = None, N
    else:
        T = closure([y for y, _ in classes])
        assert len(T) * socle_index == N, (len(T), N)
        K.T, K.Tsize = frozenset(T), len(T)
    K.cands = subgroup_candidates()
    # admissible candidates per odd class: no element of the class in H
    K.adm = []
    for cid, (y, cls) in enumerate(classes):
        cs = set(cls)
        L = [hi for hi, (hname, H) in enumerate(K.cands) if not any(h in cs for h in H)]
        K.adm.append(L)
    print(f"  [{name}] |G|={N} degree={n} odd classes={len(classes)} "
          f"involution classes={len(inv_reps)} subgroup candidates={len(K.cands)} "
          f"(largest {max(len(H) for _, H in K.cands) if K.cands else 0})  setup {time.time()-t0:.1f}s", flush=True)

def small_generating_set(H, rnd):
    Hs = set(H)
    gens = []
    cur = {K.e}
    while len(cur) < len(Hs):
        h = rnd.choice(H)
        if h in cur:
            continue
        gens.append(h)
        cur = closure(gens)
    return gens

def normaliser(H):
    rnd = random.Random(len(H))
    gens = small_generating_set(H, rnd)
    Hs = set(H)
    return [g for g in K.elems if all(mul(mul(inv(g), h), g) in Hs for h in gens)]

def stabiliser(pred):
    return [g for g in K.elems if pred(g)]

def subgroup_candidates():
    elems, N, n = K.elems, K.N, K.n
    cands = {}
    def add(name, H):
        H = list(H)
        if 1 < len(H) < N:
            fs = frozenset(H)
            if fs not in cands:
                cands[fs] = (name, H)
    # (a) centralisers and normalisers attached to odd-order classes
    for y, cls in K.odd_classes:
        o = K.ords[y]
        Cy = [g for g in elems if mul(g, y) == mul(y, g)]
        add(f"C(x{o})", Cy)
        cyc = list(closure([y]))
        add(f"Z{o}", cyc)
        add(f"N(Z{o})", normaliser(cyc))
        add(f"N(C(x{o}))", normaliser(Cy))
    # (b) centralisers of involutions
    for b, _ in K.inv_reps:
        Cb = [g for g in elems if mul(g, b) == mul(b, g)]
        add("C(inv)", Cb)
        add("N(C(inv))", normaliser(Cb))
    # (c) stabilisers of structures in the natural action
    orbits, seen = [], set()
    for p in range(n):
        if p in seen:
            continue
        orb = {p}
        frontier = [p]
        while frontier:
            new = []
            for q in frontier:
                for g in elems[:50] + [elems[i] for i in range(0, N, max(1, N // 50))]:
                    r = g[q]
                    if r not in orb:
                        orb.add(r); new.append(r)
            frontier = new
        # make sure it is the full orbit
        orb = {g[p] for g in elems}
        seen |= orb
        orbits.append(sorted(orb))
    for orb in orbits:
        p = orb[0]
        Gp = stabiliser(lambda g: g[p] == p)
        add(f"stab(pt)", Gp)
        # two-point stabilisers (pointwise and setwise), one per orbit of Gp
        sub_seen = {p}
        for q in orb:
            if q in sub_seen:
                continue
            sub_seen |= {g[q] for g in Gp}
            add("stab(2pts)", stabiliser(lambda g: g[p] == p and g[q] == q))
            add("stab(2set)", stabiliser(lambda g: {g[p], g[q]} == {p, q}))
        # k-subsets of the orbit (a few random ones) and equal-block partitions
        rnd = random.Random(0)
        for k in range(3, min(len(orb) // 2, 6) + 1):
            for _ in range(2 if len(orb) > 12 else 1):
                B = set(rnd.sample(orb, k))
                add(f"stab({k}set)", stabiliser(lambda g: {g[q] for q in B} == B))
            # consecutive k-set as well
            B = set(orb[:k])
            add(f"stab({k}set)", stabiliser(lambda g: {g[q] for q in B} == B))
        L = len(orb)
        for msize in range(2, L // 2 + 1):
            if L % msize:
                continue
            blocks = [frozenset(orb[i:i + msize]) for i in range(0, L, msize)]
            bl = frozenset(blocks)
            add(f"stab({L//msize}x{msize} partition)",
                stabiliser(lambda g: frozenset(frozenset(g[q] for q in B) for B in blocks) == bl))
    # (d) abelian subgroups of odd order (the old method's candidates)
    for H in C.odd_abelian_subgroups(elems, K.ords):
        add(f"abelian odd", H)
    # (e) conjugate subgroups give isomorphic quotients: keep at most two candidates per
    #     (order, multiset of element orders) signature
    by_sig = {}
    for fs, (name, H) in cands.items():
        sig = (len(H), tuple(sorted(K.ords[h] for h in H)))
        by_sig.setdefault(sig, []).append((name, H))
    out = []
    for sig, L in by_sig.items():
        out += L[:2]
    out.sort(key=lambda t: -len(t[1]))
    return out

# ------------------------------------------------------------------ pairs
def bfs_exceeds(gens, threshold):
    seen = {K.e}
    frontier = [K.e]
    while frontier:
        new = []
        for g in frontier:
            for s in gens:
                h = mul(g, s)
                if h not in seen:
                    seen.add(h)
                    if len(seen) > threshold:
                        return True
                    new.append(h)
        frontier = new
    return False

def generates(a, x):
    if K.T is None:
        return bfs_exceeds([a, x], K.N // K.min_index)
    if a in K.T:
        return False
    xa = mul(mul(a, x), a)
    return bfs_exceeds([x, xa], K.Tsize // K.min_index)

def enumerate_pairs():
    pairs = []
    for a, _ in K.inv_reps:
        if K.T is not None and a in K.T:
            continue
        cent = [c for c in K.elems if mul(c, a) == mul(a, c)]
        cent_inv = [inv(c) for c in cent]
        for cid, (y, cls) in enumerate(K.odd_classes):
            # X(G;a,x) = X(G;a,x^-1): if the class of x^-1 is a different class
            # with a smaller index, its pairs have already been listed
            if K.class_of[inv(y)] < cid:
                continue
            seen = set()
            for x in cls:
                if x in seen:
                    continue
                xi = inv(x)
                for c, ci in zip(cent, cent_inv):
                    seen.add(mul(mul(ci, x), c))
                    seen.add(mul(mul(ci, xi), c))
                if generates(a, x):
                    pairs.append((a, x, cid))
    return pairs

# ------------------------------------------------------------------ quotient, lift, verify
def quotient(H, a, x):
    elems, index, N = K.elems, K.index, K.N
    coset = [-1] * N
    reps = []
    for gi, g in enumerate(elems):
        if coset[gi] >= 0:
            continue
        cid = len(reps); reps.append(gi)
        for h in H:
            coset[index[mul(h, g)]] = cid
    n = len(reps)
    items, a_item, x_item = [], [None] * n, [None] * n
    for cid, gi in enumerate(reps):
        c2 = coset[index[mul(elems[gi], a)]]
        if a_item[cid] is None:
            if c2 == cid:
                items.append((cid, None)); a_item[cid] = len(items) - 1
            else:
                items.append((cid, c2)); a_item[cid] = a_item[c2] = len(items) - 1
    for cid, gi in enumerate(reps):
        c2 = coset[index[mul(elems[gi], x)]]
        if c2 == cid:
            return None                      # loop: H contains a conjugate of x
        items.append((cid, c2)); x_item[cid] = len(items) - 1
    return n, items, coset, a_item, x_item

def lift_and_verify(a, x, coset, a_item, x_item, col):
    """Independent check that the colouring induced on X(G;a,x) is proper:
    the a-edge {g,ga} gets col[a_item[coset(g)]], the x-edge {g,gx} gets
    col[x_item[coset(g)]]."""
    elems, index = K.elems, K.index
    xi = inv(x)
    for gi, g in enumerate(elems):
        c = coset[gi]
        ca = col[a_item[c]]
        if col[a_item[coset[index[mul(g, a)]]]] != ca:
            return False                     # a-edge colour must agree from both ends
        cx_out = col[x_item[c]]
        cx_in = col[x_item[coset[index[mul(g, xi)]]]]
        if len({ca, cx_out, cx_in}) != 3:
            return False
    return True

MAX_QUOTIENTS = int(os.environ.get("MAX_QUOTIENTS", 60))
FULL_TIME_LIMIT = float(os.environ.get("FULL_TIME_LIMIT", 900.0))

def quotient_time_limit(nv):
    return min(120.0, 5.0 + nv / 40.0)

def decide_pair(a, x, cid):
    tried = 0
    unsat = []
    for hi in K.adm[cid]:
        if tried >= MAX_QUOTIENTS:
            break
        hname, H = K.cands[hi]
        q = quotient(H, a, x)
        if q is None:
            continue
        tried += 1
        n, items, coset, a_item, x_item = q
        res = run_with_timeout(colour_pregraph, (n, items), quotient_time_limit(n))
        if res is None:
            continue
        ok, col = res
        if ok:
            assert lift_and_verify(a, x, coset, a_item, x_item, col), "lift failed"
            return True, f"quotient by {hname} of order {len(H)} ({n} vertices)", n
        unsat.append(n)
    # fallback: the full Cayley graph (quotient by the trivial subgroup)
    q = quotient([K.e], a, x)
    n, items, coset, a_item, x_item = q
    res = run_with_timeout(colour_pregraph, (n, items), FULL_TIME_LIMIT)
    if res is None:
        return None, f"timeout (uncolourable quotients on {unsat} vertices)", n
    ok, col = res
    if ok:
        assert lift_and_verify(a, x, coset, a_item, x_item, col)
        return True, f"full SAT ({n} vertices; uncolourable quotients on {unsat})", n
    return False, "full SAT: UNSAT -- CAYLEY SNARK", n

VERBOSE = os.environ.get("VERBOSE") == "1"

def _worker(args):
    a, x, cid = args
    t0 = time.time()
    try:
        out = decide_pair(a, x, cid)
    except Exception as ex:
        out = (None, f"error: {ex!r}", 0)
    if VERBOSE:
        print(f"    pair ord(x)={K.ords[x]} ord(ax)={K.ords[mul(a, x)]}: {out[0]} {out[1]} [{time.time()-t0:.1f}s]", flush=True)
    return out

def check_group(name, G, min_index, socle_index, pool=True):
    t0 = time.time()
    setup_group(name, G, min_index, socle_index)
    pairs = enumerate_pairs()
    print(f"  [{name}] generating pairs (a,x) up to conjugacy and inversion: {len(pairs)}  "
          f"({time.time()-t0:.1f}s)", flush=True)
    if pool and len(pairs) > 1:
        import multiprocessing as mp
        nworkers = int(os.environ.get("WORKERS", max(1, os.cpu_count() - 1)))
        with mp.Pool(nworkers) as P:
            outs = P.map(_worker, pairs, chunksize=1)
    else:
        outs = [_worker(p) for p in pairs]
    results, methods = [], {}
    for (a, x, cid), (ok, how, nv) in zip(pairs, outs):
        s, t = K.ords[x], K.ords[mul(a, x)]
        results.append(dict(ord_x=s, ord_ax=t, ok=ok, method=how, quotient_vertices=nv))
        key = how.split(" (")[0]
        methods[key] = methods.get(key, 0) + 1
        if ok is False:
            print(f"!!! NON-3-EDGE-COLOURABLE CAYLEY GRAPH on {name}: ord(x)={s}, ord(ax)={t}", flush=True)
        elif ok is None:
            print(f"??? undecided pair on {name}: ord(x)={s}, ord(ax)={t}: {how}", flush=True)
    n_bad = sum(1 for r in results if r["ok"] is False)
    n_und = sum(1 for r in results if r["ok"] is None)
    maxq = max((r["quotient_vertices"] for r in results if r["ok"]), default=0)
    print(f"{name:14s} |G|={K.N:7d}  inv.classes={len(K.inv_reps)}  pairs={len(pairs):5d}  "
          f"non-colourable={n_bad}  undecided={n_und}  largest quotient solved={maxq}  "
          f"[{time.time()-t0:.1f}s]", flush=True)
    return dict(group=name, order=K.N, degree=K.n, involution_classes=len(K.inv_reps),
                pairs_tested=len(pairs), non_colourable=n_bad, undecided=n_und,
                orders_x=sorted(set(r["ord_x"] for r in results)),
                pairs=sorted(set((r["ord_x"], r["ord_ax"]) for r in results)),
                methods=sorted(set(r["method"].split(" (")[0] for r in results)),
                max_quotient_vertices=maxq,
                subgroups_used=sorted(set(r["method"].split(" (")[0] for r in results)),
                seconds=round(time.time() - t0, 1),
                details=results)

# ------------------------------------------------------------------ catalogue
def A(n): return lambda: C.alternating(n)
def S(n): return lambda: C.symmetric(n)
def L2(q): return lambda: C.psl2(q)
def PGL2(q): return lambda: C.psl2(q, True)

GROUPS = [
    # name, constructor, min index of the socle, |G : socle|, order
    ("A5", A(5), 5, 1, 60),
    ("S5", S(5), 5, 2, 120),
    ("PSL(2,7)", L2(7), 7, 1, 168),
    ("PGL(2,7)", PGL2(7), 7, 2, 336),
    ("A6", A(6), 6, 1, 360),
    ("PSL(2,8)", L2(8), 9, 1, 504),
    ("PSL(2,11)", L2(11), 11, 1, 660),
    ("S6", S(6), 6, 2, 720),
    ("PGL(2,9)", PGL2(9), 6, 2, 720),
    ("M10", C.m10, 6, 2, 720),
    ("PSL(2,13)", L2(13), 14, 1, 1092),
    ("PGL(2,11)", PGL2(11), 11, 2, 1320),
    ("PGL(2,13)", PGL2(13), 14, 2, 2184),
    ("PSL(2,17)", L2(17), 18, 1, 2448),
    ("A7", A(7), 7, 1, 2520),
    ("PSL(2,19)", L2(19), 20, 1, 3420),
    ("PSL(2,16)", L2(16), 17, 1, 4080),
    ("PGL(2,17)", PGL2(17), 18, 2, 4896),
    ("S7", S(7), 7, 2, 5040),
    ("PSL(3,3)", C.psl3_3, 13, 1, 5616),
    ("PSU(3,3)", lambda: psu3(3, 6048), 28, 1, 6048),
    ("PSL(2,23)", L2(23), 24, 1, 6072),
    ("PGL(2,19)", PGL2(19), 20, 2, 6840),
    ("A5wrZ2", C.wreath_A5_Z2, 5, 2, 7200),
    ("PSL(2,25)", L2(25), 26, 1, 7800),
    ("M11", C.m11, 11, 1, 7920),
    ("PSigmaL(2,16)", lambda: C.psl2(16, frob=4, expected=8160), 17, 2, 8160),
    ("PSL(2,27)", L2(27), 28, 1, 9828),
    ("PSL(3,3).2", lambda: psl3(3, 11232, graph_aut=True), 13, 2, 11232),
    ("G2(2)", lambda: psu3(3, 12096, field_aut=3), 28, 2, 12096),
    ("PGL(2,23)", PGL2(23), 24, 2, 12144),
    ("PSL(2,29)", L2(29), 30, 1, 12180),
    ("PSL(2,31)", L2(31), 32, 1, 14880),
    ("PGL(2,25)", PGL2(25), 26, 2, 15600),
    ("PSigmaL(2,25)", lambda: C.psl2(25, frob=5, expected=15600), 26, 2, 15600),
    ("PSL(2,25).2_3", lambda: C.psl2(25, frob=5, frob_diag=True, expected=15600), 26, 2, 15600),
    ("PGL(2,27)", PGL2(27), 28, 2, 19656),
    ("A8", A(8), 8, 1, 20160),
    ("PSL(3,4)", lambda: psl3(4, 20160), 21, 1, 20160),
    ("PSU(4,2)", lambda: psp4_3(), 27, 1, 25920),
    ("Sz(8)", lambda: atlas_group("Sz8G1-p65B0", 29120), 65, 1, 29120),
    ("S8", S(8), 8, 2, 40320),
    ("PSL(3,4).2_gf", lambda: psl3(4, 40320, field_aut=2, graph_aut=True), 21, 2, 40320),  # graph-field automorphism (central in Out)
    ("PSL(3,4).2_f", lambda: psl3(4, 40320, field_aut=2), 21, 2, 40320),               # field automorphism
    ("PSL(3,4).2_g", lambda: psl3(4, 40320, graph_aut=True), 21, 2, 40320),              # graph automorphism
    ("PGL(2,29)", PGL2(29), 30, 2, 24360),
    ("PSL(2,37)", L2(37), 38, 1, 25308),
    ("PGL(2,31)", PGL2(31), 32, 2, 29760),
    ("PSL(2,32)", L2(32), 33, 1, 32736),
    ("PSL(2,41)", L2(41), 42, 1, 34440),
    ("PSL(2,43)", L2(43), 44, 1, 39732),
    ("PGL(2,37)", PGL2(37), 38, 2, 50616),
    ("PSL(2,47)", L2(47), 48, 1, 51888),
    ("PSL(2,49)", L2(49), 50, 1, 58800),
    ("PSU(4,2).2", lambda: psp4_3(True), 27, 2, 51840),
    ("PSL(2,7)wrZ2", lambda: wreath_Z2(C.psl2(7), 56448), 7, 2, 56448),
    ("PSU(3,4)", lambda: psu3(4, 62400), 65, 1, 62400),
    ("PGL(2,41)", PGL2(41), 42, 2, 68880),
    ("PSL(2,53)", L2(53), 54, 1, 74412),
    ("PGL(2,43)", PGL2(43), 44, 2, 79464),
    ("M12", lambda: atlas_group("M12G1-p12aB0", 95040), 12, 1, 95040),
    ("PGL(2,47)", PGL2(47), 48, 2, 103776),
    ("PSL(2,59)", L2(59), 60, 1, 102660),
    ("PSL(2,61)", L2(61), 62, 1, 113460),
    ("PGL(2,49)", PGL2(49), 50, 2, 117600),
    ("PSigmaL(2,49)", lambda: C.psl2(49, frob=7, expected=117600), 50, 2, 117600),
    ("PSL(2,49).2_3", lambda: C.psl2(49, frob=7, frob_diag=True, expected=117600), 50, 2, 117600),
    ("PSU(3,4).2", lambda: psu3(4, 124800, field_aut=4), 65, 2, 124800),
    ("PSU(3,5)", lambda: psu3(5, 126000), 50, 1, 126000),
    ("PGL(2,53)", PGL2(53), 54, 2, 148824),
    ("PSL(2,67)", L2(67), 68, 1, 150348),
    ("J1", lambda: atlas_group("J1G1-p266B0", 175560), 266, 1, 175560),
    ("PSL(2,71)", L2(71), 72, 1, 178920),
    ("A9", A(9), 9, 1, 181440),
    ("M12.2", lambda: atlas_group("M12d2G1-p24B0", 190080), 12, 2, 190080),
    ("PSL(2,73)", L2(73), 74, 1, 194472),
    ("PGL(2,59)", PGL2(59), 60, 2, 205320),
    ("PGL(2,61)", PGL2(61), 62, 2, 226920),
    ("PSL(2,79)", L2(79), 80, 1, 246480),
    ("PSU(3,5).2", lambda: psu3(5, 252000, field_aut=5), 50, 2, 252000),
    ("A6wrZ2", lambda: wreath_Z2(C.alternating(6), 259200), 6, 2, 259200),
    ("PSL(2,64)", L2(64), 65, 1, 262080),
    ("PSL(2,81)", L2(81), 82, 1, 265680),
    ("PSL(2,83)", L2(83), 84, 1, 285852),
    ("PGL(2,67)", PGL2(67), 68, 2, 300696),
    ("PSL(2,89)", L2(89), 90, 1, 352440),
    ("S9", S(9), 9, 2, 362880),
]

if __name__ == "__main__":
    max_order = int(sys.argv[1]) if len(sys.argv) > 1 else 10 ** 5
    only = set(sys.argv[2:])
    import multiprocessing as mp
    mp.set_start_method("fork")
    outfile = os.environ.get("RESULTS", "results2.json")
    summary = json.load(open(outfile)) if os.path.exists(outfile) else []
    done = {r["group"] for r in summary if r["undecided"] == 0 and r["non_colourable"] == 0}
    for name, ctor, min_index, socle_index, expected in GROUPS:
        if only and name not in only:
            continue
        if expected > max_order or (name in done and not only):
            continue
        t0 = time.time()
        G = ctor()
        assert len(G) == expected, (name, len(G), expected)
        print(f"== {name}: order {len(G)} built in {time.time()-t0:.1f}s", flush=True)
        row = check_group(name, G, min_index, socle_index)
        summary = [r for r in summary if r["group"] != name] + [row]
        summary.sort(key=lambda r: (r["order"], r["group"]))
        with open(outfile, "w") as f:
            json.dump(summary, f, indent=1)
        del G
    print("\nSUMMARY")
    for r in summary:
        print(f"{r['group']:14s} {r['order']:8d} classes={r['involution_classes']} pairs={r['pairs_tested']:5d} "
              f"bad={r['non_colourable']} undecided={r['undecided']} maxq={r['max_quotient_vertices']}")
