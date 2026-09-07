# Research review and revised attack plan

Review date: 5 September 2026. Baseline: commit `174fc2a`.
Implementation follow-up: 6 September 2026 (below).

Latest follow-up: a joint polynomial criterion combines the projective
label recurrence with the common-translation chord equations. It gives
exact realizability and simultaneous-unit-gain tests, and proves a
forbidden anchored three-chord pattern over every characteristic-two
field. A one-row polynomial certificate handles the fixed q=2048 diagram
without full rank. A negative control shows the recurrence alone cannot
replace the translation equations. The final section gives the proofs,
certificates and limitations. The best proved asymptotic full-cycle
coverage remains at least 13/16; the full conjecture remains open.

The quotient and matching-monodromy approach is still productive for family
theorems. The recent sequence of larger finite symmetric-group cases does not,
however, establish that a proof of the full Alspach–Zhang conjecture is close.
The missing ingredient remains a uniform existence theorem: a reason that a
suitable matching must exist for every relevant generating pair.

This review checks the principal reductions and barriers, confirms the stored
computational totals and the external result used for degree 23, and develops
and tests a new Sylow-subgroup target. It does not rerun the tens of millions
of cases in the earlier symmetric-group enumerations or establish priority
over the Hamiltonicity literature.

## What has actually been achieved

* The survey gives proofs for every symmetric-group generating pair with a
  transposition as the involution. This is an infinite-family result, not an
  extrapolation from small tests.
* For an odd full cycle and an involution with `r` transpositions, the proved
  ranges are `r <= 7`, `N >= 4r+1`, `(r,N)=(9,19)`, and `N=23`. The small
  boundary cases use exhaustive computations; the large sparse range has a
  direct construction.
* The degree-23 argument correctly reduces to a cycle permutation graph of
  order 46. Goedgebeur, Renders and Van Overberghe explicitly checked all
  non-Hamiltonian examples of that order for colourability, so there is no
  hidden girth or connectivity restriction preventing the application.
  See their [paper, pp. 13–14](https://arxiv.org/pdf/2411.12606v4).
* The stored catalogue reports 85 groups, 25,954 generating-pair representatives,
  and a smallest-counterexample lower bound of 352,440 vertices. These are
  verified finite ranges, not a density estimate over all possible groups.

The three-semiedge path question remains a conjecture. The exact experiment
through order 13 comprises 288,970 diagrams and 866,910 terminal-pair tests in
total; **810,810 is the number at order 13 alone**. The seeded order-21 run
tests 20,000 diagrams, or 60,000 terminal pairs. Exhausting order 21 without
symmetry reduction would mean 45,831,035,250 diagrams. This growth is a reason
to seek an induction, rather than treat another enumeration as a route to a
uniform theorem.

## Why the present results do not imply a nearby full proof

The minimal-counterexample reduction leaves all nonabelian simple groups,
appropriate extensions of a simple group by an involution, and the wreath
product cases. A full cycle in `S_N` is one special family inside the
index-two branch. Proving the three-fixed-point question would still leave
full cycles with one fixed point, generators with several cycles, and the
simple-group branch untouched.

In particular, a nonabelian simple group has no nontrivial homomorphism to
`C_2`: such a map would have a normal subgroup of index two. Thus the sign
certificate used in the recent symmetric-group work cannot be transferred
unchanged to that branch. Even the order-three generator case still needs
a general theorem about cubic graphs with an arc-regular simple-group action.

The CDC paper remains a source of methodology, but it supplies no implication
from a cycle double cover to a Tait colouring. Its affine compatibility system
constructs the cover; the extra palette-colourability requirement is a
separate global existence problem. The survey already proves that several
palette formulations are exact reformulations of Tait colourability and
records failures of dimension-only and local-filter shortcuts. See the
[CDC paper](https://cdn.openai.com/pdf/04d1d1e4-bc75-476a-97cf-49055cd98d31/cdc_proof.pdf)
and the survey's palette section.

Likewise, the explicit PSL(2,11) examples disprove universal single-circuit
repair and some minimal-block rules; other examples refute two uniform
negative-drift arguments. These results do not disprove the conjecture or
every possible exchange method. They do mean those particular proof plans
must not be revived merely because they worked on smaller groups.

There is no defensible percentage complete or completion date. A uniform
theorem for one infinite family of simple groups would be a substantive next
milestone. No current result supplies all the further steps to the full proof.

## Approach 1: use a 2-subgroup to replace the missing sign

Let `X = Cay(G,{a,x,x^-1})`, with `x` of odd order, and let `P` be a
2-subgroup of `G`, preferably a Sylow 2-subgroup containing `a`. Put
`Y = P\X`. For a perfect pregraph matching `Mbar` of `Y`, let `B(Mbar)`
be the number of odd complementary circuits whose traversal word is the
identity in `G`.

**Proved reduction.** The matching lifted to `X` has exactly

\[
                 |P|\,B(\overline M)
\]

odd complementary circuits. Consequently it gives a Tait colouring exactly
when `B(Mbar)=0`.

Proof: `P` contains no conjugate of the nonidentity odd-order element `x`, so
the quotient has no loops. A circuit based at `Pg` has word `w` satisfying
`gwg^-1 in P`, hence `ord(w)` is a power of two. A circuit of length `ell`
lifts to `|P|/ord(w)` circuits of length `ell*ord(w)`. An odd circuit remains
odd precisely when `w=1`, in which case it has `|P|` lifts. Capped paths
always lift to even circuits. This is a specialization of the survey's
matching-monodromy lemma, not a new existence theorem.

The new target is:

> For every generating pair in `PSL(2,q)`, does some Sylow 2-subgroup `P`
> containing `a` admit a quotient matching with no odd identity-word circuit?

This is stronger than colourability alone, because it asks for an invariant
matching. It is meaningful on simple groups and has a finite obstruction:
for each odd identity-word circuit `C`, require that the matching select at
least one edge of `C`. Together with the matching constraints, these are
necessary and sufficient for this particular quotient to work.

The proof task is to show these constraints always have a common integral
solution in a specified family. Solving them on additional examples does not
itself provide that proof. A failed quotient would require trying other Sylow
subgroups containing `a`, or smaller 2-subgroups, before rejecting the broader
approach. A bad word already equal to the identity cannot be repaired by
changing the subgroup while keeping exactly that same complementary circuit.

For odd `q`, matrix coordinates offer a precise arithmetic version. Choose an
`SL(2,q)` lift `W` of a circuit word. Because its projective order is a power
of two,

\[
  w=1\quad\Longleftrightarrow\quad(\operatorname{tr}W)^2=4.
\]

Indeed, trace `+2` or `-2` gives a scalar matrix or a scalar times a
nontrivial unipotent; the latter has odd projective order equal to the field
characteristic and is excluded. The converse is immediate. Thus a uniform
construction could certify the required words by finite-field identities,
without a sign character of the ambient simple group. The inequality still
has to be proved for the constructed matching.

### Diagnostic evidence collected in this review

Run `python3 code/sylow_monodromy_review.py` from the repository root.

| q | Group order | Sylow order | Quotient order | Pair orbits passing | Quotients not colourable |
|---|---:|---:|---:|---:|---:|
| 5 | 60 | 4 | 15 | 6 | 0 |
| 7 | 168 | 8 | 21 | 4 | 0 |
| 8 | 504 | 8 | 63 | 42 | 0 |
| 9 | 360 | 8 | 45 | 4 | 0 |
| 11 | 660 | 4 | 165 | 14 | 0 |
| 13 | 1092 | 4 | 273 | 45 | 0 |
| 17 | 2448 | 16 | 153 | 60 | 2 |
| 19 | 3420 | 4 | 855 | 76 | 0 |
| 23 | 6072 | 8 | 759 | 107 | 0 |

All 358 orbits pass the monodromy test. Each orbit is taken under inner
conjugation and `x -> x^-1`; the code verifies coverage of the involution
class and checks generation. A chosen Sylow subgroup works at each
representative, which proves existence of some such subgroup throughout
that orbit, not that every Sylow subgroup works. Every positive matching
is lifted and its complementary cycles traversed on the full Cayley graph.

The two non-colourable quotients at `q=17` both have `ord(x)=3`. A separate
colouring encoding also reports them unsatisfiable. Their monodromy
certificates have respectively two and one odd quotient circuits with
nonidentity words, and give even factors upstairs with nine and ten
components. Thus these examples distinguish matching monodromy from
colouring the quotient itself.

These groups were already covered by the catalogue. This experiment tests a
more specific proposed mechanism; it does not improve the counterexample
size bound or prove the infinite-family target.

## Approach 2: a counting theorem in a structurally chosen matching class

For a nonempty matching class `F` in the above quotient, choose a matching
uniformly from `F`. Let `Z=|F|`, and for every odd identity-word circuit `C`
let `Z_C` count the members whose complements contain `C`. Then

\[
 \mathbb E B = \frac{1}{Z}\sum_C Z_C.
\]

Therefore the strictly stronger sufficient bound `sum_C Z_C < Z` proves
existence of a good matching. This asks for an unsigned counting estimate,
and does not require a decreasing move from every bad matching.

For the class of all perfect pregraph matchings, each `Z_C` has a concrete
description. Force the unique item outside `C` at every vertex of `C`.
If these forced items are incompatible, `Z_C=0`. Otherwise delete all
vertices they cover, and count perfect matchings of the remaining pregraph.
The same description works with prescribed costs by retaining the residual
cost requirement. These are ordinary matching counts; the difficult part is
a useful uniform estimate on their sum.

**A necessary correction to the naive version:** uniform sampling from all
quotient matchings does not always give expectation less than one. Exact
enumeration on an `A_5` order-three example gives 16 matchings, with bad-cycle
histogram `{0:6, 1:3, 2:6, 5:1}`, so `E B=5/4`. On the PSL(2,7)
order-three example the expectation is `7/4`. Both graphs are colourable.

A distribution specified without knowing a colouring can instead restrict
to matchings minimizing the number of vertices matched by involution
items. In the ten representative quotients for `q=5,7`, this gives expectation
zero in nine cases and `2/17` in the other. The latter example also disproves
the claim that every minimum-cost Sylow-invariant matching works. All 522
unrestricted matchings in these ten quotients were checked against direct
full lifts, including the exact formula `oddness = |P| B` above.

The worthwhile target is a bound for this specified minimum-cost measure
on a whole family, possibly separating order-three generators from larger
orders. No such estimate has been proved. Selecting a measure concentrated
on an already known good matching would be circular and is not the proposal.

## Approach 3: a parity-preserving induction for three semiedges

Keep the three-fixed-point question as a secondary, more focused objective.
A proof must preserve both path connectivity and the two parity quantities
in the statement; existence of an arbitrary alternating path is insufficient.

There is an elementary reduction that does preserve the stronger,
prescribed-terminal-pair certificate. Suppose a chord joins two consecutive
nonterminal vertices `u,v` of the Hamilton cycle, with cycle segment
`p-u-v-q`. Delete `u,v` and that chord, replacing the segment by `p-q`.
Any certified alternating path in the reduced diagram extends as follows.

* If it uses `p-q`, replace that path edge by cycle edge `p-u`, chord `u-v`,
  and cycle edge `v-q`. The complement gains the even two-cycle consisting
  of the parallel chord and cycle edges between `u,v`.
* If it does not use `p-q`, retain the same path and keep chord `u-v` in the
  matching. The complementary edge `p-q` is replaced by three cycle edges,
  changing its component length by two and its chord count by zero.

In either case the prescribed terminals are unchanged, the alternating path
stays simple, and every old circuit's two parities are preserved. The base
diagram on three semiedge vertices is immediate. Thus a minimal
counterexample to the stronger path question has no chord parallel to a
cycle edge.

This is only one reduction. For a general chord insertion, its two endpoints
can affect different complementary components. Suppressing them can change
the parity of each component separately. The next useful experiment is to
retain boundary states recording matching occupancy, endpoint pairings,
length parity and chord parity, and test which reductions preserve an
adequate set of states. A successful induction needs an unavoidable
reducible configuration in every remaining diagram. A finite table of
states without that structural theorem would not complete the proof.

## Recommended priorities and evidence needed to change them

1. Make the Sylow-monodromy target for one infinite simple-group family the
   main program. Search for obstructions to this stronger target as well as
   constructions. The desired output is a uniform matching rule plus a word
   argument, not another larger positive catalogue.
2. Use the minimum-cost counting estimate as an alternative way to prove
   existence within that program. Test its proposed inequality directly and
   retain counterexamples to overly strong variants.
3. Pursue the three-semiedge induction as a contained family theorem. A proof
   of the exact three-fixed-point version would settle `(r,N)=(9,21)`, but
   would not by itself prove the version with arbitrarily many fixed points.
4. Return to the CDC palette systems only if there is a new theorem excluding
   their global inconsistency certificates using the group action. More
   affine solutions or larger nullspaces alone do not meet that milestone.

The evidence supports changing the main emphasis, while retaining the
quotient and monodromy machinery already developed. It does not support
announcing that the full conjecture is one small step away.

## Follow-up: implementing the recommendation, 6 September

The main effort now follows Approaches 1 and 2. The survey includes the
2-subgroup odd-circuit formula as a proved corollary, and a new uniform
symmetry lemma for the characteristic-two simple-group family. Neither is
an existence proof for every generating pair.

### Exact minimum-cost counting

The new `code/sylow_minimum_face_experiment.py` contracts the quotient
`x`-cycles and enumerates matchings selecting exactly one `a`-dart per cycle.
An ordinary contracted loop is excluded, since it would select two darts
on the same cycle. A semi-edge selects one dart. Once these items are
selected, each remaining even cycle-path has a unique perfect matching.
Thus the algorithm counts a class specified without using any colouring.
When nonempty it is the minimum-cost class described above; its universal
nonemptiness is not assumed.

Completion counts give both an exact enumerator and an exact uniform
sampler. The following results are **exhaustive**, not sampled:

| q | Generator orders | Pair orbits | Matchings checked | Largest exact E B |
|---|---|---:|---:|---:|
| 5 | all odd | 6 | 40 | 0 |
| 7 | all odd | 4 | 63 | 2/17 |
| 8 | all odd | 42 | 22,200 | 3/47 |
| 9 | all odd | 4 | 388 | 8/97 |
| 11 | 3 | 1 | 53,160 | 491/17720 |
| 17 | 3 | 8 | 211,012 | 254/13539 |

All 65 representatives have `E B < 1`, with 286,863 matchings checked in
total. The `q=11` histogram is `{0:51690, 1:1467, 2:3}`. These numbers
support the proposed first-moment inequality, not the false claim that
every matching in this class works. Full Cayley traversals independently
verify the odd-circuit formula for every small matching (`q=5,7,8,9`) and
for positive and negative witnesses in each larger case.

All six order-three representatives at `q=13` hit the preset 500,000-state
counting limit. They are unresolved counting tests, not obstructions and
not part of the exhaustive totals. This limitation matters: the state
space grows with the quotient, not just with the group order.

The exact data are in `code/sylow_minimum_face_small.json` and
`code/sylow_minimum_face_order3.json`. Reproduce from the repository root:

```sh
python3 code/sylow_minimum_face_experiment.py 5 7 8 9 --json-output code/sylow_minimum_face_small.json
python3 code/sylow_minimum_face_experiment.py 11 17 --x-order 3 --json-output code/sylow_minimum_face_order3.json
```

### A uniform Borel-reversal lemma, and a genuine obstruction

Let `q=2^m`, and put `a(z)=z+c`, `x(z)=1/(z+t)`, with `c,t` nonzero and
`x` a full `(q+1)`-cycle on the projective line. Every involution/full-order
pair in `PSL(2,q)` can be put in these coordinates. In the quotient by the
point stabiliser at infinity, select its unique semi-edge and alternate
cycle edges on the other `q` points.

**Proved:** this matching has precisely one odd complementary circuit,
through infinity, and its word is a translation `z -> z+beta`.

The key is `r(z)=z+t`: it commutes with `a`, reverses `x`, and preserves
the matching. It reverses the distinguished circuit, so `r w r = w^-1`.
Writing its word in the affine stabiliser as `w(z)=lambda*z+beta` forces
`lambda=lambda^-1`, hence `lambda=1` in characteristic two. All other
complementary circuits alternate cycle and chord edges and are even.
The matching therefore colours upstairs exactly when `beta != 0`.

This is a uniform word argument on an infinite simple-group family, but
the canonical nonvanishing claim is **false**. Exact tests give:

| q | Canonical successes | Identity-word failures | Repaired failures |
|---|---:|---:|---:|
| 4 | 6 | 0 | 0 |
| 8 | 21 | 0 | 0 |
| 16 | 120 | 0 | 0 |
| 32 | 310 | 0 | 0 |
| 64 | 1500 | 12 | 12 |

These count every allowed parameter pair, including nongenerating pairs,
not conjugacy orbits. At all twelve failures the odd circuit has length
49 and its matrix word is exactly the identity. For `GF(64)` represented
modulo `theta^6+theta+1`, the values `t=1+theta+theta^2` and
`c=theta+theta^3+theta^4` give one failure. A separate closure computation
verifies that this pair generates all 262,080 elements of `PSL(2,64)`.
Thus the issue cannot be dismissed as a nongenerating exception.

There is a small structured repair family: select a chord and its reflected
chord under cycle reversal, then complete with cycle edges if all residual
paths have even order. There are at most `q/4` such reflected pairs. Each
of the twelve failures has six admissible candidates, with five or six
successful. For the displayed counterexample, chords `(17,32)` and
`(33,48)` leave circuits of lengths 57 and 8, the odd one having translation
parameter `1+theta+theta^2+theta^3+theta^4+theta^5`, which is nonzero.
All circuit matrices are checked separately by composing their induced
projective permutations.

The full certificates, including failures, are retained in
`code/borel_reversal_results.json`. Reproduce with:

```sh
python3 code/borel_reversal_experiment.py --verify-generation --json-output code/borel_reversal_results.json
```

### Sharpened next proof targets

1. Prove, or refute beyond the tested range, that the canonical matching
   together with the reflected-two-chord family always contains a good
   matching for generating pairs with `q` even and `ord(x)=q+1`. This is a
   linear-size candidate family, not unrestricted matching search. A
   successful proof must control **every** odd circuit after a switch;
   symmetry of the circuit through infinity alone is insufficient.
2. Prove a bound `sum_C Z_C < Z` for the specified Sylow minimum-cost class
   on an infinite family, including its nonemptiness. The exact results
   test this actual inequality, but give no uniform bound as `q` grows.

The two-chord question covers only full-order generators in characteristic
two. It does not address all odd orders, odd characteristic, or the other
simple groups. The original conjecture remains open, with no justified
estimate of how far these targets are from a complete proof.

## Continuation: folding, scalar gains, and a proved subfamily

The next stage produces an actual infinite parameter-family theorem, as
well as a sharper version of the remaining existence problem.

### Proved: two parameter families on simple groups

For every `q=2^m >= 4`, choose `t` so that `x(z)=1/(z+t)` has order `q+1`.
If `a(z)=z+1` or `a(z)=z+t+1`, then `a,x` generate `PSL(2,q)` and their
cubic Cayley graph is 3-edge-colourable.

The canonical complementary circuit through infinity is exactly
`infinity, 0, c, c+t, t, infinity`. Its word is the nonidentity translation
`z -> z+t+1`; all other circuits alternate chord and cycle edges. This
gives an explicit proof, not an extrapolation from tests. Such `t` exist
for every `q`, by taking the trace of an element of order `q+1` in the
norm-one subgroup of `GF(q^2)^*`.

Generation follows from the standard maximal-subgroup classification:
a proper subgroup containing the full-cycle element must be its dihedral
normaliser, whose sole involution fixing infinity is `z -> z+t`.
Thus every normalized pair with `c != t` generates the simple group.
See [Devillers–Giudici–Li–Praeger, Theorem 7.17](https://api.research-repository.uwa.edu.au/ws/portalfiles/portal/1479773/11274_PID11274.pdf).
The case `c=t+1` has `ord(ax)=3` and overlaps the known Hamiltonicity
theorem; we do not claim priority for either colouring family. This result
settles specified parameters for every field, not all generating pairs on
these simple groups.

### Proved: the obstruction is a scalar path product

Assume `c != t` and fold the projective cycle by its reversal `i -> -i`.
The matching problem becomes a path on `h=q/2` vertices plus a fixed
matching of chords. The semi-edge at infinity is forced, and parity
prevents selection of the middle cycle edge. Thus perfect matchings of
this smaller graph correspond exactly to reversal-invariant matchings
of the point quotient.

Write `z_i=x^i(infinity)` and label the folded vertices by

\[
 k_i=z_i^2+t z_i+1,\qquad \Delta=c(c+t).
\]

These are distinct nonzero labels in an affine trace hyperplane. Chords
pair `k` with `k+Delta`; the path starts at label 1 and ends at label `t`.
For a matching `N`, its complement has a distinguished path `P_N` from
the infinity end to the middle-edge end. Define

\[
 \gamma(N)=\prod_{i\to j\text{ a chord traversed on }P_N} k_i/k_j.
\]

**Exact formula:** the odd circuit through infinity has word
`z -> z+t*(1+gamma(N))`. Equivalently, if the half-path matrix has
lower-left entry `d`, its translation parameter is `t+d^-2`.

If every endpoint of a selected chord is on `P_N`, every other quotient
circuit is even. Under this additional condition, `gamma(N) != 1` is
therefore necessary and sufficient for this matching to colour upstairs.
This supplies a sufficient construction without separately controlling
new odd circuits elsewhere.

The original two-reflected-chord construction selects **one folded chord**
`(i,j)` with `i<j`, `i` odd and `j` even, then completes the remaining path
segments uniquely. Its unproved existence target is now particularly
concrete: find such a chord with both endpoints on `P_N` and `gamma(N)!=1`
whenever the canonical gain is 1.

For the source labels of traversed chords, put
`f_N(T)=product(T+k)`. Then `gamma(N)=f_N(0)/f_N(Delta)`. One possible
algebraic attack is to force `f_N(0)!=f_N(Delta)` for at least one suitable
path using the recurrence `D_(i+1)=t D_i+D_(i-1)`, where `k_i=D_i^-2`.
The recurrence retains the projective ordering; arbitrary reordered labels
do not suffice.

### Exhaustive parameter tests through q=1024

The larger-field checker uses one trace representative per Frobenius orbit
and tests every nonzero `c`. Field automorphisms preserve the indexed
diagram and all circuit-word orders. The weighted counts therefore cover
all normalized parameters, including the dihedral `c=t` cases.

| q | Representatives tested | Parameters covered | Canonical failures | Minimum successful root-contained repairs |
|---|---:|---:|---:|---:|
| 64 | 252 | 1,512 | 12 | 4 |
| 128 | 762 | 5,334 | 56 | 3 |
| 256 | 4,080 | 32,640 | 96 | 7 |
| 512 | 9,198 | 82,782 | 198 | 13 |
| 1024 | 40,920 | 409,200 | 580 | 24 |

All 942 canonical failures among 531,468 covered parameters have a
successful root-contained repair. These are 55,212 direct representative
tests. They cover only full-order generators, so they do not improve the
smallest-counterexample size bound or settle the larger groups outright.

The implementation independently checks field irreducibility and the
multiplication table. Failed canonical words and a positive repair per
failure are also checked as permutations of the full projective line.
The separate folded-graph implementation exhausts all 1,130 matchings at
the representative parameters for `q=4,8,16,32`, checks the folding
bijection and scalar formula against matrix words, and tests the pentagon
family on 180 trace/parameter representatives through `q=1024`.

### Two barriers that must not be ignored

1. **Root containment alone is insufficient, even in the genuine family.**
   In `GF(1024)` modulo `theta^10+theta^3+1`, binary coordinates `t=62`,
   `c=836` and folded chord `(33,236)` keep both endpoints on the path but
   leave gain 1 and an identity-word odd circuit of length 873. All other
   circuits are even. Exact matrix Schreier steps verify generation of
   the full group of order 1,073,740,800, with the method independently
   checked against direct group closures in eight smaller controls.
   Four root-contained but identity-word candidates occur in the actual
   `q=1024` audit; each affected pair also has successful repairs.
2. **A generic path-order argument cannot prove the stronger repair target.**
   In `GF(16)`, with `t=2`, `c=6`, `Delta=11`, reorder the labels as
   `(1,14,6,10,13,9,5,2)`. The canonical gain is 1, but the only one-chord
   candidate `(1,4)` misses endpoint 4 on its distinguished path. These are
   the correct affine-hyperplane labels but in the wrong projective order.
   This is a control against an overgeneralised lemma, not a Cayley snark.

Reproduce the new checks from the repository root:

```sh
python3 code/borel_reflected_family.py 64 128 256 512 --json-output code/borel_reflected_results.json
python3 code/borel_reflected_family.py 1024 --json-output code/borel_reflected_1024.json
python3 code/borel_folded_experiment.py --reflected-results code/borel_reflected_results.json code/borel_reflected_1024.json --json-output code/borel_folded_results.json
```

The next proof problem is **existence of a root-contained, non-unit-gain
path for the projective ordering**. The formulas and the thin parameter
families are proved; that existence statement and the full conjecture
remain open.

## Further follow-up: a proved positive-density family (6 September 2026)

Changing the starting matching gives a uniform theorem, rather than just
more evidence for the canonical-repair conjecture. In the same full-cycle
coordinates, with `c` neither 0 nor `t`, the Cayley graph colours whenever

\[
 \operatorname{Tr}(1/c)=0\quad\text{or}\quad
 \operatorname{Tr}(1/(c+t))=0.
\]

For every fixed full-cycle `t`, this covers at least `q/2` of the `q-2`
generating translations. Its proportion tends uniformly to **three quarters**
as `q` grows. This is a proportion inside one specified parameter family,
not an estimate of how much of the full Alspach–Zhang conjecture is proved.
No literature-priority claim is made.

### Why the new construction works

Start with the folded matching of **all a-chords**, whose complement is the
natural spanning path but has gain 1. If a chord is parallel to a path edge,
replace that selected chord by the parallel path edge. The new complement
still spans the graph, using the chord in place of the path edge. Its gain
is `k/(k+Delta)`, which cannot be 1. Thus the point quotient has one odd
Hamiltonian circuit with a nonidentity translation word. Every circuit
upstairs has length `2(q+1)`, so the lifted factor is even.

A parallel folded chord exists exactly when either `ax` or `arx` fixes a
projective point, where `r(z)=z+t`. Their matrix traces are `c+t` and `c`.
A determinant-one matrix of nonzero trace `v` splits over the field exactly
when `Tr(v^-2)=Tr(v^-1)=0`. This proves the displayed trace test.

There is a second proved short exchange. If the folded chords contain
`(i,j)` and `(i+1,j+1)`, with `i+1<j<h`, replace those two selected chords
by path edges `(i,i+1)` and `(j,j+1)`. The complementary spanning path is

`1,...,i, j,j-1,...,i+1, j+1,...,h`.

It traverses exactly two chords. If their source labels are `k,l`, unit
gain would imply `k*l=(k+Delta)*(l+Delta)`, hence `k+l=Delta`. That would
make them the two ends of the same chord, a contradiction. This gives
another uniform sufficient configuration, but its existence is not uniform.

### Exact count and the role of the external bound

Define the **unshifted** binary Kloosterman sum

\[
 K(b)=\sum_{u\ne0}(-1)^{\operatorname{Tr}(u+b/u)}.
\]

The number of translations passing the trace test is exactly

\[
 N_t=\frac{3q-5-K(t^{-2})}{4}.
\]

This follows by expanding the two trace indicators. The individual
character sums over `c != 0,t` vanish because `Tr(1/t)=1`; the mixed sum
is `K(t^-2)-1` after substituting `s=c/(c+t)`. The standard bound
`|K(b)| <= 2*sqrt(q)` then gives

\[
 \left|N_t-\frac{3q-5}{4}\right|\le\frac{\sqrt q}{2}.
\]

The external ingredient is the elliptic-curve interpretation and Hasse
bound, as recorded in [Ahmadi–Granger, Lemma 2.1 and Section 6.3](https://arxiv.org/pdf/1104.3882).
Their Kloosterman sum includes an additional 1; ours does not. The matching
construction and the displayed parameter count are derived here, not
claimed to be results from that paper.

The elementary `N_t >= q/2` bound needs no character-sum estimate: the set
`A={c!=0: Tr(1/c)=0}` has odd size `q/2-1`, excludes `t`, and cannot be
invariant under fixed-point-free translation by `t`. Hence
`|A union (A+t)| >= |A|+1 = q/2`.

### Verification and the remaining gap

The new independent audit covers 55,224 normalized generating parameter
representatives, weighted to 531,146 parameters, for every `q=4,8,...,1024`.
These totals differ from the older run because this one includes the four
smallest fields and excludes every nongenerating `c=t` parameter.
It verifies the exact count at every full-cycle trace representative and
checks all 69,022 short-exchange candidates by their explicit spanning
paths and scalar gains. On those candidates there are 702 matrix/half-word
audits and 340 full projective-permutation audits; larger-field full-word
checks are deliberately sampled, including every applicable stored canonical
failure. All ten small-field candidates at `q=4,8` are additionally lifted
to the full Cayley graphs (60 and 504 vertices): the actual complementary
circuits and all three edge colours are checked directly. The JSON records
these separate verification levels, as well as the extra residual-case
word checks.

The old canonical failures split as follows; counts are Frobenius-weighted
normalized parameters, not graph isomorphism classes:

| q | Canonical failures | Parallel exchange | Additional square exchange | Neither short exchange |
|---|---:|---:|---:|---:|
| 64 | 12 | 12 | 0 | 0 |
| 128 | 56 | 42 | 0 | 14 |
| 256 | 96 | 48 | 0 | 48 |
| 512 | 198 | 126 | 36 | 36 |
| 1024 | 580 | 540 | 20 | 20 |
| Total | 942 | 768 | 56 | 118 |

For a concrete residual case, in `GF(128)` modulo `theta^7+theta+1`, take
`t=15`, `c=35`. Both inverse traces are 1, the canonical word is identity
on a quotient circuit of length 129, and neither new exchange exists.
The older root-contained one-chord repair `(1,56)` nevertheless has gain
25, translation parameter 11, and quotient circuit lengths `97,16,8,8`.
The full projective-permutation audit verifies that it colours upstairs.
Thus the new construction does not itself settle all canonical failures.

The reordered-label barrier is also stronger than first stated. The
`GF(16)` diagram with labels `(1,14,6,10,13,9,5,2)` and `Delta=11` has
exactly three perfect matchings. The canonical and all-chord matchings
are root-contained with gain 1; the only other matching has gain 2 but
is not root-contained. Consequently allowing **all** folded matchings,
rather than just one-chord repairs, still does not prove the sufficient
target for arbitrary label orders.

Reproduce the new audit:

```sh
python3 code/borel_spanning_exchange.py --json-output code/borel_spanning_results.json
```

The sharpened next task is to exploit the projective ordering for the
residual class `Tr(1/c)=Tr(1/(c+t))=1`, after applying the proved pentagon
and short-exchange constructions. A proof for this class would settle the
full-projective-cycle family in characteristic two; other odd generator
orders and other simple groups would still remain. The new result is a
genuine uniform family theorem, not a full proof or a basis for a completion
date.

## Further follow-up: short-cycle elimination and the existence barrier

The next step gives a proved algebraic restriction on bad three-chord
exchanges, together with verified larger exchanges for every stored
residual case. It does not settle the existence problem.

### A constant-size algebraic exception set

Start with the folded matching of all chords and exchange along one simple
alternating circuit containing three chords and three path edges. Suppose
its complementary path spans the folded graph. Put `Delta=c(c+t)` and

\[
 E_t(\Delta)=\Delta^4+t^3\Delta^2+t^4\Delta+1+t^3+t^4+t^5.
\]

**Proved:** if `c` is not one of `0,t,1,t+1` and `E_t(Delta) != 0`, the
spanning exchange has non-unit gain and colours the Cayley graph. For each
fixed `t`, at most four values of `Delta`, hence at most eight translations
`c`, can violate this nonvanishing conclusion. The excluded `c=1,t+1`
graphs are already coloured by the pentagon construction.

This is conditional on a spanning exchange existing. It must not be read
as saying that all but eight remaining translations have now been coloured.

The proof replaces a long complementary-circuit word with a six-edge
alternating-circuit word. With `R(z)=z+t`, the corrected short word is

\[
 W=R^\eta X^{\epsilon_3}A X^{\epsilon_2}A X^{\epsilon_1}A,
 \quad \eta\in\{0,1\},\quad \epsilon_i\in\{1,-1\},
\]

and fixes its starting point `z`. If its denominator there is `d`, its
oriented short-cycle gain is the derivative `d^-2`. When the three chord
orientations agree with the spanning path, the trace of this short word
already proves nonvanishing outside the pentagon families.

Otherwise orient and start the short circuit so that exactly its first
chord disagrees with the path. Writing `Q(z)=z^2+t*z+1`, unit path gain
forces `d*Q(z)=Q(z)+Delta`. Eliminate `z` between this equation and
`W(z)=z`. Exact multiplication of all sixteen short words leaves, after
excluding factors already known nonzero, only

\[
 H_0=(\Delta+1)^2+t^2(c+t+1),\qquad
 H_1=(\Delta+1)^2+t^2(c+1).
\]

Their product is exactly `E_t(Delta)`. The survey gives the elimination
identity and complete three-row factor table. The code checks these as
identities in `F_2[c,t]`, not just as evaluations in small fields.

### The exceptions and the existence obstruction are real

In `GF(16)` modulo `theta^4+theta+1`, take `t=8,c=4,Delta=5`. An exchange
selecting path edges `(1,2),(3,4),(6,7)` spans all eight folded vertices but
has gain 1; its point-quotient circuit has length 17 and identity word.
Here `E_t(Delta)=0`. A different three-chord exchange, selecting path edges
`(2,3),(4,5),(7,8)`, succeeds. Thus automatic nonidentity for arbitrary
spanning three-chord exchanges is false, even in the double-trace-one class.

There is a second non-pentagon example at `q=1024,t=95,c=533,Delta=910`,
with selected path-edge starting indices `101,321,402`. Its spanning
quotient circuit has length 1025 and identity word, again with `E=0`.
Full projective-permutation checks verify both examples. An additional
four-chord unit-gain control occurs at `q=512,t=53,c=66`, with path-edge
starting indices `27,58,179,210`.

All 118 stored canonical failures left by the parallel and crossing-square
constructions have a successful spanning exchange along a single alternating
circuit with at most five chords:

| q | Residual parameters | Minimum three chords | Minimum five chords |
|---|---:|---:|---:|
| 128 | 14 | 0 | 14 |
| 256 | 48 | 16 | 32 |
| 512 | 36 | 36 | 0 |
| 1024 | 20 | 20 | 0 |
| Total | 118 | 72 | 46 |

These are Frobenius-weighted normalized parameters, represented by 14
direct cases. Every positive witness has a spanning complement and is
checked by full projective permutations. The cycle search exhausts all
single alternating-circuit exchanges up to the stated size; it does not
claim to enumerate disconnected combinations of exchanges. For example,
`q=128,t=15,c=35` has no successful such exchange with at most four chords,
so nonidentity for three-chord exchanges cannot by itself finish the proof.

The general three-chord audit uses Frobenius representatives for `t` and
one of the folded pair `c,c+t`. It covers all generating parameters for
`q<=64`, and the residual double-inverse-trace-one class for `q=128,...,1024`.
There are 7,032 folded parameter representatives, 7,789 three-chord
alternating circuits, and 3,941 spanning exchanges in that scope. All
spanning cases satisfy the proved quartic criterion; 139 are additionally
checked by independent full projective-word composition. The 92 unit-gain
spanning exchanges all lie in the excluded pentagon families or on the
quartic, as required.
The exact scope, counts, symbolic identities, negative controls and positive
certificates are retained in `code/borel_exchange_results.json`. The cycle
enumerator is also checked independently against exhaustive perfect
matchings for `q=4,8,16`.

```sh
python3 code/borel_exchange_algebra.py
python3 code/borel_alternating_exchange.py --json-output code/borel_exchange_results.json
```

The next proof task is now more sharply separated: establish existence of
a suitable spanning exchange, then control its gain using the short-word
arithmetic. The three-chord quartic theorem makes the latter step uniform
for that size. The observed five-chord bound was only finite evidence;
the next follow-up disproves its uniform spanning version. Allowing the
exchange length to grow, or dropping spanning, may be necessary. The
full-cycle family and the full Alspach–Zhang conjecture remain open.

## Further follow-up: a certified spanning-bound obstruction and two alternatives

### The five-chord spanning bound is false

In `GF(4096)` modulo `theta^12+theta^3+1`, take binary coordinates
`t=681,c=1207,Delta=2854`. The projective map `x` has full order 4097;
`c!=t` implies generation of all `PSL(2,4096)` by the proved generation
lemma. Both inverse traces are 1. The canonical quotient circuits have
lengths `1929,1012,1140,8,8`, with identity word on the odd one.

Starting from the folded matching of all chords, the following table is
the complete list of single simple alternating cycles with at most six
chords. Cuts are one-based starting indices of the selected path edges.
The folded graph has 2048 vertices.

| Chords | Cuts | Vertices on distinguished path | Colours upstairs? |
|---:|---|---:|:---:|
| 3 | 1261,1396,1700 | 1609 | No |
| 4 | 38,252,1168,1382 | 1834 | Yes |
| 4 | 48,103,648,1311 | 1385 | No |
| 4 | 585,1594,1630,1866 | 1003 | Yes |
| 5 | 314,1388,1570,1695,1780 | 764 | No |
| 5 | 441,472,1190,1807,2010 | 1330 | Yes |
| 6 | 137,568,1555,1859,1882,1895 | 1604 | No |
| 6 | 358,725,1060,1073,1624,2041 | 2048 | Yes |

There are no one- or two-chord cycles. In particular, **no single exchange
with at most five chords spans**. The six-chord spanning witness has gain
2930 and one point-quotient circuit of length 4097 with word `U(2333)`;
the lifted circuits have length 8194.

The negative conclusion has two independent enumerations. The original
graph search is checked against all fixed points of all words
`R^eta X^eps_r A ... X^eps_1 A` for `r<=6`, `eps_i=+1/-1`, `eta=0/1`.
For a word matrix `[p,b;ell,d]`, finite fixed points are precisely roots
of `ell*z^2+(p+d)*z+b=0`. Each genuine simple folded cycle appears `4r`
times: two sheets over each of its `2r` starting vertices. The resulting
sets agree exactly. Contracted path-segment checks agree with direct
folded traversals. All eight quotient profiles, plus the failed canonical
profile, are independently checked as projective permutations (30 circuit
word checks in total). Four small-field controls also agree between the
two cycle enumerators. All data are in `borel_exchange_certificate.json`.
A small checker correction also makes identity detection accept JSON list
matrices as well as in-memory tuples; the canonical negative control is
tested in both representations.

This is not a counterexample to the Alspach–Zhang conjecture, or even to
every bounded-short-exchange strategy. The first four-chord witness has
point-quotient circuit lengths `3669,428`, with odd-circuit word `U(1199)`.
It already colours the graph, without spanning. Raising the spanning
bound from five to six would merely fit the latest example; no uniform
six-chord assertion is justified.

### Alternative 1: allow even off-path circuits

The spanning condition was a sufficient convenience. A weaker, now
explicit criterion retains a local certificate for the other components.
Give each folded path edge sheet label zero. Give a chord `{i,j}` label
0 or 1 according as `a(z_i)=z_j` or `a(z_i)=z_(-j)`. For a folded circuit
off the distinguished path, let `s` be the sum of its sheet labels modulo 2.
A circuit of length `ell` lifts to two point circuits of length `ell` when
`s=0`, and one of length `2*ell` when `s=1`.

Consequently a matching colours if its distinguished gain is non-unit and
every off-path folded circuit has even length or odd sheet sum. The survey
proves this directly from the two-sheeted cover and matching monodromy;
the new certificate independently checks the predicted lengths for all
eight exchanges above. This condition is sufficient, not asserted necessary:
an odd point circuit can also be harmless when its word has even order.
The first four-chord witness has an off-path folded circuit of length 214
with sheet sum 1, explaining its single even point circuit of length 428.

This changes the next existence question: instead of requiring all path
segments to join one component, require that every leftover component pass
the parity test. The projective sheet data matter, and no theorem yet
guarantees such an exchange for every residual parameter.

### Alternative 2: use endpoint rotations without a size bound

The residual folded graph is simple, has degree 2 at `1,h`, and degree 3
elsewhere. Classical Hamiltonian-path parity therefore guarantees another
Hamiltonian path from `1` to `h`, paired with the natural path by endpoint
rotations. This is an application of Thomason's theorem, not a new general
graph theorem; see [Cameron–Edmonds, Theorem 0′, p. 817](https://www.numdam.org/item/10.5802/aif.1694.pdf).
The survey includes the short rotation-graph proof.

For the q=4096 obstruction, anchoring at 1 gives a spanning mate after
1147 rotations, with gain 460. Its matching differs from the all-chord
matching along a 569-chord circuit. Anchoring at h gives a mate after
1379 rotations, with gain 1022 when oriented from 1 to h; its difference
has three components containing 499,85,38 chords. These are independently
replayed and certified by quotient matrices and projective permutations.
Thus the method accommodates long and disconnected exchanges naturally.
No polynomial time bound is claimed.

Parity does **not** supply non-unit gain. At q=128,t=15,c=2, the first mate
from anchor 1 has gain 1 after 54 rotations; the other anchored mate has
gain 46. This example already has canonical gain 103, so it does not
disprove a rule restricted to canonical failures. At q=32,t=7,c=1 both
first anchored mates have unit gain; that graph is already handled by the
pentagon theorem. The earlier reordered-label control also prevents a
gain theorem based only on the underlying path-plus-matching structure.

The new audit exhausts all 62 spanning paths in 11 small residual folded
diagrams through q=32, checking pairing, reversibility, and 614 exact
one-step gain updates. It checks both anchored mates at 24 specified
larger/control inputs: all 48 walks finish (55,690 rotations in total),
and all 42 mates at the 21 canonical-failure inputs have non-unit gain.
These inputs include all ten residual failures from the complete q=2048
scan and the ten encountered before the partial q=4096 scan stopped.
This is finite evidence only. The target is a field-sensitive argument
forcing non-unit gain on some suitable rotation-generated path, not an
assumption that the first mate always works.

### Coverage and reproducibility

The stress implementation uses O(q) field storage and contracts uncut path
segments for each short exchange. Logarithmic multiplication is checked
against polynomial arithmetic; the cut checker is compared with 251
direct small-field exchange traversals and 48 canonical diagrams.

At q=2048 the residual double-inverse-trace-one scan is **complete**:
15,869 representatives modulo Frobenius and `c <-> c+t` cover 349,118
normalized parameters. Ten representatives (220 parameters) have canonical
gain 1; each has a successful spanning exchange with at most five chords.
Combining this with the proved trace-test theorem certifies all 1,395,372
normalized generating parameters with full-order x at q=2048. This does
not cover other odd orders in that group.

The q=4096 scan is **partial** and stops on the spanning-bound obstruction.
Its ten encountered failures are examples, not the total for the field.
The global counterexample lower bound remains 352,440 vertices.

```sh
python3 code/borel_exchange_stress.py 2048 4096 --stop-on-failure --json-output code/borel_exchange_stress_results.json
python3 code/borel_exchange_certificate.py --json-output code/borel_exchange_certificate.json
python3 code/borel_rotation_experiment.py --stress-results code/borel_exchange_stress_results.json --json-output code/borel_rotation_results.json
```

The full-projective-cycle family and the full conjecture remain open.
The useful progress here is an independently certified obstruction to one
overstrong target, a larger complete finite check, and two better-separated
existence problems. None establishes that a full proof is close.

## Further follow-up: a uniform paired-square theorem and a 13/16 bound

Allowing harmless extra circuits now yields a further **uniform theorem**,
not merely larger finite coverage. For each fixed full-cycle parameter t,
the guaranteed proportion of generating translations has improved from
`3/4 + O(q^-1/2)` to **at least `13/16 - O(q^-1/2)`**. This is a lower
bound within the characteristic-two full-projective-cycle family, not a
percentage of the whole conjecture proved.

### Exact existence and the paired repair

Put `Delta=c(c+t)` and assume both inverse traces are 1. Then:

1. A simple alternating square in the folded graph exists **if and only if**
   `Tr(1/Delta)=0`.
2. When it exists, the square is unique.
3. Exchanging the square from the all-chord matching colours at least one
   of the two Cayley graphs with translations `c,c+t`.

The existence proof uses the four signed two-chord words, not an
extrapolation from an enumeration. For

\[
 W=R^\eta X^{\epsilon_2}A X^{\epsilon_1}A,
\]

the unreflected traces are `(c+t)^2` for equal signs and `c^2` otherwise;
all reflected traces are `Delta`. The first two trace conditions rule out
unreflected fixed points. The four reflected words have two finite fixed
points each exactly when `Tr(1/Delta)=0`. Each simple square contributes
eight fixed-point occurrences, so there is at most one. A fixed point of
`R(XA)^2` produces a square: the walk avoids infinity, repeated folded
vertices would force a parallel path/chord pair, and a folded middle-edge
loop would force a second such loop at a different vertex. All are excluded.
The survey supplies the explicit matrix and these degeneracy checks.

Let the two selected path-edge starts be `i<j`, with `i+1<j`. A crossing
square has a spanning complement and already works for both translations.
A nested square has outer chord `(i,j+1)` and inner chord `(i+1,j)`.
Its distinguished path uses only the outer chord and has gain
`k_i/k_(j+1) != 1`. The inner interval closes to a folded circuit of length
`ell=j-i`. If ell is even, both translations work. If ell is odd, the
matching colours exactly when the inner chord's sheet label is 1.
Replacing c by c+t toggles that label, so precisely one translation works
in the odd case. This explains why the earlier spanning restriction lost
useful two-chord repairs.

For example, the old q=128,t=15,c=35 residual case required five chords
to obtain a spanning path, but its two-chord non-spanning exchange already
colours. The theorem does not assume that all extra folded circuits are
absent: it controls the one created by a nested square.

### The improved count is exact before applying an error bound

Let `b=1/t` and keep the unshifted Kloosterman convention

\[
 K(v)=\sum_{u\ne0}(-1)^{\operatorname{Tr}(u+v/u)}.
\]

Define `K0=K(b^2)`, `K1=K(b^3+b^4)`, and `K2=K(b^2+b^4)`. The number of
translations in the residual class admitting the unique square is

\[
 R_t=\frac{q-3+2K_0+2K_1+K_2}{8}.
\]

This follows by expanding the three trace indicators and evaluating

\[
 S(A,B)=\sum_{c\ne0,t}(-1)^{\operatorname{Tr}(A/c+B/(c+t))}
       =(-1)^{\operatorname{Tr}(b(A+B))}K(b^2AB)-1
\]

for nonzero A,B. The substitution is `s=c/(c+t)` and then `u=b*B*s`.
Here `Tr(b)=1`, `b!=0,1`, and `K(v^2)=K(v)`; all three Kloosterman
arguments used in the bound are nonzero. The survey includes the full
indicator expansion, rather than assuming independent trace conditions.

The `R_t/2` pairs are disjoint from the original trace-test family, and
each contributes at least one additional successful translation. Therefore
the combined constructions colour at least

\[
 L_t=\frac{13q-23-2K_0+2K_1+K_2}{16}
     \ge \frac{13q-23-10\sqrt q}{16}
\]

of the q-2 generating translations for each t. The last inequality uses
the standard Hasse bound `|K(v)|<=2*sqrt(q)` in the unshifted convention;
see [Ahmadi–Granger, Lemma 2.1 and §6.3](https://arxiv.org/pdf/1104.3882v2).
After division by q-2, the lower bound tends to 13/16 uniformly in t.
The actual combined success proportion may be higher.

### A real obstruction to strengthening the theorem

At q=32,t=7,c=29, modulo `theta^5+theta^2+1`, the unique square has cut
starts `7,14`. Its nested inner circuit has length 7 and sheet sum 0;
the point-quotient circuit lengths are `19,7,7`, so this matching fails.
The partner c+t=26 succeeds. The canonical matching for c=29 already
works, but this does not make the square itself a valid repair.

More importantly, in q=2048 modulo `theta^11+theta^2+1`, take
`t=278,c=1567`. Both the canonical matching and the unique square fail.
The canonical odd circuit has length 361 and identity word. The square's
cut starts are `926,1015`, its distinguished gain is 854, and its quotient
circuits have lengths `1871,89,89`. The two 89-circuits have matrix
`[1104,108;590,109]` of nonzero trace, hence odd-order monodromy. They
remain odd upstairs despite the nontrivial distinguished word `U(91)`.

The partner c+t=1801 succeeds with the same square. The original c=1567
graph also colours, by the previously certified five-chord spanning
exchange with cuts `63,213,273,432,1019` and gain 1375. Thus neither
unconditional square success nor “canonical or square always works” is
valid. The new theorem deliberately asserts only the paired guarantee.

### Independent audits and the next target

`code/borel_square_exchange.py` checks all eight short-word identities as
polynomial identities. It constructs squares from quadratic fixed points
using an independently checked Artin–Schreier root table. Through q=4096,
837,116 parameter representatives cover 9,786,998 normalized generating
parameters. Both c partners are tested; only t is reduced by Frobenius.
Every trace representative satisfies the exact character-sum formulas.

There are 1,619 independent graph-search, folded-component and full
projective-word audits, plus 363 no-square graph checks. The graph
comparisons are exhaustive through q=128 modulo the stated symmetries and
sampled at larger fields, including every geometry/sheet-parity type at
each t and the double-failure control. All 13,100 unsuccessful square
parameters in the direct representative set are additionally checked for
canonical gain: exactly the displayed q=2048 representative fails both.
This last check says nothing about the parameters with no square.

Two further checks construct the full 32,736-vertex Cayley graphs for
q=32,t=7,c=26 and c=29, without relying on quotient-word orders. The
successful matching has 496 complementary circuits of length 38 and 32
of length 434; every vertex's three edge colours are checked. The failed
matching has 496 circuits of length 38 and 64 odd circuits of length 217.
This independently validates both outcomes of the sheet-parity test.

| q | All normalized generating parameters | Added square successes | Combined successes |
|---:|---:|---:|---:|
| 1024 | 408,800 | 45,140 | 351,780 |
| 2048 | 1,395,372 | 151,877 | 1,198,131 |
| 4096 | 7,860,480 | 859,896 | 6,756,216 |

These are Frobenius-weighted counts for the original trace test and the
new residual square construction alone. The prior complete full-cycle
verification at q=2048, using other matchings, is unchanged. The new
square-only audit at q=4096 is complete for its test; it does not turn the
earlier partial general-exchange scan into a complete colouring result.

```sh
python3 code/borel_square_exchange.py --audit-canonical-failures --json-output code/borel_square_results.json
```

The next existence problem is now a more explicit union of two classes:
parameters with `Tr(1/Delta)=1`, where no alternating square exists, and
the unsuccessful member of a nested odd-square pair. The latter cannot
be dismissed by assuming that the canonical matching succeeds. Longer
exchanges or a different matching construction are still needed. The new
positive-density theorem is genuine progress, but neither it nor the
finite tests proves the full-projective-cycle family or the full conjecture.

## Further follow-up: exact three-chord existence and a uniform short-exchange obstruction

The attempt to extend the paired-square theorem to three chords exposed a
real obstruction, and the short-word calculation now gives a uniform
theorem about it. The best proved colouring density remains **at least
13/16 asymptotically**. The new result is not a higher colouring density:
it proves that a subfamily of asymptotic density **1/32** has no alternating
exchange of one, two or three chords from the all-chord matching.

This does **not** say that these graphs are uncolourable, or that their
canonical matchings fail. It means that the proposed universal
three-chord construction is impossible even if non-spanning complements
are allowed. Longer exchanges and other starting matchings remain in play.

### An exact existence theorem, including the missing degeneracy

Write `u=c+t`, `Delta=cu`, and exclude the already-colourable pentagon
translations `c=1,t+1`. In the residual class with both inverse traces one,
put

\[
 \tau_1=u(c+1)^2,\qquad \tau_2=c(u+1)^2,\qquad
 J_t(\Delta)=\Delta^3+\Delta+t^3.
\]

The exact number of simple three-chord alternating circuits is

\[
 n_3(t,c)=2-\operatorname{Tr}(1/\tau_1)-\operatorname{Tr}(1/\tau_2)
                  -\mathbf1_{J_t(\Delta)=0}.
\]

In particular there are at most two. The sixteen signed/reflected words
have four trace types: `u(u+1)^2`, `tau1`, `tau2`, and `c(c+1)^2`, with
multiplicities 2, 6, 6, 2. The first and last are nonsplit under the two
residual trace assumptions. Each genuine six-cycle contributes twelve
fixed-point occurrences to the other types.

The raw fixed-point count is not always the cycle count. A closed walk can
use the folded middle loop at h, traverse its incident chord twice, and
go around a triangle at the other end j. There is at most one such
configuration. Its existence is exactly

\[
 \Delta=k_{j-1}+k_{j+1}
       =\frac{t^2(t+\Delta)}{(t+\Delta+1)^2},
\]

which is equivalent to `J_t(Delta)=0`. It contributes twelve degenerate
occurrences. The survey proves that these are the only degeneracies,
including checking why walks through infinity would force an excluded
pentagon translation.

For a concrete correction control, in GF(128) modulo `theta^7+theta+1`,
`t=13,c=16` gives `Delta=85`: the raw count predicts one six-cycle, but
there is none. All twelve walks are degenerate. At `t=13,c=118`, the
correction instead reduces the count from two to one.

### A proved 1/32 obstruction subfamily, not an independence heuristic

Let M_t count translations for which all five inverse traces

\[
 \operatorname{Tr}(1/c),\quad\operatorname{Tr}(1/u),\quad
 \operatorname{Tr}(1/\Delta),\quad
 \operatorname{Tr}(1/\tau_1),\quad\operatorname{Tr}(1/\tau_2)
\]

are one. These conditions exclude every alternating circuit with at most
three chords. The character-sum proof establishes

\[
 |32M_t-(q-3)|\le114\sqrt q+36,
 \qquad M_t\ge(q-39-114\sqrt q)/32.
\]

Consequently this obstruction occurs for every full-cycle t once
`q>=16384`, and its parameter proportion tends uniformly to 1/32. The
total number of non-pentagon translations with no such short circuit is
`M_t+e_t`, where `0<=e_t<=6`; the extra cases come from the cubic correction.

The five trace functions reduce to rational functions with simple poles
at `0,t,1,t+1`. Every nonzero binary combination retains a pole. A
combination with k poles defines a geometrically irreducible
Artin–Schreier curve of genus k-1; applying Hasse–Weil to all 31
combinations gives the explicit bound. The proof accounts for removed
poles and the two points over infinity, which determine the constant term.
The genus formula is supported by
[Pries–Zhu, Lemma 2.6](https://www.numdam.org/item/10.5802/aif.2692.pdf);
the standard point-count bound is recalled on
[Kresch–Wetherell–Zieve's author page](https://websites.umich.edu/~zieve/papers/kwz.html).
The resulting obstruction theorem is our application, not a claim from
those papers.

Since the symmetric difference of two perfect matchings decomposes into
simple alternating circuits, these parameters actually exclude *any*
change of at most three chords from the all-chord matching, including
disconnected changes of that total size.

### A replacement for the root-only test: exact affine circuit certificates

For any folded perfect matching, put `w_i=t+z_i` and assign to a directed
chord i→j the affine map

\[
 y\longmapsto \alpha_{ij}y+\nu_{ij},\qquad
 \alpha_{ij}=k_j/k_i,\qquad
 \nu_{ij}=w_j+\alpha_{ij}w_i.
\]

Path-edge maps are identity. Composing around each off-path circuit gives
two field elements `(alpha_C,nu_C)`. The lifted matching colours **if and
only if** its root gain is not one and every off-path circuit C satisfies
at least one of:

- C has even length;
- C has sheet sum one;
- `alpha_C=1` and `nu_C!=0`.

The last condition is precisely a nontrivial unipotent word. The proof
uses the frame `X^i` over each projective point; the survey displays the
exact upper-triangular edge matrix. This works for arbitrary folded
matchings, not only a fixed number of exchanges. It strengthens the
earlier sufficient sheet-parity test to a necessary-and-sufficient one.

The affine maps are unchanged by `c -> c+t`; the sheet sum toggles only
on circuits with an odd number of chords. For a circuit with just one or
two distinct chords, alpha cannot be one. Thus the earlier sheet-parity
test is already exact for all exchanges of at most three chords.

The new test is genuinely stronger for larger changes. In GF(64), modulo
`theta^6+theta+1`, take `t=2,c=6`, select folded chords `{19,24},{27,30}`,
and complete the matching with path-edge starts
`1,3,5,7,9,11,13,15,17,20,22,25,28,31`. Its root gain is 30 and its
point profile is `23,21,21`. The off-path folded 21-circuit has sheet sum
zero and affine data `(1,48)`: the two odd point circuits have nonidentity
involution words `[34,48;5,34]` and `[51,15;35,51]`. Every circuit therefore
becomes even upstairs, of length 46 or 42. The old sufficient parity test
rejects this matching; the exact affine test certifies it. This is checked
by full projective permutations and retained as `unipotent_rescue` in the
JSON. The canonical matching also colours this graph, so this witness is
not claimed as an additional parameter family.

Here is the precise failure of the paired three-chord shortcut. In the
same GF(128), take `t=8,c=84` and cuts `18,33,46`. The root gain is 124,
but the off-path circuits have lengths 15 and 13, each with one chord.
Their sheet sums are 0 and 1. For the partner `c+t=92` the sheet sums are
1 and 0. The two translations therefore fail on different odd circuits.
Their point profiles are respectively `73,15,15,26` and `73,30,13,13`.
The root word is `U(97)` in both cases; that alone does not colour either
matching. Full projective-permutation checks certify these assertions.

### Reproducibility, scope and revised next step

```sh
python3 code/borel_three_exchange.py --json-output code/borel_three_results.json
```

The complete trace audit through q=4096 has 836,492 representatives,
covering 9,780,236 normalized non-pentagon parameters. All 31
character-sum bounds are checked at each of 311 full-cycle trace
representatives. There are 3,905 graph-count comparisons, exhaustive
through q=512 modulo Frobenius and the folded translation pairing,
and sampled by trace/degeneracy signature thereafter. Fifty-two
parameters also have independent fixed-point-word enumerations.

For the affine criterion, all 1,130 perfect matchings on the 96 generating
parameter representatives through q=32 are independently checked against
frame matrices and 2,098 full projective-permutation component words.
All residual three-chord exchanges through q=128 are additionally audited
for both translations: 252 matching and 452 permutation-component checks.
The JSON distinguishes these scopes.

At q=4096, the complete trace count identifies 244,488 normalized
non-pentagon parameters with no exchange of at most three chords. A
single-trace control at q=16384, modulo `theta^14+theta^5+1`, has `t=3`
and `M_t=486`. Neither computation is an all-parameters colouring check
at those field sizes. The earlier complete full-cycle result at q=2048,
the partial general scan at q=4096, and the global smallest-counterexample
bound of 352,440 vertices are unchanged.

The revised target is a **global matching or rotation argument with all
circuit obstructions tracked**, rather than a universal three-chord
repair. The affine certificate makes that target exact and inexpensive
to verify: odd, sheet-zero circuits must have translation rather than
identity or semisimple monodromy, and the distinguished gain must not
be one. What is still missing is a structural reason, using the actual
projective ordering, that some matching meets all those requirements
simultaneously. No descent theorem or uniform existence result has been
proved here. The all-parameters full-cycle family and the full conjecture
remain open.

## Global recolouring follow-up (6 September 2026)

### An unrestricted-length search that automatically controls off-path circuits

Add a distinguished closing edge `e*={1,h}` to the folded graph F, giving a
cubic graph H. It has an explicit proper three-edge-colouring: all chords
have colour 0, odd-start path edges colour 1, and even-start path edges plus
the closing edge colour 2. Edge identities are retained when parallel edges
occur; H is simple in the non-pentagon residual class.

In **any** proper three-edge-colouring of H, either colour class avoiding
e* is a perfect matching N of F with every off-path circuit even. Its lift
therefore colours the Cayley graph exactly when `gamma(N)!=1`. Conversely,
every matching with even off-path circuits extends to such a colouring of H:
its root path also has an even number of vertices, since h is even.

This gives a global search space without a chord-count cap or extra circuit
tests. It is a sufficient restriction, not the full affine criterion: the
odd-folded-circuit unipotent rescue from the previous section lies outside it.

Swapping two colours along an entire two-colour circuit is a Kempe change.
Two simple subfamilies independently switch subsets of the natural (0,1)
or (0,2) circuits that avoid e*, then use colour class 0. The selected path
edges have a common starting-index parity. Testing these unrestricted-length
families repairs 100 of the 102 stored canonical-failure representatives
through q=1024, accounting for **926 of the 942 weighted failures**.
It also repairs the specified q=2048 square-plus-canonical failure and
q=4096 five-chord-spanning obstruction. These graphs were already known to
be colourable; this does not increase exhaustive coverage or prove a density.

### An entire three-colour Kempe class can fail

In GF(256), modulo `theta^8+theta^4+theta^3+theta+1`, take
`t=15,c=151`, with partner `c+t=152`. H is simple with 128 vertices and
192 edges. In its natural colouring **each of the three two-colour
subgraphs is a single Hamiltonian circuit**. Therefore every Kempe change
just globally swaps two colour names. The entire reachable class contains
exactly six labelled colourings, all global permutations of the start.

There are only two eligible matchings throughout this class: all chords
and the canonical matching. Both have gain one and a single point circuit
of length 257 with identity word. Thus no sequence of three-colour Kempe
changes from this start can give a successful lift, regardless of length.
The checker exhausts every move from all six states and certifies closure.
This is a finite exact proof of an unbounded obstruction, not a timeout.

Other Kempe classes do work. Select the folded chord `{3,122}` and complete
the matching with path-edge starts `1,4,6,...,120,123,125,127`. The root
spans all 128 vertices, with gain **28** and word **U(187)** on its point
circuit of length 257. The lifted complementary circuits have length 514.
The same matching works for both translations.

Its difference from the canonical matching is one 120-edge circuit using
1, 60 and 59 edges of the three natural colours. Switching this matching
circuit is not a Kempe change: it involves all three colours. The repair's
completed H-colouring belongs to another three-colour Kempe class. Its
difference from the all-chord matching has two alternating circuits with
38 and 25 chords. The root words and both translations are checked by
matrix multiplication and full projective permutations.

### A temporary fourth colour gives a concrete escape, but not a proof

McDonald, Mohar and Scheide prove connectivity of four-edge-colourings under
Kempe changes for simple subcubic graphs. It connects existing colourings;
it does not produce one satisfying our gain condition.
See [their paper, Theorem 3 and Lemma 7](https://arxiv.org/pdf/1005.2248).

The new script implements the paper's vertex-deletion method on F, then
separately attempts to extend each change across e*. For the rigid example,
121 changes on F extend to **122 changes on H**, temporarily using at most
28 edges of the fourth colour. The final state is the known successful
three-colouring. An independent replay checks that each switched set is
connected and maximal in its two-colour subgraph, and that every state is
proper. No intermediate boundary relaxation is needed in this H certificate.
The extension routine can fail elsewhere; it is not a general constructive
implementation of the cubic theorem.

Crucially, the target was supplied in advance. The sequence proves escape
from the rigid class but supplies neither gain descent nor a way to discover
a suitable target in every graph. Four-colour connectivity is not a substitute
for arithmetic nonvanishing on the three-colour states.

### Reproducibility and the remaining proof target

```sh
python3 code/borel_kempe_exchange.py --json-output code/borel_kempe_results.json
```

The small audit covers every even-off-path folded matching through q=32:
850 matchings on 96 generating parameter representatives, with independent
matrix/permutation checks. For the simple folded graphs among them, 168
buffer-colour sequences comprising 1,292 changes are also replayed.
Stored-failure tests stop at the first successful subset; the two failures
exhaust both circuit-subset families. The JSON retains the entire closed
six-state-class certificate, the repaired matching and both buffer sequences.

The new sufficient target is: **in the projectively ordered folded graph,
some proper three-edge-colouring of H has an eligible class with non-unit
root gain.** Equivalently, find a folded matching with even off-path circuits
and non-unit gain. Arbitrarily long three-colour Kempe search from a fixed
start cannot prove this. A fourth colour makes different classes reachable,
so a next approach can study returns to three colours and their gain changes,
or work directly with changes involving all three colour classes. Any such
argument must exploit the actual projective ordering, not just cubicity or
reconfiguration connectivity.

This is a sharper global formulation and a verified escape mechanism, not
a new infinite colouring family. The **13/16 asymptotic guarantee**, complete
full-cycle check at q=2048, partial general scan at q=4096, and **352,440-vertex
smallest-counterexample lower bound** are unchanged. A proof of the full
conjecture is not yet in hand.

## Two-buffer follow-up: seven changes without a supplied target (6 September 2026)

The preceding 122-step escape used a known successful target and up to
28 temporary-colour edges. The new construction needs no target and uses
only two temporary edges. Its endpoint-pairing return condition is proved;
its uniform arithmetic success is not.

### One temporary edge provably cannot change the Kempe class

In a loopless cubic graph with ordinary colours 0,1,2 and precisely one
colour-3 edge uv, the endpoints miss the same ordinary colour. Indeed,
each ordinary colour is missing at an even number of vertices, and only
u,v can miss one. There is therefore a unique completion to three colours.

Every allowed Kempe move that keeps at most one temporary edge induces
either no change or one ordinary Kempe change on completions. An ordinary
alternating path closes through uv; an ordinary circuit remains a circuit.
Introducing or deleting the temporary edge leaves the completion fixed.
A swap involving its colour and another present colour would create two
temporary edges, except for a harmless parallel two-edge circuit.

Thus an arbitrarily long one-buffer sequence cannot change the three-colour
Kempe class. For the rigid q=256 example, the reachable one-buffer class has
exactly `6*(192+1)=1158` states. All 5,778 allowed moves are checked to stay
in this set. This is a structural lower bound, not a search cutoff.

### A seven-change rule, with an exact return criterion

After renaming ordinary colours, choose a colour-1 pivot uv whose incident
colour-0 edges uw,vz have four distinct endpoints. Recolour uw and vz to 3
(two changes), recolour uv to 0, and swap the unique (0,2)-path from w to z.
Now u,v miss colour 1 and w,z miss colour 2. The (1,2)-subgraph has two
paths joining those four endpoints.

| Endpoint pairs of those paths | Can one path swap restore both temporary edges? |
|---|---|
| uw and vz | No |
| uv and wz | Yes |
| uz and vw | Yes |

In a successful row, swapping either path makes the missing colours agree
at the ends of each temporary edge. Restoring those two edges gives seven
Kempe changes in total, with temporary-edge counts `1,2,2,2,2,1,0`.
Every step is a swap on a maximal two-colour component. The proof is only
an endpoint-pairing argument; it applies to any properly three-edge-coloured
cubic graph with the stated seed. It does not assert a different Kempe
class or non-unit root gain on return.

In the folded graph we take temporary colour 0 to be the chord colour and
try each natural path edge, including the closing edge, as pivot. No target
matching or random seed is supplied. Both path choices must be tested.

### The rigid example now has an optimal two-buffer escape

At q=256, t=15,c=151, start with pivot `{1,2}` and temporary chords
`{1,41},{2,65}`. The (0,2)-path has length 69. The two (1,2)-paths have
lengths 17 and 39, with endpoint pairs `{1,2}` and `{41,65}`. Swapping the
17-edge path and restoring both chords to colour 2 gives eligible gains
**255 and 31**. The first matching has root word **U(125)** and point
profile `117,54,16,54,16`, so all lifted complementary circuits are even.

This is a different successful target from the earlier supplied one.
Two temporary-colour edges are necessary by the completion lemma and
sufficient by this explicit sequence. Seven is the certified number of
changes, not a proved minimum number of changes. The swaps can traverse
long paths; this is not another bounded-chord claim.

Exhausting all 128 pivots in this example gives 14 blocked pairings,
66 pivots returning only global permutations of the natural colouring,
and 48 with successful returns. Pivot `{2,3}` is blocked; pivot `{3,4}`
returns but has both gains one. The return condition alone therefore
cannot be substituted for the gain test. At q=128,t=15,c=55 the first
pivot is blocked, and only the second branch of the second pivot succeeds.

### Scope of the target-free audit

```sh
python3 code/borel_two_buffer.py --json-output code/borel_two_buffer_results.json
```

The rule repairs all 102 stored representative canonical failures through
q=1024, accounting for all 942 Frobenius-weighted failures, with at most
three pivots tested on these inputs. It also repairs all 20 stored larger
cases at q=2048 and q=4096. These successes are finite evidence, not a
uniform three-pivot or seven-change theorem.

Small controls cover 96 generating parameter representatives through q=32:
3,636 pivot/temporary-colour choices, comprising 2,970 returning choices,
446 blocked pairings and 220 parallel seeds. The two branches of each
return and their two eligible matchings receive 11,880 independent
matrix/permutation audits. Another 1,914 one-buffer states and 12,568
moves check the completion lemma. The rigid seven-change certificate and
its one-buffer closed class are separately audited.

A further complete single-trace control at q=8192, modulo
`theta^13+theta^4+theta^3+theta+1`, takes t=13. It tests all 4,095 folded
translation pairs, covering 8,190 generating translations for this fixed
trace. The only canonical-failure pair is c=2468,2473; the first pivot
repairs it with an eligible gain 4513 and root word U(6778). Both partners
are checked by full projective permutations. This is **not** a complete
all-traces q=8192 result and does not alter the global counterexample bound.

An additional exhaustive pivot audit on the 102 stored representatives
finds 100, 772, 2,304, 9,512 and 50,780 nontrivial returns at q=64,128,256,
512,1024 respectively. None of these 63,468 returns has both eligible
gains one. These are counts with pivot/branch multiplicity, not distinct
colourings. The optional `--exhaustive-stored` flag reproduces this larger
audit, retaining per-case return counts; the default stops at first success.
The observation suggests a nonvanishing question but is not its proof.

### Why a graph-theoretic return proof is not enough

There is a stronger abstract warning than a return to the same colouring.
Over GF(32), modulo `theta^5+theta^2+1`, take t=6,c=3,Delta=15 and reorder
the correct label set as

`1,25,27,9,22,30,4,3,11,12,14,20,28,19,17,6`.

The endpoints remain 1 and t and the canonical gain is one. At pivot
`{1,2}`, the second seven-change return changes the colouring even modulo
global colour permutations, but both eligible gains remain one. The
checker certifies the entire proper recolouring sequence and the scalar
path gains. This is **not** the projective order, has no claimed Cayley
lift, and does not refute the proposed rule for actual generating pairs.
It shows that the arithmetic conclusion cannot follow just from the
affine label set, the return pairing and a nontrivial colour change.

The concrete sufficient target is now: for every remaining projectively
ordered canonical failure, **some natural path/closing-edge pivot admits
a seven-change return with at least one non-unit eligible root gain**.
The search is target-free, and its combinatorial return condition is exact.
What is missing is a proof that a successful pivot exists and that the
actual projective recurrence prevents the relevant gain cancellations.
The full conjecture and full-cycle all-parameters family remain open;
the 13/16 asymptotic guarantee and 352,440-vertex global lower bound remain
unchanged.

## 6 September follow-up: collective nonvanishing via signed root-path lattices

The new direction is to force **some** good gain by an integer relation
among candidate paths, rather than prove noncancellation for each return.
This supplies a proved sufficient criterion and a sharper open target,
not a proof of the full-cycle family or the full conjecture.

Orient the m=q/4 folded chords and put g_j=k_u/k_v on chord u to v.
Every g_j is non-unit. For each eligible matching N, let r_N be the signed
chord-incidence vector of its root path from 1 to h. The root formula is

`gamma(N) = product_j g_j^(r_N[j])`.

Let L be the integer lattice generated by these vectors. If
`sum_N a_N r_N = d e_j` with gcd(d,q-1)=1, simultaneous unit root gains
would force g_j^d=1 and hence g_j=1, a contradiction. All the matchings
here have even off-root circuits, so some matching colours its Cayley
lift. In particular, one doubled coordinate `2 e_j in L` suffices. A
finite lattice index coprime to q-1 also suffices, but is stronger.
This elementary lemma is now proved in the survey. It is an integer
lattice statement, not a binary cycle-space argument. Rational full
rank alone is insufficient: odd torsion can matter.

The specified family comprises the natural eligible matchings and both
eligible matchings from both branches of every returning seven-change
natural path/closing-edge pivot, using natural colour-0 chords as the
two temporary edges. There is no supplied target colouring. On a
canonical failure the two natural gains are one, so a lattice certificate
forces a successful seven-change return.

### Exact audit and genuine limitations

`code/borel_gain_lattice.py --stored --json-output code/borel_gain_lattice_results.json`
reproduces the lattice audit. It checks all 51 distinct folded diagrams
from the stored canonical failures: 1,4,6,11,29 at q=64,128,256,512,1024.
Restoring partners and Frobenius weights gives 102 direct cases and the
same 942 weighted failures already known to colour. **Every one has a
doubled coordinate in L.** This is not additional exhaustive coverage.

The seven displayed exact lattice controls are:

| q | t | c | Nonzero rows up to sign | Rank | Index |
|---|---|---|---|---|---|
| 64 | 7 | 26 | 27 | 16 | 2 |
| 128 | 8 | 20 | 52 | 32 | 4 |
| 128 | 15 | 35 | 63 | 32 | 2 |
| 256 | 15 | 151 | 124 | 64 | 4 |
| 256 | 25 | 13 | 82 | 64 | 4 |
| 512 | 47 | 10 | 151 | 128 | 4 |
| 512 | 53 | 66 | 217 | 128 | 4 |

All contain twice the entire coordinate lattice. Exact integer
combinations are independently checked by direct summation. The ten
designated controls, including the three q=1024 exceptions below, also
check both lattice inclusions: every input row reduces in the basis,
and every basis row is verified as an integer combination of input rows.
The broader stored-case audit computes exact echelon bases and checks
input-row membership, without retaining all change-of-basis matrices.
Pivot/branch source seeds make every saved row reproducible.

The larger checks rule out stronger guesses for *all canonical failures*:

- q=1024,t=49,c=80: rank 255 out of 256; exactly three doubled coordinates
  belong to L (zero-based coordinates 120,227,235).
- q=1024,t=49,c=910: full rank, index 6, and 253 doubled coordinates.
  The index shares the factor 3 with q-1=1023.
- q=1024,t=62,c=836: full rank, a 477-bit index whose gcd with 1023 is 11,
  but five doubled coordinates still belong to L.
- q=1024,t=287,c=683: rank 253 out of 256, with doubled coordinates
  8 and 219. This last case is in the broader stored audit.

These exceptions are already covered by the earlier trace theorem and
do **not** refute the stronger inclusion restricted to the hard class.
In fact, all ten stored folded diagrams with both inverse traces equal
to one contain `2 Z^(q/4)` and have power-of-two index. Five of these
also have inverse trace of Delta equal to one, excluding squares.
No uniform conclusion is inferred from these ten examples.

The collective certificate can be visibly small even when the paths
are not. At q=256,t=25,c=13, every nonzero candidate row has support at
least five, but two returned rows differ by exactly the coordinate of
chord {17,18}. They come from the second branch at pivot {1,2}, excluding
colour 2, and the first branch at pivot {16,17}, excluding colour 0.
Their gains are 137 and 39; their ratio is k_17/k_18 != 1. Both receive
matrix/permutation checks. This graph has a parallel path/chord pair and
is already covered by the trace theorem; it illustrates the mechanism.
The two displayed q=512 cases instead have all three inverse traces one.
Neither admits the simplest certificate `+/-r +/-s = d e_j`, d=1 or 2,
even allowing the zero row and repeated rows. Their verified lattice
relations go beyond that two-row shortcut.

### A genuine both-unit colouring and an abstract obstruction

There is now an actual projective counterexample to the broad shortcut
“every non-natural three-colouring has some non-unit eligible root gain”.
Over GF(128) modulo theta^7+theta+1, take t=8,c=20. The saved proper
colouring is not a global permutation of the natural colouring, but
both eligible gains are one. Its point profiles are 65,48,16 and 129;
both root words are identity. At the partner c=28 the profiles are
65,24,24,8,8 and 129, again with identity root words. All four matching
certificates are checked by matrices and full projective permutations.
The complete 96-edge colour assignment is in the JSON file.

This was found by enumerating pinned folded matchings and complement
colour orientations, stopping at the first example. It is not an
exhaustive classification of q=128 colourings. More importantly, it is
**not asserted to be a direct seven-change return** and does not refute
the proposed local-seed rule on actual projective diagrams.

For the earlier abstract GF(16) reordered labels
`1,14,6,10,13,9,5,2`, Delta=11, the entire seven-change family has just
one nonzero signed row up to sign, `(1,1,-1,1)`, with unit gain. Its
rank-one lattice contains no nonzero coordinate multiple. The ordering
is not projective, and this is not a Cayley counterexample. It shows
that the label set and return rule alone do not guarantee the inclusion.

### Next proof target and assessment

**Superseded by the next follow-up:** the fixed-family coordinate-inclusion
target in this paragraph is false at q=2048. The conditional lemma and
all finite results above remain valid.

Prove that, for every *remaining projectively ordered canonical failure*,
the specified root lattice contains `2^s e_j` for some j and s>=0, or a
coordinate multiple coprime to q-1. The stronger target `2 Z^m subset L`
remains plausible specifically in the residual class, but is not needed.
An eventual proof must derive a relation from the projective ordering
or show that a hypothetical obstructing character contradicts that
ordering. Additional finite-rank checks alone cannot provide this step.

This direction is worth pursuing because it gives a precise collective
obstruction and avoids demanding success from an arbitrary returned
colouring. It is not evidence that the full conjecture is almost proved.
Uniform inclusion remains the central unsolved step even in the
full-cycle PSL(2,2^m) family, and other groups and generator orders remain
outside this route. The 13/16 asymptotic guarantee, full all-parameter
coverage through q=2048, partial q=4096 general scan, single-trace q=8192
control, and 352,440-vertex global lower bound are unchanged.

## 6 September follow-up: a fixed-family obstruction and a one-prime criterion

The previous turn's sufficient lattice lemma is correct, but its proposed
uniform inclusion is now **disproved**, even in the no-square residual
class. This is a limitation of that proof target, not a Cayley snark or
a failure of the original seven-change colouring rule.

### A proved arithmetic improvement

The q/4 chord gains generate the whole group GF(q)^*. Indeed, their two
orientations give q/2 distinct nonidentity elements: k -> k/(k+Delta) is
injective, and swapping the two labels of a chord inverts the gain.
Together with 1, these cannot fit in a proper subgroup, whose order is
at most (q-1)/3 because q-1 is odd.

Consequently **full column rank of the root-row matrix modulo any one
prime p dividing q-1 suffices**. If every candidate gain were a pth power,
the row matrix would annihilate the chord-logarithm vector modulo p.
That vector is nonzero because the chord gains generate GF(q)^*.
Full rank is a contradiction. Thus some candidate is not a pth power,
in particular not one, and its even off-root complement gives a colouring.

There is also an index refinement. If L has finite index D and H is the
subgroup generated by the candidate root gains, the gain map induces a
surjection Z^(q/4)/L -> GF(q)^*/H. Hence the index of H divides
gcd(D,q-1). In particular `(q-1) does not divide D` is sufficient;
coprimality of D and q-1 is unnecessary. This handles the old q=1024
index-6 and index-divisible-by-11 controls using a different prime.
It does not handle rank-deficient matrices by itself.

### Exact counterexample to every coordinate-isolation target

Over GF(2048), modulo theta^11+theta^2+1, take

`t=343, c=1766, partner c=1969`.

The generator has order 2049, the canonical gain is one, and all inverse
traces of c,c+t,Delta are one. There are no parallel path/chord edges or
alternating squares. Of the 1024 original path/closing-edge pivots,
956 return and 68 are blocked. The 3826 eligible root paths, including
the natural two, yield exactly **425 nonzero rows up to sign on 512
coordinates**. Their ranks modulo 23 and 89 are both 425; 2047=23*89.

Three saved vectors in GF(23)^512 are annihilated by all 425 rows, and
every coordinate is nonzero in at least one vector. All 1275 dot
products are independently verified. This proves more than failure of
full rank: **L contains no nonzero integer multiple of any coordinate**.

For the proof, the 425 rows are independent modulo 23, hence have
rational rank 425. For each j, the kernel witness nonzero at j shows
that appending e_j raises rank modulo 23 to 426. It therefore also raises
rational rank. Thus e_j is not in the rational row span, ruling out
every nonzero multiple d*e_j. This is an exact finite counterexample
to the previously proposed uniform coordinate-inclusion theorem,
including any choice of j, any power of two, or any other nonzero d.

Nevertheless, the actual candidate root gains generate all of GF(2048)^*.
The seven-change family has good colourings, just no coordinate-isolating
proof of their existence. The sufficient condition was too strong.

### An ordinary prefix restores rank

In this same example, the natural (0,1) circuits have lengths
290,156,300,238,40; the (0,2) circuits have lengths 552,402,70.
Swap 0 and 1 on the 290-edge component containing path edge {1,2}.
From the resulting colouring, run the seven-change rule using its
current colour-0 edges as temporary seeds and all nonzero-colour edges
as pivots. These seeds need not be the original physical chord edges.

The new family has 920 returning pivots, 104 blocked pivots, and 666
distinct nonzero root rows. Its rank modulo 23 is **512**. A retained
512-by-512 minor has determinant **10 modulo 23**, independently checked
using ordinary scalar elimination. The generation lemma now forces
a non-unit candidate gain. The stored audit also replays an explicit
eight-change sequence and verifies the resulting lifts at both partners
by matrices and full projective permutations.

This is a rank-restoring prefix, not a minimal colouring route. Even
the prefix alone gives root gain 1733, root word U(633), and point
profile 2037,6,6 at both translations. The original seven-change family
already colours as well. No uniform eight-change theorem is proposed.

### Audit scope and next proof attempt

`python3 code/borel_gain_rank.py --larger --json-output code/borel_gain_rank_results.json`
reproduces the new scan and certificates. The one-prime test passes nine
of the ten stored q=2048 residual cases modulo 23; the example above is
the sole failure. All ten stored q=4096 cases pass modulo 3. These are
stored cases only, not an all-parameter rank audit. Seven smaller
controls include the two q=1024 rank-deficient cases, where this new
test fails but the earlier coordinate certificates work.

The bit-packed modular arithmetic is checked against scalar operations
in 1920 tests and against scalar ranks/determinants on 64 small matrices.
Seven saved root matrices through q=512 receive additional independent
scalar rank/minor checks. The counterexample's kernel vectors are checked
by direct integer dot products modulo 23, independently of elimination.
Its complete signed rows and replayable source seeds are retained.
The prefixed 512-row minor is independently checked without packed
arithmetic. Packed arithmetic is a computational optimisation, not a
mathematical assumption in the proof.

The revised direction is a **move-closed family**: retain root rows from
successive colourings reached by ordinary Kempe changes and the
seven-change rule. The accumulated span can only grow. Full rank at
one prime dividing q-1 suffices; alternatively, one can try to exclude
an annihilating character using the actual projective chord labels.
The key missing theorem is that this larger family must supply such
a certificate. Monotonicity of the accumulated span does not prove
strict growth or eventual full rank. Ordinary moves alone fail in the
old rigid q=256 class; one batch of seven-change moves has the new
q=2048 rank obstruction. Allowing both successively addresses those
two concrete limitations without inventing another universal move bound.

This is a corrected proof direction, not a proof that we are close to
settling the conjecture. No new uniform colouring family or larger
exhaustive group catalogue has been obtained. The 13/16 guarantee,
q=2048 full-cycle coverage, partial q=4096 general scan, single-trace
q=8192 control, and 352,440-vertex global lower bound remain unchanged.

## Follow-up: unrestricted temporary colours, minimum support, and exact closures

The most useful correction in this continuation is that the natural chord
colour is not privileged in the seven-change theorem. Allowing other
temporary ordinary colours removes the three recorded rank deficiencies
without an ordinary prefix. This is a finite improvement of the candidate
family, not a uniform full-rank theorem or a proof of the conjecture.

### Remove the temporary-colour restriction

For each ordinary colour a, use every pivot whose colour is not a, and
both return branches, temporarily removing its two incident colour-a
edges. Equivalently, permute a to colour 0 and apply the existing rule.
Deduplicate the union of all resulting root rows up to sign.

| q | t | c | Union rows | Prime | Rank | Selected minor determinant |
|---|---|---|---|---|---|---|
| 1024 | 49 | 80 | 831 | 3 | 256 | 1 |
| 1024 | 287 | 683 | 780 | 3 | 256 | 1 |
| 2048 | 343 | 1766 | 1403 | 23 | 512 | 21 |

All three minors receive independent scalar-elimination checks. At
q=2048 the separate temporary-colour families have 425, 722, and 886
rows; their successive unions have 425, 1097, and 1403 rows, of ranks
425, 512, and 512 modulo 23. Thus colours 0 and 1 already suffice for
full rank in this example. The previous ordinary prefix remains valid
but is unnecessary for this enlarged rank test. The earlier negative
certificate concerns **only the explicitly chord-seeded family**.

### The restricted root span has minimum support exactly three

For the 425-by-512 matrix R at q=2048,t=343,c=1766, eight saved kernel
vectors over GF(23) give 512 nonzero, pairwise projectively distinct
column signatures. All 3400 row/vector dot products are checked directly.

This rules out every nonzero vector of support one or two in the
modular row span. It also rules them out in the **rational** row span:
scale any hypothetical such vector to a primitive integer vector w.
Its nonzero reduction modulo 23 has support at most two. Projective
distinctness gives a kernel vector not orthogonal to w, so adjoining w
raises the modular rank to 426. This is incompatible with rational
rank 425, which equals the number of input rows. A recorded row has
support three, proving the minimum is exactly three over both fields.

Consequently neither coordinate isolation nor any two-coordinate
relation can work in this restricted rational span. This strengthens
the previous coordinate obstruction; it does not obstruct colouring.
Indeed the first return branch at path pivot {6,7}, excluding colour 0,
has signed row e_6-e_7+e_142 (one-based chord coordinates). Its oriented
endpoint labels are

```
starts: 802, 928, 1838
ends:   1463, 1333, 443
```

In the recorded binary field, their products are 1583 and 106, so the
root gain is 178, not one. The root has 748 folded vertices and root
word U(1971). The partner point profiles are 1497,276,276 at c=1766 and
1497,552 at c=1969. Both lifts and all seven recolouring steps are
independently audited, with at most two temporary-colour edges. Three
counts root-path chords, not path length or spanning-exchange size.

### Rank-neutral states must not be discarded

Let W be the restricted natural row span modulo 23. Each of the nine
ordinary Kempe neighbours of the natural colouring has both root rows
in W. The certificate retains all eighteen coefficient vectors and
verifies their linear combinations directly, not merely against a
sample of null vectors. Nevertheless the 290-edge ordinary prefix,
which is rank-neutral relative to W, has a following colour-0
seven-change neighbourhood of full rank.

Thus a colouring's current root rows do not determine what its next
moves can add. An exhaustive search should deduplicate colouring states
up to global names, not discard them because their rows add no rank.
This is **not** a plateau for the unrestricted all-temporary-colour
family: that family already has full rank at the natural start here.

### Exhaustive small state-space comparison

A new breadth-first search permits every ordinary Kempe component
change and every temporary-colour/pivot/return-branch choice of the
seven-change rule. It never prunes rank-neutral states. An independent
edge-domain backtracking enumeration, pinning the three colours at one
vertex, supplies every three-colouring up to global colour names.

| q | Normalized parameters | All three-colourings | Parameters missed by ordinary-only closure | Parameters failing the all-state rank test |
|---|---|---|---|---|
| 4 | 2 | 2 | 0 | 0 |
| 8 | 6 | 10 | 0 | 2 |
| 16 | 28 | 86 | 6 | 24 |
| 32 | 60 | 472 | 42 | 48 |

The combined moves reach every one of these 570 states on 96 parameters.
These counts retain both c,c+t, use Frobenius-representative full-cycle
traces, and sum states across parameters rather than isomorphism types.
Every state has a non-unit eligible gain; 1140 matrix/permutation lift
checks independently verify its two eligible matchings. None of these
small parameters is a canonical failure, limiting the relevance to the
hard residual case.

In 74 parameters the root-row matrix of **all** three-colourings is
rank-deficient at every prime dividing q-1, despite those successful
gains. Therefore full rank is sufficient, not necessary, even for the
actual projective ordering. This does not disprove a rank theorem
restricted to canonical failures.

The earlier reordered-label q=16 control has exactly one three-colouring
up to global names, with both gains one. The combined moves reach it,
but there is no other three-colouring to help. This strengthens the
warning against replacing the projective ordering by an arbitrary order
of the correct label set. It is not a Cayley counterexample.

### Reproduction and the actual next proof obligation

```
python3 code/borel_gain_closure.py --json-output code/borel_gain_closure_results.json
python3 code/borel_kempe_closure.py --json-output code/borel_kempe_closure_results.json
```

The next argument must address arithmetic nonvanishing, not just
reachability. Under the contrary assumption that all reachable root
gains are one, the actual chord-logarithm vector annihilates every
reachable row modulo q-1. The target is to contradict that assumption
using the projective label recurrence and the unrestricted move rules.
Proving connectivity alone would leave this step unresolved. Full rank
at one prime is one sufficient route, not a condition that must be
forced in every already-colourable parameter.

There is no new uniform colouring family or improved asymptotic/global
bound in this continuation. The full conjecture remains open, and no
reliable estimate of distance to a proof follows from these finite tests.

## Follow-up: joint polynomial elimination and an anchored forbidden pattern

This continuation addresses the arithmetic issue directly. It produces
an exact polynomial reformulation, a uniform forbidden-pattern lemma,
and independent certificates on six stored diagrams. It does not prove
that every obstruction contains a forbidden pattern, and does not extend
the verified parameter range or any global/asymptotic bound.

### A recurrence-based nonvanishing certificate

Let n=q+1 and h=q/2. Regard D_0=0, D_1=1 and
D_(i+1)=T D_i+D_(i-1) as polynomials over GF(2). Let Psi_n be the monic
polynomial whose roots are lambda+lambda^(-1), one per inverse pair of
primitive nth roots of unity. For n=q+1 its roots are exactly the
full-cycle traces in GF(q), all simple.

For a fixed ordered chord diagram and an eligible root path, set

```
A_N(T) = product of D_i(T) over the source ends of its oriented chords
B_N(T) = product of D_j(T) over the destination ends
p_N(T) = A_N(T) + B_N(T).
```

All D_i(t), 1 <= i <= h, are nonzero at full-cycle traces, and
gamma(N)=(B_N(t)/A_N(t))^2. Hence unit gain is exactly p_N(t)=0.
If gcd(Psi_n,p_N : N in the chosen family)=1, some candidate has non-unit
gain at every full-cycle trace on this **fixed ordered diagram**.
This is a sufficient certificate, not a uniform existence theorem.

The trace polynomial has two independently checked constructions:
multiply the minimal polynomials of the full-cycle Frobenius orbits, or
use the exact identity

```
sqrt(D_n(T)) = product of Psi_d(T) for d|n, d>1    (n odd).
```

For the latter, substitute T=lambda+lambda^(-1) into
D_n=(lambda^n+lambda^(-n))/(lambda+lambda^(-1)). Nonidentity nth roots
give all roots, each trace with multiplicity two; monicity and degree
give the displayed polynomial identity.

### Translation equations give an exact realizability test

Fix an arbitrary perfect pairing M of the h ordered vertices, with
reference chord {a,b}. For every other chord {i,j}, put

```
C_ij = (D_i+D_j) D_a D_b + (D_a+D_b) D_i D_j
J_M  = gcd(Psi_n, all C_ij).
```

The roots of J_M are exactly the full-cycle traces realising this
chord diagram, with two translations per trace. Thus 2 deg(J_M) is the
exact number of normalized generating pairs (t,c) with this diagram.

The converse avoids spurious algebraic solutions. At a root of Psi_n all
denominators are nonzero. The equations C_ij=0 say that
(D_i+D_j)/(D_i D_j) has one common value. Squaring gives one nonzero
difference Delta=k_i+k_j for every paired label. Since all labels lie in
Tr(k/t^2)=1, their difference satisfies Tr(Delta/t^2)=0. Consequently
c^2+tc+Delta=0 has exactly two roots in GF(q), neither 0 nor t, and the
resulting actual folded diagram has precisely the pairing M.

For any combinatorially specified eligible matching family N, define

```
K_(M,N) = gcd(J_M, all p_N).
```

Its roots are exactly the realizable full-cycle traces at which all
candidate gains are one. This is an exact elimination criterion. A gcd
of one is verified by polynomial Bezout identities; a nonconstant gcd
identifies failures of the chosen family, not necessarily failures of
Cayley colourability by other methods.

### A one-row certificate despite the q=2048 rank obstruction

For the previously audited path at q=2048,t=343,c=1766, the oriented
chords have index pairs 6 -> 573, 162 -> 7, and 161 -> 850. Its polynomial is

```
p_N = D_6 D_162 D_161 + D_573 D_7 D_850.
```

This has degree 1427 and is coprime to Psi_2049, of degree 682. The
retained Bezout identity is independently checked with coefficient-array
arithmetic. Thus this one row has non-unit gain at all 682 full-cycle
traces on the fixed diagram, despite the restricted lattice containing
no nonzero vector on just one or two coordinates. This uses the actual
recurrence, not full row rank or coordinate isolation.

Only eleven traces actually realise the diagram: the Frobenius orbit of
343. They yield the twenty-two already-covered parameter pairs. Thus the
broader symbolic statement does not add actual Cayley parameter coverage.
Direct field-recurrence checks cover all 62 full-cycle trace orbits, and
the actual two partner lifts are audited.

All six stored controls admit a one-row polynomial certificate:

| q | t | c | Root chords | Actual root gain | deg(Psi_(q+1)) | deg(J_M) |
|---|---|---|---|---|---|---|
| 64 | 7 | 26 | 1 | 55 | 24 | 6 |
| 128 | 8 | 20 | 1 | 45 | 42 | 7 |
| 256 | 15 | 151 | 1 | 105 | 128 | 8 |
| 512 | 47 | 10 | 3 | 5 | 162 | 9 |
| 512 | 53 | 66 | 2 | 294 | 162 | 9 |
| 2048 | 343 | 1766 | 3 | 178 | 682 | 11 |

Each J_M has just the original trace's Frobenius orbit as its roots.
Independent direct-field checks reconstruct every represented parameter
orbit and verify that no other full-cycle trace orbit realises the same
diagram. These successes do not justify a universal one-row bound.

### A forbidden pattern over every characteristic-two field

The anchored chords **{1,3}, {2,7}, {4,6} cannot coexist** with the
projective recurrence and a common additive label difference. This is
not an extrapolation from field-size tests.

Using {1,3} as the reference chord, the other two translation equations
are P(t)=Q(t)=0, where

```
P(T) = T^9 + T^8 + T^7 + T^4 + T^2 + T + 1
Q(T) = T^10 + T^7 + T^6 + T.
```

The exact identity over GF(2) is

```
(T^6+T+1) P(T) + (T^5+T^4+T^3+T^2) Q(T) = 1.
```

It excludes a common root over every characteristic-two extension field.
Clearing the D_i denominators is valid in every full-cycle diagram.
The positions are anchored at D_1=1; an arbitrary shift of this pattern
has not been justified.

### Why both arithmetic ingredients are necessary

Over GF(16) modulo theta^4+theta+1, the actual t=8 projective labels are

```
1, 10, 3, 15, 13, 4, 6, 8.
```

Keep that order but impose the arbitrary pairing
{1,3}, {2,7}, {4,6}, {5,8}. The closed graph has exactly one three-colouring
up to global names, and its two eligible gains are both one. Its only
nonzero root row is (1,1,-1,1). Modulo Psi_17, its unit-gain polynomial is
T^4+T^3+T^2+T+1, which vanishes at t=8.

This graph violates the translation condition: its chord-label differences
are 2,12,11,5. The forbidden-pattern identity proves it cannot arise from
the required recurrence and a common translation in any characteristic-two
field. Thus recurrence alone is insufficient. Conversely, the earlier
reordered-label control retained a common translation while dropping the
recurrence and also failed. Neither control is a Cayley counterexample.

An independent enumeration exhausts all arbitrary pairings at the
full-cycle Frobenius representatives: three diagrams at q=8 have no
recurrence-only failure, while 210 diagrams at q=16 have four. All four
are eliminated when the translation equations are included.

### Reproduction and next proof obligation

```
python3 code/borel_gain_polynomial.py --json-output code/borel_gain_polynomial_results.json
```

The checker includes 128 independent multiplication/division/Bezout test
triples, both trace-polynomial constructions, coefficient-array checks of
the stored gain identities, independent field-recurrence comparisons at
all full-cycle trace representatives, and twelve matrix/permutation lift
audits, one for each selected path at each partner translation. The JSON
retains the fixed-diagram and certificate data.

The useful new target is an algebraic exclusion of realizable gain-one
obstructions: show that failure of the chosen move family forces chord
relations whose translation/recurrence equations are inconsistent. The
displayed pattern is one rigorously excluded configuration, but there is
no theorem making such configurations unavoidable. Equivalently, the
missing statement is K_(M,N)=1 uniformly for the relevant realizable
diagrams. This would settle the corresponding full-cycle family, still
not all simple groups or the full Alspach–Zhang conjecture. The conjecture
remains open; all coverage and counterexample-size bounds are unchanged.

## Strategy diagnostic: gain distribution over all three-colourings (6 September 2026, late)

Before choosing between a counting proof and a growing-size construction
for the residual class, this experiment measures how the eligible root
gains are distributed over **every** proper three-edge-colouring of the
closed folded graph H, not just those reached by a move rule. It is a
diagnostic, not a theorem, and it adds no coverage.

Each colouring has two eligible matchings (the classes avoiding the
closing edge), and each eligible matching has a root path from folded
vertex 1 to h whose signed chord row determines the gain. Many colourings
share one root path: they differ only in the even off-path circuits. So
the statistics are taken at two levels, per eligible matching and per
distinct root path.

```sh
python3 code/borel_gain_distribution.py --fields 64 128 --controls 64 --workers 6 --json-output code/borel_gain_distribution_results.json
```

Scope: the stored canonical-failure diagrams (one at q=64, four at q=128,
each covering the partner pair c, c+t), plus every folded diagram for the
first stored trace as controls (t=7 at q=64, t=8 at q=128). The
enumeration is the independent edge-domain backtracker of the closure
audit; four matchings per diagram are additionally lifted and checked
by matrices and full projective permutations (396 checks). Runtime was
three minutes for 99 diagrams.

| q | Diagrams | Colourings (mod names) | Distinct root paths | Unit-gain root paths | Unit fraction times (q-1) | chi-square per class, root paths |
|---|---|---:|---:|---:|---|---|
| 64 | failure (7,26) | 170 | 141 | 4 | 1.79 | 1.0 |
| 64 | 30 controls, non-pentagon | 21 to 148 each | 37 to 148 each | 68 total | 0.59 to 4.50, median 1.56 | median 1.0, max 1.4 |
| 128 | 4 failures | 6,862 to 20,212 each | 5,331 to 12,938 each | 42 to 97 each | 0.95 to 1.17 | 0.9 to 1.0 |
| 128 | 62 controls, non-pentagon | 3,068 to 20,164 each | 3,134 to 13,884 each | 3,371 total | 0.64 to 1.34, median 0.99 | median 1.0, max 1.3 |

The chi-square column is the chi-square statistic of the gain histogram
against the uniform distribution on GF(q)^*, divided by q-1; a uniformly
random sample gives 1.0. Over distinct root paths the gains are
statistically indistinguishable from uniform in every non-pentagon diagram,
including all five canonical failures. In particular the unit value is
not privileged: about one root path in q-1 has gain one, and no diagram
has zero unit-gain root paths. The pentagon diagram c=1 is the one genuine
outlier (127,816 colourings at q=128, all sharing one root path in one
class), and it is already covered by the pentagon theorem.

Per eligible matching the histograms are far from uniform (chi-square per
class up to 59 at q=128), but only because of hub root paths shared by up
to 3,008 colourings. That multiplicity is combinatorial (the number of
completions of the off-path structure) and carries no arithmetic
information; it disappears at the root-path level.

Two further facts. First, in every failure diagram the natural colouring
is the only colouring known a priori, and both of its root paths have
unit gain; there are just one (q=64) or two (q=128, two of the four
diagrams) other both-unit colourings among thousands. Second, the number
of root paths grows roughly by a factor 90 when h doubles from 32 to 64,
so at q=128 more than 99% of eligible matchings colour the Cayley graph.

**Conclusion for strategy.** The supply of good matchings is not scarce;
it is overwhelming, and the gains look like independent uniform values on
the multiplicative group. The obstruction to a proof is therefore not
existence of a good matching but our inability to *exhibit* one without
computing in the field. This favours a counting formulation over further
move rules: the sufficient statement is a nontrivial cancellation bound
for the character sums

```
S(chi) = sum over root paths P of chi(gamma(P)),   chi a nontrivial character of GF(q)^*,
```

since the number of unit-gain root paths is (1/(q-1)) sum over all chi of
S(chi), and S(1) is the total. The observed chi-square of 1.0 means the
average of |S(chi)|^2 over nontrivial chi equals the number of root paths:
square-root cancellation on average. Even a much weaker bound,
|S(chi)| < S(1) for every nontrivial chi, would prove the residual class,
and bounds of that type are what the recurrence and translation structure
should be aimed at. Growing-size constructions such as the seven-change
rule remain valid search tools, but the data give no reason to expect
that any specific rule is forced to succeed; the two-buffer successes are
consistent with picking essentially random root paths from a set in which
99% succeed.

Limitations: exhaustive enumeration is limited to q<=128 by the colouring
count (an attempted q=256 run was stopped because storing the colouring
set would have exhausted memory; a streaming counter or sampler would be
needed there). The controls cover a single trace per field. No estimate
for S(chi) has been proved; the conclusion is that this is the right
target, not that it is within reach. The full conjecture and all bounds
are unchanged.
