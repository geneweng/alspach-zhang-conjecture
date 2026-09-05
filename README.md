# The Alspach–Zhang conjecture: a survey

**Conjecture (Alspach–Zhang, 1992).** Every connected cubic Cayley graph is 3-edge-colourable.
Equivalently (Mader + Jaeger + Tutte): every Cayley graph of valency at least 2 admits a
nowhere-zero 4-flow. In the terminology of Nedela and Škoviera, there are no *Cayley snarks*.

This repository contains a survey of the conjecture and a small computational verification.
The project's aim is to settle the conjecture: prove it, or construct a Cayley snark.

**Project page:** https://geneweng.github.io/alspach-zhang-conjecture/ (survey as a web page and PDF, status, plan).

## Contents

| Path | What it is |
|---|---|
| `survey.tex`, `survey.pdf` | The survey (LaTeX source and compiled PDF). |
| `docs/` | GitHub Pages site: project page (`index.html`), the survey as HTML and PDF. |
| `results_table.tex` | Table of computational results, generated from the JSON files below. |
| `code/cayley_snark_check2.py` | **Current checker.** Exhaustive 3-edge-colourability check for all cubic Cayley graphs on a catalogue of non-solvable permutation groups (simple groups and their index-2 extensions up to order about 260 000: PSL/PGL(2,q), A_n, S_n, PSL(3,4) and its three index-2 extensions, PSU(4,2) = PSp(4,3) and PSU(4,2).2, Sz(8), PSU(3,4), PSU(3,4).2, PSU(3,5), M11, M12, M12.2, J1, PSL(2,7) ≀ Z2, A6 ≀ Z2, ...). Reduces each Cayley graph to a quotient pregraph by a large subgroup avoiding the conjugacy class of x (semi-edges allowed), solves the quotient with CaDiCaL, lifts and verifies the colouring on the full graph. |
| `code/cayley_snark_check.py` | First checker (quotients by odd-order abelian subgroups only); also provides the group constructors (finite fields, PSL/PGL(2,q), A_n, S_n, M10, M11, PSL(3,3), A5 ≀ Z2) used by the current one. |
| `code/cdc_palette_experiment.py` | Reproduces the eight-point palette analysis of the 2026 cycle-double-cover proof: exact affine-system enumeration for K4, Petersen, and a cubic Cayley graph on S5, with exact palette chromatic numbers, a direct SAT check of the full-rank tetrahedral reformulation, a separated transvection-matching witness beyond all parallel matching certificates, weighted-cut-code and restriction ranks, and an unbiased sample of generator-separated flows. |
| `code/layer_potential_experiment.py` | Tests the symmetric quotient bilinear form against the original CDC equations and searches its lift family for a large layer-potential space. |
| `code/transvection_switch_experiment.py` | Checks the affine palette-layer equation, the order-five quotient-switch table, an all-lifts obstruction, and a direct colouring construction from filtered quotient flows. |
| `code/colour_support_experiment.py` | Exact signed-cycle test for a prescribed quotient colour support: all 243 local words, a comparison with restricted switches on S5, translated alternating-cycle repairs, independent quotient-matching starts, and all 32 Petersen supports. |
| `code/translated_repair_audit.py` | Strict descent audit with independently checked flip certificates, complete SAT/backtracking agreement on all 10 500 order-five A5 supports, subgroup-invariant and successive-repair tests on larger groups, stabiliser divisibility checks, and a transitive non-Cayley stall control. |
| `code/involution_repair_experiment.py` | Checks the involution-exchange lemma, verifies paired-circuit descent on all obstructed A5 supports, and certifies an explicit matching whose every decreasing flip loses stabiliser symmetry. |
| `code/translated_block_repair.py` | Two explicit PSL(2,11) counterexamples to restricted repair rules, an independent all-translates stall certificate, the exact implication criterion for valid simultaneous flips, and successive closed-set repairs on larger samples. |
| `code/quotient_face_experiment.py` | Checks the minimum-a-edge matching face and quotient/lift exchanges, independently exhausts restricted invariant matching spaces, tests lexicographic minima and unrestricted PSL(2,11) descent samples, rules out a two-consecutive-matching shortcut on A5, and certifies positive mean oddness drift for two natural random repair rules. |
| `code/loop_congruence_experiment.py` | Checks the unsigned/Tait count congruence, exhausts negative circuit-subset orbits, verifies sharp moduli and the nonfree Petersen control, and tests nonzero residues on complete small faces and invariant subfamilies. |
| `code/transposition_quotient_experiment.py` | Checks the infinite transposition–odd-cycle construction on two-subset quotients, enumerates all invariant minimum-face matchings through S12, verifies cycle monodromies, and independently checks all 21 lifts on S6. |
| `code/two_cycle_monodromy_experiment.py` | Checks the canonical circuit-free quotient matching for transpositions joining two coprime odd cycles, and exhausts small quotient matching spaces including genuine bad-word controls. |
| `code/odd_transposition_monodromy_experiment.py` | Checks the explicit one-cycle chord construction for every nonadjacent coprime chord, including the exceptional distances 2 and 3. |
| `code/sparse_involution_cycle_experiment.py` | Exhaustively certifies the odd full-cycle/sparse-involution theorem: point- and two-set-stabilizer quotients for involutions with three or five transpositions and for seven transpositions on fifteen points, with exact Schreier checks of the proper-group exceptions. |
| `code/parity_interlacement.py` | Computes oddness both as the number of odd-degree vertices in the defect-incidence multigraph and as a difference of two binary interlacement nullities, using a two-sheeted cover of the four-regular quotient minus a matching. |
| `code/extra_groups.py` | The same check for PSU(3,3), G2(2) = PSU(3,3).2 and PSL(3,3).2 (first-checker method). |
| `code/atlas/` | Permutation generators (MeatAxe text format) for Sz(8), M12, M12.2, J1, J2, M22, PSL(3,5), PSL(3,7), copied from the ATLAS of Finite Group Representations (brauer.maths.qmul.ac.uk/Atlas). |
| `code/make_table.py` | Builds `results_table.tex` from the JSON result files. |
| `code/results2.json`, `code/results.json`, `code/results_extra.json` | Raw results (one record per group; `results2.json` is the current run and includes a per-pair record of the subgroup used). |
| `code/run2.log`, `code/run.log`, `code/run_extra.log` | Logs of the runs reported in the survey. |

## What the survey covers

1. Origin and equivalent formulations (flows vs. colourings; reduction to the cubic case).
2. The trivial cases, and the reduction to `Cay(G, {a, x, x⁻¹})` with `a` an involution and `x` of odd order.
3. Elementary reformulations: normal subgroups containing `a` (always colourable), small girth,
   generators of order 3 ⇔ arc-regular cubic graphs, index-2 normal subgroups ⇔ bi-Cayley graphs.
   The index-2 case is proved when `x` has order 3; for general odd order it is reduced to finding
   a consecutive 2-factor in a cyclically ordered bipartite coset-incidence graph.
4. Solvable groups (Alspach–Liu–Zhang 1996) and solvable vertex-transitive groups (Potočnik 2004).
5. The structure of a minimal counterexample (Nedela–Škoviera 2001), with a self-contained proof.
6. The hamiltonicity route: (2,s,3)-Cayley graphs (Glover, Kutnar, Malnič, Marušič) and Cayley maps.
7. Computational evidence: the census of cubic vertex-transitive graphs up to 1280 vertices
   (Potočnik–Spiga–Verret) and the new computations in `code/`.
8. The 2026 cycle-double-cover proof, an exact eight-point palette reformulation, and its present obstruction.
9. Related flow problems and open questions.

## Status and next steps (5 September 2026)

**Latest result: full odd cycles paired with sparse involutions.** Let N be odd,
let x be an N-cycle in S_N, and let a be a product of r disjoint transpositions,
with <a,x> = S_N. The survey proves that Cay(S_N,{a,x,x^-1}) is 3-edge-colourable
whenever N >= 4r+1, for every possible N when r <= 5, and also for (r,N)=(7,15). The large range follows
because a then has two consecutive fixed points. Exact, directly verified quotient
certificates handle the remaining cases r=3 and r=5. They use point stabilizers where
possible and switch to two-set stabilizers for the parity-blocked and Petersen-like
cases. The only exceptional types not covered by those certificates are a degree-9 pair
generating a solvable group of order 162 and cycle reflections generating the dihedral
groups of orders 14 and 22, so none generates the required symmetric group.

For r=3 the checker verifies 11, 85 and 350 dihedral orbit types in degrees 7, 9
and 11. For r=5 it verifies the nonconsecutive-fixed-point types in degrees
11,13,15,17,19: respectively 513, 5832, 12152, 5832 and 513 types. Of these,
46 degree-13 types and 512 degree-11 types require the two-set quotient.
For (r,N)=(7,15), the two-set quotient colours 68144 of 68219 dihedral types;
exact Schreier recursion shows that all 75 failures generate proper subgroups, of
order at most 466560, rather than S15.

**Previous result: the transposition case for every symmetric group.** If N >= 3,
a is a transposition, x in S_N has odd order, and <a,x> = S_N, the survey now proves
that Cay(S_N,{a,x,x^-1}) is 3-edge-colourable. In even degree, x necessarily has two
coprime odd cycles; in odd degree it is one N-cycle and the endpoints of a have coprime
cyclic distance. Both cases use the quotient by C_G(a), whose vertices are two-subsets.
Explicit matchings make the complementary quotient a union of paths capped by semi-edges,
apart from one harmless even quotient circuit in the smallest exceptional case.

The two-cycle checker verifies all 33 coprime pairs through total degree 30. The one-cycle
checker verifies all 246 nonadjacent coprime chord cases through degree 51 (270 including
the adjacent cases covered by the Hamiltonicity theorem); every constructed quotient
is circuit-free except (N,k)=(5,2), where its only circuit has even length two. Complete
full lifts independently check the three exceptional S5 and S7 matchings. Small
two-cycle enumerations also expose the selection issue: for cycle lengths (3,5),
8 of 111 invariant quotient matchings have a bad 9-cycle of odd monodromy order 5, while
the explicit construction avoids all quotient circuits.

**Previous result: a uniform monodromy construction.** For every odd s >= 5, take
G = S_(s+1), x an s-cycle, and a the transposition joining its fixed point to a point
of the cycle. The survey proves that every C_G(a)-invariant matching in the
minimum-a-edge face has an even complement. This class is nonempty and has exactly
(s-2)^((s-1)/2) + (s-1)(s-2)^((s-3)/2) members. The proof controls the quotient cycle
words: possible odd quotient pentagons have monodromy of even order s-1, so they lift
to even cycles. Some of these matchings admit no colouring with all colours invariant
under C_G(a), so this goes beyond colouring the subgroup quotient directly.

The independent quotient checker exhausts 21, 275, 5145 and 124659 matchings for
s = 5, 7, 9 and 11, and checks all 21 S6 lifts directly against a separate exact-cover
enumerator. That earlier theorem covers one generating pair per symmetric group of even
degree, whereas the new result above covers every pair whose involution is a
transposition; no claim of priority over the Hamiltonicity literature is made.
The Alspach–Zhang conjecture remains open. The next symmetric-group target is a full
odd cycle with seven transpositions in degrees 17 through 27; the broader target is to control
odd quotient-circuit words for other almost simple groups.

**Previous result: a congruence using the regular action.** For any invariant family of
perfect matchings, let U sum 2 to the number of complementary circuits, ignoring their
parities, and let T count the corresponding labelled Tait colourings. The survey proves
U = T modulo 4|G|_2 when 4 divides |G| (modulo 4 when |G|_2 = 2).
When 4 divides |G|, negative terms in the circuit-subset expansion have odd-order stabilisers
under translation and complementation, which forces the divisibility. A nonzero unsigned residue therefore
certifies a colouring. Exact checks on both A5 and S5 faces give residues 12 modulo 16
and 24 modulo 32. The full order-80 affine face has zero residue, but each of its six
40-matching orbits with two complementary circuits gives 160 = 32 modulo 64 and hence
a certificate with no Hamiltonian complement. Equivalently, if the largest power of two
dividing a matching's stabiliser order is 2^r, at most r+1 complementary circuits certify
a colouring when 4 divides |G|. Each S5 case has 63 certifying non-Hamiltonian matching
orbits. The order-50 wreath-product face has no certifying orbit, so even selecting an
invariant subfamily cannot make this test work there. Existence of a matching satisfying
the bound in the remaining non-solvable cases is unproved; the original conjecture remains open.
The extended PSL(2,11) check finds 240 certificates among 4336 centraliser-invariant
matchings, including 224 non-Hamiltonian ones. In contrast, all 21 such S6 matchings
have at least 16 complementary circuits and fail the bound of five; the monodromy
construction proves their evenness instead.

**Latest direction: an exchange-closed minimum-cost face.** Restrict to matchings selecting
exactly one a-edge at each x-cycle. Whenever the contracted graph has a perfect matching
(proved here in the non-solvable order-five case), this class is nonempty. The survey now proves
that it is the minimum-a-edge face of the perfect-matching polytope, linearly equivalent to the
quotient's matching polytope. All translated circuit flips stay in this face and correspond
exactly to quotient circuit exchanges. Thus no local-validity obstruction occurs there.
A two-sheeted construction and Traldi's extended Cohn–Lempel equality give an exact
parity-refined binary-nullity formula, independently checked against direct cycle counts.
The missing step is to force an even complement somewhere in this face; the formula alone
does not do so. An explicit restricted PSL(2,11) matching also refutes uniform negative drift
for both random single flips and random subsets of difference circuits, although it repairs.

**Earlier obstruction: the universal single-circuit repair rule is false.** On a 660-vertex
PSL(2,11) Cayley graph, an explicit locally valid matching has four odd complementary cycles,
but all 168 distinct single-circuit flips create an all-a pentagon. An independent graph
reconstruction and traversal verify the obstruction. Flipping two circuits together colours
the graph, so this is **not a counterexample to Alspach–Zhang**.
A second matching rules out requiring an inclusion-minimal valid block to improve:
two individually valid, non-improving flips jointly give a colouring.
The survey proves that valid simultaneous flips are exactly implication-closed sets of
difference circuits. Existence of a parity-improving closed set remains unproved.

**Done.** 85 non-solvable groups verified (25 954 generating pairs, all 3-edge-colourable): every group
admissible for a smallest Cayley snark of order below 352 440, i.e. a smallest Cayley snark has more than
352 439 vertices. The survey has a new Section 4 (quotients by arbitrary subgroups with semi-edges, a
cyclewise transition criterion, a Z_3 reformulation, and the Petersen-minus-a-vertex obstruction), proves
the index-2 case for generators of order 3, and gives a consecutive-2-factor criterion for the remaining
index-2 cases. It also analyses the 2026 proof of the cycle double cover theorem: a 4-colourable graph on
its eight flow labels would pull back to a Tait colouring, but random full-rank flows produce dense palette
graphs and a fully Cayley-invariant characteristic-two flow is impossible. Full rank alone is still an
exact reformulation of Tait colourability. A generator-separated flow reduces the missing step to a
four-by-four Hall condition; twelve of its twenty-four matching cases can be imposed by augmenting the
paper's affine system with further linear equations. For a fixed flow, all solution freedom beyond global
translation is exactly a weighted cut code. The successful S5 transvection witness has code dimension five,
and its 40 dangerous-edge equations reduce to a rank-four syndrome after global translations are removed.
Its smallest layer potential comes from a consecutive monochromatic 4-cycle in the contracted graph; the
survey proves that every such quotient cycle can be converted into a nonzero code direction by choosing the lifts.
More generally, it gives an exact interlacement-parity obstruction to extending any quotient subflow and
classifies all 480 local configurations when the odd generator has order five.
The layer-potential space is the kernel of an explicit symmetric quotient form, modulo constants.
Searching this matrix family produces six directions on S5, but none of the twelve transvection targets
succeeds for the resulting flow: dimension alone is insufficient. Local triangle orders are only changes
of coordinates and cannot alter feasibility.
The attainable palette layers are the characteristic vectors of the same symmetric form.
For order five, transvection switches give a direct affine colouring test on the quotient.
A local repeated-value obstruction excludes all 2^48 lifts of the previously unsuccessful S5 quotient flow.
Retaining only one colour's quotient-edge support gives a broader test: a disjoint union of signed cycles
must have even sign sum on each component. On twelve locally filtered S5 quotient flows, the restricted
switch test succeeds once and the support test eight times; all constructed colourings are verified.
The twelve flows have eight distinct supports. All three failing distinct supports can be repaired by
one alternating-cycle flip against a group translate of the forced matching. More importantly, a locally
valid starting support always exists in the non-solvable order-five case: take a perfect matching of the
quotient, whose existence follows from edge connectivity and Tutte's theorem. Of twelve sampled quotient
matchings on S5, nine colour directly and the other three are repaired by one translated flip.
An exhaustive check of all 32 Petersen supports fails, as expected. These finite single-flip
successes cannot be extrapolated: the PSL(2,11) certificate above refutes the universal rule.
The repair rule is now verified on every locally valid support for both order-five A5 cases:
two independent enumerators agree on all 10 500 supports, and every one of the 4090 obstructed states
has a decreasing flip. Thus any such A5 start colours in at most two flips. Successive repairs also
pass the documented A6 and S6 tests. A structural lemma uses the free Cayley action: if a subgroup H
preserves a matching, its complementary odd-circuit count is divisible by the largest power of two
dividing |H|. This restricts possible obstructions but does not yet force the count to vanish.
An involution-exchange lemma now identifies exactly when a matching and its translate can be
combined into an invariant matching: no difference circuit may be fixed by the involution;
any such fixed circuit has length 2 modulo 4. Flipping a paired circuit increases agreement
with that translate. All 4090 obstructed A5 supports admit a decreasing flip of this restricted
kind, as do the obstructed states in new 512-matching samples on S5, A6 and S6.
An explicit A5 example rules out monotone stabiliser growth as a universal repair strategy:
every valid decreasing single flip changes four odd cycles to two and stabiliser order 12 to 3.
For every involution outside its stabiliser, the union with its translate has no invariant
perfect matching. This is a barrier to stronger symmetry strategies, not a stalled descent.
Solvable-group stress tests also pass: all 1480 locally valid supports on C5 wreath C2
are enumerated, and a 30 000-matching sample is checked on the 80-element affine group.
The more flexible closed-set repair reaches colourings from all 561 obstructed starts among
6144 sampled centralizer-invariant matchings on twelve specified PSL(2,11) generating pairs,
using at most nine decreasing steps. These are bounded samples, not an exhaustive theorem.
Inside the smaller minimum-a-edge face, two independent enumerators exhaust all 4336
centralizer-invariant matchings on those twelve pairs: 3906 colour directly and all 430
obstructed starts colour in at most three single flips. Intermediate states may lose symmetry
but remain in the face. The face has 125 matchings in each A5 case, 120 on the wreath-product
example and 705 on the affine example; all 1075 are checked against the parity-nullity formula.
On 6144 additional unrestricted PSL(2,11) face samples, all 3563 obstructed starts have an
immediate strict single-circuit repair. Every minimum-circuit factor is even in each of the
twelve exact invariant spaces, although the wreath-product face shows that minimum circuit
count does not force evenness in general. The quotient is naturally an orientably regular map,
and the exact flow target is a nowhere-zero F2^2 flow whose multiplicity-three value is
consecutive at every pentavalent vertex; the face target additionally fixes one singleton value
globally. Requiring both singleton values globally would reduce the problem to two consecutive
disjoint quotient perfect matchings, but exact enumeration shows that this stronger shortcut
already fails on A5. A proper three-colouring of the dual map is sufficient, and a local
four-colouring version succeeds on all four completely enumerated small examples. This
face-potential route is not universal: in four PSL(2,11) cases with ord(ax)=11, the 60-vertex
11-valent dual has chromatic number five, checked by a SAT refutation for four colours and a
verified five-colouring. Every quotient flow is instead a face-boundary flow plus a class in
the first homology of the regular map; those four genus-70 examples force a nonzero class.
The fixed-colour minimum-face count is now an explicit signed-loop state sum: a quotient
perfect matching selects one gap at every pentavalent vertex, the other four darts form
signed smoothed loops, and a loop contributes 2 or 0 according as its sign is even or odd.
A Fourier expansion makes this a fixed three-state local vertex model and grouping states
gives an exact G-orbit formula. The values on the two A5 maps, the order-50 wreath product,
and the order-80 affine map are respectively 540, 540, 440 and 8320; optional exhaustive
S5 runs give 115960 on each of two maps. The new congruence compares this coefficient with
the unsigned count modulo a power of two. It supplies a sufficient non-vanishing test,
but a universal nonzero-residue theorem is missing, and the Petersen control has value zero.
No Hamiltonian complement exists in the face for the latter two graphs, so simply minimising
the number of complementary circuits is not an adequate general substitute for controlling parity.
No complete proof and no Cayley snark.

**Next computational steps**, in order of the bound they give (each is about an hour on 12 cores):

```
cd code
WORKERS=11 VERBOSE=1 RESULTS=results2g.json python3 cayley_snark_check2.py 400000 "PSL(2,89)"   # order 352 440 -> bound 362 879 (S9 next)
python3 summarize_results.py          # prints the new bound and the next unverified admissible groups
python3 make_table.py results2g.json results2f.json results2e.json results2d.json results2c.json results2.json results2b.json results2_j1.json results.json results_extra.json
```

Then update the literal numbers in `survey.tex` (abstract, Table 1 row, Section 7.2 text, Theorem 7.x, Corollary 7.x)
and in `docs/index.html` (the "Verified" row, "What a counterexample must look like" item 3, Track B), rebuild
(`pdflatex survey` twice, `python3 code/build_html.py`, `cp survey.pdf docs/`), commit and push. After PSL(2,89)
the next admissible groups are S9 (362 880, in the catalogue), PGL(2,71), PGL(2,73), PSL(2,97), PSL(2,101), M22 (443 520;
ATLAS generators are in `code/atlas/`), PGL(2,79), PSL(2,103), PSL(2,107), PSL(2,109), PSL(2,113), PSL(2,121), PSL(2,125);
groups above about 400 000 elements need more memory per worker (use `WORKERS=6`).

**Theory.** Open directions recorded in the survey (Sections 9–10): the involution-class quotient, the
Z_3 reformulation for ord(x) = 5, adaptive choice of the quotient subgroup (the natural point
quotient can be the Petersen snark even when the Cayley graph is colourable), and forcing a
4-colourable eight-point palette in the cycle-double-cover construction without assuming a rank-two flow.
For ord(x) = 5, the most direct current target is to choose a quotient colour support whose signed
circuits all have even parity. A quotient perfect matching supplies a start in the nonempty,
exchange-closed minimum-a-edge face. The new target is to force zero defect in the two-sheeted
interlacement-nullity formula while staying in that face. Equivalently, prove non-vanishing of
the signed-loop state sum's lowest a-edge coefficient using the regular G-action. The
unsigned-count congruence now gives a sufficient arithmetic criterion: choose an invariant
family whose unsigned residue is nonzero. This is impossible within the order-50 face;
existence in the remaining non-solvable cases is an open sufficient target.
Ordinary Penrose/transition identities do not supply it.
Uniform averaging does not guarantee
negative oddness drift, and forcing a Hamiltonian complement is too strong within this class.
For more general supports, odd-sign circuits come in pairs. A single translated alternating-cycle flip works on every obstructed
locally valid A5 support, but the universal rule fails on PSL(2,11). The current repair target
allows an implication-closed union of difference circuits. Local validity is characterised
exactly; proving that some closed union decreases the odd-circuit count remains open.
Even restricting to inclusion-minimal valid blocks fails, so interactions between distinct
valid flips matter. The stabiliser-divisibility lemma provides a constraint that genuinely uses
regularity; it fails for the transitive, nonregular Petersen action. This is the even-2-factor
obstruction in quotient coordinates, not a proof by itself. The involution-exchange lemma describes
paired difference circuits, but they cannot suffice universally. Increased local agreement must not be
confused with growth of the full stabiliser: the explicit A5 example forces symmetry loss at
every first decreasing flip. Within the CDC approach, the remaining
problem is to hit a transvection syndrome after first avoiding the new quotient-level obstructions;
lift changes preserving the symmetric form preserve the entire affine set of layer assignments.

## Reproducing the computation

```
cd code
python3 -m pip install python-sat
RESULTS=results2.json python3 cayley_snark_check2.py 270000   # all catalogue groups of order <= 270000
python3 cayley_snark_check2.py 270000 J1 A9                    # or just the named groups
python3 extra_groups.py                  # PSU(3,3), PSL(3,3).2, G2(2), writes results_extra.json
python3 cdc_palette_experiment.py         # exact eight-point palette diagnostic
python3 layer_potential_experiment.py     # quotient form, direct checks, and lift search
python3 transvection_switch_experiment.py # affine layers and direct order-five switches
python3 colour_support_experiment.py      # exact signed-cycle support criterion
python3 quotient_face_experiment.py --s5-state-sum # signed-loop sum and optional exhaustive S5 values
python3 loop_congruence_experiment.py --s5 # unsigned residues and every negative-state orbit
python3 loop_congruence_experiment.py --psl11-centralizers # 240 certificates in 12 invariant classes
python3 transposition_quotient_experiment.py --max-s 11 # quotient monodromies and independent S6 lifts
python3 two_cycle_monodromy_experiment.py --max-degree 10 # all small matchings and canonical constructions
python3 odd_transposition_monodromy_experiment.py --max-N 51 # every coprime one-cycle chord through S51
python3 sparse_involution_cycle_experiment.py # sparse involutions paired with full odd cycles
python3 translated_repair_audit.py --samples 100000 --groups A5 A5_alt --independent
python3 translated_repair_audit.py --samples 1024 --groups A6 S6 --invariant three --descend
python3 translated_repair_audit.py --samples 2000 --groups A6 --invariant centralizer --descend
python3 translated_repair_audit.py --samples 30000 --groups F80
python3 translated_repair_audit.py --samples 2000 --groups W50 --independent --descend
python3 involution_repair_experiment.py --exhaustive-a5
python3 involution_repair_experiment.py --groups S5 A6 S6 --samples 512
python3 translated_block_repair.py        # exact local check and both PSL(2,11) barriers
python3 translated_block_repair.py --stress-psl11 512 --descend
python3 quotient_face_experiment.py        # closure, two enumerators, parity matrices, drift barrier
python3 quotient_face_experiment.py --psl11-centralizers
python3 quotient_face_experiment.py --psl11-unrestricted 512
python3 make_table.py results2.json results.json results_extra.json
cd .. && pdflatex survey && pdflatex survey && python3 code/build_html.py
```

For each group `G` the program enumerates all pairs `(a, x)` with `a` an involution, `x` of odd
order and `<a, x> = G`, up to conjugation and inversion of `x` (the other cubic Cayley graphs on
`G` are trivially 3-edge-colourable). For each pair it takes the quotient of the Cayley graph by
the left action of a subgroup `H` containing no conjugate of `x` (point and set stabilisers,
centralisers, normalisers, odd-order abelian subgroups, ... — largest first). The quotient is a
cubic pregraph (semi-edges where `H` contains conjugates of `a`) covered by the Cayley graph, so a
3-edge-colouring of the quotient (found by CaDiCaL through python-sat) lifts to the Cayley graph;
the lifted colouring is verified independently on the full graph. If no quotient is colourable the
full graph is passed to the SAT solver with a time limit; this fallback was never needed. Group
orders are checked against the known values and the Petersen graph serves as a negative control.
