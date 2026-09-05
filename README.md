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
An exhaustive check of all 32 Petersen supports fails, as expected. A general decreasing repair move
has not been proved; the finite repair successes do not establish the conjecture.
The repair rule is now verified on every locally valid support for both order-five A5 cases:
two independent enumerators agree on all 10 500 supports, and every one of the 4090 obstructed states
has a decreasing flip. Thus any such A5 start colours in at most two flips. Successive repairs also
pass the documented A6 and S6 tests. A structural lemma uses the free Cayley action: if a subgroup H
preserves a matching, its complementary odd-circuit count is divisible by the largest power of two
dividing |H|. This restricts possible obstructions but does not yet force the count to vanish.
No complete proof and no counterexample.

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
circuits all have even parity. A quotient perfect matching guarantees a locally valid start. Odd-sign
circuits come in pairs; flipping an alternating cycle between the forced matching and a group translate
is now verified to give a decreasing move on every obstructed locally valid A5 support, but no general decreasing
move has been proved. The new stabiliser-divisibility lemma provides a constraint that genuinely uses
regularity; it fails for the transitive, nonregular Petersen action. This is the even-2-factor
obstruction in quotient coordinates, not a proof by itself. Within the CDC approach, the remaining
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
python3 translated_repair_audit.py --samples 100000 --groups A5 A5_alt --independent
python3 translated_repair_audit.py --samples 1024 --groups A6 S6 --invariant three --descend
python3 translated_repair_audit.py --samples 2000 --groups A6 --invariant centralizer --descend
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
