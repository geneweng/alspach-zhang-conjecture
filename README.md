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
4. Solvable groups (Alspach–Liu–Zhang 1996) and solvable vertex-transitive groups (Potočnik 2004).
5. The structure of a minimal counterexample (Nedela–Škoviera 2001), with a self-contained proof.
6. The hamiltonicity route: (2,s,3)-Cayley graphs (Glover, Kutnar, Malnič, Marušič) and Cayley maps.
7. Computational evidence: the census of cubic vertex-transitive graphs up to 1280 vertices
   (Potočnik–Spiga–Verret) and the new computations in `code/`.
8. Related flow problems and open questions.

## Status and next steps (4 September 2026)

**Done.** 85 non-solvable groups verified (25 954 generating pairs, all 3-edge-colourable): every group
admissible for a smallest Cayley snark of order below 352 440, i.e. a smallest Cayley snark has more than
352 439 vertices. The survey has a new Section 4 (quotients by arbitrary subgroups with semi-edges, parity
lemma, Z_3 reformulation, Hamiltonian quotients, the Petersen-minus-a-vertex obstruction). No proof, no
counterexample.

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

**Theory.** Open directions recorded in the survey (Section 9): the involution-class quotient, the
Z_3 reformulation for ord(x) = 5, and adaptive choice of the quotient subgroup (the natural point
quotient can be the Petersen snark even when the Cayley graph is colourable).

## Reproducing the computation

```
cd code
python3 -m pip install python-sat
RESULTS=results2.json python3 cayley_snark_check2.py 270000   # all catalogue groups of order <= 270000
python3 cayley_snark_check2.py 270000 J1 A9                    # or just the named groups
python3 extra_groups.py                  # PSU(3,3), PSL(3,3).2, G2(2), writes results_extra.json
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
