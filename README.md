# The Alspach–Zhang conjecture: a survey

**Conjecture (Alspach–Zhang, 1992).** Every connected cubic Cayley graph is 3-edge-colourable.
Equivalently (Mader + Jaeger + Tutte): every Cayley graph of valency at least 2 admits a
nowhere-zero 4-flow. In the terminology of Nedela and Škoviera, there are no *Cayley snarks*.

This repository contains a survey of the conjecture and a small computational verification.

## Contents

| Path | What it is |
|---|---|
| `survey.tex`, `survey.pdf` | The survey (LaTeX source and compiled PDF). |
| `results_table.tex` | Table of computational results, generated from the JSON files below. |
| `code/cayley_snark_check.py` | Exhaustive 3-edge-colourability check for all cubic Cayley graphs on a list of non-solvable permutation groups (PSL/PGL(2,q), alternating and symmetric groups, M10, M11, PSL(3,3), PΣL(2,16), the index-2 extensions of PSL(2,25), A5 ≀ Z2). |
| `code/extra_groups.py` | The same check for PSU(3,3), G2(2) = PSU(3,3).2 and PSL(3,3).2. |
| `code/make_table.py` | Builds `results_table.tex` from `code/results.json` and `code/results_extra.json`. |
| `code/results.json`, `code/results_extra.json` | Raw results (one record per group). |
| `code/run.log`, `code/run_extra.log` | Logs of the runs reported in the survey. |

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

## Reproducing the computation

```
cd code
python3 -m pip install python-sat
python3 cayley_snark_check.py 60000      # groups of order up to 60000, writes results.json
python3 extra_groups.py                  # PSU(3,3), PSL(3,3).2, G2(2), writes results_extra.json
python3 make_table.py results.json results_extra.json
cd .. && pdflatex survey && pdflatex survey
```

For each group `G` the program enumerates all pairs `(a, x)` with `a` an involution, `x` of odd
order and `<a, x> = G`, up to conjugation and inversion of `x` (the other cubic Cayley graphs on
`G` are trivially 3-edge-colourable). For each pair it takes the quotient of the Cayley graph by
the left action of an abelian subgroup `H` of odd order containing no element of order `ord(x)`
(largest such `H` first); the quotient is a cubic multigraph covered by the Cayley graph, so a
3-edge-colouring of the quotient (found by CaDiCaL through python-sat, using the "perfect matching
whose complement is bipartite" formulation) lifts to the Cayley graph, and the lifted colouring is
verified independently. If no quotient is colourable the full graph is passed to the SAT solver
with a time limit; that fallback was never needed in the reported runs. Group orders are checked
against the known formulas and the Petersen graph serves as a negative control.
