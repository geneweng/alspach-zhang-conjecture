#!/usr/bin/env python3
"""Turn results.json (written by cayley_snark_check.py) into results_table.tex."""
import json, sys, os

NAMES = {
    "A5": r"$A_5$", "S5": r"$S_5$", "A6": r"$A_6$", "S6": r"$S_6$", "A7": r"$A_7$",
    "S7": r"$S_7$", "A8": r"$A_8$", "S8": r"$S_8$", "M10": r"$M_{10}$", "M11": r"$M_{11}$",
    "A5wrZ2": r"$A_5\wr\mathbb{Z}_2$", "PSL(3,3)": r"$\mathrm{PSL}(3,3)$",
    "PSL(3,3).2": r"$\mathrm{PSL}(3,3).2$", "PSU(3,3)": r"$\mathrm{PSU}(3,3)$",
    "G2(2)": r"$G_2(2)=\mathrm{PSU}(3,3).2$",
    "PSigmaL(2,16)": r"$\mathrm{P\Sigma L}(2,16)$", "PSigmaL(2,25)": r"$\mathrm{P\Sigma L}(2,25)$",
    "PSL(2,25).2_3": r"$\mathrm{PSL}(2,25).2_3$",
}
def name(n):
    if n in NAMES: return NAMES[n]
    if n.startswith("PSL(2,"): return r"$\mathrm{PSL}(2,%s)$" % n[6:-1]
    if n.startswith("PGL(2,"): return r"$\mathrm{PGL}(2,%s)$" % n[6:-1]
    return n

files = sys.argv[1:] or ["results.json"]
rows = []
for f in files:
    if os.path.exists(f):
        rows += json.load(open(f))
rows.sort(key=lambda r: (r["order"], r["group"]))
seen = set()
out = [r"\begin{tabular}{@{}lrrrrl@{}}", r"\toprule",
       r"group $G$ & $|G|$ & inv.\ classes & pairs $(a,x)$ & orders of $x$ & result\\", r"\midrule"]
for r in rows:
    if r["group"] in seen: continue
    seen.add(r["group"])
    if r["pairs_tested"] == 0:
        res = "no such pairs"
    elif r["non_colourable"] == 0 and r["undecided"] == 0:
        res = "all colourable"
    else:
        res = f"{r['non_colourable']} non-col., {r['undecided']} undecided"
    orders = ", ".join(str(o) for o in r["orders_x"]) or "--"
    out.append(r"%s & %d & %d & %d & %s & %s\\" % (name(r["group"]), r["order"],
               r["involution_classes"], r["pairs_tested"], orders, res))
out += [r"\bottomrule", r"\end{tabular}"]
open("../results_table.tex", "w").write("\n".join(out) + "\n")
print("\n".join(out))
