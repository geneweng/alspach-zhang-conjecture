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
    "PSL(2,49).2_3": r"$\mathrm{PSL}(2,49).2_3$", "PSigmaL(2,49)": r"$\mathrm{P\Sigma L}(2,49)$",
    "PSL(3,4)": r"$\mathrm{PSL}(3,4)$", "PSL(3,4).2_f": r"$\mathrm{PSL}(3,4).2$ (field)",
    "PSL(3,4).2_g": r"$\mathrm{PSL}(3,4).2$ (graph)", "PSL(3,4).2_gf": r"$\mathrm{PSL}(3,4).2$ (graph-field)",
    "PSU(4,2)": r"$\mathrm{PSU}(4,2)=\mathrm{PSp}(4,3)$", "PSU(4,2).2": r"$\mathrm{PSU}(4,2).2=W(E_6)$",
    "Sz(8)": r"$Sz(8)$", "PSU(3,4)": r"$\mathrm{PSU}(3,4)$", "PSU(3,4).2": r"$\mathrm{PSU}(3,4).2$",
    "PSU(3,5)": r"$\mathrm{PSU}(3,5)$", "PSU(3,5).2": r"$\mathrm{PSU}(3,5).2$",
    "M12": r"$M_{12}$", "M12.2": r"$M_{12}.2$", "J1": r"$J_1$", "J2": r"$J_2$", "J2.2": r"$J_2.2$",
    "PSL(2,7)wrZ2": r"$\mathrm{PSL}(2,7)\wr\mathbb{Z}_2$", "A6wrZ2": r"$A_6\wr\mathbb{Z}_2$",
    "A9": r"$A_9$", "S9": r"$S_9$", "A10": r"$A_{10}$", "M22": r"$M_{22}$", "M22.2": r"$M_{22}.2$",
    "PSL(3,5)": r"$\mathrm{PSL}(3,5)$", "PSL(2,64)": r"$\mathrm{PSL}(2,64)$",
}
def name(n):
    if n in NAMES: return NAMES[n]
    if n.startswith("PSL(2,"): return r"$\mathrm{PSL}(2,%s)$" % n[6:-1]
    if n.startswith("PGL(2,"): return r"$\mathrm{PGL}(2,%s)$" % n[6:-1]
    return n

def main():
  files = sys.argv[1:] or ["results.json"]
  rows = []
  for f in files:
      if os.path.exists(f):
          rows += json.load(open(f))
  # a complete record (no undecided pairs) wins over an incomplete one for the same group
  rows.sort(key=lambda r: (r["order"], r["group"], r.get("undecided", 0)))
  seen = set()
  out = [r"\begin{tabular}{@{}lrrrrrl@{}}", r"\toprule",
         r"group $G$ & $|G|$ & inv.\ classes & pairs $(a,x)$ & orders of $x$ & $\max|\bar V|$ & result\\", r"\midrule"]
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
      maxq = r.get("max_quotient_vertices")
      if maxq is None:
          maxq = max((r["order"] // o for o in r.get("quotient_orders", [])), default=r["order"]) if r["pairs_tested"] else 0
      out.append(r"%s & %d & %d & %d & %s & %s & %s\\" % (name(r["group"]), r["order"],
                 r["involution_classes"], r["pairs_tested"], orders, (str(maxq) if maxq else "--"), res))
  out += [r"\bottomrule", r"\end{tabular}"]
  open("../results_table.tex", "w").write("\n".join(out) + "\n")
  print("\n".join(out))

if __name__ == "__main__":
  main()
