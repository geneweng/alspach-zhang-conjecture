#!/usr/bin/env python3
"""Fill the @PLACEHOLDERS@ in survey.tex and docs/index.html from the result files.
Run from the code/ directory after make_table.py."""
import json, os, re, sys
sys.path.insert(0, os.path.dirname(__file__))
import summarize_results as S   # runs the summary (prints), exposes rows/admissible

rows = S.rows
adm = S.admissible
missing = [(n, o) for n, o in adm if n not in rows]
bound = missing[0][1]
ngroups = len(rows)
maxorder = max(r["order"] for r in rows.values())
npairs = sum(r["pairs_tested"] for r in rows.values())
maxq = max(r.get("max_quotient_vertices", 0) for r in rows.values())

def tex_num(n):
    s = f"{n:,}".replace(",", "\\,")
    return s

# --- theorem list
simple_orders = [(n, o) for n, o in S.SIMPLE]
verified_simple_prefix = 0
for n, o in simple_orders:
    if n in rows:
        verified_simple_prefix = o
    else:
        break
qs_psl = sorted(int(m.group(1)) for m in (re.match(r"PSL\(2,(\d+)\)$", n) for n in rows) if m)
qs_pgl = sorted(int(m.group(1)) for m in (re.match(r"PGL\(2,(\d+)\)$", n) for n in rows) if m)
an = sorted(int(m.group(1)) for m in (re.match(r"A(\d+)$", n) for n in rows) if m)
sn = sorted(int(m.group(1)) for m in (re.match(r"S(\d+)$", n) for n in rows) if m)
def contiguous_psl(qs, start, odd_only=False):
    # report "q <= Q" for the largest prime power Q such that all prime powers start <= q <= Q are present
    pp = [q for q in range(start, 200) if all(q % p for p in range(2, q)) or q in (4, 8, 9, 16, 25, 27, 32, 49, 64, 81, 121, 125, 128)]
    if odd_only:
        pp = [q for q in pp if q % 2 == 1]
    Q = 0
    for q in pp:
        if q in qs: Q = q
        else: break
    extra = [q for q in qs if q > Q]
    return Q, extra
# PSL(2,4) = PSL(2,5) = A5, PSL(2,9) = A6, PGL(2,5) = S5 appear under their alternating/symmetric names
if "A5" in rows: qs_psl += [4, 5]
if "A6" in rows: qs_psl += [9]
if "S5" in rows: qs_pgl += [5]
Qpsl, extra_psl = contiguous_psl(sorted(qs_psl), 4)
Qpgl, extra_pgl = contiguous_psl(sorted(q for q in qs_pgl if q % 2 == 1), 5, odd_only=True)
others = [n for n in rows if not re.match(r"(PSL\(2,\d+\)|PGL\(2,\d+\)|A\d+|S\d+)$", n)]
others_sorted = sorted(others, key=lambda n: rows[n]["order"])

from make_table import name as tname
def fmt_names(names):
    return ", ".join(tname(n) for n in names)

theoremlist = (f"every non-abelian simple group of order at most ${tex_num(verified_simple_prefix)}$, "
               f"for $\\PSL(2,q)$ with $q\\le {Qpsl}$" + (f" and $q\\in\\{{{', '.join(map(str, extra_psl))}\\}}$" if extra_psl else "") +
               f", for $\\PGL(2,q)$ with $q\\le {Qpgl}$" + (f" and $q\\in\\{{{', '.join(map(str, extra_pgl))}\\}}$" if extra_pgl else "") +
               f", for $A_n$ with $n\\le {max(an)}$ and $S_n$ with $n\\le {max(sn)}$, and for the groups " + fmt_names(others_sorted))

firstunverified = ", ".join(f"{tname(n)} (order ${tex_num(o)}$)" for n, o in missing[:5])

subs = {"@NGROUPS@": str(ngroups), "@MAXORDER@": f"${tex_num(maxorder)}$", "@BOUND@": tex_num(bound),
        "@BOUNDM1@": f"${tex_num(bound - 1)}$", "@NPAIRS@": f"${tex_num(npairs)}$", "@MAXQ@": f"${tex_num(maxq)}$",
        "@FIRSTUNVERIFIED@": firstunverified, "@THEOREMLIST@": theoremlist}
root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
tex = open(os.path.join(root, "survey.tex")).read()
# @BOUND@ is used both inside and outside math mode; handle "$<@BOUND@$" and "less than @BOUND@"
tex = tex.replace("$<@BOUND@$", f"$<{tex_num(bound)}$").replace("less than @BOUND@", f"less than ${tex_num(bound)}$").replace("at least @BOUND@", f"at least ${tex_num(bound)}$")
for k, v in subs.items():
    tex = tex.replace(k, v)
assert "@" not in re.sub(r"\\url\{[^}]*\}|@\{\}", "", tex), [l for l in tex.split("\n") if "@" in l and "\\url" not in l and "@{}" not in l][:5]
open(os.path.join(root, "survey.tex"), "w").write(tex)

# --- project page
html_path = os.path.join(root, "docs", "index.html")
html = open(html_path).read()
def html_name(n):
    t = tname(n)
    t = t.replace("$", "").replace("\\mathrm{", "").replace("\\wr", " ≀ ").replace("\\mathbb{Z}_2", "Z<sub>2</sub>")
    t = re.sub(r"_\{(\d+)\}", r"<sub>\1</sub>", t)
    t = re.sub(r"_(\d)", r"<sub>\1</sub>", t)
    t = t.replace("\\Sigma ", "Σ").replace("}", "").replace("\\", "")
    return t
hsubs = {"@NGROUPS@": str(ngroups), "@MAXORDER@": f"{maxorder:,}".replace(",", " "), "@BOUND@": f"{bound:,}".replace(",", " "),
         "@BOUNDM1@": f"{bound-1:,}".replace(",", " "), "@NPAIRS@": f"{npairs:,}".replace(",", " "), "@MAXQ@": str(maxq),
         "@FIRSTUNVERIFIED@": ", ".join(html_name(n) for n, o in missing[:6]),
         "@QPSL@": str(Qpsl), "@QPGL@": str(Qpgl), "@AN@": str(max(an)), "@SN@": str(max(sn)),
         "@OTHERS@": ", ".join(html_name(n) for n in others_sorted)}
for k, v in hsubs.items():
    html = html.replace(k, v)
assert "@" not in re.sub(r"https?://[^\s\"<]*|mailto:[^\"<]*|[\w.]+@[\w.]+|@(media|import|font-face|keyframes|supports)", "", html), "unfilled placeholder in index.html"
open(html_path, "w").write(html)
print("filled:", subs)
