#!/usr/bin/env python3
"""Summarise the JSON result files: which groups are verified, and the smallest
group admissible for a minimal Cayley snark (Nedela-Skoviera) that is NOT yet
verified.  The list of admissible groups (non-abelian simple groups, and T.2 /
S wr Z_2 with T, S non-abelian simple) is enumerated from the known orders of the
simple groups up to 300 000."""
import json, sys, os

files = sys.argv[1:] or ["results2f.json", "results2e.json", "results2d.json", "results2c.json", "results2.json", "results2b.json", "results2_j1.json", "results.json", "results_extra.json"]
rows = {}
for f in files:
    if os.path.exists(f):
        for r in json.load(open(f)):
            if r["undecided"] == 0 and r["non_colourable"] == 0:
                rows.setdefault(r["group"], r)

# non-abelian simple groups of order <= 300000: (name, order, number of index-2 extensions
# with involutions outside the socle is not needed -- we list the index-2 extensions by name)
SIMPLE = [("A5", 60), ("PSL(2,7)", 168), ("A6", 360), ("PSL(2,8)", 504), ("PSL(2,11)", 660),
          ("PSL(2,13)", 1092), ("PSL(2,17)", 2448), ("A7", 2520), ("PSL(2,19)", 3420),
          ("PSL(2,16)", 4080), ("PSL(3,3)", 5616), ("PSU(3,3)", 6048), ("PSL(2,23)", 6072),
          ("PSL(2,25)", 7800), ("M11", 7920), ("PSL(2,27)", 9828), ("PSL(2,29)", 12180),
          ("PSL(2,31)", 14880), ("A8", 20160), ("PSL(3,4)", 20160), ("PSL(2,37)", 25308),
          ("PSU(4,2)", 25920), ("Sz(8)", 29120), ("PSL(2,32)", 32736), ("PSL(2,41)", 34440),
          ("PSL(2,43)", 39732), ("PSL(2,47)", 51888), ("PSL(2,49)", 58800), ("PSU(3,4)", 62400),
          ("PSL(2,53)", 74412), ("M12", 95040), ("PSL(2,59)", 102660), ("PSL(2,61)", 113460),
          ("PSU(3,5)", 126000), ("PSL(2,67)", 150348), ("J1", 175560), ("PSL(2,71)", 178920),
          ("A9", 181440), ("PSL(2,73)", 194472), ("PSL(2,79)", 246480), ("PSL(2,64)", 262080),
          ("PSL(2,81)", 265680), ("PSL(2,83)", 285852)]
# index-2 extensions T.2 (all classes of such extensions; groups with no involution outside
# the socle, e.g. M10, PSL(2,25).2_3, still have to be listed since they are admissible a priori)
EXT2 = [("S5", 120), ("PGL(2,7)", 336), ("S6", 720), ("PGL(2,9)", 720), ("M10", 720),
        ("PGL(2,11)", 1320), ("PGL(2,13)", 2184), ("PGL(2,17)", 4896), ("S7", 5040),
        ("PGL(2,19)", 6840), ("PSigmaL(2,16)", 8160), ("PSL(3,3).2", 11232), ("G2(2)", 12096),
        ("PGL(2,23)", 12144), ("PGL(2,25)", 15600), ("PSigmaL(2,25)", 15600), ("PSL(2,25).2_3", 15600),
        ("PGL(2,27)", 19656), ("PGL(2,29)", 24360), ("PGL(2,31)", 29760), ("S8", 40320),
        ("PSL(3,4).2_gf", 40320), ("PSL(3,4).2_f", 40320), ("PSL(3,4).2_g", 40320),
        ("PGL(2,37)", 50616), ("PSU(4,2).2", 51840), ("PGL(2,41)", 68880), ("PGL(2,43)", 79464),
        ("PGL(2,47)", 103776), ("PGL(2,49)", 117600), ("PSigmaL(2,49)", 117600), ("PSL(2,49).2_3", 117600),
        ("PSU(3,4).2", 124800), ("PGL(2,53)", 148824), ("M12.2", 190080), ("PGL(2,59)", 205320),
        ("PGL(2,61)", 226920), ("PSU(3,5).2", 252000), ("PGL(2,67)", 300696)]
# S wr Z2
WR = [("A5wrZ2", 7200), ("PSL(2,7)wrZ2", 56448), ("A6wrZ2", 259200)]
# (PSL(2,8) has no index-2 extension: Out = 3; Sz(8): Out = 3; PSL(2,32): Out = 5; PSL(2,64): Out = 6 gives
#  PSL(2,64).2 of order 524160, beyond the range; PSL(2,81): Out = 2x4, three extensions of order 531360.)

admissible = sorted(SIMPLE + EXT2 + WR, key=lambda t: (t[1], t[0]))
missing = [(n, o) for n, o in admissible if n not in rows]
print("verified groups:", len(rows), "  largest order:", max(r["order"] for r in rows.values()))
print("admissible groups listed:", len(admissible))
print("first unverified admissible groups:", missing[:6])
if missing:
    bound = missing[0][1]
    print(f"=> a smallest Cayley snark has more than {bound - 1} vertices "
          f"(every admissible group of order < {bound} is verified)")
print("verified but not in the admissible list (extra coverage):",
      sorted((n, r["order"]) for n, r in rows.items() if n not in dict(admissible)))
print("total generating pairs checked:", sum(r["pairs_tested"] for r in rows.values()))
print("largest quotient solved:", max(r.get("max_quotient_vertices", 0) for r in rows.values()))
