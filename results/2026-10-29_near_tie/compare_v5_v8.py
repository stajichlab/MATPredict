"""Before (v5: per-locus rule with a 3-point near-tie margin) and after (v8: the family decided once per genome by mean identity,
`family_ambiguous` on every unscoped cross-call) on the 354 yeast genomes. usage: python3 compare_v5_v8.py <yeast_routing dir>"""
import collections, glob, os, sys
import yaml

R = sys.argv[1]

def load(run):
    out = {}
    for d in sorted(glob.glob(f"{R}/{run}/runs/*")):
        f = d + "/detection_report.yaml"
        if os.path.exists(f):
            out[os.path.basename(d)] = yaml.safe_load(open(f)).get("detected") or []
    return out

v5, v8 = load("out_v5"), load("out_v8")
fam = lambda loci: sorted({l["family"].split(":")[1] for l in loci})
print("genomes", len(v5), len(v8))
mixed = lambda d: sum(1 for loci in d.values() if len(fam([l for l in loci if l.get("cross_call")])) > 1)
print("genomes whose cross-called loci keep two different families: v5", mixed(v5), "v8", mixed(v8))
amb = [(g, l) for g, loci in v8.items() for l in loci if l.get("family_ambiguous")]
print("v8 loci with family_ambiguous:", len(amb), "in", len({g for g, _ in amb}), "genomes; kept family", collections.Counter(l["family"].split(":")[1] for g, l in amb))
lower = sum(1 for g, l in amb if (l["family_ambiguous"]["gap"] or 0) < 0)
print("kept call has the lower identity:", lower)
gaps = [l["family_ambiguous"]["gap"] for g, l in amb if l["family_ambiguous"]["gap"] is not None]
print("gap (kept minus other), min/median/max:", min(gaps), sorted(gaps)[len(gaps) // 2], max(gaps))
st = collections.Counter((l.get("family_status") or {}).get("status", "") + ":" + (l.get("family_status") or {}).get("reason", "") for loci in v8.values() for l in loci)
print("family_status v8:", dict(st))
st5 = collections.Counter((l.get("family_status") or {}).get("status", "") + ":" + (l.get("family_status") or {}).get("reason", "") for loci in v5.values() for l in loci)
print("family_status v5:", dict(st5))
changed = []
for g in v5:
    a = sorted((l["contig"], l["start"], l["family"].split(":")[1], l.get("idiomorph")) for l in v5[g])
    b = sorted((l["contig"], l["start"], l["family"].split(":")[1], l.get("idiomorph")) for l in v8.get(g, []))
    if a != b: changed.append((g, a, b))
print("genomes whose calls differ:", len(changed))
for g, a, b in changed[:40]:
    print(" ", g, [x[2:] for x in a], "=>", [x[2:] for x in b])
