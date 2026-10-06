"""Audit flank-carried calls (no modelled MTL core gene) in the Serinales-wide scan.
Splits them by (a) whether the genome also has a core-modelled call and (b) whether
the core hit lies inside the PAP1/OBP1/PIK1 span +-3 kb. usage: flank_carried_audit.py RUN_DIR"""
import collections, csv, os, sys, yaml
from yaml import CSafeLoader as L
T = sys.argv[1]
meta = {r["ASMID"]: r for r in csv.DictReader(open("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"))}
FL = {"PAP1", "OBP1", "PIK1"}; CORE = {"MTLA1", "MTLA2", "MTLalpha1", "MTLalpha2"}
def modelled(e): return (e.get("status") or "").startswith("polished") or e.get("method") == "diamond_proteome"
cat = collections.Counter(); sp = collections.defaultdict(collections.Counter); total = 0
for g in sorted(os.listdir(T)):
    p = f"{T}/{g}/detection_report.yaml"
    if not os.path.exists(p): continue
    det = (yaml.load(open(p), Loader=L) or {}).get("detected") or []; total += len(det)
    has_core = any(any(e["gene"] in CORE and modelled(e) for e in x.get("gene_evidence") or []) for x in det)
    for x in det:
        ev = x.get("gene_evidence") or []
        if any(e["gene"] in CORE and modelled(e) for e in ev) or not any(e["gene"] in FL and modelled(e) for e in ev): continue
        f = [e for e in ev if e["gene"] in FL]; lo, hi = min(e["start"] for e in f), max(e["end"] for e in f)
        inside = all(lo - 3000 <= e["start"] and e["end"] <= hi + 3000 for e in ev if e["gene"] in CORE)
        k = ("also core-modelled call" if has_core else "only locus", "inside flank span" if inside else "outside flank span")
        cat[k] += 1; sp[k][meta.get(g, {}).get("SPECIES", "?")] += 1
print(f"calls {total}; flank-carried {sum(cat.values())}")
for k, v in cat.most_common():
    print(f"{v:4d}  {k[0]:24s} {k[1]:20s} {dict(sp[k].most_common(8))}")
