#!/usr/bin/env python3
"""Identity-first tie-break (run-supported6) against (a) the order-dependent code (default33 runs) and (b) the lower-id-only tie-break
(tiebreak_fix runs). Calls must be identical everywhere; every other difference is listed with what it is.
usage: verify_identity_tiebreak.py (in this folder)"""
import glob, re, yaml
from collections import Counter
R = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/"
def load(d): return {f.split("/")[-2]: yaml.safe_load(open(f)) for f in glob.glob(f"{d}/runs/*/detection_report.yaml")}
def strip(r): r = dict(r); r.pop("run", None); return r
def diff(a, b, path=""):
    out = []
    if isinstance(a, dict) and isinstance(b, dict):
        for k in set(a) | set(b):
            if k not in a or k not in b: out.append(path + "/" + k + " (one side only)")
            else: out += diff(a[k], b[k], path + "/" + k)
    elif isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
        for i, (x, y) in enumerate(zip(a, b)): out += diff(x, y, path + f"[{i}]")
    elif a != b: out.append(path)
    return out
def callkey(x): return (x["family"], x["contig"], x["start"], x["end"], x["confidence"], x["idiomorph"], x.get("locus_class"),
                        tuple(sorted(x.get("genes_found") or [])), tuple(sorted(str(m) for m in (x.get("merged_from") or []))))
def sig(r): return ([callkey(x) for x in r.get("detected") or []],
                    sorted((s["family"], s["contig"], s["start"], s["end"], s.get("withheld_reason"), tuple(s.get("genes_found") or [])) for s in r.get("suppressed_loci") or []))
SETS = (("Basidiomycota panel (34)", "basidio", "2026-10-08_supported_span_default33/basidio", "2026-10-08_tiebreak_fix/basidio"),
        ("Ascomycota + Mucoromycota (110)", "other", "2026-10-08_supported_span_default33/other", "2026-10-08_tiebreak_fix/other"))
for name, new, ref_order, ref_id in SETS:
    N = load(new)
    print(f"\n{name}: {len(N)} reports")
    for lab, ref in (("order-dependent code (default33 run)", ref_order), ("lower-id-only tie-break", ref_id)):
        O = load(R + ref); common = [g for g in N if g in O]
        bad = [g for g in common if sig(N[g]) != sig(O[g])]
        c = Counter(); chg = []
        for g in common:
            d = diff(strip(N[g]), strip(O[g]))
            for p in d: c[re.sub(r"\[\d+\]", "[]", p)] += 1
            if d: chg.append(g)
        print(f"   vs {lab}: calls/spans/genes/merged/withheld differ in {len(bad)} genomes; any difference in {len(chg)} genomes {[g[:28] for g in chg]}")
        if c: print("      fields:", dict(c.most_common(6)))
print("\nLeppa1 repeats (gene[2] record):")
t = Counter()
for d in sorted(glob.glob("leppa_run*")):
    f = glob.glob(d + "/runs/*/detection_report.yaml")
    if not f: t["no report"] += 1; continue
    x = (yaml.safe_load(open(f[0])).get("detected") or [None])[0]
    t[x["gene_evidence"][2]["reference_record"].split("_")[1] if x else "no call"] += 1
print("  ", dict(t), "of", sum(t.values()))
