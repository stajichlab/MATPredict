#!/usr/bin/env python3
"""The polish tie-break fix: (a) both panels against the earlier default-33 runs (only reference-record fields may differ; calls identical),
(b) 16 Leppa1 repeats (record chosen). usage: verify_tiebreak.py (in this folder)"""
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
for name, new, old in (("Basidiomycota panel (34)", "basidio", R + "2026-10-08_supported_span_default33/basidio"),
                       ("Ascomycota + Mucoromycota (110)", "other", R + "2026-10-08_supported_span_default33/other")):
    N, O = load(new), load(old)
    common = [g for g in N if g in O]
    bad_calls = [g for g in common if sig(N[g]) != sig(O[g])]
    c = Counter(); changed = []
    for g in common:
        d = diff(strip(N[g]), strip(O[g]))
        for p in d: c[re.sub(r"\[\d+\]", "[]", p)] += 1
        if d: changed.append(g)
    print(f"\n{name}: {len(N)} new reports, {len(common)} compared with the earlier default-33 run")
    print(f"   calls, cluster spans, genes, merged_from, withheld loci differ in {len(bad_calls)} genomes {bad_calls[:3]}")
    print(f"   genomes with any difference in the report: {len(changed)}; fields: {dict(c.most_common(8))}")
    for g in changed[:6]:
        print("      ", g[:36], [(re.sub(r"\[\d+\]", "[]", p)) for p in diff(strip(N[g]), strip(O[g]))][:3])
print("\nLeppa1 repeats (gene[2] record):")
tally = Counter()
for d in sorted(glob.glob("leppa_run*")):
    f = glob.glob(d + "/runs/*/detection_report.yaml")
    if not f: tally["no report"] += 1; continue
    x = (yaml.safe_load(open(f[0])).get("detected") or [None])[0]
    tally[x["gene_evidence"][2]["reference_record"].split("_")[1] if x else "no call"] += 1
print("  ", dict(tally), "of", sum(tally.values()))
