#!/usr/bin/env python3
"""All genomes: v0.7.0 definitive run vs frozen cleanrun-2026-10b. Writes diff lists per campaign."""
import glob, os, sys, collections, multiprocessing as mp
import yaml
L = getattr(yaml, "CSafeLoader", yaml.SafeLoader)
R = "/bigdata/stajichlab/jstajich/projects/MATPredict/results"
OUT = f"{R}/2026-10-26_v070_submit/compare"
os.makedirs(OUT, exist_ok=True)
PAT = {
 "asco": (f"{R}/2026-10-26_v070_ascomycota/wave_*/small_*/runs/*/detection_report.yaml", f"{R}/2026-10-23_ascomycota_cleanrun/wave_*/small_*/runs/*/detection_report.yaml"),
 "basidio": (f"{R}/2026-10-26_v070_basidiomycota/runs/*/runs/*/detection_report.yaml", f"{R}/2026-10-23_basidiomycota_cleanrun/runs/*/runs/*/detection_report.yaml"),
 "mucoro": (f"{R}/2026-10-26_v070_mucoromycotina/runs/*/detection_report.yaml", f"{R}/2026-10-24_mucoromycotina_cleanrun/runs/*/detection_report.yaml"),
}
def call_sig(d):
    return sorted((l.get("family"), l.get("contig"), l.get("start"), l.get("end"), l.get("idiomorph"), l.get("confidence"),
                   tuple(sorted(l.get("genes_found") or [])), l.get("polished_genes"), l.get("locus_class"))
                  for l in (d.get("detected") or []))
def model_sig(d):
    return sorted((l.get("contig"), l.get("start"), e.get("gene"), e.get("status"), e.get("start"), e.get("end"), e.get("strand"),
                   tuple((x["start"], x["end"]) for x in (e.get("exons") or [])))
                  for l in (d.get("detected") or []) for e in (l.get("gene_evidence") or []))
def work(a):
    g, new, old = a
    n = yaml.load(open(new), Loader=L) or {}; o = yaml.load(open(old), Loader=L) or {}
    if call_sig(n) != call_sig(o): return g, "call"
    if model_sig(n) != model_sig(o): return g, "model"
    return g, "same"
if __name__ == "__main__":
    for camp, (pn, po) in PAT.items():
        new = {p.split("/")[-2]: p for p in glob.glob(pn)}; old = {p.split("/")[-2]: p for p in glob.glob(po)}
        both = sorted(set(new) & set(old))
        with mp.Pool(14) as pool:
            res = pool.map(work, [(g, new[g], old[g]) for g in both], chunksize=50)
        c = collections.Counter(r for _, r in res)
        print(f"[{camp}] new {len(new)} old {len(old)} both {len(both)} only_new {len(set(new)-set(old))} only_old {len(set(old)-set(new))}: {dict(c)}", flush=True)
        for kind in ("call", "model"):
            open(f"{OUT}/diff_{camp}_{kind}.txt", "w").write("".join(g + "\n" for g, r in res if r == kind))
