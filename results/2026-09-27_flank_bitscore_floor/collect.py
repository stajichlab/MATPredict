"""Collect every flank-carried locus (no modelled core gene; flanks modelled) per run.
Pre-rule runs: detected calls whose gene_evidence has no modelled core gene but a modelled flank.
Rule runs: detected with idiomorph_unmodelled true + suppressed_loci withheld as flank_carried*.
Writes loci.tsv."""
import csv, glob, os, sys, yaml, tarfile, io
from yaml import CSafeLoader as L
R = "/bigdata/stajichlab/jstajich/projects/MATPredict/results"
FL = {"SLA2","APN2","COX13","sla2","apn2","cox13","PAP1","OBP1","PIK1","tptA","rnhA","glrA","algA","btbA",
      "MIP1","beta_fg","RibL6","BAP31","CAF1","STE20"}
def modelled(e): return (e.get("status") or "").startswith("polished") or e.get("method") == "diamond_proteome"
RUNS = [  # (label, glob of report dirs, has_rule)
 ("mucoro_4457c3a", f"{R}/2026-09-27_btbA_homothallic/Mucoromycota_4457c3a/runs/*", True),
 ("serinales_882aa01", f"{R}/2026-09-26_serinales_all_882aa01/runs/*", False),
 ("basidio_full_ad1f865", f"{R}/2026-09-26_basidiomycota_full/*/runs/*", True),
 ("dothideo_recheck_9932c7c", f"{R}/2026-09-27_dothideo_recheck/run*/runs/*", True),
 ("early_div_634dda4", f"{R}/2026-09-26_early_diverging/*/runs/*", False),
]
meta = {r["ASMID"]: r for r in csv.DictReader(open("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"))}
out = csv.writer(open("loci.tsv","w"), delimiter="\t")
out.writerow(["run","genome","phylum","order","species","family","contig","start","end","outcome","core_genes","rundir"])
n = 0
for lab, pat, rule in RUNS:
    dirs = sorted(glob.glob(pat))
    if lab.startswith("dothideo") and not dirs: dirs = sorted(glob.glob(f"{R}/2026-09-27_dothideo_recheck/**/runs/*", recursive=True))
    for d in dirs:
        p = f"{d}/detection_report.yaml"
        if not os.path.exists(p): continue
        g = os.path.basename(d); m = meta.get(g, {})
        r = yaml.load(open(p), Loader=L) or {}
        rows = []
        for x in r.get("detected") or []:
            ev = x.get("gene_evidence") or []
            core = [e for e in ev if e["gene"] not in FL]; fl = [e for e in ev if e["gene"] in FL]
            fc = (x.get("idiomorph_unmodelled") is True) if rule else (not any(modelled(e) for e in core) and any(modelled(e) for e in fl))
            if fc: rows.append((x["family"], x["contig"], x["start"], x["end"], "kept_low" if rule else "prerule_call",
                                ",".join(sorted({e["gene"] for e in core}))))
        if rule:
            for x in r.get("suppressed_loci") or []:
                if "flank_carried" in str(x.get("withheld_reason")):
                    rows.append((x["family"], x["contig"], x["start"], x["end"], "withheld",
                                 ",".join(sorted(set(x.get("genes_found", [])) - FL))))
        for fam, c, s, e, o, cg in rows:
            out.writerow([lab, g, m.get("PHYLUM",""), m.get("ORDER",""), m.get("SPECIES",""), fam, c, s, e, o, cg, d]); n += 1
print("loci", n)
