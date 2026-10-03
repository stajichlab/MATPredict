"""Collect every reported call with its routing, confidence and best core identity.

Writes calls.tsv (one row per detected locus) from the runs listed in RUNS.
best_core_any      = best identity among gene_evidence rows with role core_MAT
best_core_modelled = same, restricted to modelled rows (status polished_* / exonerate / miniprot)
"""
import csv, glob, os, sys
import yaml

R = "/bigdata/stajichlab/jstajich/projects/MATPredict/results"
RUNS = {
    "basidio_full": f"{R}/2026-09-26_basidiomycota_full/*/runs/*/detection_report.yaml",
    "basidio_capoff": f"{R}/2026-09-26_basidio_capoff_test/capoff_*/runs/*/detection_report.yaml",
    "pilots_0924": f"{R}/2026-09-24_pilots/*/runs/*/detection_report.yaml",
    "dothideo_f81dad1": f"{R}/2026-09-26_dothideo_curation/f81dad1/runs/*/detection_report.yaml",
    "dothideo_dca6ccc": f"{R}/2026-09-26_dothideo_curation/dca6ccc/runs/*/detection_report.yaml",
    "dothideo_d18ec5a": f"{R}/2026-09-26_dothideo_curation2/d18ec5a/runs/*/detection_report.yaml",
    "dothideo_52cd292": f"{R}/2026-09-26_dothideo_curation2/52cd292/runs/*/detection_report.yaml",
    "serinales_882aa01": f"{R}/2026-09-26_serinales_all_882aa01/runs/*/detection_report.yaml",
    "polishcap_cap6": f"{R}/2026-09-26_polish_cap/cap6/*/runs/*/detection_report.yaml",
    "early_diverging": f"{R}/2026-09-26_early_diverging/*/runs/*/detection_report.yaml",
}
MODELLED = ("polished", "exonerate", "miniprot", "rescued")

samples = {}
for r in csv.DictReader(open("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv")):
    samples[r["ASMID"]] = r

out = []
for run, pat in RUNS.items():
    files = glob.glob(pat)
    if not files:  # some runs nest one level less
        files = glob.glob(pat.replace("/*/runs/", "/runs/"))
    for f in files:
        genome = f.split("/")[-2]
        try:
            rep = yaml.safe_load(open(f))
        except Exception:
            continue
        s = samples.get(genome, {})
        for d in rep.get("detected") or []:
            core = [g for g in d.get("gene_evidence") or [] if g.get("role") == "core_MAT"]
            ids_any = [g["identity"] for g in core if g.get("identity") is not None]
            ids_mod = [g["identity"] for g in core if g.get("identity") is not None
                       and any(k in str(g.get("status", "")) for k in MODELLED)]
            bc = max((g for g in core if g.get("identity") is not None), key=lambda g: g["identity"], default=None)
            out.append(dict(best_core_gene=bc["gene"] if bc else "", core_contig=bc.get("contig", "") if bc else "",
                            core_start=bc.get("start", "") if bc else "", core_end=bc.get("end", "") if bc else "",
                            run=run, genome=genome, order=s.get("ORDER", ""), cls=s.get("CLASS", ""),
                            phylum=s.get("PHYLUM", ""), species=s.get("SPECIES", ""),
                            routing=rep.get("routing_mode", ""), family=d.get("family", ""),
                            contig=d.get("contig", ""), start=d.get("start", ""), end=d.get("end", ""),
                            confidence=d.get("confidence", ""), locus_class=d.get("locus_class", ""),
                            idiomorph=d.get("idiomorph", ""),
                            best_core_any=max(ids_any) if ids_any else "",
                            best_core_modelled=max(ids_mod) if ids_mod else "",
                            core_genes=";".join(f"{g['gene']}:{g.get('identity')}:{g.get('status')}" for g in core),
                            genes_found=";".join(d.get("genes_found") or [])))
    print(run, len(files), "reports", file=sys.stderr)
with open("calls.tsv", "w", newline="") as fo:
    w = csv.DictWriter(fo, fieldnames=list(out[0]), delimiter="\t")
    w.writeheader(); w.writerows(out)
print(len(out), "calls", file=sys.stderr)
