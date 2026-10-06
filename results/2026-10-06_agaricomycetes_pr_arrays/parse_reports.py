#!/usr/bin/env python3
"""Parse v0.6.0 detection reports of the Agaricomycetes genomes (reports_all.tar.zst, extracted to REPDIR).
Writes pr_calls.tsv (every detected or withheld PR locus, one row), hd_loci.tsv (detected and withheld HD loci)
and per-call receptor evidence (n receptor gene rows, identity).
Usage: parse_reports.py REPDIR"""
import glob, os, sys
import pandas as pd
import yaml

rep = sys.argv[1]
ag = set(pd.read_csv("agari_genomes.tsv", sep="\t").genome)
pr, hd = [], []
for f in glob.glob(f"{rep}/*/runs/*/detection_report.yaml"):
    asm = f.split("/")[-2]
    if asm not in ag:
        continue
    d = yaml.load(open(f), Loader=yaml.CSafeLoader)
    for x in d.get("detected") or []:
        fam = x["family"].split(":")[1]
        ev = x.get("gene_evidence") or []
        rec = [e for e in ev if e["gene"] == "pheromone_receptor"]
        if fam in ("PR", "Balpha", "Bbeta"):
            ver = x.get("verification")
            pr.append(dict(genome=asm, family=fam, status="called", contig=x["contig"], start=x["start"], end=x["end"],
                           confidence=x["confidence"], verification=(ver or {}).get("status", "verified_other") if ver else "none",
                           detection_pass=x["detection_pass"], genes_found="|".join(x["genes_found"]),
                           polished_genes=x["polished_genes"], n_receptor_rows=len(rec),
                           rec_identity=max([e.get("identity") or 0 for e in rec], default=0),
                           caax_motif="|".join(str(e.get("caax_motif")) for e in ev if e["gene"] == "caax_precursor"),
                           n_merged=len(x.get("merged_from") or []), n_subloci=len(x.get("subloci") or [])))
        elif fam in ("HD", "bLocus"):
            hd.append(dict(genome=asm, status="called", family=fam, contig=x["contig"], start=x["start"], end=x["end"],
                           genes_found="|".join(x["genes_found"])))
    for x in d.get("suppressed_loci") or []:
        fam = x["family"].split(":")[1]
        if fam in ("PR", "Balpha", "Bbeta"):
            pr.append(dict(genome=asm, family=fam, status="withheld:" + str(x.get("withheld_reason")), contig=x["contig"], start=x["start"], end=x["end"],
                           confidence="", verification="", detection_pass="", genes_found="|".join(x["genes_found"]),
                           polished_genes=x["polished_genes"], n_receptor_rows=0, rec_identity=x.get("best_identity") or 0, caax_motif="",
                           n_merged=0, n_subloci=0))
        elif fam in ("HD", "bLocus"):
            hd.append(dict(genome=asm, status="withheld:" + str(x.get("withheld_reason")), family=fam, contig=x["contig"],
                           start=x["start"], end=x["end"], genes_found="|".join(x["genes_found"])))
pd.DataFrame(pr).to_csv("pr_calls.tsv.gz", sep="\t", index=False)
pd.DataFrame(hd).to_csv("hd_loci.tsv.gz", sep="\t", index=False)
print(len(pr), "PR rows", len(hd), "HD rows")
