"""Legacy (flag off) against patched (`--homothallic-policy both_genes`) for every genome of the patched run.
python compare_patch.py RUNS_DIR AUDIT_DIR OUT_DIR   (AUDIT_DIR = results/2026-10-18_dothideomycetes_audit, from branch dothideomycetes-mat-audit)
A re-classed locus carries `homothallic_note.previous_locus_class`, so the legacy label is recovered from the patched report;
equality of everything else with the audit's own legacy run is checked locus by locus (dothideo set).
Writes OUT_DIR/loci_patched.tsv, genomes_patched.tsv and summary tables."""
import csv, sys
from pathlib import Path
import pandas as pd
import yaml

RUNS, AUDIT, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
OUT.mkdir(parents=True, exist_ok=True)
L = getattr(yaml, "CSafeLoader", yaml.SafeLoader)
sets = pd.read_csv(OUT / "sets.tsv", sep="\t", dtype=str)
rows, grows = [], []
for rep in sorted(RUNS.glob("*/detection_report.yaml")):
    g = rep.parent.name
    d = yaml.load(open(rep), Loader=L)
    par = (d.get("run") or {}).get("parameters") or {}
    det = d.get("detected") or []
    grows.append(dict(genome=g, n_loci=len(det), homothallic_policy=par.get("homothallic_policy"),
                      partial_core_policy=par.get("partial_core_policy")))
    for x in det:
        note = x.get("homothallic_note")
        rows.append(dict(genome=g, family=x["family"], contig=x["contig"], start=int(x["start"]), end=int(x["end"]),
                         idiomorph=x["idiomorph"], confidence=x["confidence"], class_after=x.get("locus_class"),
                         class_before=(note or {}).get("previous_locus_class", x.get("locus_class")),
                         reclassed=bool(note), genes_found="|".join(x.get("genes_found") or []),
                         alpha_gene=(note or {}).get("alpha_box_gene", {}).get("gene"),
                         alpha_bits=(note or {}).get("alpha_box_gene", {}).get("alpha_bits"),
                         alpha_aa=(note or {}).get("alpha_box_gene", {}).get("protein_aa"),
                         hmg_gene=(note or {}).get("hmg_box_gene", {}).get("gene"),
                         hmg_bits=(note or {}).get("hmg_box_gene", {}).get("hmg_bits"),
                         hmg_aa=(note or {}).get("hmg_box_gene", {}).get("protein_aa"),
                         hmg_lost_to=(note or {}).get("hmg_box_gene", {}).get("lost_to"),
                         separation_bp=(note or {}).get("separation_bp")))
P = pd.DataFrame(rows)
P.to_csv(OUT / "loci_patched.tsv", sep="\t", index=False)
pd.DataFrame(grows).to_csv(OUT / "genomes_patched.tsv", sep="\t", index=False)
print(len(grows), "genomes,", len(P), "loci,", int(P.reclassed.sum()), "re-classed;", "policies", pd.DataFrame(grows).homothallic_policy.value_counts().to_dict())

# equality with the audit's legacy run (Dothideomycetes set)
A = pd.read_csv(AUDIT / "loci_audit_flagged.tsv", sep="\t", low_memory=False)
A = A[["genome", "contig", "start", "idiomorph", "confidence", "locus_class", "genes_found"]].rename(
    columns={"idiomorph": "idiomorph_audit", "confidence": "confidence_audit", "locus_class": "class_audit", "genes_found": "genes_audit"})
dg = set(sets[sets.set == "dothideo"].genome)
M = P[P.genome.isin(dg)].merge(A, on=["genome", "contig", "start"], how="outer", indicator=True)
M = M[M.genome.isin(dg)]
both = M[M._merge == "both"]
diff = both[(both.idiomorph != both.idiomorph_audit) | (both.confidence != both.confidence_audit) |
            (both.class_before != both.class_audit) | (both.genes_found != both.genes_audit)]
print("Dothideomycetes loci in both runs:", len(both), "; only in patched run:", int((M._merge == "left_only").sum()),
      "; only in audit run:", int((M._merge == "right_only").sum()), "; differ in idiomorph/confidence/legacy class/genes:", len(diff))
diff.to_csv(OUT / "differences_vs_audit.tsv", sep="\t", index=False)
M[M._merge != "both"].to_csv(OUT / "loci_not_in_both.tsv", sep="\t", index=False)
