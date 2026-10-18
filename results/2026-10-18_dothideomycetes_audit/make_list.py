"""Genome list: BFD genomes whose CLASS is Dothideomycetes or whose taxid lies under 147541 (taxonkit, taxonkit_check.slurm)."""
import csv, os, sys
CSV = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/ncbi_accessions_taxonomy.csv"
LIB = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes"
out = sys.argv[1]
under = set(open(os.path.join(out, "taxids_under_147541.txt")).read().split())
rows = []
for r in csv.DictReader(open(CSV)):
    if r["CLASS"] == "Dothideomycetes" or r["NCBI_TAXID"] in under:
        r["has_fasta"] = os.path.exists(f"{LIB}/{r['ASM_ACCESSION']}.fa.gz")
        r["class_by_csv"] = r["CLASS"] == "Dothideomycetes"
        r["taxid_under_147541"] = r["NCBI_TAXID"] in under
        rows.append(r)
cols = ["ASM_ACCESSION", "NCBI_TAXID", "SPECIES_IN", "STRAIN", "CLASS", "ORDER", "FAMILY", "GENUS", "SPECIES",
        "has_fasta", "class_by_csv", "taxid_under_147541"]
with open(os.path.join(out, "genome_list.tsv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, cols, delimiter="\t", extrasaction="ignore")
    w.writeheader(); w.writerows(rows)
run = [r for r in rows if r["has_fasta"] and r["NCBI_TAXID"]]
os.makedirs(os.path.join(out, "chunks"), exist_ok=True)
CH = 25
for i in range(0, len(run), CH):
    with open(os.path.join(out, "chunks", f"chunk_{i // CH:03d}.tsv"), "w") as fh:
        for r in run[i:i + CH]:
            fh.write(f"{r['ASM_ACCESSION']}\t{r['NCBI_TAXID']}\n")
print(len(rows), "listed;", len(run), "to run;", (len(run) + CH - 1) // CH, "chunks")
print("class!=Doth but under 147541:", sum(1 for r in rows if not r["class_by_csv"]))
print("class==Doth but not under 147541:", sum(1 for r in rows if not r["taxid_under_147541"]))
