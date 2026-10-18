"""Genome lists for the patched runs (run on the HPCC): python build_lists.py OUTDIR REPO FLIP_TSV AUDIT_LIST
sets: dothideo (the 2,730 of the audit), truth (BioSample-stated mating types, all classes), flip (flip-review loci),
holdout (genomes of the curated-record holdout manifest, all phyla), extra (A. fumigatus 40, Neurospora, Sordaria).
Writes OUTDIR/sets.tsv (set, genome, taxid) and OUTDIR/chunks/chunk_NNN.tsv (genome<TAB>taxid, 25 per chunk, unique genomes)."""
import csv, os, sys
from pathlib import Path

out, repo, flip_tsv, audit_list = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3], sys.argv[4]
BFD = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/ncbi_accessions_taxonomy.csv"
LIB = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes"
bfd = list(csv.DictReader(open(BFD)))
by_id = {r["ASM_ACCESSION"]: r for r in bfd}
by_acc = {}
for r in bfd:
    by_acc.setdefault("_".join(r["ASM_ACCESSION"].split("_")[:2]), r)
sets = []


def add(name, gid):
    r = by_id.get(gid) or by_acc.get(gid)
    if r and os.path.exists(f"{LIB}/{r['ASM_ACCESSION']}.fa.gz") and r["NCBI_TAXID"]:
        sets.append((name, r["ASM_ACCESSION"], r["NCBI_TAXID"]))


for r in csv.DictReader(open(audit_list), delimiter="\t"):
    if r["has_fasta"] == "True":
        add("dothideo", r["ASM_ACCESSION"])
for r in csv.DictReader(open(repo / "results/2026-10-17_biosample_truth/truth_set.tsv"), delimiter="\t"):
    add("truth", r["genome"])
for r in csv.DictReader(open(flip_tsv), delimiter="\t"):
    add("flip", r["genome"])
for r in csv.DictReader(open(repo / "testset/holdout/holdout_manifest.tsv"), delimiter="\t"):
    if r["assembly"]:
        add("holdout", r["assembly"])
afum = [r for r in bfd if r["SPECIES"] == "Aspergillus fumigatus"][:40]
for r in afum:
    add("extra_afum", r["ASM_ACCESSION"])
for r in bfd:
    if r["SPECIES"] in ("Neurospora crassa", "Sordaria macrospora", "Neurospora tetrasperma", "Podospora anserina"):
        add("extra_neurospora_sordaria", r["ASM_ACCESSION"])
out.mkdir(parents=True, exist_ok=True)
with open(out / "sets.tsv", "w") as fh:
    fh.write("set\tgenome\ttaxid\n")
    for s in sets:
        fh.write("\t".join(s) + "\n")
seen, uniq = set(), []
for _, g, t in sets:
    if g not in seen:
        seen.add(g); uniq.append((g, t))
(out / "chunks").mkdir(exist_ok=True)
for i in range(0, len(uniq), 25):
    with open(out / "chunks" / f"chunk_{i // 25:03d}.tsv", "w") as fh:
        for g, t in uniq[i:i + 25]:
            fh.write(f"{g}\t{t}\n")
from collections import Counter
print(Counter(s for s, _, _ in sets), len(uniq), "unique genomes,", (len(uniq) + 24) // 25, "chunks")
