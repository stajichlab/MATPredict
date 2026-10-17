#!/usr/bin/env python3
"""Join the Coccidioides raw-read inventory to the assembly freeze and apply the species corrections.

Inputs (copied read-only from hpcc, see README.md for the commands):
  samples.csv      .../2025_All_Cocci/Assembly/samples.csv           one row per read set (RunAcc, Strain, Organism, Center, FileBase, LocusTag)
  species.csv      .../NovInvenio_Investigations/studies/fungi/coccidioides_pangenome/species.csv   corrected species + freeze assembly paths
  input_files.txt  `stat -c '%n %s' *` in .../Assembly/input         file names and sizes
Corrections (coccidioides_pangenome/CORRECTIONS.md, 2026-09-21): 485B-1_L_OLD_CPA0023 immitis -> posadasii; B3476 posadasii -> immitis.
  aaftf.txt        `ls Assembly/asm/AAFTF/*.sorted.fasta | xargs -n1 basename | sed 's/.sorted.fasta//'`   strains with an AAFTF assembly
The freeze's DNA/*.scaffolds.fa links all dangle (target Assembly/unscaffolded_annotation/ is gone), so `assembly` is the AAFTF
`<strain>.sorted.fasta` (unscaffolded assembly that went into the annotation); `freeze_assembly_link` keeps the dangling path for reference.
Usage: build_manifest.py samples.csv species.csv input_files.txt aaftf.txt > manifest.tsv
"""
import csv, re, sys

CORRECTED = {"485B-1_L_OLD_CPA0023": "Coccidioides posadasii", "B3476": "Coccidioides immitis"}


AAFTF = "/bigdata/stajichlab/shared/projects/Coccidioides/PopGenomics/2025_All_Cocci/Assembly/asm/AAFTF/{}.sorted.fasta"


def main(samples, species, files, aaftf):
    have = {l.strip() for l in open(aaftf)}
    sp = {r["Short"]: r for r in csv.DictReader(open(species))}
    size = {}
    for line in open(files):
        n, b = line.split()
        size[n] = int(b)
    out = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    out.writerow(["strain", "read_set", "species_samples_csv", "species", "relabelled", "center", "locus_tag", "read1", "read2", "bytes",
                  "assembly_in_freeze", "assembly", "freeze_assembly_link", "notes"])
    for r in csv.DictReader(open(samples)):
        strain, rs, fb = r["Strain"], r["RunAcc"], r["FileBase"]
        suf = (".fastq.gz", ".fastq.gz")
        pat = fb.replace("[12]", "{}")
        r1, r2 = rs + pat.format(1), rs + pat.format(2)
        assert r1 in size and r2 in size, (strain, r1, r2)
        orig = r["Organism"]
        s = sp.get(strain)
        species_ = s["Species"] if s else orig  # species.csv already holds the corrected label
        relabelled = "yes" if strain in CORRECTED else ""
        if strain in CORRECTED:
            assert species_ == CORRECTED[strain], strain
        if not orig.startswith("Coccidioides"):
            continue  # TX13 is labelled Fusarium solani; dropped at the user's request
        notes = []
        if not s:
            notes.append("no assembly in the pangenome freeze (species.csv)")
        if r["Notes"]:
            notes.append(r["Notes"])
        out.writerow([strain, rs, orig, species_, relabelled, r["Center"], r["LocusTag"], r1, r2, size[r1] + size[r2],
                      "yes" if s else "no", AAFTF.format(strain) if strain in have else "", s["Genome_Accession"] if s else "", "; ".join(notes)])


main(*sys.argv[1:5])
