#!/usr/bin/env python3
"""Pilot step 1: receptor proteins for the labelled STE3-like loci of 6 Agaricomycete genomes.

miniprot (--gff) of ste3_all.faa against each genome; for every labelled locus
(panel_loci.tsv) keep the highest-scoring mRNA that overlaps it and its
translated protein (##STA line). Writes pilot_proteins.faa and pilot_loci.tsv.
Run on SLURM: miniprot takes about 30 s per genome with 8 threads.
"""
import csv
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(sys.argv[0]))
BIN = "/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin"
LIB = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes"
PILOT = {  # genome -> record-id prefix of its own curated REF receptors (held out), short name
    "GCA_016772295.1_ASM1677229v1": ("5346_", "Coprinopsis_cinerea"),
    "GCF_000143185.2_Schco3": ("5334_", "Schizophyllum_commune"),
    "GCF_000271585.1_Trametes_versicolor_v1.0": ("", "Trametes_versicolor"),
    "GCA_001683735.1_ASM168373v1": ("", "Grifola_frondosa"),
    "GCA_984573805.1_gfRusNobi1.hap1.1": ("", "Russula_nobilis"),
    "GCF_000320585.1_Heterobasidion_irregulare_v2.0": ("", "Heterobasidion_irregulare"),
}
THREADS = os.environ.get("SLURM_CPUS_PER_TASK", "8")


def parse_gff(path):
    """[(contig, start, end, score, protein)] for each mRNA."""
    out, cur = [], None
    for line in open(path):
        if line.startswith("##STA"):
            if cur is not None:
                cur["prot"] = line.split("\t", 1)[1].strip().replace("*", "").replace("-", "")
            continue
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) >= 9 and f[2] == "mRNA":
            attrs = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
            cur = {"contig": f[0], "start": int(f[3]), "end": int(f[4]),
                   "score": int(f[5]) if f[5].isdigit() else 0, "prot": ""}
            out.append(cur)
    return out


rows = list(csv.DictReader(open(os.path.join(HERE, "panel_loci.tsv")), delimiter="\t"))
fa = open(os.path.join(HERE, "pilot_proteins.faa"), "w")
tsv = open(os.path.join(HERE, "pilot_loci.tsv"), "w")
tsv.write("species\tasm\tcontig\tstart\tend\tmating\tT_10kb\tHx_10kb\tprot_len\tprot_id\n")
os.makedirs(os.path.join(HERE, "pilot_tmp"), exist_ok=True)
for asm, (prefix, species) in PILOT.items():
    gff = os.path.join(HERE, "pilot_tmp", asm + ".gff")
    with open(gff, "w") as fo:
        subprocess.run([f"{BIN}/miniprot", "-t", THREADS, "-I", "--gff", "--trans", "--outn=5",
                        f"{LIB}/{asm}.fa.gz", os.path.join(HERE, "ste3_all.faa")],
                       stdout=fo, stderr=subprocess.DEVNULL, check=True)
    mrnas = parse_gff(gff)
    if not any(m["prot"] for m in mrnas):
        sys.exit(f"no ##STA protein sequences in {gff}")
    n_missing = 0
    for r in (x for x in rows if x["asm"] == asm):
        s, e = int(r["start"]), int(r["end"])
        hits = [m for m in mrnas if m["contig"] == r["contig"] and m["start"] <= e and m["end"] >= s and m["prot"]]
        if not hits:
            n_missing += 1
            continue
        best = max(hits, key=lambda m: (m["score"], len(m["prot"])))
        pid = f"{species}|{r['contig']}:{s}-{e}|{r['mating']}"
        fa.write(f">{pid}\n{best['prot']}\n")
        tsv.write("\t".join([species, asm, r["contig"], str(s), str(e), r["mating"],
                             r["T_10kb"], r["Hx_10kb"], str(len(best["prot"])), pid]) + "\n")
    print(asm, "loci", sum(x["asm"] == asm for x in rows), "no protein", n_missing)
fa.close()
tsv.close()
