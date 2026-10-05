#!/usr/bin/env python3
"""sexM alignments that include the Mycotypha africana record (curator request 2026-10-04).

A. training: the 12 sexM training proteins exactly as the classifier build uses them
   (Mycotypha trimmed to aa 26-218 by classifier_training_region), MAFFT L-INS-i with the
   build's options. B. full length: the same 12 records with Mycotypha's full 531-aa ORF.
Writes sexM_training.afa, sexM_fulllength.afa and PNG/SVG views (pyMSAviz).
Run from the repo root with the locked env; PYTHONPATH must include src and pyMSAviz.
"""
import subprocess, sys
from pathlib import Path
from MATPredict.detect.classifier_build import training_set, read_fasta
from MATPredict.detect.family_registry import load_all_families
from pymsaviz import MsaViz

DB = Path("db"); OUT = Path("results/2026-10-04_mycotypha_alignment")
fam = next(f for f in load_all_families(DB) if (f.key.phylum, f.key.locus_name) == ("Mucoromycota", "MAT"))
rows = [r for r in training_set(DB, fam, DB / "Mucoromycota/classifiers/MAT") if r["gene"] == "sexM"]

def short(i):
    rec = i.split("|")[1] if i.startswith("REF|") else i.split("|")[0]
    return rec.replace("_MAT_Minus", "").replace("_MAT_combined", "_combined")

def mafft(seqs, path):
    tmp = path.with_suffix(".in.faa")
    tmp.write_text("".join(f">{k}\n{v}\n" for k, v in sorted(seqs.items())))
    out = subprocess.run(["mafft", "--localpair", "--maxiterate", "1000", "--quiet", "--thread", "1", str(tmp)],
                         capture_output=True, text=True, check=True).stdout
    path.write_text(out); tmp.unlink()

train = {short(r["id"]): r["sequence"] for r in rows}
mafft(train, OUT / "sexM_training.afa")
full = dict(train)
myc = [k for k in full if "64632" in k][0]
prot = read_fasta(DB / "Mucoromycota/Mucorales/64632_nrrl-2978_MAT_combined/proteins.faa")
full[myc] = next(v for k, v in prot.items() if "name=sexM" in k)
mafft(full, OUT / "sexM_fulllength.afa")
for name in ("sexM_training", "sexM_fulllength"):
    mv = MsaViz(str(OUT / f"{name}.afa"), wrap_length=100, show_count=True, show_consensus=True,
                color_scheme="Clustal")
    fig = mv.plotfig()
    fig.savefig(OUT / f"{name}.png", dpi=150, bbox_inches="tight")
    fig.savefig(OUT / f"{name}.svg", bbox_inches="tight")
    print(name, len(read_fasta(OUT / f"{name}.afa")), "seqs")
