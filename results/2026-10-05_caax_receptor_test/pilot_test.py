#!/usr/bin/env python3
"""Pilot step 2: does protein similarity to known mating receptors separate
mating from non-mating STE3-like copies in a held-out species?

Leave-one-species-out over the 6 pilot species. For each held-out species s:
  training mating set  = labelled-mating pilot receptors of the other species
                         + curated REF receptors (variant A: Agaricomycete REFs
                           only, 5346_/5334_, minus s's own; variant B: all 34 REFs
                           minus s's own)
  training other set   = labelled-other pilot copies of the other species
Scores for each held-out locus:
  sim   max normalised local-alignment score to the mating set minus the same to
        the other set (BLOSUM62, open -11, extend -1)
  hmm   bit score of a profile HMM built (mafft + hmmbuild) from the mating set
        (hmmsearch --max, 0 if no hit)
  T     the strict-CAAX flag (reference, binary)
AUC pooled over all held-out loci with a bootstrap 95% interval over loci.
Pilot only: 6 species, 14 mating receptors.
Needs modules hmmer/3.4 and mafft on PATH.
"""
import csv
import os
import random
import subprocess
import sys
import tempfile

from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices

HERE = os.path.dirname(os.path.abspath(sys.argv[0]))
AGARICO_REF = ("5346_", "5334_")
OWN_PREFIX = {"Coprinopsis_cinerea": "5346_", "Schizophyllum_commune": "5334_"}

pilot = list(csv.DictReader(open(os.path.join(HERE, "pilot_loci.tsv")), delimiter="\t"))
seq = {r.id: str(r.seq) for r in SeqIO.parse(os.path.join(HERE, "pilot_proteins.faa"), "fasta")}
ref = {r.id: str(r.seq) for r in SeqIO.parse(os.path.join(HERE, "ste3_all.faa"), "fasta") if r.id.startswith("REF|")}
for r in pilot:
    r["seq"] = seq[r["prot_id"]]

aligner = PairwiseAligner()
aligner.mode = "local"
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score, aligner.extend_gap_score = -11, -1
_self = {}


def selfscore(s):
    if s not in _self:
        _self[s] = aligner.score(s, s)
    return _self[s]


def norm(a, b):
    return aligner.score(a, b) / (selfscore(a) * selfscore(b)) ** 0.5


def auc(pos, neg):
    if not pos or not neg:
        return float("nan")
    w = sum((p > n) + 0.5 * (p == n) for p in pos for n in neg)
    return w / (len(pos) * len(neg))


def hmm_scores(train_pos, targets):
    """{target_id: bit score} from an HMM of train_pos (dict id->seq)."""
    with tempfile.TemporaryDirectory(dir=os.environ.get("SCRATCH", None)) as d:
        for name, items in (("train", train_pos), ("tgt", targets)):
            with open(f"{d}/{name}.faa", "w") as fo:
                for k, v in items.items():
                    fo.write(f">{k}\n{v}\n")
        with open(f"{d}/train.aln", "w") as fo:
            subprocess.run(["mafft", "--auto", "--quiet", f"{d}/train.faa"], stdout=fo, check=True)
        subprocess.run(["hmmbuild", "--amino", "-n", "m", f"{d}/m.hmm", f"{d}/train.aln"],
                       stdout=subprocess.DEVNULL, check=True)
        subprocess.run(["hmmsearch", "--max", "-E", "1000", "--tblout", f"{d}/t.tbl", f"{d}/m.hmm", f"{d}/tgt.faa"],
                       stdout=subprocess.DEVNULL, check=True)
        out = {k: 0.0 for k in targets}
        for line in open(f"{d}/t.tbl"):
            if not line.startswith("#"):
                f = line.split()
                out[f[0]] = float(f[5])
        return out


def boot_auc(scored, key, n=2000):
    rng = random.Random(11)
    xs = []
    for _ in range(n):
        s = [rng.choice(scored) for _ in scored]
        a = auc([r[key] for r in s if r["mating"] == "mating"], [r[key] for r in s if r["mating"] != "mating"])
        if a == a:
            xs.append(a)
    xs.sort()
    return xs[int(0.025 * len(xs))], xs[int(0.975 * len(xs)) - 1]


out_lines = []
for variant in ("A", "B"):
    scored = []
    for sp in sorted({r["species"] for r in pilot}):
        held = [r for r in pilot if r["species"] == sp]
        rest = [r for r in pilot if r["species"] != sp]
        own = OWN_PREFIX.get(sp, "\0")
        refs = {k: v for k, v in ref.items()
                if not k.split("|")[1].startswith(own)
                and (variant == "B" or k.split("|")[1].startswith(AGARICO_REF))}
        pos = {r["prot_id"]: r["seq"] for r in rest if r["mating"] == "mating"}
        pos.update(refs)
        neg = [r["seq"] for r in rest if r["mating"] != "mating"]
        if len(pos) < 3:
            continue
        hs = hmm_scores(pos, {r["prot_id"]: r["seq"] for r in held})
        for r in held:
            sp_pos = max(norm(r["seq"], p) for p in pos.values())
            sp_neg = max(norm(r["seq"], n) for n in neg)
            scored.append({"species": sp, "id": r["prot_id"], "mating": r["mating"],
                           "sim": sp_pos - sp_neg, "hmm": hs[r["prot_id"]],
                           "T": float(int(r["T_10kb"]) > 0), "train_pos": len(pos)})
    with open(os.path.join(HERE, f"pilot_scores_{variant}.tsv"), "w") as fo:
        fo.write("species\tid\tmating\tsim\thmm\tT\ttrain_pos\n")
        for r in scored:
            fo.write(f"{r['species']}\t{r['id']}\t{r['mating']}\t{r['sim']:.4f}\t{r['hmm']:.1f}\t{r['T']:.0f}\t{r['train_pos']}\n")
    out_lines.append(f"\nVariant {variant} ({'Agaricomycete REFs only' if variant == 'A' else 'all 34 REFs'}); "
                     f"{len(scored)} held-out loci, {sum(r['mating'] == 'mating' for r in scored)} mating")
    out_lines.append("score\tpooled AUC\tboot95\tper-species AUC (mean over species with both classes)")
    for key in ("sim", "hmm", "T"):
        a = auc([r[key] for r in scored if r["mating"] == "mating"], [r[key] for r in scored if r["mating"] != "mating"])
        lo, hi = boot_auc(scored, key)
        per = []
        for sp in {r["species"] for r in scored}:
            g = [r for r in scored if r["species"] == sp]
            x = auc([r[key] for r in g if r["mating"] == "mating"], [r[key] for r in g if r["mating"] != "mating"])
            if x == x:
                per.append(x)
        out_lines.append(f"{key}\t{a:.2f}\t{lo:.2f}-{hi:.2f}\t{sum(per)/len(per):.2f} (n={len(per)})")
    # top-ranked locus per species is mating?
    for key in ("sim", "hmm"):
        top = []
        for sp in sorted({r["species"] for r in scored}):
            g = [r for r in scored if r["species"] == sp]
            if any(r["mating"] == "mating" for r in g):
                best = max(g, key=lambda r: r[key])
                top.append(best["mating"] == "mating")
        out_lines.append(f"top-ranked locus is mating, by {key}: {sum(top)}/{len(top)} species")

text = "\n".join(out_lines)
print(text)
open(os.path.join(HERE, "pilot_summary.txt"), "w").write(text + "\n")
