"""Calibration of the alpha-box / HMG-box profile thresholds from the CURATED records only (no BFD genome):
leave-one-genus-out score of every curated MAT1-1-1 on the alpha profile and every MAT1-2-1 on the HMG profile, plus the
scores of the curated MAT1-1-3 (an HMG-containing protein of the MAT1-1 idiomorph) and MAT1-2-10 on both profiles."""
import sys, tempfile, json
from pathlib import Path
from MATPredict.detect import classifier_build as cb
from MATPredict.detect.family_registry import FamilyKey, load_all_families
import pyhmmer

db = Path(sys.argv[1])
family = next(f for f in load_all_families(db) if f.key == FamilyKey("Ascomycota", "MAT"))
work = Path(tempfile.mkdtemp())
GENES = {"MAT1-1-1": "MAT1-1", "MAT1-2-1": "MAT1-2"}
rows = [r for r in cb.training_set(db, family, work) if r["gene"] in GENES]
loo = cb.leave_one_genus_out(rows, GENES, work)
print("LOGO rows", len(loo))
for r in sorted(loo, key=lambda r: (r["gene"], r["own_score"])):
    print("LOGO", r["gene"], r["genus"], r["id"].split("|")[1], "own", round(r["own_score"], 1), "other", round(r["other_score"], 1), "len",
          len(next(x["sequence"] for x in rows if x["id"] == r["id"])))
# other genes on the full profiles
hm = {}
for g in GENES:
    with pyhmmer.plan7.HMMFile(str(db / "Ascomycota/classifiers/MAT_core" / f"{g}.hmm")) as fh:
        hm[g] = fh.read()
ref = Path(tempfile.mkdtemp()) / "ref.faa"
from MATPredict.detect.reference_fasta import build_reference_fasta
build_reference_fasta(db, ref, family_keys={family.key})
seqs = cb.read_fasta(ref)
for h, s in seqs.items():
    gene = h.split("|")[-1]
    if gene in ("MAT1-1-3", "MAT1-2-10", "MAT1-1-2", "MAT1-1-5", "MAT1-2-4", "MAT1-2-3", "MAT1-2-7", "MAT1-1-4"):
        sc = {g: cb.score(hm[g], {"p": s}).get("p", 0.0) for g in GENES}
        print("OTHER", gene, h.split("|")[0], len(s), {k: round(v, 1) for k, v in sc.items()})
