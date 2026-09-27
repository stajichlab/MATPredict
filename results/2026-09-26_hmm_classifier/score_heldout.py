"""Score held-out proteins with the SHIPPED Mucoromycota:MAT HMMs.

Sets (none of them in the training set):
  zygo      the HMG-box protein annotated inside each Zygo 23 truth locus
            (results/2026-09-26_sexMP_hmm/zygo_locus_proteins.faa, filtered to
            the one PF00505-bearing protein per locus as that experiment did)
  disputed  the 49 disputed proteins (results/2026-09-26_sexMP_hmm/disputed_scores.tsv)
Writes heldout_scores.tsv and prints accuracy / margins.
"""
import csv, sys
sys.path.insert(0, "/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/polish-scope-cuts/src")
from pathlib import Path
from MATPredict.detect.classifier_build import read_fasta, score
import pyhmmer

CL = Path("/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/polish-scope-cuts/db/Mucoromycota/classifiers/MAT")
HM = Path("../2026-09-26_sexMP_hmm")
PH = Path("../2026-09-26_sexMP_phylogeny")
hmms = {}
for g in ("sexM", "sexP"):
    with pyhmmer.plan7.HMMFile(str(CL / f"{g}.hmm")) as fh:
        hmms[g] = fh.read()

zygo = read_fasta(HM / "zygo_locus_proteins.faa")
zdom = {l.split()[0] for l in open(HM / "zygo_pf.domtbl") if not l.startswith("#")}
zygo = {k: v for k, v in zygo.items() if k in zdom}
allp = read_fasta(PH / "all_proteins.faa")
disp = list(csv.DictReader(open(HM / "disputed_scores.tsv"), delimiter="\t"))

rows = []
zs = {g: score(h, zygo) for g, h in hmms.items()}
for k in zygo:
    truth = "sexP" if k.split("|")[2] == "Plus" else "sexM"
    m, p = zs["sexM"][k], zs["sexP"][k]
    call = "sexP" if p > m else "sexM"
    rows.append(dict(set="zygo", id=k, group="", truth=truth, sexM=round(m, 1), sexP=round(p, 1),
                     margin=round(p - m, 1), call=call, correct=call == truth))
did = {r[list(r)[1]] if "candidate_id" not in r else r["candidate_id"]: r for r in disp}
dseq = {k: allp[k] for k in did if k in allp}
ds = {g: score(h, dseq) for g, h in hmms.items()}
for k in dseq:
    m, p = ds["sexM"][k], ds["sexP"][k]
    r = did[k]
    rows.append(dict(set="disputed", id=k, group=r.get("set", r.get("group", "")), truth="",
                     sexM=round(m, 1), sexP=round(p, 1), margin=round(p - m, 1),
                     call="sexP" if p > m else "sexM", correct=""))
with open("heldout_scores.tsv", "w", newline="") as fo:
    w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
z = [r for r in rows if r["set"] == "zygo"]
print(f"zygo held-out: {sum(r['correct'] for r in z)}/{len(z)} correct; "
      f"|margin| min {min(abs(r['margin']) for r in z)} (correct ones)")
for r in z:
    if not r["correct"]:
        print("  WRONG", r)
import collections
d = [r for r in rows if r["set"] == "disputed"]
print("disputed:", len(d), "of", len(did), "scored")
for grp, rs in collections.defaultdict(list, {g: [r for r in d if r["group"] == g] for g in {r['group'] for r in d}}).items():
    c = collections.Counter(r["call"] for r in rs)
    print(f"  {grp:32s} n={len(rs)} {dict(c)} margins {sorted(r['margin'] for r in rs)[:3]}..{sorted(r['margin'] for r in rs)[-3:]}")
