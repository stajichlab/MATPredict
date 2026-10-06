"""F3: does the Mucoromycota classifier's min_margin (25 bits) separate correct
calls from HMG paralogs when the input is a FRAGMENT, not a full model?

Positives (truth known): the classifier's training proteins (curated records +
training_extra, via classifier_build.training_set) and the Zygo 23 locus HMG
proteins (held out of training). Each genus is held out in turn: sexM/sexP HMMs
are rebuilt without it (classifier_build.build_hmm) and its proteins scored as
  full     the whole protein
  hmgbox   the PF00505 envelope +-5 aa
  win      3 random 50-90 aa windows that contain the HMG-box midpoint
Positive margin = own-gene score - other-gene score (negative = wrong call).

Negatives: non-locus HMG copies (FastTree tips, status nonlocus, clade
other_HMG), scored with the SHIPPED HMMs (they are not training data) as the
same three input types. Negative "margin" = |sexP - sexM|: the confidence with
which the classifier would (wrongly) assign a paralog.

Read-only: imports the pipeline code from the polish-scope-cuts worktree.
"""
import csv
import random
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

WT = Path("/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/polish-scope-cuts")
sys.path.insert(0, str(WT / "src"))
import pyhmmer  # noqa: E402
from MATPredict.detect.classifier_build import (  # noqa: E402
    build_hmm, read_fasta, score, training_set)
from MATPredict.detect.family_registry import load_all_families  # noqa: E402

R = Path("/bigdata/stajichlab/jstajich/projects/MATPredict/results")
CL = WT / "db/Mucoromycota/classifiers/MAT"
PF = R / "2026-09-26_sexMP_phylogeny/PF00505.hmm"
rng = random.Random(11)
AA = pyhmmer.easel.Alphabet.amino()

with pyhmmer.plan7.HMMFile(str(PF)) as fh:
    PFHMM = fh.read()


def hmg_env(seqs):
    """PF00505 envelope (start, end), 0-based half-open, per sequence."""
    block = pyhmmer.easel.DigitalSequenceBlock(AA, [
        pyhmmer.easel.TextSequence(name=k.encode(), sequence=v).digitize(AA)
        for k, v in seqs.items()])
    out = {}
    pipe = pyhmmer.plan7.Pipeline(AA, E=10, domE=10)
    for hit in pipe.search_hmm(PFHMM, block):
        name = hit.name.decode() if isinstance(hit.name, bytes) else hit.name
        best = max(hit.domains, key=lambda d: d.score)
        out[name] = (best.env_from - 1, best.env_to)
    return out


def fragments(seqs):
    env = hmg_env(seqs)
    out = {}
    for k, s in seqs.items():
        out[(k, "full")] = s
        if k not in env:
            continue
        a, b = env[k]
        out[(k, "hmgbox")] = s[max(0, a - 5):min(len(s), b + 5)]
        mid = (a + b) // 2
        for i in range(3):
            L = rng.randint(50, 90)
            lo = max(0, min(mid - rng.randint(5, L - 5), len(s) - L))
            if len(s) >= 50:
                out[(k, f"win{i}")] = s[lo:lo + L]
    return out


def scoreset(hmms, frags):
    flat = {f"{k}||{t}": v for (k, t), v in frags.items() if len(v) >= 20}
    sc = {g: score(h, flat) for g, h in hmms.items()}
    return {tuple(fk.split("||")): (sc["sexM"][fk], sc["sexP"][fk]) for fk in flat}


def kind(t):
    return "win" if t.startswith("win") else t


def main():
    fam = next(f for f in load_all_families(WT / "db") if f.key.phylum == "Mucoromycota" and f.key.locus_name == "MAT")
    rows = training_set(WT / "db", fam, CL)
    # Zygo 23 held-out locus HMG proteins
    zy = read_fasta(R / "2026-09-26_sexMP_hmm/zygo_locus_proteins.faa")
    zdom = {l.split()[0] for l in open(R / "2026-09-26_sexMP_hmm/zygo_pf.domtbl") if not l.startswith("#")}
    zyrows = []
    for k, v in zy.items():
        if k not in zdom:
            continue
        org = k.split("|")[1]
        zyrows.append(dict(id=k, gene="sexP" if k.split("|")[2] == "Plus" else "sexM",
                           genus=org.split("_")[0], source="zygo", sequence=v))
    pos_rows = []
    work = Path(tempfile.mkdtemp(prefix="f3_"))
    genera = sorted({r["genus"] for r in rows + zyrows})
    for g in genera:
        train = [r for r in rows if r["genus"] != g]
        hmms = {}
        for gene in ("sexM", "sexP"):
            seqs = {r["id"]: r["sequence"] for r in train if r["gene"] == gene}
            hmms[gene], _ = build_hmm(seqs, f"loo_{gene}", work)
        held = [r for r in rows + zyrows if r["genus"] == g]
        frs = fragments({r["id"]: r["sequence"] for r in held})
        sc = scoreset(hmms, frs)
        truth = {r["id"]: (r["gene"], r["source"]) for r in held}
        for (k, t), (m, p) in sc.items():
            gene, src = truth[k]
            own, oth = (p, m) if gene == "sexP" else (m, p)
            pos_rows.append(dict(set="positive", source=src, genus=g, id=k, frag=kind(t),
                                 own=round(own, 1), other=round(oth, 1), margin=round(own - oth, 1)))
    # negatives: non-locus HMG paralogs with the shipped HMMs
    shipped = {}
    for gene in ("sexM", "sexP"):
        with pyhmmer.plan7.HMMFile(str(CL / f"{gene}.hmm")) as fh:
            shipped[gene] = fh.read()
    tips = [r for r in csv.DictReader(open(R / "2026-09-27_sexMP_fasttree/tip_names.tsv"), delimiter="\t")
            if r["status"] == "nonlocus" and r["tree_clade"] == "other_HMG"]
    allp = read_fasta(R / "2026-09-27_sexMP_fasttree/all_proteins.faa")
    neg = {r["source"]: allp[r["source"]] for r in tips if r["source"] in allp}
    sc = scoreset(shipped, fragments(neg))
    neg_rows = [dict(set="paralog", source="nonlocus_other_HMG", genus="", id=k, frag=kind(t),
                     own="", other="", margin=round(abs(p - m), 1), sexM=round(m, 1), sexP=round(p, 1))
                for (k, t), (m, p) in sc.items()]
    out = Path(__file__).parent
    with open(out / "f3_scores.tsv", "w", newline="") as fo:
        keys = ["set", "source", "genus", "id", "frag", "own", "other", "margin", "sexM", "sexP"]
        w = csv.DictWriter(fo, fieldnames=keys, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(pos_rows + neg_rows)
    # summary
    lines = [f"tips negatives available: {len(tips)}; with sequence: {len(neg)}",
             f"positive proteins: {len({r['id'] for r in pos_rows})} over {len(genera)} genera"]
    by = defaultdict(list)
    for r in pos_rows:
        by[("pos", r["frag"])].append(r["margin"])
    for r in neg_rows:
        by[("neg", r["frag"])].append(r["margin"])
    for (s, f), v in sorted(by.items()):
        v = sorted(v)
        q = lambda p: v[min(len(v) - 1, int(p * len(v)))]
        extra = ""
        if s == "pos":
            extra = f" wrong(<=0)={sum(x <= 0 for x in v)}"
        lines.append(f"{s} {f:6s} n={len(v):4d} min={v[0]:7.1f} p05={q(.05):7.1f} p50={q(.5):7.1f} max={v[-1]:7.1f}{extra}")
    for floor in (10, 15, 20, 25, 30, 40, 50):
        for f in ("full", "hmgbox", "win"):
            P = by[("pos", f)]
            N = by[("neg", f)]
            lines.append(f"floor {floor:3d} {f:6s}: positives called correctly {sum(x >= floor for x in P)}/{len(P)}"
                         f"  wrong-but-above-floor {sum(x <= -floor for x in P)}"
                         f"  paralogs above floor {sum(x >= floor for x in N)}/{len(N)}")
    (out / "f3_summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
