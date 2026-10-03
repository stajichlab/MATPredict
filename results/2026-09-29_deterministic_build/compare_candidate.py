"""Compare the deterministic candidate classifier with the shipped one.

1. Score deltas (candidate - shipped) of the final HMMs on: the classifier's
   training proteins, the Zygo 23 locus proteins (held out), and the 189 HMG
   paralog negatives.
2. Replay the latest PR #9 Mucoromycota reports (Mucoromycota_41bd471, code
   41bd471 = shipped HMMs): rebuild each model-typed call's core proteins from
   the report's exon coordinates, score them with both classifiers, check the
   shipped scores reproduce the report, and list calls whose verdict would
   change (typing at min_margin 25, paralog class, MAT-gene gate threshold).

Run from the deterministic-build worktree with PYTHONPATH=src.
"""
import csv
import gzip
import sys
import tempfile
import types
from pathlib import Path

import yaml

from MATPredict.detect import classifier as C
from MATPredict.detect.classifier_build import read_fasta, training_set
from MATPredict.detect.family_registry import FamilyKey, load_all_families
from MATPredict.detect.pipeline import _translate_model

HERE = Path(__file__).resolve().parent
DB = Path("db")
SHIPPED = DB / "Mucoromycota/classifiers/MAT"
CAND = HERE / "candidate"
RUNS = HERE.parent / "2026-09-29_r4_paralog/Mucoromycota_41bd471/runs"
GENOMES = Path("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes")
ZYGO = HERE.parent / "2026-09-26_sexMP_hmm/zygo_locus_proteins.faa"
MIN_FLANK_IDENTITY, MIN_FLANKS = 40.0, 2

fam = next(f for f in load_all_families(DB) if f.key == FamilyKey("Mucoromycota", "MAT"))
spec = dict(fam.idiomorph_classifier)
C._CACHE.clear()
ship = C.load_classifier({**spec, "dir": str(SHIPPED.resolve())}, fam)
cand = C.load_classifier({**spec, "dir": str(CAND.resolve())}, fam)
gate = {n: (yaml.safe_load((d / "manifest.yaml").read_text()).get("mat_gene_gate") or {}).get("min_score")
        for n, d in (("shipped", SHIPPED), ("candidate", CAND))}
roster_flanks = {g["name"] for g in fam.genes if str(g.get("role", "")).startswith("flanking")}

out = open(HERE / "score_deltas.tsv", "w")
out.write("set\tid\tshipped_Plus\tshipped_Minus\tcand_Plus\tcand_Minus\tdelta_best\tdelta_margin\n")
summary = []


def per_protein(clf, seqs):
    return {k: C.score_proteins(clf, [v]) for k, v in seqs.items()}


def deltas(name, seqs):
    s, c = per_protein(ship, seqs), per_protein(cand, seqs)
    db, dm, flips = [], [], 0
    for k in seqs:
        sb, cb = max(s[k].values()), max(c[k].values())
        sm = s[k]["Plus"] - s[k]["Minus"]
        cm = c[k]["Plus"] - c[k]["Minus"]
        db.append(cb - sb)
        dm.append(cm - sm)
        flips += (sm > 0) != (cm > 0)
        out.write(f"{name}\t{k}\t{s[k]['Plus']:.1f}\t{s[k]['Minus']:.1f}\t{c[k]['Plus']:.1f}\t"
                  f"{c[k]['Minus']:.1f}\t{cb - sb:.1f}\t{cm - sm:.1f}\n")
    gate_s = sum(max(s[k].values()) >= gate["shipped"] for k in seqs)
    gate_c = sum(max(c[k].values()) >= gate["candidate"] for k in seqs)
    summary.append(f"{name}: n={len(seqs)} best-score delta min {min(db):.1f} max {max(db):.1f} "
                   f"mean|d| {sum(map(abs, db)) / len(db):.2f}; margin delta min {min(dm):.1f} "
                   f"max {max(dm):.1f}; Plus/Minus sign flips {flips}; at/above gate "
                   f"shipped {gate_s} ({gate['shipped']}) vs candidate {gate_c} ({gate['candidate']})")


rows = training_set(DB, fam, SHIPPED)
deltas("training", {r["id"]: r["sequence"] for r in rows})
deltas("zygo23_heldout", read_fasta(ZYGO))
deltas("paralog_negatives", read_fasta(SHIPPED / "paralog_negatives.faa"))
out.close()

# --- replay the reports ---------------------------------------------------------
changes = []
n_calls = n_repro = 0
repro_fail = []
for rep in sorted(RUNS.glob("*/detection_report.yaml")):
    r = yaml.safe_load(rep.read_text())
    calls = [x for x in (r.get("detected") or [])
             if x.get("family", "").startswith("Mucoromycota") and x.get("idiomorph_classifier")
             and x["idiomorph_classifier"].get("classifier_input") == "model"]
    if not calls:
        continue
    asm = rep.parent.name
    gz = GENOMES / f"{asm}.fa.gz"
    if not gz.exists():
        continue
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=True) as tmp:
        with gzip.open(gz, "rt") as fh:
            for line in fh:
                tmp.write(line)
        tmp.flush()
        cache = {}
        for x in calls:
            n_calls += 1
            clf_rep = x["idiomorph_classifier"]
            prots = []
            for g in x.get("gene_evidence") or []:
                if g.get("role") != "core_MAT" or not g.get("exons"):
                    continue
                if str(g.get("status", "")).startswith(("unpolished", "not_polish")):
                    continue
                m = types.SimpleNamespace(contig=g["contig"], start=g["start"], end=g["end"],
                                          strand=g["strand"], exons=[types.SimpleNamespace(**e)
                                                                     for e in g["exons"]])
                p = _translate_model(Path(tmp.name), m, r.get("genetic_code") or 1, cache)
                if p:
                    prots.append(p)
            if not prots:
                continue
            vs = C.classify(ship, prots)
            vc = C.classify(cand, prots)
            rep_scores = clf_rep["scores"]
            ok = all(abs(vs.scores[k] - rep_scores[k]) <= 0.15 for k in rep_scores)
            n_repro += ok
            if not ok:
                repro_fail.append((asm, x["contig"], rep_scores, {k: round(v, 1) for k, v in vs.scores.items()}))
            flanks = {g["gene"] for g in x.get("gene_evidence") or []
                      if g.get("gene") in roster_flanks and (g.get("identity") or 0) >= MIN_FLANK_IDENTITY
                      and not str(g.get("status", "")).startswith(("unpolished", "not_polish"))}

            def verdict(v, thr):
                if v.paralog_class:
                    return "paralog"
                if v.idiomorph == C.UNDETERMINED:
                    return "undetermined"
                ok_gate = max(v.scores.values()) >= thr or len(flanks) >= MIN_FLANKS
                return v.idiomorph if ok_gate else "gate_withheld"

            a, b = verdict(vs, gate["shipped"]), verdict(vc, gate["candidate"])
            if a != b:
                changes.append(dict(reproduced=ok, genome=asm, contig=x["contig"], start=x["start"],
                                    reported=x["idiomorph"], shipped=a, candidate=b,
                                    shipped_scores={k: round(v, 1) for k, v in vs.scores.items()},
                                    candidate_scores={k: round(v, 1) for k, v in vc.scores.items()},
                                    flanks=sorted(flanks)))

with open(HERE / "replay_changes.tsv", "w") as fh:
    w = csv.DictWriter(fh, fieldnames=["reproduced", "genome", "contig", "start", "reported", "shipped",
                                       "candidate", "shipped_scores", "candidate_scores", "flanks"],
                       delimiter="\t")
    w.writeheader()
    for c in changes:
        w.writerow(c)
lines = summary + [
    f"replay: model-typed Mucoromycota calls {n_calls}; shipped scores reproduced (<=0.15 bits) "
    f"{n_repro}; verdict changes candidate vs shipped: {sum(c['reproduced'] for c in changes)} among "
    f"reproduced, {sum(not c['reproduced'] for c in changes)} among the {n_calls - n_repro} "
    f"not reproduced (proxy: rebuilt report model only)",
    f"gate thresholds: shipped {gate['shipped']} candidate {gate['candidate']}",
]
for c in changes:
    lines.append(f"  CHANGE (reproduced={c['reproduced']}) {c['genome']} {c['contig']}:{c['start']} {c['shipped']} -> {c['candidate']} "
                 f"ship {c['shipped_scores']} cand {c['candidate_scores']} flanks {c['flanks']}")
for f in repro_fail[:15]:
    lines.append(f"  NOT_REPRODUCED {f}")
(HERE / "compare_output.txt").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
