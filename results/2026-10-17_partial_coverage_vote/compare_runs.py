#!/usr/bin/env python
"""Compare partial-core-policy variants against stated mating type.

    compare_runs.py --runs RUNSDIR --jobs jobs_all.tsv --off off --on on --out DIR

jobs.tsv columns: id, taxid, stated, organism, fasta. A genome's call is its best Ascomycota:MAT locus
(confidence high > medium > low, then more polished genes); `undetermined` is a call; no such locus is `none`.
Writes per_genome.tsv (labels for both variants, vote note) and prints the outcome tables.
"""
import argparse, csv, sys
from collections import Counter
from pathlib import Path
import yaml

L = getattr(yaml, "CSafeLoader", yaml.SafeLoader)
RANK = {"high": 3, "medium": 2, "low": 1}


def best_call(report: Path):
    if not report.exists():
        return None
    doc = yaml.load(report.open(), Loader=L)
    loci = [c for c in doc.get("detected") or [] if c.get("family") == "Ascomycota:MAT"]
    if not loci:
        return dict(label="none", confidence="", note=None, labels="", cand="")
    loci.sort(key=lambda c: (-RANK.get(c["confidence"], 0), -(c.get("polished_genes") or 0)))
    c = loci[0]
    note = c.get("idiomorph_vote_note")
    return dict(label=c["idiomorph"], confidence=c["confidence"], note=note,
                labels=",".join(sorted({x["idiomorph"] for x in loci})),
                cand=";".join(f"{x['idiomorph']}={x['score']}" for x in c.get("idiomorph_candidates") or []))


def outcome(stated, label):
    if label is None:
        return "missing"
    if label == "none":
        return "no_locus"
    if label == "undetermined":
        return "undetermined"
    return "correct" if label == stated else "wrong"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", required=True); ap.add_argument("--jobs", required=True)
    ap.add_argument("--off", default="off"); ap.add_argument("--on", default="on"); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rows = []
    for line in open(a.jobs):
        gid, taxid, stated, org, _fa = line.rstrip("\n").split("\t")
        off = best_call(Path(a.runs) / a.off / gid / "detection_report.yaml")
        on = best_call(Path(a.runs) / a.on / gid / "detection_report.yaml")
        rows.append(dict(genome=gid, organism=org, stated=stated,
                         off=off["label"] if off else "missing", off_conf=off["confidence"] if off else "",
                         on=on["label"] if on else "missing", on_conf=on["confidence"] if on else "",
                         off_outcome=outcome(stated, off["label"] if off else None),
                         on_outcome=outcome(stated, on["label"] if on else None),
                         note_decision=(on["note"] or {}).get("decision", "") if on else "",
                         note_reason=(on["note"] or {}).get("reason", "") if on else "",
                         off_candidates=off["cand"] if off else "", on_loci=on["labels"] if on else ""))
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    with open(out / "per_genome.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
    def table(sel, title):
        print(f"\n## {title} (n={len(sel)})")
        print("variant\tcorrect\twrong\tundetermined\tno_locus\tmissing")
        for v in ("off", "on"):
            c = Counter(r[v + "_outcome"] for r in sel)
            print(f"{v}\t{c['correct']}\t{c['wrong']}\t{c['undetermined']}\t{c['no_locus']}\t{c['missing']}")
        ch = [r for r in sel if (r["off"], r["off_conf"]) != (r["on"], r["on_conf"])]
        print(f"label or confidence changed: {len(ch)}")
        for r in ch:
            print(f"  {r['genome']} {r['organism']} stated={r['stated']} {r['off']}({r['off_conf']}) -> {r['on']}({r['on_conf']}) {r['note_decision']}")
    bot = [r for r in rows if "Botrytis" in r["organism"]]
    table(bot, "Botrytis")
    table([r for r in rows if "Botrytis" not in r["organism"]], "non-Botrytis truth genomes")
    table(rows, "all")


if __name__ == "__main__":
    main()
