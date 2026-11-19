"""Re-measure the weak-fragment cap (#275) under the old and the redefined (#341) 'strongest core gene'.

Reads the `detection_report.yaml` files of a frozen run from a tar stream (`zstd -dc reports_all.tar.zst | python
weak_fragment_cap_remeasure.py - OUTDIR`, or a directory of reports instead of `-`), applies `weak_fragment_cap.evaluate` twice to every `Ascomycota:MAT` locus of the reports
and writes:

* `summary.tsv`  one row per rule (old: every MAT1-* gene ranked; new: only the genes of the called idiomorph): loci in scope,
  matched, already low, would change tier (by tier);
* `changes.tsv`  every locus matched by either rule, with genome, tier, idiomorph, and what each rule says (the strongest gene,
  its status and bitscore). Rows are sorted by genome so a curator can go genome by genome.

Nothing is written to the reports. Run on the HPCC; the test is `tests/test_weak_fragment_cap_remeasure.py`.
"""
from __future__ import annotations

import csv
import sys
import tarfile
from collections import Counter
from pathlib import Path

import yaml

from MATPredict.detect import weak_fragment_cap as wfc

FAMILY = "Ascomycota:MAT"


def _stats(locus: dict) -> list[wfc.CoreStat]:
    return [wfc.CoreStat(e["gene"], e.get("status"), e.get("bitscore")) for e in locus.get("gene_evidence") or [] if not e.get("lost_to")]


def evaluate_both(locus: dict):
    """(old decision, new decision) for one report locus; each is a note dict or None."""
    args = (locus["confidence"], locus.get("idiomorph"), ("Ascomycota", "MAT"), _stats(locus))
    new = wfc.evaluate(*args)
    original = wfc.idiomorph_gene_prefix
    wfc.idiomorph_gene_prefix = lambda idiomorph: wfc.CORE_GENE_PREFIX  # the #275 definition: every MAT1-* gene
    try:
        old = wfc.evaluate(*args)
    finally:
        wfc.idiomorph_gene_prefix = original
    return (old.note if old else None), (new.note if new else None)


def run(docs, outdir: Path) -> dict:
    """`docs` yields (genome, report dict). Returns the counts of the summary."""
    outdir.mkdir(parents=True, exist_ok=True)
    counts = {"old": Counter(), "new": Counter()}
    rows = []
    for genome, doc in docs:
        for locus in doc.get("detected") or []:
            if locus.get("family") != FAMILY or locus.get("idiomorph") not in wfc.CALLED_IDIOMORPHS:
                continue
            old, new = evaluate_both(locus)
            for name, note in (("old", old), ("new", new)):
                counts[name]["loci"] += 1
                if note is not None:
                    counts[name]["matched"] += 1
                    counts[name]["change_" + locus["confidence"]] += 1
            if locus["confidence"] == "low":
                continue
            counts["old"]["loci_medium_high"] += 1
            counts["new"]["loci_medium_high"] += 1
            if old or new:
                def f(n, k):
                    return "" if n is None else n[k]
                rows.append([genome, locus["contig"], locus["start"], locus["end"], locus["idiomorph"], locus["confidence"],
                             "capped" if old else "kept", f(old, "strongest_core_gene"), f(old, "strongest_core_status"), f(old, "strongest_core_bitscore"),
                             "capped" if new else "kept", f(new, "strongest_core_gene"), f(new, "strongest_core_status"), f(new, "strongest_core_bitscore")])
    rows.sort(key=lambda r: (r[0], r[1], r[2]))
    with open(outdir / "changes.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["genome", "contig", "start", "end", "idiomorph", "tier", "old_rule", "old_gene", "old_status", "old_bits",
                    "new_rule", "new_gene", "new_status", "new_bits"])
        w.writerows(rows)
    with open(outdir / "summary.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["rule", "loci_called", "loci_medium_high", "would_change_tier", "from_medium", "from_high"])
        for name in ("old", "new"):
            c = counts[name]
            w.writerow([name, c["loci"], c["loci_medium_high"], c["matched"], c["change_medium"], c["change_high"]])
    return counts


def _iter_tar(stream):
    with tarfile.open(fileobj=stream, mode="r|") as tf:
        for m in tf:
            if m.isfile() and m.name.endswith("detection_report.yaml"):
                yield Path(m.name).parent.name, yaml.safe_load(tf.extractfile(m).read()) or {}


def _iter_dir(root: Path):
    for report in sorted(root.rglob("detection_report.yaml")):
        yield report.parent.name, yaml.safe_load(report.read_text()) or {}


def main(argv: list[str]) -> int:
    src, out = argv
    if src != "-" and Path(src).is_dir():
        docs = _iter_dir(Path(src))
    else:
        docs = _iter_tar(sys.stdin.buffer if src == "-" else open(src, "rb"))
    run(docs, Path(out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
