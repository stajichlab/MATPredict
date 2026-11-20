"""Census of every `polished_disagree` gene model in a frozen run: how often would the other tool's model be the cleaner one (#290, #349).

Reads the `detection_report.yaml` files of a run from a tar stream (`zstd -dc reports_all.tar.zst | python model_disagreement_census.py - ...`),
and for every gene of every detected locus with `status: polished_disagree` and an `alternate_model` applies the splice-motif and
model-quality check of `model_disagreement_check.py` (the same CLEAN definition) to both models on the genome sequence.

Per gene (`census.tsv`): the family, gene, role, tier, both tools, both models' aa length, introns (class), internal stops, the identity of
each, and the class
  both_clean            neither model has an issue
  primary_issue_alt_clean   the primary has an issue and the alternate is clean (the case a "prefer the alternate" rule would switch)
  primary_clean_alt_issue   the reverse
  both_issue            neither is clean
  not_checked           genome, contig or genetic code not available
plus `trigger_noncanonical` (the primary has a noncanonical intron and the alternate has none), the change a switch would make to the gene
span (`d_start`, `d_end`, `d_aa`, `alt_outside_locus`: the alternate model lies partly outside the reported locus span).
`summary.tsv` counts the classes by run, family, gene and tier. Nothing here changes a report.

Usage: model_disagreement_census.py REPORTS_TAR|- GENOME_DIR DB_ROOT OUT_PREFIX
"""
from __future__ import annotations

import csv
import sys
import tarfile
from collections import Counter
from functools import lru_cache
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import model_disagreement_check as mdc  # noqa: E402


def find_genome(genome_dir: Path, asm: str) -> Path | None:
    for p in sorted(genome_dir.glob(asm + "*")):
        if p.name.endswith((".fasta.gz", ".fa.gz", ".fna.gz", ".fasta", ".fa", ".fna")) and ".masked" not in p.name:
            return p
    return None


def classify(p_issues: list[str], a_issues: list[str]) -> str:
    if not p_issues and not a_issues:
        return "both_clean"
    if p_issues and not a_issues:
        return "primary_issue_alt_clean"
    if a_issues and not p_issues:
        return "primary_clean_alt_issue"
    return "both_issue"


def noncanonical(m: dict) -> bool:
    return any(i["class"] == "noncanonical" for i in m["introns"])


def row_for(sample: str, locus: dict, e: dict, seq: str | None, table: int, ref_aa: int | None) -> dict:
    base = {"sample": sample, "family": locus["family"], "tier": locus["confidence"], "gene": e["gene"], "role": e.get("role"), "contig": e["contig"],
            "strand": e["strand"], "primary_tool": e.get("method"), "alternate_tool": e["alternate_model"].get("method"),
            "primary_identity": e.get("identity"), "alternate_identity": e["alternate_model"].get("identity")}
    if seq is None or len(seq) < max(e["end"], e["alternate_model"]["end"]):
        return {**base, "class": "not_checked", "note": "genome or contig not available"}
    try:
        r = mdc.check_gene(e, seq, table, ref_aa)
    except ValueError as exc:
        return {**base, "class": "not_checked", "note": str(exc)[:80]}
    p, a = r["primary"], r["alternate"]
    pi, ai = mdc.issues(p), mdc.issues(a)
    alt = e["alternate_model"]
    return {**base, "class": classify(pi, ai), "primary_aa": p["aa"], "alternate_aa": a["aa"], "ref_aa": ref_aa,
            "primary_introns": ",".join(i["class"] for i in p["introns"]), "alternate_introns": ",".join(i["class"] for i in a["introns"]),
            "primary_stops": p["internal_stops"], "alternate_stops": a["internal_stops"],
            "primary_issues": "; ".join(pi), "alternate_issues": "; ".join(ai),
            "trigger_noncanonical": noncanonical(p) and not noncanonical(a),
            "d_start": alt["start"] - e["start"], "d_end": alt["end"] - e["end"], "d_aa": a["aa"] - p["aa"],
            "alt_outside_locus": alt["start"] < locus["start"] or alt["end"] > locus["end"], "note": ""}


def _iter_tar(stream):
    with tarfile.open(fileobj=stream, mode="r|") as tf:
        for m in tf:
            if m.isfile() and m.name.endswith("detection_report.yaml"):
                yield Path(m.name).parent.name, yaml.safe_load(tf.extractfile(m).read()) or {}


COLS = ["sample", "family", "tier", "gene", "role", "contig", "strand", "primary_tool", "alternate_tool", "class", "trigger_noncanonical", "primary_aa",
        "alternate_aa", "ref_aa", "primary_introns", "alternate_introns", "primary_stops", "alternate_stops", "primary_identity", "alternate_identity",
        "d_start", "d_end", "d_aa", "alt_outside_locus", "primary_issues", "alternate_issues", "note"]


def run(docs, genome_dir: Path, db_root: Path, out: Path, label: str = "run") -> list[dict]:
    @lru_cache(maxsize=None)
    def ref_len(record, gene):
        return mdc.reference_aa(db_root, record, gene)

    rows = []
    for sample, doc in docs:
        table = int(doc.get("genetic_code") or 1)
        genes = [(locus, e) for locus in doc.get("detected") or [] for e in locus.get("gene_evidence") or []
                 if e.get("status") == "polished_disagree" and e.get("alternate_model")]
        if not genes:
            continue
        fa = find_genome(genome_dir, sample)
        contigs: dict[str, str | None] = {}
        for locus, e in genes:
            if fa is not None and e["contig"] not in contigs:
                contigs[e["contig"]] = mdc.read_contig(fa, e["contig"])
            rows.append({**row_for(sample, locus, e, contigs.get(e["contig"]), table, ref_len(e.get("reference_record"), e["gene"])), "label": label})
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(f"{out}.census.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", restval="", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    summ = []
    def add(group, key, rs):
        for k, v in sorted(Counter(r["class"] for r in rs).items()):
            summ.append([label, group, key, k, v])
    add("all", "all", rows)
    for fam in sorted({r["family"] for r in rows}):
        add("family", fam, [r for r in rows if r["family"] == fam])
    for tier in sorted({r["tier"] for r in rows}):
        add("tier", tier, [r for r in rows if r["tier"] == tier])
    for gene, _ in Counter(r["gene"] for r in rows).most_common(25):
        add("gene", gene, [r for r in rows if r["gene"] == gene])
    summ.append([label, "trigger", "noncanonical_primary_alt_none", "genes", sum(1 for r in rows if r.get("trigger_noncanonical"))])
    sw = [r for r in rows if r["class"] == "primary_issue_alt_clean"]
    summ.append([label, "switch", "alt_outside_locus", "genes", sum(1 for r in sw if r.get("alt_outside_locus"))])
    summ.append([label, "switch", "span_changes", "genes", sum(1 for r in sw if r.get("d_start") or r.get("d_end"))])
    with open(f"{out}.summary.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["run", "group", "key", "class", "genes"])
        w.writerows(summ)
    return rows


def main(argv: list[str]) -> int:
    src, gdir, db, out = argv[:4]
    label = argv[4] if len(argv) > 4 else Path(out).name
    stream = sys.stdin.buffer if src == "-" else open(src, "rb")
    rows = run(_iter_tar(stream), Path(gdir), Path(db), Path(out), label)
    print(label, dict(Counter(r["class"] for r in rows)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
