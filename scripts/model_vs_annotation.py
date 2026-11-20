"""Which of the two models of a `polished_disagree` gene reproduces the annotated intron structure (#290, #349).

For every `polished_disagree` gene (with an `alternate_model`) of the detection reports of a run on annotated genomes (the FungiDB release-71
run), the gene's two models are compared with the annotated mRNA that overlaps them most: for each model, the annotated introns that lie
inside the model's span are collected and compared with the model's own introns. A model `matches` when its intron set equals the annotated
intron set inside its span (a single-exon model matches when the annotation has no intron there). Per gene:

  primary_matches / alternate_matches / both_match / neither_matches / no_annotation

and the splice-motif class of the same gene from the census (`model_disagreement_census.classify`), so the rule "prefer the alternate when
the primary has an issue and the alternate is clean" can be scored against the annotation. The annotation is a cross-check, not the truth:
a model that differs from it may be the better one. Nothing here changes a report.

Usage: model_vs_annotation.py RUNS_DIR GENOME_TEMPLATE GFF_TEMPLATE DB_ROOT OUT_PREFIX
  GENOME_TEMPLATE / GFF_TEMPLATE  with `{sample}`, e.g. /.../genome/FungiDB-71_{sample}_Genome.fasta.gz  /.../gff/FungiDB-71_{sample}.gff.gz
"""
from __future__ import annotations

import csv
import gzip
import sys
from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import model_disagreement_census as census  # noqa: E402
import model_disagreement_check as mdc  # noqa: E402


def load_mrnas(gff: Path) -> dict[str, list[dict]]:
    """{contig: [{id, strand, cds: [(s, e)], introns: [(s, e)]}]} from the CDS features of a GFF3 (grouped by Parent mRNA)."""
    cds: dict[str, dict] = defaultdict(lambda: {"cds": []})
    op = gzip.open if str(gff).endswith(".gz") else open
    with op(gff, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            parent = next((a.split("=", 1)[1] for a in f[8].split(";") if a.startswith("Parent=")), None)
            if parent is None:
                continue
            m = cds[parent]
            m.update(contig=f[0], strand=f[6], id=parent)
            m["cds"].append((int(f[3]), int(f[4])))
    out: dict[str, list[dict]] = defaultdict(list)
    for m in cds.values():
        m["cds"].sort()
        m["introns"] = [(e1 + 1, s2 - 1) for (_, e1), (s2, _) in zip(m["cds"], m["cds"][1:])]
        out[m["contig"]].append(m)
    return out


def overlap(a: tuple[int, int], b: tuple[int, int]) -> int:
    return max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1)


def best_mrna(mrnas: dict, contig: str, strand: str, exons: list[tuple[int, int]]) -> dict | None:
    """The annotated mRNA on `strand` whose CDS overlaps the model's exons most (None when none overlaps)."""
    best, best_n = None, 0
    for m in mrnas.get(contig, []):
        if m["strand"] != strand:
            continue
        n = sum(overlap(c, x) for c in m["cds"] for x in exons)
        if n > best_n:
            best, best_n = m, n
    return best


def matches(model_exons: list[tuple[int, int]], mrna: dict) -> bool:
    """True when the model's introns equal the annotated introns that lie inside the model's span."""
    exons = sorted(model_exons)
    lo, hi = exons[0][0], exons[-1][1]
    mine = {(e1 + 1, s2 - 1) for (_, e1), (s2, _) in zip(exons, exons[1:])}
    annotated = {i for i in mrna["introns"] if i[0] >= lo and i[1] <= hi}
    return mine == annotated


def verdict(primary_matches: bool, alternate_matches: bool) -> str:
    return {(True, False): "primary_matches", (False, True): "alternate_matches", (True, True): "both_match", (False, False): "neither_matches"}[(primary_matches, alternate_matches)]


COLS = ["sample", "family", "tier", "gene", "role", "class", "trigger_noncanonical", "annotation", "annotated_mrna", "primary_introns", "alternate_introns", "annotated_introns_in_span"]


def run(docs, genome_tmpl: str, gff_tmpl: str, db_root: Path, out: Path) -> list[dict]:
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
        genome, gff = Path(genome_tmpl.format(sample=sample)), Path(gff_tmpl.format(sample=sample))
        mrnas = load_mrnas(gff) if gff.exists() else None
        contigs: dict[str, str | None] = {}
        for locus, e in genes:
            if genome.exists() and e["contig"] not in contigs:
                contigs[e["contig"]] = mdc.read_contig(genome, e["contig"])
            r = census.row_for(sample, locus, e, contigs.get(e["contig"]), table, ref_len(e.get("reference_record"), e["gene"]))
            row = {k: r.get(k, "") for k in ("sample", "family", "tier", "gene", "role", "trigger_noncanonical")}
            row["class"] = r["class"]
            if r["class"] == "not_checked" or mrnas is None:
                row["annotation"] = "no_annotation" if mrnas is None else "not_checked"
                rows.append(row)
                continue
            alt = e["alternate_model"]
            pe, ae = mdc.model_exons(e), mdc.model_exons(alt)
            m = best_mrna(mrnas, e["contig"], e["strand"], pe + ae)
            if m is None:
                row["annotation"] = "no_annotation"
            else:
                pm, am = matches(pe, m), matches(ae, m)
                row["annotation"] = verdict(pm, am)
                row["annotated_mrna"] = m["id"]
                p = mdc.analyse_model(contigs[e["contig"]], pe, e["strand"], table)
                a = mdc.analyse_model(contigs[e["contig"]], ae, alt.get("strand", e["strand"]), table)
                row["primary_introns"], row["alternate_introns"] = len(p["introns"]), len(a["introns"])
                lo, hi = min(pe[0][0], ae[0][0]), max(pe[-1][1], ae[-1][1])
                row["annotated_introns_in_span"] = sum(1 for i in m["introns"] if i[0] >= lo and i[1] <= hi)
            rows.append(row)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(f"{out}.genes.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", restval="", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    with open(f"{out}.summary.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["group", "class", "annotation", "genes"])
        for group, rs in (("all", rows), ("core_MAT", [r for r in rows if r["role"] == "core_MAT"]),
                          ("trigger_noncanonical", [r for r in rows if str(r["trigger_noncanonical"]) == "True"])):
            for (c, a), n in sorted(Counter((r["class"], r["annotation"]) for r in rs).items()):
                w.writerow([group, c, a, n])
    return rows


def _iter_dir(runs: Path):
    for report in sorted(runs.glob("*/detection_report.yaml")):
        yield report.parent.name, yaml.safe_load(report.read_text()) or {}


def main(argv: list[str]) -> int:
    runs, gt, ft, db, out = argv
    rows = run(_iter_dir(Path(runs)), gt, ft, Path(db), Path(out))
    print(dict(Counter((r["class"], r["annotation"]) for r in rows)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
