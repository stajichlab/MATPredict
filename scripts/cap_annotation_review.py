"""Check the loci the weak-fragment cap would lower against an independent gene annotation (#341).

Inputs
  changes.tsv   from `scripts/weak_fragment_cap_remeasure.py` (rows with `new_rule == capped` are reviewed)
  reports       a tar stream of detection reports (`zstd -dc reports_all.tar.zst | python cap_annotation_review.py - ...`) or a folder
  samples.csv   BFD `Fungi_BFD/samples.csv` (ASMID, SPECIES_IN, STRAIN), to find the genome's annotation folder
  annotation    `Fungi_BFD_runs/genome_annotation` (funannotate `predict_results/*.gff3.gz` per genome; contig names equal the assembly's)

For each capped locus: the core genes of the called idiomorph (the genes the cap ranks) and of the other idiomorph, each with its
bitscore and polish status, and the verdict of the annotation over the gene model: `annotated` (an annotated gene covers at least 50%
of the model), `partial` (overlaps, under 50%), `novel` (no annotated gene), or `no_annotation` (no annotation folder for the genome).
The annotation is a cross-check, not the truth: `novel` does not make the call wrong (the predictor may miss the gene), and `annotated`
does not make it right. In particular an annotated gene over a 40-bit model only says the region is inside a predicted gene: an HMG-box
paralog would be annotated too, so it cannot tell a real MAT1-2-1 from a paralog. The `review` column says only what was measured, for
the called idiomorph's strongest core gene: `annotated_gene_present`, `no_annotated_gene`, `annotation_partial`, or `no_annotation`.
"""
from __future__ import annotations

import csv
import gzip
import re
import sys
import tarfile
from collections import Counter, defaultdict
from pathlib import Path

import yaml

from MATPredict.detect import weak_fragment_cap as wfc


def folder_candidates(species_in: str, strain: str) -> list[str]:
    """Names the BFD annotation folder may have for one sample: the species with the strain, then the species alone."""
    def clean(s: str) -> str:
        return re.sub(r"[^A-Za-z0-9._-]+", "_", s.strip()).strip("_")
    base = clean(species_in)
    out = []
    if strain and clean(strain) not in base:
        out.append(f"{base}_{clean(strain)}")
    out.append(base)
    return out


def load_genes(gff_gz: Path) -> dict[str, list[tuple[int, int]]]:
    genes: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with gzip.open(gff_gz, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > 8 and f[2] == "gene":
                genes[f[0]].append((int(f[3]), int(f[4])))
    return genes


def verdict(genes: dict, contig: str, start: int, end: int) -> str:
    """annotated / partial / novel for one model span against the annotated genes of its contig."""
    span = end - start + 1
    best = max((min(end, e) - max(start, s) + 1 for s, e in genes.get(contig, []) if s <= end and e >= start), default=0)
    if best <= 0:
        return "novel"
    return "annotated" if best / span >= 0.5 else "partial"


def classify(core_verdict: str | None) -> str:
    return {"annotated": "annotated_gene_present", "novel": "no_annotated_gene", None: "no_annotation"}.get(core_verdict, "annotation_partial")


def review_locus(locus: dict, genes: dict | None) -> dict:
    prefix = wfc.idiomorph_gene_prefix(locus["idiomorph"])
    core = [e for e in locus.get("gene_evidence") or [] if not e.get("lost_to") and e["gene"].startswith(wfc.CORE_GENE_PREFIX)]
    own = [e for e in core if e["gene"].startswith(prefix)]
    other = [e for e in core if not e["gene"].startswith(prefix)]
    stats = [wfc.CoreStat(e["gene"], e.get("status"), e.get("bitscore")) for e in own]
    strongest = max(stats, key=wfc._strength).gene if stats else None
    def fmt(es):
        return ";".join(f'{e["gene"]}:{e.get("status")}:{e.get("bitscore")}' + ("" if genes is None else ":" + verdict(genes, e["contig"], e["start"], e["end"])) for e in es)
    sv = None
    if genes is not None and strongest:
        e = next(e for e in own if e["gene"] == strongest)
        sv = verdict(genes, e["contig"], e["start"], e["end"])
    return {"called_idiomorph_genes": fmt(own), "other_idiomorph_genes": fmt(other), "strongest_own_gene": strongest or "",
            "strongest_own_annotation": sv or "", "review": classify(sv) if genes is not None else "no_annotation",
            "genes_found": "|".join(locus.get("genes_found") or [])}


def _iter_tar(stream):
    with tarfile.open(fileobj=stream, mode="r|") as tf:
        for m in tf:
            if m.isfile() and m.name.endswith("detection_report.yaml"):
                yield Path(m.name).parent.name, yaml.safe_load(tf.extractfile(m).read()) or {}


def _iter_dir(root: Path):
    for report in sorted(root.rglob("detection_report.yaml")):
        yield report.parent.name, yaml.safe_load(report.read_text()) or {}


def run(docs, changes: Path, samples: Path, annotation: Path, out: Path) -> Counter:
    capped = {(r["genome"], r["contig"], int(r["start"]), int(r["end"])): r for r in csv.DictReader(changes.open(), delimiter="\t") if r["new_rule"] == "capped"}
    sample = {r["ASMID"]: r for r in csv.DictReader(samples.open())}
    rows, counts, cache = [], Counter(), {}
    for genome, doc in docs:
        for locus in doc.get("detected") or []:
            key = (genome, locus.get("contig"), locus.get("start"), locus.get("end"))
            if key not in capped or locus.get("family") != "Ascomycota:MAT":
                continue
            s = sample.get(genome)
            gff = None
            if s:
                for name in folder_candidates(s["SPECIES_IN"], s["STRAIN"]):
                    hits = sorted((annotation / name / "predict_results").glob("*.gff3.gz"))
                    if hits:
                        gff = hits[0]
                        break
            if gff is not None and gff not in cache:
                cache[gff] = load_genes(gff)
            r = review_locus(locus, cache.get(gff) if gff else None)
            counts[r["review"]] += 1
            rows.append([genome, *key[1:], locus["idiomorph"], locus["confidence"], r["strongest_own_gene"], r["strongest_own_annotation"], r["review"],
                         r["called_idiomorph_genes"], r["other_idiomorph_genes"], r["genes_found"], gff.parent.parent.name if gff else ""])
    out.mkdir(parents=True, exist_ok=True)
    rows.sort(key=lambda r: (r[0], r[1], r[2]))
    with open(out / "review.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["genome", "contig", "start", "end", "idiomorph", "tier", "strongest_own_gene", "annotation", "review",
                    "called_idiomorph_genes(gene:status:bits:annotation)", "other_idiomorph_genes", "genes_found", "annotation_folder"])
        w.writerows(rows)
    with open(out / "summary.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["review", "loci"])
        for k, v in sorted(counts.items()):
            w.writerow([k, v])
    return counts


def main(argv: list[str]) -> int:
    src, changes, samples, annotation, out = argv
    docs = _iter_dir(Path(src)) if src != "-" and Path(src).is_dir() else _iter_tar(sys.stdin.buffer if src == "-" else open(src, "rb"))
    print(dict(run(docs, Path(changes), Path(samples), Path(annotation), Path(out))))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
