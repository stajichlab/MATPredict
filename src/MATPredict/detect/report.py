"""GFF3 and validation-style YAML report writers for detection results.

Both writers consume a `DetectionOutcome` so that families which were
attempted but fell short of the ambiguity floor are reported as "not
detected" with a reason, rather than silently omitted (spec section 3).

The GFF3 follows `db/gff_export.write_gff3`'s conventions so predictions and
curated records are directly comparable: 1-based fully-closed coordinates, a
`##sequence-region` line per segment, and one `gene` feature per detected
gene carrying `Name`, `role` and `present`.
"""
from __future__ import annotations

from pathlib import Path

import yaml

from MATPredict.detect.pipeline import DetectionOutcome, DetectionResult


def _family_label(key) -> str:
    return f"{key.phylum}:{key.locus_name}"


def _locus_id(result: DetectionResult, index: int) -> str:
    return f"{result.family_key.phylum}_{result.family_key.locus_name}_{index}"


def write_detection_gff3(outcome: DetectionOutcome, out_path: Path) -> None:
    """Write one locus-region feature plus one gene feature per detected gene.

    Genes the family expects but which were not found are emitted with
    `present=false` at the locus region's own coordinates (they have no
    coordinates of their own), mirroring sub-project 1's `present` semantics
    so an absence is explicit rather than an empty space. Genes classified
    "not searchable by this method" (short-ORF case) carry
    `not_searchable=true` so an absence caused by a tool limitation is never
    read as a curated negative.
    """
    lines = ["##gff-version 3"]
    for index, r in enumerate(outcome.results):
        segments = r.segments or [None]
        for segment in segments:
            contig = segment.contig if segment else r.contig
            start = segment.start if segment else r.start
            end = segment.end if segment else r.end
            lines.append(f"##sequence-region {contig} {start} {end}")

        locus_id = _locus_id(r, index)
        attrs = (
            f"ID={locus_id};family={_family_label(r.family_key)};confidence={r.confidence}"
            f";idiomorph={r.idiomorph};fragmented={str(r.fragmented).lower()}"
        )
        if r.ambiguous_with:
            attrs += ";ambiguous_with=" + ",".join(_family_label(k) for k in r.ambiguous_with)
        if r.reference_records:
            attrs += ";reference_records=" + ",".join(r.reference_records)
        lines.append("\t".join([
            r.contig, "MATPredict", "MAT_locus", str(r.start), str(r.end), ".", ".", ".", attrs,
        ]))

        for gene_index, evidence in enumerate(r.gene_evidence):
            gene_attrs = (
                f"ID={locus_id}.gene{gene_index};Parent={locus_id};Name={evidence.gene_name}"
                f";role={evidence.role};present=true;identity={evidence.identity}"
                f";reference_record={evidence.reference_record_id};method={evidence.method}"
            )
            if evidence.coverage is not None:
                gene_attrs += f";coverage={evidence.coverage}"
            lines.append("\t".join([
                evidence.contig, "MATPredict", "gene", str(evidence.start), str(evidence.end),
                ".", evidence.strand or ".", ".", gene_attrs,
            ]))

        absent_index = len(r.gene_evidence)
        for gene_name in r.genes_missing:
            lines.append("\t".join([
                r.contig, "MATPredict", "gene", str(r.start), str(r.end), ".", ".", ".",
                f"ID={locus_id}.gene{absent_index};Parent={locus_id};Name={gene_name};present=false",
            ]))
            absent_index += 1
        for gene_name in r.genes_not_searchable:
            lines.append("\t".join([
                r.contig, "MATPredict", "gene", str(r.start), str(r.end), ".", ".", ".",
                f"ID={locus_id}.gene{absent_index};Parent={locus_id};Name={gene_name}"
                f";present=false;not_searchable=true",
            ]))
            absent_index += 1

    out_path.write_text("\n".join(lines) + "\n")


def _result_doc(r: DetectionResult) -> dict:
    return {
        "family": _family_label(r.family_key),
        "contig": r.contig,
        "start": r.start,
        "end": r.end,
        "confidence": r.confidence,
        "idiomorph": r.idiomorph,
        "ambiguous_with": [_family_label(k) for k in r.ambiguous_with],
        "genes_found": r.genes_found,
        "genes_missing": r.genes_missing,
        "genes_not_searchable": r.genes_not_searchable,
        "fragmented": r.fragmented,
        "reference_records": r.reference_records,
        "segments": [
            {
                "contig": s.contig,
                "start": s.start,
                "end": s.end,
                "contig_edge_distance": s.contig_edge_distance,
            }
            for s in r.segments
        ],
        "gene_evidence": [
            {
                "gene": e.gene_name,
                "role": e.role,
                "contig": e.contig,
                "start": e.start,
                "end": e.end,
                "strand": e.strand,
                "identity": e.identity,
                "coverage": e.coverage,
                "reference_record": e.reference_record_id,
                "method": e.method,
            }
            for e in r.gene_evidence
        ],
    }


def write_detection_report(outcome: DetectionOutcome, out_path: Path) -> None:
    doc = {
        "families_attempted": [_family_label(k) for k in outcome.families_attempted],
        "detected": [_result_doc(r) for r in outcome.results],
        "not_detected": [
            {
                "family": _family_label(n.family_key),
                "reason": n.reason,
                "best_fraction_found": n.best_fraction_found,
                "genes_found": n.genes_found,
                "genes_missing": n.genes_missing,
                "genes_not_searchable": n.genes_not_searchable,
            }
            for n in outcome.not_detected
        ],
    }
    out_path.write_text(yaml.safe_dump(doc, sort_keys=False))
