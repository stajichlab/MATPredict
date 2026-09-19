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

import logging
from pathlib import Path

import yaml
from Bio import SeqIO

from MATPredict.detect.benchmark import _extract_translated_gene, _open_fasta_text
from MATPredict.detect.pipeline import DetectionOutcome, DetectionResult

logger = logging.getLogger(__name__)


def _family_label(key) -> str:
    return f"{key.phylum}:{key.locus_name}"


def _locus_id(result: DetectionResult, index: int) -> str:
    return f"{result.family_key.phylum}_{result.family_key.locus_name}_{index}"


def write_detection_gff3(
    outcome: DetectionOutcome, out_path: Path, genome_fasta: Path | None = None,
) -> None:
    """Write one locus-region feature plus one gene feature per detected gene.

    When `genome_fasta` is given (the rollout genome's own FASTA, already on
    local disk at detection time -- no network fetch needed), two additions
    are made, both purely additive so a caller that omits `genome_fasta`
    (the default) gets byte-identical output to before this parameter
    existed:

    * A companion FASTA is written at `out_path.with_suffix(".fasta")`,
      using the SAME contig/seqid names the GFF3 itself uses, so clinker's
      GFF3+FASTA input convention (same base filename) is satisfied. Each
      contig referenced anywhere in this outcome gets one FASTA record
      holding its REAL, full sequence sliced directly out of `genome_fasta`
      -- never a placeholder -- so the GFF3's own absolute contig
      coordinates line up against it exactly as written. A contig the GFF3
      references but `genome_fasta` doesn't contain (e.g. a report/genome
      contig-name mismatch) is logged and simply omitted from the companion
      FASTA -- never fabricated.
    * Each gene feature gains a sibling `CDS` feature (same coordinates,
      `Parent` pointing at the `gene` feature's own ID) carrying a
      `translation=` attribute, re-derived via `benchmark.py`'s
      `_extract_translated_gene` -- reused directly rather than
      re-implemented here, since it already fetches from a real genome
      FASTA and splices `GeneEvidence.exons` (when present) via
      `_splice_transcript` with the minus-strand exon-order handling
      already tested there. A gene whose sequence can't be extracted (same
      contig-mismatch case, or any other extraction failure) simply gets no
      `CDS` feature -- logged, never crashing the rest of the write.

    Genes the family expects but which were not found are emitted with
    `present=false` at the locus region's own coordinates (they have no
    coordinates of their own), mirroring sub-project 1's `present` semantics
    so an absence is explicit rather than an empty space. Genes classified
    "not searchable by this method" (short-ORF case) carry
    `not_searchable=true` so an absence caused by a tool limitation is never
    read as a curated negative.

    A fragmented locus's segments can land on different contigs. GFF3's
    Parent/child model assumes a child feature's parent is declared on the
    SAME seqid, so one shared `MAT_locus` feature spanning two contigs would
    make a second-contig gene's `Parent` point at a feature GFF3 tooling
    would consider undeclared on that seqid. Rather than the larger
    redesign of a true spec-correct multi-contig feature model, each segment
    gets its OWN `MAT_locus` feature (scoped to its own contig, `ID=
    <locus>.seg<N>` when there is more than one segment), and every gene
    feature's `Parent` points at whichever segment shares its contig -- so
    every Parent reference is same-contig and valid. A shared
    `locus_group=<locus_id>` attribute on each segment's feature keeps the
    multi-segment grouping visible without a cross-contig Parent claim.

    `##sequence-region` pragmas are collected across the WHOLE outcome
    first and deduplicated per contig (widened to the min start/max end
    seen for that contig across all results), so two separate
    `DetectionResult`s that happen to reference the same contig do not each
    emit their own, duplicate pragma for it.
    """
    contig_extent: dict[str, tuple[int, int]] = {}
    for r in outcome.results:
        segments = r.segments or [None]
        for segment in segments:
            contig = segment.contig if segment else r.contig
            start = segment.start if segment else r.start
            end = segment.end if segment else r.end
            if contig in contig_extent:
                prev_start, prev_end = contig_extent[contig]
                contig_extent[contig] = (min(prev_start, start), max(prev_end, end))
            else:
                contig_extent[contig] = (start, end)

    lines = ["##gff-version 3"]
    for contig, (start, end) in contig_extent.items():
        lines.append(f"##sequence-region {contig} {start} {end}")

    for index, r in enumerate(outcome.results):
        segments = r.segments or [None]
        multi_segment = len(segments) > 1
        locus_id = _locus_id(r, index)
        # Per-contig id for each segment's own MAT_locus feature, so every
        # gene's Parent can point at a feature declared on its own contig.
        segment_ids: dict[str, str] = {}

        for seg_index, segment in enumerate(segments):
            contig = segment.contig if segment else r.contig
            start = segment.start if segment else r.start
            end = segment.end if segment else r.end
            seg_id = f"{locus_id}.seg{seg_index + 1}" if multi_segment else locus_id
            segment_ids[contig] = seg_id

            attrs = (
                f"ID={seg_id};family={_family_label(r.family_key)};confidence={r.confidence}"
                f";idiomorph={r.idiomorph};fragmented={str(r.fragmented).lower()}"
            )
            if multi_segment:
                attrs += f";locus_group={locus_id}"
            if r.ambiguous_with:
                attrs += ";ambiguous_with=" + ",".join(_family_label(k) for k in r.ambiguous_with)
            if r.reference_records:
                attrs += ";reference_records=" + ",".join(r.reference_records)
            lines.append("\t".join([
                contig, "MATPredict", "MAT_locus", str(start), str(end), ".", ".", ".", attrs,
            ]))

        primary_id = segment_ids.get(r.contig, locus_id)

        for gene_index, evidence in enumerate(r.gene_evidence):
            parent_id = segment_ids.get(evidence.contig, primary_id)
            gene_attrs = (
                f"ID={locus_id}.gene{gene_index};Parent={parent_id};Name={evidence.gene_name}"
                f";role={evidence.role};present=true;identity={evidence.identity}"
                f";reference_record={evidence.reference_record_id};method={evidence.method}"
                f";status={evidence.status}"
            )
            if evidence.coverage is not None:
                gene_attrs += f";coverage={evidence.coverage}"
            if evidence.alternate_model is not None:
                alt = evidence.alternate_model
                gene_attrs += (
                    f";alt_method={alt['method']};alt_contig={alt['contig']}"
                    f";alt_start={alt['start']};alt_end={alt['end']};alt_identity={alt['identity']}"
                )
            gene_id = f"{locus_id}.gene{gene_index}"
            lines.append("\t".join([
                evidence.contig, "MATPredict", "gene", str(evidence.start), str(evidence.end),
                ".", evidence.strand or ".", ".", gene_attrs,
            ]))

            if genome_fasta is not None:
                exons = list(evidence.exons) if evidence.exons else None
                protein = _extract_translated_gene(
                    genome_fasta, evidence.contig, evidence.start, evidence.end,
                    evidence.strand, exons,
                )
                if protein:
                    cds_attrs = f"ID={locus_id}.cds{gene_index};Parent={gene_id};translation={protein}"
                    lines.append("\t".join([
                        evidence.contig, "MATPredict", "CDS", str(evidence.start), str(evidence.end),
                        ".", evidence.strand or ".", ".", cds_attrs,
                    ]))
                else:
                    logger.warning(
                        "write_detection_gff3: could not extract/translate gene %r on "
                        "contig %r from genome FASTA %s (contig not found, or extraction "
                        "failed) -- omitting its CDS feature",
                        evidence.gene_name, evidence.contig, genome_fasta,
                    )

        absent_index = len(r.gene_evidence)
        for gene_name in r.genes_missing:
            lines.append("\t".join([
                r.contig, "MATPredict", "gene", str(r.start), str(r.end), ".", ".", ".",
                f"ID={locus_id}.gene{absent_index};Parent={primary_id};Name={gene_name};present=false",
            ]))
            absent_index += 1
        for gene_name in r.genes_not_searchable:
            lines.append("\t".join([
                r.contig, "MATPredict", "gene", str(r.start), str(r.end), ".", ".", ".",
                f"ID={locus_id}.gene{absent_index};Parent={primary_id};Name={gene_name}"
                f";present=false;not_searchable=true",
            ]))
            absent_index += 1

    out_path.write_text("\n".join(lines) + "\n")

    if genome_fasta is not None:
        needed_contigs = set(contig_extent)
        contig_sequences: dict[str, str] = {}
        try:
            with _open_fasta_text(genome_fasta) as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    if record.id in needed_contigs:
                        contig_sequences[record.id] = str(record.seq)
        except OSError:
            logger.warning(
                "write_detection_gff3: could not open genome FASTA %s for the companion "
                "FASTA -- no companion sequence written", genome_fasta,
            )
            contig_sequences = {}

        missing_contigs = needed_contigs - contig_sequences.keys()
        for contig in sorted(missing_contigs):
            logger.warning(
                "write_detection_gff3: contig %r referenced in detection results not "
                "found in genome FASTA %s -- omitting it from the companion FASTA",
                contig, genome_fasta,
            )

        if contig_sequences:
            fasta_lines = []
            for contig in contig_extent:
                seq = contig_sequences.get(contig)
                if seq is None:
                    continue
                fasta_lines.append(f">{contig}")
                fasta_lines.append(seq)
            out_path.with_suffix(".fasta").write_text("\n".join(fasta_lines) + "\n")


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
                "status": e.status,
                "alternate_model": e.alternate_model,
                "exons": [{"start": s, "end": end} for s, end in e.exons] if e.exons else None,
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
