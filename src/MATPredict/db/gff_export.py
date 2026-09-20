"""Build locus.gff3, locus.gbk, and proteins.faa from an accepted record's metadata."""
from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

if TYPE_CHECKING:
    from MATPredict.db.ncbi_client import NcbiClient


def write_gff3(record: dict, out_path: Path) -> None:
    """Write a minimal GFF3 for the core locus and its present genes (1-based, fully-closed)."""
    segments = record["locus"]["core"]["segments"]
    lines = ["##gff-version 3"]
    for segment in segments:
        seq_region = segment["sequence_source"]["seq_region"]
        lines.append(f"##sequence-region {seq_region} {segment['start']} {segment['end']}")

    for gene in record["genes"]:
        if not gene.get("present", True):
            continue
        segment = segments[gene["segment_index"]]
        seq_region = segment["sequence_source"]["seq_region"]
        attrs = f"ID={record['record_id']}.gene{gene['gene_index']};Name={gene['name']};role={gene['role']}"
        lines.append(
            "\t".join([
                seq_region, "MATPredict", "gene", str(gene["start"]), str(gene["end"]),
                ".", gene["strand"] or ".", ".", attrs,
            ])
        )
    out_path.write_text("\n".join(lines) + "\n")


def write_proteins_fasta(record: dict, sequences: dict[int, str], out_path: Path) -> None:
    """Write one FASTA entry per present gene with a sequence available.

    Header convention: >{record_id}|gene_index={gene_index}|name={name}|role={role}
    """
    lines = []
    for gene in record["genes"]:
        if not gene.get("present", True) or gene["gene_index"] not in sequences:
            continue
        header = f">{record['record_id']}|gene_index={gene['gene_index']}|name={gene['name']}|role={gene['role']}"
        lines.append(header)
        lines.append(sequences[gene["gene_index"]])
    out_path.write_text("\n".join(lines) + "\n")


def _cds_location(gene: dict, segment: dict, strand: int) -> FeatureLocation | CompoundLocation:
    """Build the CDS feature's location: a real multi-part `CompoundLocation` from
    `gene["exons"]` when present, one `FeatureLocation` per exon, falling back to a
    single-span `FeatureLocation` covering the gene's outer bounds only when `exons`
    is absent/empty (a gene with no known intron structure).

    `gene["exons"]` entries are `{"start": ..., "end": ...}` in this project's 1-based,
    fully-closed genomic coordinates.

    Exon order (verified against real data, not assumed): this project's CURATED-RECORD
    schema stores `gene["exons"]` ALREADY IN TRANSCRIPT ORDER (5'->3') -- ascending
    genomic order for a plus-strand gene, DESCENDING genomic order for a minus-strand
    gene, matching GenBank's own `complement(join(...))` listing order. This is
    documented and empirically verified in two other places in this codebase:
    `db/validate.py`'s `_assemble_transcript` docstring ("exons entries must already be
    listed in transcript order ... for a minus-strand gene this is descending genomic
    coordinate order ... this function does not reorder them") and `db/ncbi_client.py`'s
    real-fixture-verified `CdsStructure` parsing notes ("Biopython's `location.parts`
    already iterates in the REVERSE of the GenBank text listing order -- i.e. already in
    descending-genomic-coordinate (transcript, 5'->3') order ... matching this project's
    exon-list convention"). A real curated record (COX13,
    `db/Ascomycota/Onygenales/199306_rmscc1040_MAT_MAT1-1/metadata.yaml`, minus strand)
    confirms this directly: its `exons` list runs 2057->1362 (descending).

    `Bio.SeqFeature.CompoundLocation.extract()` reverse-complements each part
    individually (when `strand=-1`) and then concatenates the parts IN THE ORDER THEY
    ARE GIVEN in the `parts` list -- verified directly in this task's own test
    (`test_cds_location_minus_strand_compound_location_part_order`) with a synthetic
    two-block minus-strand location. Since `gene["exons"]` is already in the order
    `CompoundLocation` needs (descending for minus strand), the parts below are built
    and used in the SAME order the list stores them -- NO reversal is performed.

    This is a DIFFERENT field, with a DIFFERENT (and incompatible) convention, from
    `detect/pipeline.py`'s `GeneEvidence.exons` (used by `detect/benchmark.py`'s
    `_splice_transcript` for rollout/detection results), which is documented as always
    ASCENDING genomic order regardless of strand and requires an explicit reversal for a
    minus-strand gene. The two conventions must not be conflated -- this function only
    ever operates on the curated-record schema's `gene["exons"]`.
    """
    exons = gene.get("exons")
    if not exons:
        local_start = gene["start"] - segment["start"]
        local_end = gene["end"] - segment["start"] + 1
        return FeatureLocation(local_start, local_end, strand=strand)

    parts = []
    for exon in exons:
        local_start = exon["start"] - segment["start"]
        local_end = exon["end"] - segment["start"] + 1
        parts.append(FeatureLocation(local_start, local_end, strand=strand))
    if len(parts) == 1:
        return parts[0]
    return CompoundLocation(parts)


def _fetch_accession(source: dict) -> str | None:
    """The accession to subrange-fetch this segment's nucleotide sequence with, or
    None when the segment cites nothing fetchable (caller then keeps the all-"N"
    placeholder for that segment -- never a fabricated sequence).

    WHY the two branches differ, and why both are correct:

    A segment's `start`/`end` are, in every case, coordinates ON `seq_region` -- the
    contig/scaffold/chromosome record the curator read them off. What differs between
    the two `sequence_source.type` values is only whether `accession` names that same
    record or something coarser.

    - `insdc_nucleotide`: `accession` IS the nucleotide record the coordinates refer
      to (in the curated DB these segments carry `accession` == `seq_region`, e.g.
      `ACC1.1`/`ACC1.1`). Fetching by `accession` is kept EXACTLY as before, so this
      branch's behavior is unchanged -- it is the one the 56 `insdc_nucleotide`
      segments in the DB already rely on.
    - `assembly`: `accession` is an assembly-level `GCA_`/`GCF_` identifier, which is
      NOT a nucleotide record at all -- `efetch -db nuccore` cannot return sequence
      for it, and `NcbiClient` cannot resolve it. It was never a usable fetch target,
      which is why these segments previously fell through to the all-"N" placeholder.
      Their `seq_region` however IS a real, individually fetchable contig accession,
      and is the record the coordinates are relative to. Verified live on 2026-09-19
      against NCBI efetch: the two `seq_region` values the 5 `assembly` segments in
      the curated DB carry, `NW_026089539.1` and `JAAGWA010000001.1`, both return
      real sequence.

    So `seq_region` is the semantically right fetch target in BOTH branches; the
    `insdc_nucleotide` branch keeps using `accession` only because that is its
    existing, already-verified behavior and the two values agree there.
    """
    source_type = source.get("type")
    if source_type == "insdc_nucleotide":
        return source.get("accession") or None
    if source_type == "assembly":
        return source.get("seq_region") or None
    return None


def write_genbank(
    record: dict, sequences: dict[int, str], out_path: Path, ncbi: "NcbiClient | None" = None
) -> None:
    """Write a GenBank record for the core locus, from the same segments/genes data as write_gff3.

    Builds one Bio.SeqRecord per segment. When `ncbi` is given, each segment's REAL
    nucleotide sequence is fetched via NcbiClient.fetch_nucleotide_sequence -- the
    same mechanism db/validate.py's _independent_translation already uses -- and
    falls back to an all-"N" placeholder ONLY for that segment, on a fetch failure
    or a sequence_source that names nothing fetchable, never fabricating a sequence.
    Which accession each segment is fetched with (and why an `assembly`-typed
    segment is fetched by its `seq_region` contig, not its GCA_/GCF_ accession) is
    documented on `_fetch_accession`. Passing
    no `ncbi` (the default) preserves the prior all-placeholder behavior exactly,
    for any caller/test that doesn't need real sequence.

    Each present gene with a sequence available in `sequences` gets a real `CDS`
    feature (not just `gene`) carrying a `translation` qualifier and `role`/
    `gene_class`/`present_in_idiomorphs` qualifiers copied from the gene's own
    schema fields -- the per-gene attributes clinker's `--colour_map`/
    `--gene_functions` (or pyGenomeViz's `--feature_type2color`) need to color/label
    by MAT-domain biology directly, without a separate manual mapping step.
    """
    segments = record["locus"]["core"]["segments"]
    genes_by_segment: dict[int, list[dict]] = {}
    for gene in record["genes"]:
        if not gene.get("present", True):
            continue
        genes_by_segment.setdefault(gene["segment_index"], []).append(gene)

    seq_records = []
    for segment in segments:
        segment_index = segment["segment_index"]
        segment_length = segment["end"] - segment["start"] + 1
        seq_region = segment["sequence_source"]["seq_region"]
        source = segment.get("sequence_source", {})

        nucleotide_sequence = None
        fetch_accession = _fetch_accession(source)
        if ncbi is not None and fetch_accession:
            try:
                fetched = ncbi.fetch_nucleotide_sequence(
                    fetch_accession, segment["start"], segment["end"], None
                )
            except Exception:
                fetched = None
            # Only trust an actual string sequence -- a caller-supplied test double
            # (e.g. an unconfigured MagicMock in unrelated tests) returning a non-string,
            # non-exception value is not a real fetch result, so it falls back to the
            # placeholder just like a genuine fetch failure would.
            if isinstance(fetched, str) and fetched:
                nucleotide_sequence = fetched
        if not nucleotide_sequence:
            nucleotide_sequence = "N" * segment_length

        seq_record = SeqRecord(
            Seq(nucleotide_sequence),
            id=f"{record['record_id']}.segment{segment_index}",
            name=seq_region[:16] if seq_region else f"segment{segment_index}",
            description=f"{record['record_id']} core locus segment {segment_index} ({seq_region})",
        )
        seq_record.annotations["molecule_type"] = "DNA"

        for gene in genes_by_segment.get(segment_index, []):
            # GFF3/schema coordinates are 1-based fully-closed absolute coordinates.
            # Biopython's FeatureLocation is 0-based half-open, and coordinates here must
            # be relative to the start of this segment's SeqRecord sequence.
            # Conversion: relative_start_0based = gene.start - segment.start
            #             relative_end_halfopen = gene.end - segment.start + 1
            local_start = gene["start"] - segment["start"]
            local_end = gene["end"] - segment["start"] + 1
            strand = 1 if gene.get("strand") != "-" else -1
            location = FeatureLocation(local_start, local_end, strand=strand)

            gene_feature = SeqFeature(location, type="gene", qualifiers={
                "gene": [gene["name"]], "role": [gene["role"]],
            })
            seq_record.features.append(gene_feature)

            translation = sequences.get(gene["gene_index"])
            if translation:
                qualifiers = {
                    "gene": [gene["name"]], "role": [gene["role"]],
                    "translation": [translation],
                }
                if gene.get("gene_class"):
                    qualifiers["gene_class"] = [gene["gene_class"]]
                if gene.get("present_in_idiomorphs"):
                    qualifiers["present_in_idiomorphs"] = [",".join(gene["present_in_idiomorphs"])]
                cds_location = _cds_location(gene, segment, strand)
                cds_feature = SeqFeature(cds_location, type="CDS", qualifiers=qualifiers)
                seq_record.features.append(cds_feature)

        seq_records.append(seq_record)

    SeqIO.write(seq_records, out_path, "genbank")
