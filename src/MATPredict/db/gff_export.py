"""Build locus.gff3, locus.gbk, and proteins.faa from an accepted record's metadata."""
from __future__ import annotations

from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord


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


def write_genbank(record: dict, sequences: dict[int, str], out_path: Path) -> None:
    """Write a GenBank record for the core locus, from the same segments/genes data as write_gff3.

    Builds one Bio.SeqRecord per segment, with one "gene" Bio.SeqFeature per present gene at
    that segment (matching write_gff3's present-gene exclusion behavior), and writes it via
    Bio.SeqIO.write in GenBank format.

    `sequences` is accepted for interface symmetry with write_proteins_fasta and to leave room
    for a future CDS/translation feature, but is not currently used: we only have protein
    sequences per gene index here, not the segment's real nucleotide sequence, so there is
    nothing correct to translate a CDS feature against yet.

    LIMITATION: we do not have the actual nucleotide sequence for the segment. The segment's
    nucleotide sequence is therefore a placeholder of "N" characters at the segment's length.
    The gene feature coordinates and qualifiers are accurate; the nucleotide sequence itself is
    a placeholder until sub-project 2's tooling can fetch the real assembly sequence for the
    region.
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

        seq_record = SeqRecord(
            Seq("N" * segment_length),
            id=f"{record['record_id']}.segment{segment_index}",
            name=seq_region[:16] if seq_region else f"segment{segment_index}",
            description=f"{record['record_id']} core locus segment {segment_index} ({seq_region})",
        )
        seq_record.annotations["molecule_type"] = "DNA"

        for gene in genes_by_segment.get(segment_index, []):
            # GFF3/schema coordinates are 1-based fully-closed absolute coordinates.
            # Biopython's FeatureLocation is 0-based half-open, and coordinates here must
            # be relative to the start of this segment's (placeholder) SeqRecord sequence.
            # Conversion: relative_start_0based = gene.start - segment.start
            #             relative_end_halfopen = gene.end - segment.start + 1
            location = FeatureLocation(
                gene["start"] - segment["start"],
                gene["end"] - segment["start"] + 1,
                strand=1 if gene["strand"] == "+" else (-1 if gene["strand"] == "-" else None),
            )
            feature = SeqFeature(location, type="gene", qualifiers={"gene": [gene["name"]]})
            seq_record.features.append(feature)

        seq_records.append(seq_record)

    SeqIO.write(seq_records, out_path, "genbank")
