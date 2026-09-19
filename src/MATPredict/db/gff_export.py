"""Build locus.gff3, locus.gbk, and proteins.faa from an accepted record's metadata."""
from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
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


def write_genbank(
    record: dict, sequences: dict[int, str], out_path: Path, ncbi: "NcbiClient | None" = None
) -> None:
    """Write a GenBank record for the core locus, from the same segments/genes data as write_gff3.

    Builds one Bio.SeqRecord per segment. When `ncbi` is given, each segment's REAL
    nucleotide sequence is fetched via NcbiClient.fetch_nucleotide_sequence -- the
    same mechanism db/validate.py's _independent_translation already uses -- and
    falls back to an all-"N" placeholder ONLY for that segment, on a fetch failure
    or an unfetchable sequence_source.type (e.g. an assembly-level GCA_/GCF_
    accession NcbiClient can't resolve yet), never fabricating a sequence. Passing
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
        if ncbi is not None and source.get("type") == "insdc_nucleotide" and source.get("accession"):
            try:
                fetched = ncbi.fetch_nucleotide_sequence(
                    source["accession"], segment["start"], segment["end"], None
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
                cds_feature = SeqFeature(location, type="CDS", qualifiers=qualifiers)
                seq_record.features.append(cds_feature)

        seq_records.append(seq_record)

    SeqIO.write(seq_records, out_path, "genbank")
