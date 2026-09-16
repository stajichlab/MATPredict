"""Walk db/, rebuild the DuckDB query cache from metadata.yaml files."""
from __future__ import annotations

from pathlib import Path

import duckdb
import yaml

_SCHEMA_SQL_PATH = Path(__file__).resolve().parents[3] / "db" / "_schema" / "duckdb_schema.sql"


def _find_metadata_files(db_root: Path) -> list[Path]:
    return sorted(db_root.glob("**/metadata.yaml"))


def _phylum_from_path(db_root: Path, metadata_path: Path) -> str:
    relative = metadata_path.relative_to(db_root)
    parts = relative.parts
    return parts[1] if parts[0] == "candidates" else parts[0]


def _generated_file_paths(record: dict, record_dir: Path) -> tuple[str | None, str | None, str | None]:
    """Compute gff3/gbk/proteins-fasta paths for a record, or (None, None, None).

    Only accepted records can have these generated files (candidates never
    do -- accept happens before generation, per the curation workflow). Even
    for an accepted record, each path is only reported if the file actually
    exists on disk.
    """
    if record["validation"]["status"] != "accepted":
        return (None, None, None)

    gff3_path = record_dir / "locus.gff3"
    gbk_path = record_dir / "locus.gbk"
    proteins_fasta_path = record_dir / "proteins.faa"

    return (
        str(gff3_path) if gff3_path.exists() else None,
        str(gbk_path) if gbk_path.exists() else None,
        str(proteins_fasta_path) if proteins_fasta_path.exists() else None,
    )


def build(db_root: Path, out_path: Path) -> None:
    """Rebuild the DuckDB cache at out_path from every metadata.yaml under db_root."""
    if out_path.exists():
        out_path.unlink()
    con = duckdb.connect(str(out_path))
    con.execute(_SCHEMA_SQL_PATH.read_text())

    for metadata_path in _find_metadata_files(db_root):
        record = yaml.safe_load(metadata_path.read_text())
        phylum = _phylum_from_path(db_root, metadata_path)
        order_or_family = metadata_path.relative_to(db_root).parts[-2]
        record_dir = metadata_path.parent

        con.execute(
            """
            INSERT OR IGNORE INTO organism (taxid, species, lineage, lineage_resolved_date)
            VALUES (?, ?, ?, ?)
            """,
            [record["taxonomy"]["taxid"], record["organism"]["species"], record["taxonomy"]["lineage"],
             record["taxonomy"]["lineage_resolved_date"]],
        )

        strain = record["organism"]["strain"]
        locus = record["locus"]
        core = locus.get("core", {})
        segments = core.get("segments", [])
        validation = record["validation"]
        gff3_path, gbk_path, proteins_fasta_path = _generated_file_paths(record, record_dir)
        con.execute(
            """
            INSERT INTO locus_record (
                record_id, taxid, strain_name, strain_known, locus_name, idiomorph_key,
                mating_system, phylum, order_or_family, record_version, coordinate_provenance,
                excluded_from_coordinate_benchmark, completeness, reference_orientation, definition_note,
                validation_status, rejection_reason, accession_resolved, accession_resolved_version,
                sequence_match_status, taxonomy_current, proposed_by, proposal_dedupe_key, reviewed_by,
                reviewed_date, gff3_path, gbk_path, proteins_fasta_path, metadata_path
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            [
                record["record_id"], record["taxonomy"]["taxid"], strain["name"], strain["known"],
                record["mating_type"]["locus_name"], "+".join(record["mating_type"]["idiomorphs"]),
                record["mating_type"]["system"], phylum, order_or_family, record["record_version"],
                locus["coordinate_provenance"], locus["excluded_from_coordinate_benchmark"],
                core.get("completeness"), core.get("reference_orientation"), core.get("definition_note"),
                validation["status"], validation.get("rejection_reason"),
                validation.get("accession_resolved"), validation.get("accession_resolved_version"),
                (validation.get("sequence_match") or {}).get("status"), validation.get("taxonomy_current"),
                record["curation"]["proposed_by"], record["curation"].get("proposal_dedupe_key"),
                record["curation"].get("reviewed_by"), record["curation"].get("reviewed_date"),
                gff3_path, gbk_path, proteins_fasta_path, str(metadata_path),
            ],
        )

        for segment in segments:
            source = segment["sequence_source"]
            con.execute(
                """
                INSERT INTO locus_segment (record_id, segment_index, sequence_source_type, accession,
                    seq_region, start_pos, end_pos, contig_edge_distance, sequence_checksum)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                [record["record_id"], segment["segment_index"], source["type"], source.get("accession"),
                 source.get("seq_region"), segment["start"], segment["end"],
                 segment.get("contig_edge_distance"), segment.get("sequence_checksum")],
            )

        for gene in record["genes"]:
            con.execute(
                """
                INSERT INTO locus_gene (record_id, gene_index, segment_index, name, protein_accession,
                    role, present, locus_tag, start_pos, end_pos, strand, order_in_locus)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                [record["record_id"], gene["gene_index"], gene.get("segment_index"), gene["name"],
                 gene.get("protein_accession"), gene["role"], gene["present"], gene.get("locus_tag"),
                 gene.get("start"), gene.get("end"), gene.get("strand"), gene.get("order_in_locus")],
            )

        for claim, claim_data in record["evidence"].items():
            con.execute(
                "INSERT INTO evidence_claim (record_id, claim, tier, experimental_method) VALUES (?, ?, ?, ?)",
                [record["record_id"], claim, claim_data["tier"], claim_data.get("experimental_method")],
            )
            for i, citation in enumerate(claim_data.get("citations", [])):
                citation_id = f"{record['record_id']}|{claim}|{i}"
                con.execute(
                    "INSERT INTO citation (citation_id, record_id, claim, pmid, doi) VALUES (?, ?, ?, ?, ?)",
                    [citation_id, record["record_id"], claim, citation.get("pmid"), citation.get("doi")],
                )

    con.close()
