"""Genomes whose deposited name is likely wrong (curator ruling 2026-10-01).

The shared BFD samples.csv and the LCG metadata are never edited. Identities
that MATPredict finds (for example a "Rhizopus" genome whose locus genes are
100% identical to Circinella) are kept in `taxon_overrides.tsv` in the curated
database root, beside `suppress.txt`, so a frozen worktree carries the list it
was launched with.

The file is tab-separated with `#` comments and a header row. Every row must
give a basis, a status (`unconfirmed` = locus genes only, no rDNA/ITS or
genome-wide test; `confirmed`) and the use (for example "exclude from
species-level scoring"). Overrides change how results are scored and
reported; they do not change routing.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

TAXON_OVERRIDES_FILENAME = "taxon_overrides.tsv"
COLUMNS = ("genome_id", "name_in_source", "likely_identity", "identity_rank",
           "basis", "status", "use", "evidence", "ruled")
STATUSES = ("unconfirmed", "confirmed")
RANKS = ("species", "genus", "lineage")


@dataclass(frozen=True)
class TaxonOverride:
    genome_id: str
    name_in_source: str
    likely_identity: str
    identity_rank: str
    basis: str
    status: str
    use: str
    evidence: str
    ruled: str


def load_taxon_overrides(path: Path) -> dict[str, TaxonOverride]:
    """{genome_id: TaxonOverride}. A missing file gives {}; a malformed row
    raises ValueError, because a silent skip would score a genome under the
    wrong name."""
    path = Path(path)
    if path.is_dir():
        path = path / TAXON_OVERRIDES_FILENAME
    if not path.is_file():
        return {}
    lines = [l for l in path.read_text().splitlines() if l.strip() and not l.lstrip().startswith("#")]
    reader = csv.DictReader(lines, delimiter="\t")
    if tuple(reader.fieldnames or ()) != COLUMNS:
        raise ValueError(f"{path}: header must be {COLUMNS}, got {reader.fieldnames}")
    out: dict[str, TaxonOverride] = {}
    for n, row in enumerate(reader, start=2):
        if None in row or any(row[c] is None for c in COLUMNS):
            raise ValueError(f"{path}: row {n} does not have {len(COLUMNS)} fields")
        row = {c: row[c].strip() for c in COLUMNS}
        for c in ("genome_id", "likely_identity", "basis", "use", "ruled"):
            if not row[c]:
                raise ValueError(f"{path}: row {n} has an empty {c}")
        if row["status"] not in STATUSES:
            raise ValueError(f"{path}: row {n} status {row['status']!r} not in {STATUSES}")
        if row["identity_rank"] not in RANKS:
            raise ValueError(f"{path}: row {n} identity_rank {row['identity_rank']!r} not in {RANKS}")
        if row["genome_id"] in out:
            raise ValueError(f"{path}: duplicate genome_id {row['genome_id']}")
        out[row["genome_id"]] = TaxonOverride(**row)
    return out
