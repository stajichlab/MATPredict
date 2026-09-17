"""Concatenate every accepted record's proteins.faa for search.py to use as a diamond/exonerate database."""
from __future__ import annotations

import re
from pathlib import Path

from MATPredict import logger

_HEADER_RE = re.compile(
    r"^>(?P<record_id>[^|]+)\|gene_index=(?P<gene_index>\d+)\|name=(?P<name>[^|]+)\|role=(?P<role>.+)$"
)


def build_reference_fasta(db_root: Path, out_path: Path) -> Path:
    """Rewrite every db/**/proteins.faa header to `record_id|geneN|name` and concatenate.

    Transforms gff_export.write_proteins_fasta's header format
    `>{record_id}|gene_index={gene_index}|name={name}|role={role}`
    into search.py's `_parse_reference_header` expected format
    `>{record_id}|gene{gene_index}|{name}`.
    """
    lines: list[str] = []
    for faa in sorted(db_root.glob("*/*/*/proteins.faa")):
        text = faa.read_text()
        for chunk in text.split(">")[1:]:
            header, _, seq = chunk.partition("\n")
            m = _HEADER_RE.match(">" + header)
            if not m:
                logger.warning(f"Skipping malformed header in {faa}: >{header}")
                continue
            lines.append(f">{m['record_id']}|gene{m['gene_index']}|{m['name']}")
            lines.append(seq.rstrip("\n"))
    out_path.write_text("\n".join(lines) + "\n")
    return out_path
