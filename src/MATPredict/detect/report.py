"""GFF3 and validation-style YAML report writers for detection results."""
from __future__ import annotations

from pathlib import Path

import yaml

from MATPredict.detect.pipeline import DetectionResult


def write_detection_gff3(results: list[DetectionResult], out_path: Path) -> None:
    lines = ["##gff-version 3"]
    for r in results:
        attrs = f"ID={r.family_key.phylum}_{r.family_key.locus_name};confidence={r.confidence};idiomorph={r.idiomorph}"
        lines.append("\t".join([r.contig, "MATPredict", "MAT_locus", str(r.start), str(r.end), ".", ".", ".", attrs]))
    out_path.write_text("\n".join(lines) + "\n")


def write_detection_report(results: list[DetectionResult], out_path: Path) -> None:
    doc = [
        {
            "family": f"{r.family_key.phylum}:{r.family_key.locus_name}",
            "contig": r.contig,
            "start": r.start,
            "end": r.end,
            "confidence": r.confidence,
            "idiomorph": r.idiomorph,
            "ambiguous_with": [f"{k.phylum}:{k.locus_name}" for k in r.ambiguous_with],
            "genes_found": r.genes_found,
            "genes_missing": r.genes_missing,
            "genes_not_searchable": r.genes_not_searchable,
            "fragmented": r.fragmented,
        }
        for r in results
    ]
    out_path.write_text(yaml.safe_dump(doc, sort_keys=False))
