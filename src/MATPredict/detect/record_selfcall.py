"""Does each curated record's own source genome get called at the record?

Review finding F2 (2026-09-28, results/2026-09-28_fable_review/): a record
that cannot call the genome it was built from signals a record or code defect
(the Syncephalastrum racemosum NRRL 2496 record's locus was dropped at the
fraction floor). This module turns a record plus a detection report into a
verdict; `scripts/check_record_selfcall.py` runs it over real genomes.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import yaml

_ASM = re.compile(r"\b(GC[AF]_\d{9}\.\d+)\b")


@dataclass(frozen=True)
class RecordLocation:
    record_id: str
    locus_name: str
    idiomorphs: tuple[str, ...]
    assembly: str | None
    segments: tuple[tuple[str, int, int], ...]  # (contig accession, start, end)


def record_location(metadata_path: Path) -> RecordLocation:
    """Where a record says its locus is. `assembly` is the record's
    `locus.assembly_accession` (added 2026-09-28); for records without it, the
    first assembly accession in its curation or locus text; None if absent."""
    doc = yaml.safe_load(Path(metadata_path).read_text())
    segs = []
    for s in (doc.get("locus", {}).get("core", {}) or {}).get("segments", []) or []:
        acc = (s.get("sequence_source") or {}).get("accession")
        if acc and s.get("start") and s.get("end"):
            segs.append((acc, int(s["start"]), int(s["end"])))
    assembly = (doc.get("locus") or {}).get("assembly_accession")
    if not assembly:
        text = yaml.safe_dump(doc.get("curation", {})) + yaml.safe_dump(doc.get("locus", {}))
        m = _ASM.search(text)
        assembly = m.group(1) if m else None
    mt = doc.get("mating_type", {}) or {}
    return RecordLocation(
        record_id=doc["record_id"], locus_name=mt.get("locus_name", ""),
        idiomorphs=tuple(mt.get("idiomorphs") or ()), assembly=assembly,
        segments=tuple(segs),
    )


def _same_contig(a: str, b: str) -> bool:
    strip = lambda x: x.split(".")[0]
    return a == b or strip(a) == strip(b)


def selfcall_verdict(loc: RecordLocation, report: dict) -> dict:
    """`called` (a detected locus of the record's family overlaps a record
    segment), `withheld` (only a suppressed locus overlaps), or `missed`."""
    def overlaps(entry):
        fam = str(entry.get("family", ""))
        if not fam.endswith(":" + loc.locus_name):
            return False
        return any(_same_contig(entry.get("contig", ""), c) and entry.get("start", 0) <= e
                   and entry.get("end", 0) >= s for c, s, e in loc.segments)

    hits = [d for d in report.get("detected") or [] if overlaps(d)]
    if hits:
        return {"record_id": loc.record_id, "verdict": "called",
                "idiomorph": hits[0].get("idiomorph"), "confidence": hits[0].get("confidence"),
                "expected": list(loc.idiomorphs)}
    held = [d for d in report.get("suppressed_loci") or [] if overlaps(d)]
    if held:
        return {"record_id": loc.record_id, "verdict": "withheld",
                "withheld_reason": held[0].get("withheld_reason"), "expected": list(loc.idiomorphs)}
    return {"record_id": loc.record_id, "verdict": "missed", "expected": list(loc.idiomorphs)}
