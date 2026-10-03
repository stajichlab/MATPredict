#!/usr/bin/env python3
"""Check that each curated record's source genome is called at the record.

Review finding F2 (2026-09-28). For every record whose assembly accession is
in the record and whose genome is in the BFD library, find an existing
detection report for that genome (a results tree) and report the verdict:
called / withheld (with reason) / missed / no_report, or not_possible when
the record has no source assembly and no verified same-strain assembly.

Usage:
  check_record_selfcall.py --db DB_ROOT --reports DIR [DIR ...] [--out TSV]
Reports are found as <DIR>/**/<ASMID>*/detection_report.yaml, or inside
<DIR>/**/reports*.tar.zst archives. It does not run detection; point it at a
fresh run of the code under test.
"""
import argparse
import csv
import io
import subprocess
import sys
import tarfile
from pathlib import Path

import yaml

from MATPredict.detect.record_selfcall import record_location, selfcall_verdict


def _find(asm, dirs):
    for d in dirs:
        for p in Path(d).rglob(f"{asm}*/detection_report.yaml"):
            return yaml.safe_load(p.read_text())
        for tz in Path(d).rglob("reports*.tar.zst"):
            data = subprocess.run(["zstd", "-dc", str(tz)], capture_output=True).stdout
            with tarfile.open(fileobj=io.BytesIO(data)) as t:
                for m in t.getmembers():
                    if f"/{asm}" in "/" + m.name and m.name.endswith("detection_report.yaml"):
                        return yaml.safe_load(t.extractfile(m).read())
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True)
    ap.add_argument("--reports", nargs="+", required=True)
    ap.add_argument("--out")
    a = ap.parse_args()
    rows = []
    for md in sorted(Path(a.db).glob("*/*/*/metadata.yaml")):
        loc = record_location(md)
        if not loc.assembly:
            # Curator ruling 2026-10-01: a record with no source or verified
            # same-strain assembly cannot be self-checked. That is not a failure.
            rows.append({"record_id": loc.record_id, "assembly": "", "verdict": "not_possible"})
            continue
        rep = _find(loc.assembly, a.reports)
        if rep is None:
            rows.append({"record_id": loc.record_id, "assembly": loc.assembly, "verdict": "no_report"})
            continue
        v = selfcall_verdict(loc, rep)
        v["assembly"] = loc.assembly
        rows.append(v)
    keys = ["record_id", "assembly", "verdict", "expected", "idiomorph", "confidence", "withheld_reason"]
    w = csv.DictWriter(open(a.out, "w") if a.out else sys.stdout, fieldnames=keys,
                       delimiter="\t", extrasaction="ignore")
    w.writeheader()
    w.writerows(rows)


if __name__ == "__main__":
    main()
