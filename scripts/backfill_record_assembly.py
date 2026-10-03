#!/usr/bin/env python3
"""Backfill `locus.assembly_accession` in curated records from NCBI.

Added 2026-09-28 so `scripts/check_record_selfcall.py` can place every record
(review finding F2). For each record's core-segment sequence accession, NCBI
E-utilities `elink` (nuccore -> assembly) finds the assemblies that contain
it, and `esummary` on the ASSEMBLY database gives each one's accession. (Do
not use nuccore esummary's `assemblyacc`: it returns the INSDC sequence
accession, not the assembly -- measured 0/18 vs 17/18.)

A value is written only when it is unambiguous: every segment links to the
same single GenBank (GCA) or RefSeq (GCF) base accession. A GCA/GCF pair for
the same assembly is resolved by the segment's own type (an NC_/NW_/NZ_
RefSeq sequence gives the GCF, an INSDC one the GCA). A locus-specific
deposit that belongs to no assembly stays null. A record that already has the
field is left alone. Only the one line is inserted; nothing else in the
record changes.

Usage: backfill_record_assembly.py --db DB_ROOT [--write] [--out TSV]
Without --write it reports what it would do.
"""
import argparse
import csv
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

import yaml

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
REFSEQ_PREFIX = re.compile(r"^(NC|NW|NZ|NT)_")


def _get(url, params, tries=4):
    q = urllib.parse.urlencode({**params, "tool": "MATPredict"})
    for i in range(tries):
        try:
            with urllib.request.urlopen(f"{EUTILS}/{url}?{q}", timeout=60) as r:
                return r.read()
        except Exception:  # noqa: BLE001 - retried, then raised
            if i == tries - 1:
                raise
            time.sleep(2 * (i + 1))
    return b""


def assemblies_for(seq_acc):
    """Assembly accessions that contain `seq_acc`, via elink + assembly esummary."""
    root = ET.fromstring(_get("elink.fcgi", {"dbfrom": "nuccore", "db": "assembly",
                                             "id": seq_acc}))
    uids = [e.text for e in root.findall(".//LinkSetDb/Link/Id")]
    time.sleep(0.4)
    if not uids:
        return []
    root = ET.fromstring(_get("esummary.fcgi", {"db": "assembly", "id": ",".join(uids)}))
    time.sleep(0.4)
    return sorted({e.text for e in root.iter("AssemblyAccession") if e.text})


def choose(seq_acc, found):
    """One accession for this segment, or None if ambiguous or absent."""
    want = "GCF_" if REFSEQ_PREFIX.match(seq_acc) else "GCA_"
    pick = [a for a in found if a.startswith(want)] or found
    bases = {a.split(".")[0] for a in pick}
    if len(bases) != 1:
        return None
    return max(pick, key=lambda a: int(a.split(".")[1]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True)
    ap.add_argument("--write", action="store_true")
    ap.add_argument("--out")
    a = ap.parse_args()
    rows = []
    for md in sorted(Path(a.db).glob("*/*/*/metadata.yaml")):
        doc = yaml.safe_load(md.read_text())
        locus = doc.get("locus") or {}
        rid = doc.get("record_id")
        if locus.get("assembly_accession"):
            rows.append({"record_id": rid, "status": "already_set",
                         "assembly": locus["assembly_accession"]})
            continue
        segs = (locus.get("core") or {}).get("segments") or []
        accs = sorted({(s.get("sequence_source") or {}).get("accession") for s in segs} - {None})
        if not accs:
            rows.append({"record_id": rid, "status": "no_segment_accession"})
            continue
        picks, detail = set(), []
        for acc in accs:
            found = assemblies_for(acc)
            c = choose(acc, found)
            detail.append(f"{acc}->{','.join(found) or '-'}")
            picks.add(c)
        if picks == {None}:
            status, asm = "no_assembly", None
        elif None in picks or len(picks) != 1:
            status, asm = "ambiguous", None
        else:
            status, asm = "resolved", picks.pop()
        rows.append({"record_id": rid, "status": status, "assembly": asm or "",
                     "links": "; ".join(detail)})
        if a.write and asm:
            text = md.read_text()
            lines = text.splitlines(keepends=True)
            idx = next(i for i, ln in enumerate(lines) if ln.rstrip("\n") == "locus:")
            lines.insert(idx + 1, f"  assembly_accession: {asm}\n")
            md.write_text("".join(lines))
    w = csv.DictWriter(open(a.out, "w") if a.out else sys.stdout,
                       fieldnames=["record_id", "status", "assembly", "links"],
                       delimiter="\t", extrasaction="ignore")
    w.writeheader()
    w.writerows(rows)


if __name__ == "__main__":
    main()
