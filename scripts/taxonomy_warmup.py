#!/usr/bin/env python3
"""Taxonomy warm-up and snapshot record for the panel scripts (#324, #342).

`scripts/run_clade_panel.slurm` runs this before the genomes start:

    cut -f2 panel.tsv | sort -u | python scripts/taxonomy_warmup.py            # check every taxid
    python scripts/taxonomy_warmup.py --describe > TAXONOMY.txt                 # what answers lookups

Check: for each taxid on stdin the lineage, phylum and genetic code are looked up once. A taxid fails when the lookup raises, or
returns an empty answer (no lineage and no phylum), which is what a cache miss gives on a node without network (#324: taxid 5158,
seven *C. fimbriata* assemblies silently not searched). With `$MATPREDICT_OFFLINE=1` and `$MATPREDICT_TAXONOMY` set a taxid that
the snapshot does not hold raises `TaxidNotInSnapshot`, so it fails here instead of being skipped later. Every failing taxid is
listed before the exit code (3 when any failed, else 0), so one run shows them all.

`--describe` prints one line per fact (source, table path, snapshot date, sha256) for the run folder, so a campaign records which
snapshot it used. Nothing here changes a call.
"""
from __future__ import annotations

import hashlib
import os
import sys
import time


def check(taxids, sleep: float = 0.0, err=sys.stderr) -> list[int]:
    """The taxids that fail the lookup (raise, or no lineage and no phylum); each is reported on `err`."""
    from MATPredict.db.taxonomy import default_genetic_code, default_lineage_phylum_name, default_lineage_taxids

    bad: list[int] = []
    for t in taxids:
        try:
            lineage = default_lineage_taxids(t)
            phylum = default_lineage_phylum_name(t)
            default_genetic_code(t)
            if not lineage and not phylum:
                bad.append(t)
                print(f"taxonomy warm-up: taxid {t} returned no lineage and no phylum (empty answer)", file=err)
        except Exception as exc:  # noqa: BLE001 -- reported, never swallowed
            bad.append(t)
            print(f"taxonomy warm-up: taxid {t} failed: {exc}", file=err)
        if sleep:
            time.sleep(sleep)
    return bad


def describe() -> str:
    """Lines for TAXONOMY.txt: the lookup source, the table path and snapshot date, the table's sha256."""
    from MATPredict.db import local_taxonomy as lt

    lines = [f"source: {lt.source_label()}", f"offline: {lt.offline()}"]
    path = lt.configured_path()
    if path is not None:
        h = hashlib.sha256()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        lines += [f"table: {path}", f"snapshot: {lt.get_local().snapshot}", f"sha256: {h.hexdigest()}"]
    else:
        lines.append("table: none ($MATPREDICT_TAXONOMY unset or missing)")
    return "\n".join(lines) + "\n"


def main(argv: list[str]) -> int:
    if "--describe" in argv:
        sys.stdout.write(describe())
        return 0
    taxids = [int(line) for line in sys.stdin if line.strip()]
    bad = check(taxids, sleep=float(os.environ.get("TAXONOMY_WARMUP_SLEEP", "0.35")))
    print(f"taxonomy warm-up done, {len(bad)} failed" + (f": {' '.join(map(str, bad))}" if bad else ""), file=sys.stderr)
    return 3 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
