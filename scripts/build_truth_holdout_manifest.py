#!/usr/bin/env python3
"""Held-out manifest for genomes with a stated mating type (BioSample truth set), in the runner's format.

    build_truth_holdout_manifest.py --truth truth_set.tsv --summary assembly_summary_GCA.txt [assembly_summary_GCF.txt] \
        --db db --out manifest.tsv

One row per Pezizomycotina genome whose stated mating type is MAT1-1 or MAT1-2. The genome is held out already (it is
not a curated record); the radius removes its relatives from the references: `genus_drop` lists the curated records of
the genome's genus and `order_drop` those of its order (the order is the campaign's `order` column; a genome without
one has no order run). `run_asco_holdout.py run` uses a row's drop list when it has one. The reference counts say how
many curated Ascomycota `MAT` records, and how many of the stated idiomorph, remain at each radius.
"""
import argparse
import csv
import json
import re
import sys
from pathlib import Path

import yaml

RADII = ("genus", "order")


def rank_of(lineage: str, prefix: str) -> str | None:
    m = re.search(rf"(?:^|;){re.escape(prefix)}([^;]+)", lineage or "")
    return m.group(1) if m else None


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--summary", nargs="+", required=True, help="NCBI assembly_summary files (taxid lookup)")
    ap.add_argument("--db", type=Path, required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    taxid = {}
    for path in a.summary:
        for line in open(path):
            if line.startswith("#") or not line.strip():
                continue
            f = line.split("\t")
            taxid[f[0]] = f[5]
    recs = []
    for p in sorted(a.db.glob("*/*/*/metadata.yaml")):
        d = yaml.safe_load(p.read_text())
        lin = (d.get("taxonomy") or {}).get("lineage") or ""
        mt = d.get("mating_type") or {}
        recs.append(dict(id=d["record_id"], phylum=p.parts[-4], family=mt.get("locus_name"),
                         idiomorphs=set(mt.get("idiomorphs") or []), genus=rank_of(lin, "g__"),
                         order=rank_of(lin, "o__")))
    mat = [r for r in recs if r["phylum"] == "Ascomycota" and r["family"] == "MAT"]
    rows = []
    for t in csv.DictReader(open(a.truth), delimiter="\t"):
        if t["stated"] not in ("MAT1-1", "MAT1-2") or t["class_"] in ("Saccharomycetes", "Pichiomycetes", "Dipodascomycetes"):
            continue
        genus = (t["organism"].split() or [""])[0]
        row = dict(record_id=t["genome"], phylum="Ascomycota", subphylum="Pezizomycotina", order=t["order"],
                   species=t["organism"], family="MAT", system="heterothallic", idiomorphs=t["stated"],
                   taxid=taxid.get(t["accession"], ""), assembly=t["accession"], truth_spans="{}", testable=True)
        for radius, key, value in (("genus", "genus", genus), ("order", "order", t["order"])):
            drop = [r["id"] for r in recs if value and r[key] == value]
            left = [r for r in mat if r["id"] not in drop]
            row[f"{radius}_withheld"] = len(drop)
            row[f"{radius}_drop"] = ",".join(sorted(drop)) if value else ""
            row[f"{radius}_same_family_left"] = len(left)
            row[f"{radius}_same_idiomorph_left"] = sum(t["stated"] in r["idiomorphs"] for r in left)
        row["order_run"] = bool(t["order"])
        rows.append(row)
    fields = list(rows[0].keys())
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader(); w.writerows(rows)
    print(f"{len(rows)} genomes ({sum(1 for r in rows if not r['taxid'])} without a taxid, "
          f"{sum(1 for r in rows if not r['order_run'])} without an order) -> {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
