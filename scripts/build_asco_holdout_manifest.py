#!/usr/bin/env python3
"""Manifest of curated Ascomycota records usable as held-out genomes, by subphylum and mating system.

    python scripts/build_asco_holdout_manifest.py --db db --out testset/holdout/ascomycota_holdout.tsv

One row per curated record under db/Ascomycota. The truth is the record's own `mating_type` (idiomorphs and
`system`, heterothallic or homothallic), and the genome is its `locus.assembly_accession`; a record without an
assembly cannot be a held-out genome and is listed as such. For each radius (genus, family, order) the row says how
many curated records the holdout withholds and how many references of the same family and of the same idiomorph
remain: the idiomorphs of Pezizomycotina MAT1-1 and MAT1-2 are not homologous, so a radius that leaves none of the
record's own idiomorph is `no_reference_idiomorph`, a gap in the database, not a miss.
"""
import argparse
import csv
import re
import sys
from collections import Counter
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from MATPredict.detect.holdout import Radius, records_to_withhold  # noqa: E402

RADII = ("genus", "family", "order")


def subphylum(lineage: str) -> str:
    m = re.search(r"subphylum__([^;]+)", lineage or "")
    return m.group(1) if m else "?"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", type=Path, default=Path("db"))
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    recs = {}
    for p in sorted(a.db.glob("Ascomycota/*/*/metadata.yaml")):
        d = yaml.safe_load(p.read_text())
        mt = d.get("mating_type") or {}
        recs[d["record_id"]] = dict(
            record_id=d["record_id"], subphylum=subphylum((d.get("taxonomy") or {}).get("lineage")),
            order=p.parts[-3], species=(d.get("organism") or {}).get("species"),
            family=mt.get("locus_name"), system=mt.get("system"),
            idiomorphs="|".join(mt.get("idiomorphs") or []), taxid=(d.get("taxonomy") or {}).get("taxid"),
            assembly=(d.get("locus") or {}).get("assembly_accession") or "")
    rows = []
    for rid, r in recs.items():
        for radius in RADII:
            drop = records_to_withhold(a.db, rid, Radius(radius))
            rest = [x for k, x in recs.items() if k not in drop]
            r[f"{radius}_withheld"] = len(drop)
            r[f"{radius}_same_family_left"] = sum(x["family"] == r["family"] for x in rest)
            mine = set(r["idiomorphs"].split("|"))
            r[f"{radius}_same_idiomorph_left"] = sum(
                x["family"] == r["family"] and bool(mine & set(x["idiomorphs"].split("|"))) for x in rest)
        r["testable"] = bool(r["assembly"])
        rows.append(r)
    fields = ["record_id", "subphylum", "order", "species", "family", "system", "idiomorphs", "taxid", "assembly",
              "testable"] + [f"{rd}_{k}" for rd in RADII for k in ("withheld", "same_family_left", "same_idiomorph_left")]
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: (r["subphylum"], r["family"] or "", r["order"], r["record_id"])))
    c = Counter((r["subphylum"], r["family"], r["system"], r["testable"]) for r in rows)
    for k, v in sorted(c.items(), key=lambda kv: str(kv[0])):
        print(v, k)
    print(f"{len(rows)} records, {sum(r['testable'] for r in rows)} testable -> {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
