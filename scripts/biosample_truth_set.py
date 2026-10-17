#!/usr/bin/env python3
"""Turn harvested BioSample mating-type attributes into a truth set and compare it with a campaign.

    biosample_truth_set.py --harvest DIR --genomes genomes.tsv --loci loci.tsv --db db --out-dir OUT

Reads biosample_mating.tsv and biosample_map.tsv from `harvest_biosample_mating_types.py`. Normalises the
submitter's value (see `normalise`), drops genomes whose BioSample value is not a mating type, and writes
`truth_set.tsv` (one row per genome: normalised class, raw values, organism, strain, whether the genome or its
strain is already a curated record) and `agreement.tsv` (the campaign's idiomorph for the same genome against the
stated one). The stated mating type is the submitter's, so disagreements are to be read, not assumed to be errors.
"""
import argparse
import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import yaml

#: attribute names that are really "sample material type" (MIxS samp_mat_type) and say nothing about mating
NOT_MATING = re.compile(r"^samp_", re.I)


def normalise(value: str) -> str | None:
    """MAT1-1, MAT1-2, MATa, MATalpha, homothallic (both idiomorphs or self-fertile), or None."""
    v = value.strip().lower().replace("_", "-")
    v = re.sub(r"\s+", "", v)
    if re.fullmatch(r"(mat)?-?1-1(-1)?|mat-?1|mat1-1|matα?1", v) and v not in ("mat-a",):
        return "MAT1-1" if v not in ("mata", ) else None
    if re.fullmatch(r"(mat)?-?1-2(-1)?|mat-?2|mat1-2", v):
        return "MAT1-2"
    if v in ("mata", "mat-a", "a", "mat(a)"):
        return "MATa"
    if v in ("matalpha", "mat-alpha", "alpha", "matα", "α"):
        return "MATalpha"
    if v in ("mat1+2", "mat1-1/mat1-2", "homothallic", "self-fertile", "a/alpha", "a/a", "alpha/alpha",
             "a/alphadiploid", "diploid"):
        return "homothallic_or_diploid"
    return None


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--harvest", type=Path, required=True)
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--loci", required=True)
    ap.add_argument("--db", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    a = ap.parse_args(argv)
    a.out_dir.mkdir(parents=True, exist_ok=True)
    amap = {r["genome"]: r for r in csv.DictReader(open(a.harvest / "biosample_map.tsv"), delimiter="\t")}
    raw = defaultdict(list)
    for r in csv.DictReader(open(a.harvest / "biosample_mating.tsv"), delimiter="\t"):
        if not NOT_MATING.match(r["attribute"]):
            raw[r["genome"]].append(r["value"])
    info = {r["genome"]: r for r in csv.DictReader(open(a.genomes), delimiter="\t")}
    curated = set()
    for p in a.db.glob("*/*/*/metadata.yaml"):
        d = yaml.safe_load(p.read_text())
        acc = (d.get("locus") or {}).get("assembly_accession")
        if acc:
            curated.add(acc)
        curated.add(((d.get("organism") or {}).get("species"), ((d.get("organism") or {}).get("strain") or {}).get("name")))
    loci = defaultdict(list)
    for r in csv.DictReader(open(a.loci), delimiter="\t"):
        loci[r["genome"]].append(r)
    rows, dropped = [], Counter()
    for g, vals in raw.items():
        classes = {normalise(v) for v in vals} - {None}
        if not classes:
            dropped["value is not a mating type"] += 1
            continue
        if len(classes) > 1:
            dropped["conflicting values"] += 1
            continue
        m, i = amap[g], info.get(g, {})
        acc = m["accession"]
        rows.append(dict(genome=g, accession=acc, biosample=m["biosample"], organism=m["organism"], strain=m["strain"],
                         class_=i.get("class_", ""), order=i.get("order", ""), stated=classes.pop(),
                         raw="|".join(sorted(set(vals))),
                         curated_assembly=acc in curated,
                         curated_strain=(m["organism"], m["strain"]) in curated))
    fields = list(rows[0].keys())
    with open(a.out_dir / "truth_set.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader(); w.writerows(sorted(rows, key=lambda r: (r["class_"], r["organism"], r["genome"])))
    # agreement with the campaign
    ag, c = [], Counter()
    for r in rows:
        calls = [l for l in loci.get(r["genome"], []) if l["family_called"].startswith("Ascomycota:")]
        labels = sorted({l["idiomorph"] for l in calls}) or ["none"]
        fams = sorted({l["family_called"] for l in calls})
        stated = r["stated"]
        if stated in ("MAT1-1", "MAT1-2"):
            verdict = ("agree" if labels == [stated] else "none" if labels == ["none"]
                       else "both" if set(labels) >= {"MAT1-1", "MAT1-2"} else "undetermined" if "undetermined" in labels
                       and len(labels) == 1 else "DISAGREE")
        else:
            verdict = "other_system"
        ag.append(dict(r, labels="|".join(labels), families="|".join(fams), verdict=verdict))
        c[(r["class_"] or "?", stated, verdict)] += 1
    with open(a.out_dir / "agreement.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(ag[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(ag)
    print(f"{len(raw)} genomes with a mating attribute; dropped {dict(dropped)}; {len(rows)} in the truth set")
    print(Counter(r["stated"] for r in rows))
    for k, v in sorted(c.items(), key=lambda kv: (kv[0][0], -kv[1])):
        print(v, k)
    return 0


if __name__ == "__main__":
    sys.exit(main())
