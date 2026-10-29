#!/usr/bin/env python3
"""Choose the species of a wave from the survey: reference tiers given, not in earlier waves, and (with --genomes) with at least --min-genomes BFD assemblies of the species in a genomes table (the panel is cut from them).

    wave_species.py --summary ../species_summary.tsv --tiers A A1 B --exclude wave1_species.tsv --genomes asco_genomes.tsv --out species.tsv
Output columns as wave1_species.tsv: species, slug, strains_pooled, reference_tier. Standard library only.
"""
import argparse, csv, re
from collections import Counter


def slug(s):
    return re.sub(r"[^a-z0-9]+", "_", s.lower()).strip("_")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0]); ap.add_argument("--summary", required=True); ap.add_argument("--tiers", nargs="+", default=["A", "A1", "B"])
    ap.add_argument("--exclude", nargs="*", default=[]); ap.add_argument("--genomes"); ap.add_argument("--min-genomes", type=int, default=3); ap.add_argument("--min-strains", type=int, default=4); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rd = lambda p: list(csv.DictReader(open(p, newline=""), delimiter="\t"))
    skip = {r["species"] for p in a.exclude for r in rd(p)}
    have = Counter(r["species"] for r in rd(a.genomes)) if a.genomes else None
    rows = [r for r in rd(a.summary) if r["reference_tier"] in a.tiers and r["species"] not in skip and int(r["strains_pooled"]) >= a.min_strains and (have is None or have[r["species"]] >= a.min_genomes)]
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n"); w.writerow(["species", "slug", "strains_pooled", "reference_tier"])
        for r in rows: w.writerow([r["species"], slug(r["species"]), r["strains_pooled"], r["reference_tier"]])
    print(len(rows), "species;", dict(Counter(r["reference_tier"] for r in rows)))


if __name__ == "__main__":
    main()
