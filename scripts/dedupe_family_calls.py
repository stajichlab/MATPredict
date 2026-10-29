#!/usr/bin/env python3
"""Mark loci that two MAT families called at identical coordinates (one locus, two panels), keeping the better-scoring family. Output only: no pipeline or database change.

    dedupe_family_calls.py --loci results/<run>/loci.tsv --out loci_deduped.tsv

A genome routed by phylum fallback is searched with every Ascomycota family, so one yeast locus can be called by, say, MATsc and MTL with the same contig, start and end (v0.7.0: 101 of the 181
genomes with loci from two systems have only such duplicates; the same rule also catches 69 S. pombe loci called by both mat1 and mat2). The family that keeps the locus is the better scoring one: locus class (mat_locus, then partial_locus, then idiomorph_gene_only), then
confidence tier (high, medium, low), then polished genes, then genes found, then the family with more loci in the genome, then the family name. The others get `duplicate_family_call = true`,
`kept_family` and `kept_idiomorph` (idiomorph names differ between families, MATa against A, so a disagreement is not computed; compare them by hand); no row is removed, so counts can be made with or without them. Standard library only.
"""
import argparse, collections, csv

CLASS_RANK = {"mat_locus": 0, "partial_locus": 1, "idiomorph_gene_only": 2}
TIER_RANK = {"high": 0, "medium": 1, "low": 2}


def _score(row, family_loci):
    return (CLASS_RANK.get(row.get("locus_class", ""), 9), TIER_RANK.get(row.get("confidence", ""), 9), -int(row.get("polished_genes") or 0),
            -len([g for g in (row.get("genes_found") or "").split("|") if g]), -family_loci[row["family_called"]], row["family_called"])


def mark_duplicates(rows):
    """Add `duplicate_family_call` ('true' or '') and `kept_family` to every row (in place) and return the rows."""
    fam = collections.defaultdict(collections.Counter)
    for r in rows: fam[r["genome"]][r["family_called"]] += 1
    groups = collections.defaultdict(list)
    for r in rows: groups[(r["genome"], r["contig"], r["start"], r["end"])].append(r)
    for r in rows: r["duplicate_family_call"], r["kept_family"], r["kept_idiomorph"] = "", "", ""
    for g in groups.values():
        if len({r["family_called"] for r in g}) < 2: continue
        best = min(g, key=lambda r: _score(r, fam[r["genome"]]))
        for r in g:
            r["kept_family"], r["kept_idiomorph"] = best["family_called"], best.get("idiomorph", "")
            if r is not best and r["family_called"] != best["family_called"]: r["duplicate_family_call"] = "true"
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0]); ap.add_argument("--loci", required=True); ap.add_argument("--out", required=True); a = ap.parse_args()
    rows = list(csv.DictReader(open(a.loci, newline=""), delimiter="\t")); cols = list(rows[0]) + ["duplicate_family_call", "kept_family", "kept_idiomorph"]
    mark_duplicates(rows)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
    dup = [r for r in rows if r["duplicate_family_call"]]
    print(f"{len(rows)} loci; {len(dup)} duplicate family calls in {len({r['genome'] for r in dup})} genomes; by family dropped:", dict(collections.Counter(r["family_called"] for r in dup)))


if __name__ == "__main__":
    main()
