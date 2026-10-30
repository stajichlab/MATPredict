#!/usr/bin/env python3
"""Mark loci that two MAT families called at identical coordinates (one locus, two panels), keeping the better-scoring family. Output only: no pipeline or database change.

    dedupe_family_calls.py --loci results/<run>/loci.tsv --out loci_deduped.tsv [--match identical|overlap] [--min-overlap 0.9]

`--match overlap` (opt-in; the default stays `identical`) also groups loci of different families on one contig whose overlap covers at least `--min-overlap` (default 0.9) of the shorter locus, which catches
the same locus called with different boundaries (in v0.7.0 73 of the 80 yeast genomes that `idiomorph_crosscall_flags.py` sends to curation have one locus nested in the other). Loci are grouped by overlap
transitively within a genome and contig; the best-scoring family keeps the group as below.

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


def _overlap_groups(rows, min_overlap):
    """Groups of rows of one genome and contig linked by an overlap of at least `min_overlap` of the shorter locus (union-find; only loci of different families are linked)."""
    by = collections.defaultdict(list)
    for r in rows: by[(r["genome"], r["contig"])].append(r)
    groups = []
    for rs in by.values():
        par = list(range(len(rs)))
        def f(i):
            while par[i] != i: par[i] = par[par[i]]; i = par[i]
            return i
        for i in range(len(rs)):
            for j in range(i + 1, len(rs)):
                a, b = rs[i], rs[j]
                if a["family_called"] == b["family_called"]: continue
                a0, a1, b0, b1 = int(a["start"]), int(a["end"]), int(b["start"]), int(b["end"])
                ov = min(a1, b1) - max(a0, b0) + 1
                if ov > 0 and ov >= min_overlap * min(a1 - a0 + 1, b1 - b0 + 1): par[f(i)] = f(j)
        comp = collections.defaultdict(list)
        for i, r in enumerate(rs): comp[f(i)].append(r)
        groups += comp.values()
    return groups


def mark_duplicates(rows, match="identical", min_overlap=0.9):
    """Add `duplicate_family_call` ('true' or '') and `kept_family` to every row (in place) and return the rows. `match` is `identical` (same contig, start and end) or `overlap` (see the module docstring)."""
    fam = collections.defaultdict(collections.Counter)
    for r in rows: fam[r["genome"]][r["family_called"]] += 1
    if match == "overlap": groups = _overlap_groups(rows, min_overlap)
    else:
        by = collections.defaultdict(list)
        for r in rows: by[(r["genome"], r["contig"], r["start"], r["end"])].append(r)
        groups = by.values()
    for r in rows: r["duplicate_family_call"], r["kept_family"], r["kept_idiomorph"] = "", "", ""
    for g in groups:
        if len({r["family_called"] for r in g}) < 2: continue
        best = min(g, key=lambda r: _score(r, fam[r["genome"]]))
        for r in g:
            r["kept_family"], r["kept_idiomorph"] = best["family_called"], best.get("idiomorph", "")
            if r is not best and r["family_called"] != best["family_called"]: r["duplicate_family_call"] = "true"
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0]); ap.add_argument("--loci", required=True); ap.add_argument("--out", required=True); ap.add_argument("--match", choices=("identical", "overlap"), default="identical"); ap.add_argument("--min-overlap", type=float, default=0.9); a = ap.parse_args()
    rows = list(csv.DictReader(open(a.loci, newline=""), delimiter="\t")); cols = list(rows[0]) + ["duplicate_family_call", "kept_family", "kept_idiomorph"]
    mark_duplicates(rows, a.match, a.min_overlap)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
    dup = [r for r in rows if r["duplicate_family_call"]]
    print(f"{len(rows)} loci; {len(dup)} duplicate family calls in {len({r['genome'] for r in dup})} genomes; by family dropped:", dict(collections.Counter(r["family_called"] for r in dup)))


if __name__ == "__main__":
    main()
