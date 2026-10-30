#!/usr/bin/env python3
"""Add idiomorph and gene labels to a clean-run loci table (output only; the raw columns are kept).

    label_loci_table.py --loci results/<run>/loci.tsv --out loci_labelled.tsv [--summary summary.tsv]
Columns added: idiomorph_label, lineage, system, label_mapping (labelled / pending / unmapped: did the registry map the family and name), label_status (curator, 2026-10-29, #211: ok / family_unconfirmed /
cross_call / no_reference / similarity_only / undetermined, one value per locus, first that applies), label_reasons (every reason that applies, `;`-separated), registry_version (`<date>+<hash>` of the registry
that made the labels, #210), genes_labelled (the labels of `genes_found`, in order). The class column `class_` picks the Pezizomycotina system.
Cross-calls are read from the table itself: a locus overlapping (0.9 of the shorter) a locus of another family on one contig, or in a genome with loci of two Saccharomycotina systems, or a non-empty `cross_call`
column; `family_unconfirmed` from a `family_status` column when the table has one (reports of `--cross-family-policy on`). `label_status` was the registry-mapping column before 2026-10-29; that is now `label_mapping`.
Labels come from scripts/idiomorph_label_registry.json through scripts/idiomorph_labels.py. Standard library only.
"""
import argparse, collections, csv, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import idiomorph_labels as L


MIN_OVERLAP = 0.9      # of the shorter locus, the rule of detect/family_status.py
# The yeast systems a genome can be called by twice (scripts/idiomorph_crosscall_flags.py FAMS, plus MATpk): 181 genomes in v0.7.0.
YEAST_SYSTEMS = ("Ascomycota:MATsc", "Ascomycota:MTL", "Ascomycota:MATog", "Ascomycota:MATkp", "Ascomycota:MATpk")


def _overlap(a, b):
    if a["contig"] != b["contig"]: return 0.0
    a0, a1, b0, b1 = int(a["start"]), int(a["end"]), int(b["start"]), int(b["end"])
    shorter = min(a1 - a0, b1 - b0) + 1
    return max(0, min(a1, b1) - max(a0, b0) + 1) / shorter if shorter > 0 else 0.0


def cross_called(rows):
    """Indexes of the loci of `rows` called by two families: overlapping a locus of another family, or in a genome whose loci come from two or more Saccharomycotina systems, or a `cross_call` value."""
    sacch = set(YEAST_SYSTEMS)
    by = collections.defaultdict(list)
    for i, r in enumerate(rows): by[r.get("genome", "")].append(i)
    out = set()
    for idx in by.values():
        yeast = {rows[i]["family_called"] for i in idx} & sacch
        if len(yeast) > 1: out.update(i for i in idx if rows[i]["family_called"] in yeast)
        for i in idx:
            if rows[i].get("cross_call"): out.add(i)
            for j in idx:
                if j > i and rows[i]["family_called"] != rows[j]["family_called"] and _overlap(rows[i], rows[j]) >= MIN_OVERLAP: out.update((i, j))
    return out


def label_rows(rows):
    crossed = cross_called(rows)
    version = L.registry_version()
    for n, r in enumerate(rows):
        fam, cls = r["family_called"], r.get("class_", "")
        lab = L.label_locus(fam, r["idiomorph"], cls)
        st = L.locus_status(fam, r["idiomorph"], cls, family_unconfirmed=(r.get("family_status") == "family_unconfirmed"), cross_call=n in crossed)
        r.update(idiomorph_label=lab["label"], lineage=lab["lineage"], system=lab["system"], label_mapping=lab["status"], label_status=st["label_status"],
                 label_reasons=st["label_reasons"], registry_version=version,
                 genes_labelled="|".join(L.gene_label(g, family=fam, class_=cls) for g in (r.get("genes_found") or "").split("|") if g))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0]); ap.add_argument("--loci", required=True); ap.add_argument("--out", required=True); ap.add_argument("--summary")
    a = ap.parse_args()
    rows = list(csv.DictReader(open(a.loci, newline=""), delimiter="\t")); cols = list(rows[0]) + ["idiomorph_label", "lineage", "system", "label_mapping", "label_status", "label_reasons", "registry_version", "genes_labelled"]
    label_rows(rows)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
    c = collections.Counter((r["family_called"], r["idiomorph"], r["idiomorph_label"], r["label_mapping"]) for r in rows)
    if a.summary:
        with open(a.summary, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n"); w.writerow(["family_called", "idiomorph", "label", "status", "loci"])
            for k, n in sorted(c.items()): w.writerow([*k, n])
    st = collections.Counter(r["label_mapping"] for r in rows)
    print(f"{len(rows)} loci:", dict(st), "label_status:", dict(collections.Counter(r["label_status"] for r in rows)), "registry", L.registry_version(), "unmapped:", sorted({(k[0], k[1]) for k in c if k[3] == "unmapped"}))


if __name__ == "__main__":
    main()
