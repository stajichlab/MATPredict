#!/usr/bin/env python3
"""Add idiomorph and gene labels to a clean-run loci table (output only; the raw columns are kept).

    label_loci_table.py --loci results/<run>/loci.tsv --out loci_labelled.tsv [--summary summary.tsv]
Columns added: idiomorph_label, lineage, system, label_status (labelled / pending / unmapped), genes_labelled (the labels of `genes_found`, in order). The class column `class_` picks the Pezizomycotina system.
Labels come from scripts/idiomorph_label_registry.json through scripts/idiomorph_labels.py. Standard library only.
"""
import argparse, collections, csv, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import idiomorph_labels as L


def label_rows(rows):
    for r in rows:
        fam, cls = r["family_called"], r.get("class_", "")
        lab = L.label_locus(fam, r["idiomorph"], cls)
        r.update(idiomorph_label=lab["label"], lineage=lab["lineage"], system=lab["system"], label_status=lab["status"],
                 genes_labelled="|".join(L.gene_label(g, family=fam, class_=cls) for g in (r.get("genes_found") or "").split("|") if g))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0]); ap.add_argument("--loci", required=True); ap.add_argument("--out", required=True); ap.add_argument("--summary")
    a = ap.parse_args()
    rows = list(csv.DictReader(open(a.loci, newline=""), delimiter="\t")); cols = list(rows[0]) + ["idiomorph_label", "lineage", "system", "label_status", "genes_labelled"]
    label_rows(rows)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
    c = collections.Counter((r["family_called"], r["idiomorph"], r["idiomorph_label"], r["label_status"]) for r in rows)
    if a.summary:
        with open(a.summary, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n"); w.writerow(["family_called", "idiomorph", "label", "status", "loci"])
            for k, n in sorted(c.items()): w.writerow([*k, n])
    st = collections.Counter(r["label_status"] for r in rows)
    print(f"{len(rows)} loci:", dict(st), "unmapped:", sorted({(k[0], k[1]) for k in c if k[3] == "unmapped"}))


if __name__ == "__main__":
    main()
