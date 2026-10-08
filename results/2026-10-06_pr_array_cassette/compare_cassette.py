#!/usr/bin/env python3
"""Before/after check for the report-only receptor_cassette_* fields.

Baseline = pr-array-support 0bdde27 (receptor_array_* fields, no cassette); candidate = pr-array-cassette.
Everything must be identical (receptor_array_* fields included) once the receptor_cassette_* keys and
receptor_arrays_note are dropped. Writes summary.txt, loci.tsv (candidate calls with the columns of
LOCI_ARRAY_COLUMNS) and arrays.tsv; exits 1 on any difference.

Usage: compare_cassette.py OUTDIR BASE_RUNS_DIR CAND_RUNS_DIR
"""
import collections, csv, glob, os, sys

import yaml

from MATPredict.detect.receptor_arrays import LOCI_ARRAY_COLUMNS, loci_columns

CAS = {"receptor_cassette_loci", "receptor_cassette_class", "receptor_cassette_members", "receptor_cassette_max_caax_orfs"}
SKIP_TOP = {"receptor_arrays_note"}


def load(d):
    return {os.path.basename(os.path.dirname(f)): yaml.load(open(f), Loader=yaml.CSafeLoader)
            for f in glob.glob(os.path.join(d, "**", "detection_report.yaml"), recursive=True)}


def strip(e):
    return {k: v for k, v in e.items() if k not in CAS}


def main():
    out, bd, cd = sys.argv[1:4]
    os.makedirs(out, exist_ok=True)
    base, cand = load(bd), load(cd)
    diffs, rows, arr_rows = [], [], []
    n_calls = n_pr = n_pr_cas = 0
    cls_calls, cls_arrays = collections.Counter(), collections.Counter()
    for g in sorted(set(base) | set(cand)):
        b, c = base.get(g), cand.get(g)
        if b is None or c is None:
            diffs.append((g, "report missing")); continue
        for k in sorted((set(b) | set(c)) - SKIP_TOP - {"detected", "receptor_arrays"}):
            if b.get(k) != c.get(k):
                diffs.append((g, f"top-level {k} differs"))
        if [strip(x) for x in b["detected"]] != [strip(x) for x in c["detected"]]:
            diffs.append((g, "detected differs outside receptor_cassette_*"))
        if [strip(x) for x in b["receptor_arrays"]] != [strip(x) for x in c["receptor_arrays"]]:
            diffs.append((g, "receptor_arrays differ outside receptor_cassette_*"))
        if any(CAS & set(x) for x in b["detected"] + b["receptor_arrays"]):
            diffs.append((g, "baseline already has cassette fields"))
        for x in c["detected"]:
            n_calls += 1
            if "receptor_array_id" in x:
                n_pr += 1
                cl = x.get("receptor_cassette_class") or "no_array"
                cls_calls[cl] += 1
                n_pr_cas += cl in ("B", "C")
                rows.append({"genome": g, "family": x["family"], "contig": x["contig"], "start": x["start"], "end": x["end"],
                             "verification": (x.get("verification") or {}).get("status", ""), **loci_columns(x)})
            else:
                if CAS & set(x):
                    diffs.append((g, "cassette field on a non-PR call"))
        for a in c["receptor_arrays"]:
            cls_arrays[a["receptor_cassette_class"]] += 1
            arr_rows.append({"genome": g, **{k: ("|".join(v) if isinstance(v, list) else v) for k, v in a.items()}})
    cols = ["genome", "family", "contig", "start", "end", "verification", *LOCI_ARRAY_COLUMNS]
    with open(os.path.join(out, "loci.tsv"), "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=cols, delimiter="\t"); w.writeheader(); w.writerows(rows)
    if arr_rows:
        with open(os.path.join(out, "arrays.tsv"), "w", newline="") as fo:
            w = csv.DictWriter(fo, fieldnames=list(arr_rows[0]), delimiter="\t"); w.writeheader(); w.writerows(arr_rows)
    with open(os.path.join(out, "summary.txt"), "w") as fo:
        fo.write(f"genomes: base {len(base)}, cand {len(cand)}\n")
        fo.write(f"detected entries (cand): {n_calls}; PR calls with array fields: {n_pr}\n")
        fo.write(f"differences outside receptor_cassette_* (receptor_array_* compared too): {len(diffs)}\n")
        for g, d in diffs:
            fo.write(f"  {g}\t{d}\n")
        fo.write(f"arrays: {len(arr_rows)}; by cassette class: {dict(cls_arrays)}\n")
        fo.write(f"arrays with a cassette: {sum(1 for a in arr_rows if int(a['receptor_cassette_loci']) > 0)}; "
                 f"cassette loci: {sum(int(a['receptor_cassette_loci']) for a in arr_rows)} of "
                 f"{sum(int(a['receptor_array_size']) for a in arr_rows)} loci\n")
        fo.write(f"PR calls by cassette class: {dict(cls_calls)}; PR calls with a cassette: {n_pr_cas}\n")
    print(open(os.path.join(out, "summary.txt")).read())
    sys.exit(1 if diffs else 0)


main()
