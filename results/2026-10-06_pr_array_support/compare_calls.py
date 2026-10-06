#!/usr/bin/env python3
"""Before/after check for the report-only receptor arrays (branch pr-array-support).

Reads two sets of detection reports (baseline = origin/main, candidate = the branch),
checks that every `detected` entry and every `suppressed_loci` entry is identical once the
array fields are dropped, and writes loci.tsv (candidate calls with the five array columns
added through `receptor_arrays.loci_columns`) plus an array_support summary.

Usage: compare_calls.py OUTDIR BASE_RUNS_DIR CAND_RUNS_DIR [BASE_RUNS_DIR CAND_RUNS_DIR ...]
Each RUNS_DIR holds <genome>/detection_report.yaml (or <chunk>/runs/<genome>/...).
"""
import collections, csv, glob, os, sys

import yaml

from MATPredict.detect.receptor_arrays import LOCI_ARRAY_COLUMNS, loci_columns

ARRAY_KEYS = {"array_id", "array_size", "array_members", "array_support", "array_support_reasons"}
TOP_NEW = {"receptor_arrays", "receptor_arrays_note"}


def load(runs_dir):
    out = {}
    for f in glob.glob(os.path.join(runs_dir, "**", "detection_report.yaml"), recursive=True):
        out[os.path.basename(os.path.dirname(f))] = yaml.load(open(f), Loader=yaml.CSafeLoader)
    return out


def strip(entry):
    return {k: v for k, v in entry.items() if k not in ARRAY_KEYS}


def main():
    outdir, pairs = sys.argv[1], sys.argv[2:]
    os.makedirs(outdir, exist_ok=True)
    base, cand = {}, {}
    for b, c in zip(pairs[0::2], pairs[1::2]):
        base.update(load(b)); cand.update(load(c))
    genomes = sorted(set(base) | set(cand))
    diffs, rows, nb, nc = [], [], 0, 0
    sup = collections.Counter()
    arr_rows = []
    for g in genomes:
        b, c = base.get(g), cand.get(g)
        if b is None or c is None:
            diffs.append((g, "report missing on " + ("base" if b is None else "cand")))
            continue
        for key in sorted((set(b) | set(c)) - TOP_NEW):
            if key in ("detected",):
                continue
            if b.get(key) != c.get(key):
                diffs.append((g, f"top-level {key} differs"))
        db, dc = b.get("detected") or [], c.get("detected") or []
        nb += len(db); nc += len(dc)
        if [strip(x) for x in db] != [strip(x) for x in dc]:
            diffs.append((g, "detected entries differ outside the array fields"))
        if any(ARRAY_KEYS & set(x) for x in db):
            diffs.append((g, "baseline already has array fields"))
        for x in dc:
            row = {"genome": g, "family": x["family"], "contig": x["contig"], "start": x["start"], "end": x["end"],
                   "confidence": x["confidence"], "locus_class": x.get("locus_class"),
                   "detection_pass": x.get("detection_pass"),
                   "verification": (x.get("verification") or {}).get("status", ""),
                   "genes_found": "|".join(x.get("genes_found") or [])}
            row.update(loci_columns(x))
            rows.append(row)
            if "array_id" in x:
                sup[(row["verification"] or "verified_other", row["array_support"] or "no_array")] += 1
        for a in c.get("receptor_arrays") or []:
            arr_rows.append({"genome": g, **{k: (",".join(v) if isinstance(v, list) else v) for k, v in a.items()}})
    cols = ["genome", "family", "contig", "start", "end", "confidence", "locus_class", "detection_pass",
            "verification", "genes_found", *LOCI_ARRAY_COLUMNS]
    with open(os.path.join(outdir, "loci.tsv"), "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=cols, delimiter="\t"); w.writeheader(); w.writerows(rows)
    if arr_rows:
        with open(os.path.join(outdir, "arrays.tsv"), "w", newline="") as fo:
            w = csv.DictWriter(fo, fieldnames=list(arr_rows[0]), delimiter="\t"); w.writeheader(); w.writerows(arr_rows)
    with open(os.path.join(outdir, "summary.txt"), "w") as fo:
        fo.write(f"genomes compared: {len(genomes)} (base {len(base)}, cand {len(cand)})\n")
        fo.write(f"detected entries: base {nb}, cand {nc}\n")
        fo.write(f"differences outside the array fields: {len(diffs)}\n")
        for g, d in diffs:
            fo.write(f"  {g}\t{d}\n")
        fo.write(f"PR calls with array fields: {sum(sup.values())}\n")
        for (v, s), n in sorted(sup.items()):
            fo.write(f"  verification={v}\tarray_support={s}\t{n}\n")
        fo.write(f"arrays: {len(arr_rows)}, by size: "
                 f"{dict(sorted(collections.Counter(int(a['array_size']) for a in arr_rows).items()))}\n")
    print(open(os.path.join(outdir, "summary.txt")).read())
    sys.exit(1 if diffs else 0)


main()
