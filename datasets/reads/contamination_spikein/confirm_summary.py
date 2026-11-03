#!/usr/bin/env python3
"""Summarise confirm_flagged.slurm: one row per read set with the unmapped fraction and what GmG calls the unmapped reads.

  confirm_summary.py --out DIR --ids ids.tsv [--table results.tsv]
A set is `confirmed` when GmG names the expected second genus as the primary genus of the unmapped reads (or as a mixture member) and
the unmapped fraction is above the largest unmapped fraction of the negative controls; `not_confirmed` when it does not; `no_unmapped`
when too few reads were left to classify.
"""
import argparse, csv, json, os, sys


def gmg_summary(path):
    try:
        rec = json.loads(open(path).readline())
    except (OSError, ValueError):
        return {}
    q = (rec.get("qc") or {}); m = q.get("mixture") or {}; s = m.get("second") or {}
    calls = (rec.get("guess") or {}).get("calls") or {}
    return {"genus": (calls.get("genus") or {}).get("taxon"), "species": (calls.get("species") or {}).get("taxon") or (calls.get("species") or {}).get("vote"),
            "status": rec.get("status"), "primary": m.get("primary"), "second": s.get("taxon"), "second_hashes": s.get("hashes"),
            "n_hashes": q.get("n_hashes"), "coverage_est": q.get("coverage_est"), "mixed": bool(m.get("suspected"))}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", required=True); ap.add_argument("--ids", required=True); ap.add_argument("--table")
    a = ap.parse_args(argv)
    rows = []
    for line in open(a.ids):
        parts = line.rstrip("\n").split("\t")
        if not parts[0]:
            continue
        sid, exp = parts[0], (parts[1] if len(parts) > 1 else "")
        rp = os.path.join(a.out, sid, "row.tsv")
        if not os.path.exists(rp):
            rows.append({"id": sid, "expected": exp, "state": "not_run"}); continue
        f = open(rp).read().rstrip("\n").split("\t")
        pairs, unm = int(f[2]), int(f[3])
        r = {"id": sid, "expected": exp, "state": "done", "pairs": pairs, "unmapped_pairs": unm, "unmapped_pct": round(100 * unm / pairs, 2) if pairs else None,
             "aligned_pct": f[4], **{("gmg_" + k): v for k, v in gmg_summary(os.path.join(a.out, sid, "gmg", "results.jsonl")).items()}}
        rows.append(r)
    neg = [r["unmapped_pct"] for r in rows if r.get("state") == "done" and not r["expected"] and r.get("unmapped_pct") is not None]
    ceiling = max(neg) if neg else None
    for r in rows:
        if r.get("state") != "done" or not r["expected"]:
            r["verdict"] = "control" if r.get("state") == "done" else r.get("state"); continue
        if r.get("gmg_genus") is None:
            r["verdict"] = "no_unmapped"
        else:
            named = r["expected"] in (r.get("gmg_genus"), r.get("gmg_primary"), r.get("gmg_second"))
            r["verdict"] = "confirmed" if named and (ceiling is None or r["unmapped_pct"] > ceiling) else "not_confirmed"
    cols = ["id", "expected", "state", "pairs", "unmapped_pairs", "unmapped_pct", "aligned_pct", "gmg_genus", "gmg_species", "gmg_status", "gmg_primary", "gmg_second", "gmg_second_hashes", "gmg_n_hashes", "gmg_coverage_est", "gmg_mixed", "verdict"]
    w = csv.DictWriter(open(a.table, "w") if a.table else sys.stdout, cols, delimiter="\t", extrasaction="ignore", lineterminator="\n"); w.writeheader(); w.writerows(rows)
    print(f"negative-control unmapped ceiling: {ceiling}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
