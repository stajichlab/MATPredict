#!/usr/bin/env python3
"""Leave-one-out recall against assembly-mapped coordinates, at several radii.

For each record in the coordinate-truth table: withhold it (and, by radius, its
species / genus / family / order) from the reference set, run detection on the
assembly its locus really sits on, and ask whether any reported locus overlaps
the true coordinates.

WHAT A NUMBER HERE MEANS DEPENDS ENTIRELY ON THE RADIUS. At RECORD radius the
remaining references usually include the same species, so a hit shows
robustness to strain and assembly and little else. At ORDER radius nothing from
the clade remains, which is the real question -- can a genome from an uncurated
group be called at all. Recall is therefore reported per radius, never pooled.

A radius that leaves NO reference for the family is reported as `no_reference`,
not as a miss. Scoring it as failure would blame the pipeline for a gap in the
curated database, which is a different problem with a different fix.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from collections import Counter
from pathlib import Path

import yaml

RADII = ["record", "species", "genus", "family", "order"]


def truth_rows(wt: Path) -> dict[str, dict]:
    """record_id -> {accession: (start, end)} spanning that record's genes."""
    out: dict[str, dict] = {}
    doc = yaml.safe_load((wt / "testset/coordinate_truth/ascomycota.yaml").read_text())
    for r in doc["records"]:
        spans: dict[str, list[int]] = {}
        for g in r["genes"]:
            for p in g["placements"]:
                a = p["sequence_accession"]
                lo, hi = spans.get(a, [p["start"], p["end"]])
                spans[a] = [min(lo, p["start"]), max(hi, p["end"])]
        if spans:
            out[r["record_id"]] = {"species": r["species"], "spans": spans}
    sup = yaml.safe_load((wt / "testset/coordinate_truth/supplemental_locus_tags.yaml").read_text())
    for rid, r in sup["records"].items():
        lo = min(g["start"] for g in r["genes"].values())
        hi = max(g["end"] for g in r["genes"].values())
        out[rid] = {"species": r["species"], "spans": {r["sequence_accession"]: [lo, hi]}}
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--worktree", type=Path, required=True)
    ap.add_argument("--assemblies", type=Path, required=True)
    ap.add_argument("--library", type=Path,
                    default=Path("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes"))
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--radii", nargs="*", default=RADII)
    args = ap.parse_args()

    sys.path.insert(0, str(args.worktree / "src"))
    from MATPredict.detect.holdout import Radius, records_to_withhold   # noqa: E402

    wt, db = args.worktree, args.worktree / "db"
    truth = truth_rows(wt)
    asm = json.loads(args.assemblies.read_text())
    args.work.mkdir(parents=True, exist_ok=True)
    py = "/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin/python"

    results = []
    for rid, t in sorted(truth.items()):
        a = asm.get(rid) or {}
        if not a.get("bfd_asmid"):
            results.append({"record_id": rid, "species": t["species"], "status": "no_assembly"})
            print(f"{rid[:33]:35}SKIP no assembly", file=sys.stderr)
            continue
        # Decompressed ONCE and shared between the per-radius jobs, which all
        # walk this record list in lockstep and so reach each new genome at the
        # same moment. `zcat > path` truncates before it writes, so a naive
        # shared write lets one job blank the file while another's detect is
        # mid-read -- silently scoring against a partial genome. Write to a
        # per-process temp and os.replace, which is atomic within a filesystem:
        # a reader sees either the old complete file or the new one, never a
        # half-written one.
        fna = args.work / f"{a['bfd_asmid']}.fna"
        if not fna.exists():
            tmp = args.work / f".{a['bfd_asmid']}.{os.getpid()}.tmp"
            subprocess.run(f"zcat {args.library}/{a['bfd_asmid']}.fa.gz > {tmp}",
                           shell=True, check=True)
            os.replace(tmp, fna)
        for radius in args.radii:
            drop = records_to_withhold(db, rid, Radius(radius))
            out = args.work / "runs" / rid / radius
            rep = out / "detection_report.yaml"
            if not rep.exists():
                out.mkdir(parents=True, exist_ok=True)
                cmd = [py, "-m", "MATPredict", "detect", "--genome", str(fna),
                       "--out-dir", str(out), "--exclude-records", ",".join(sorted(drop))]
                env = {"PYTHONPATH": str(wt / "src"),
                       "PATH": "/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin:/usr/bin:/bin"}
                r = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=7200)
                if r.returncode != 0:
                    results.append({"record_id": rid, "radius": radius, "status": "error",
                                    "withheld": len(drop), "stderr": r.stderr.strip()[-300:]})
                    print(f"{rid[:28]:30}{radius:8}ERROR {r.stderr.strip()[-90:]}", file=sys.stderr)
                    continue
            doc = yaml.safe_load(rep.read_text()) or {}
            det = doc.get("detected") or []
            hit = None
            for L in det:
                span = t["spans"].get(L.get("contig"))
                if span and not (L["end"] < span[0] or L["start"] > span[1]):
                    hit = L
                    break
            results.append({
                "record_id": rid, "species": t["species"], "radius": radius,
                "withheld": len(drop), "loci": len(det),
                "status": "hit" if hit else "miss",
                "confidence": (hit or {}).get("confidence"),
                "idiomorph": (hit or {}).get("idiomorph"),
                "families_attempted": doc.get("families_attempted") or [],
            })
            print(f"{rid[:28]:30}{radius:8}withheld={len(drop):3} loci={len(det):3} "
                  f"{'HIT ' + str((hit or {}).get('confidence')) if hit else 'MISS'}", file=sys.stderr)
    args.out.write_text(yaml.safe_dump({"results": results}, sort_keys=False))
    print(f"\nwrote {args.out}", file=sys.stderr)
    for radius in args.radii:
        rows = [r for r in results if r.get("radius") == radius and r["status"] in ("hit", "miss")]
        if rows:
            h = sum(1 for r in rows if r["status"] == "hit")
            print(f"  {radius:8} {h}/{len(rows)} recovered", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
