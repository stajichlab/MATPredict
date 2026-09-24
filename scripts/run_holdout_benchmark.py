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

A radius that leaves NO reference is reported as `no_reference_family` or
`no_reference_idiomorph`, not as a miss. Scoring it as failure would blame the
pipeline for a gap in the curated database, which is a different problem with
a different fix. Found loci are checked for the right idiomorph, and a locus
the modelled-gene bar withheld at the right place is `suppressed`. See
`holdout.score_holdout` for every status.

REPORTS ARE REUSED ONLY WHEN THE SAME COMMIT WROTE THEM. Each run directory
gets a `source_commit` stamp; a report from any other commit is re-run. The
first measurement scored reports written before the routing fix it reported
on, because the runner reused anything already on disk.
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
    from MATPredict.detect.family_registry import load_record_families  # noqa: E402
    from MATPredict.detect.holdout import (  # noqa: E402
        FOUND, NO_REFERENCE, Radius, load_record_idiomorphs, records_to_withhold,
        score_holdout,
    )

    wt, db = args.worktree, args.worktree / "db"
    truth = truth_rows(wt)
    asm = json.loads(args.assemblies.read_text())
    args.work.mkdir(parents=True, exist_ok=True)
    # The interpreter running this script, and its env's bin/ for tblastn,
    # exonerate and miniprot -- not a hard-coded checkout.
    py = sys.executable
    env_bin = str(Path(sys.executable).parent)
    commit = subprocess.run(["git", "-C", str(wt), "rev-parse", "HEAD"],
                            capture_output=True, text=True, check=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", str(wt), "status", "--porcelain", "--", "src", "db"],
                           capture_output=True, text=True, check=True).stdout.strip()
    if dirty:
        sys.exit(f"{wt} has uncommitted changes under src/ or db/; a stamp could not "
                 f"say which code wrote a report. Commit first.")
    record_families = load_record_families(db)
    record_idiomorphs = load_record_idiomorphs(db)

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
            stamp = out / "source_commit"
            fresh = rep.exists() and stamp.exists() and stamp.read_text().strip() == commit
            if not fresh:
                out.mkdir(parents=True, exist_ok=True)
                # --taxid is NOT optional here. Without it routing falls back
                # to `exhaustive` and every run searches all 19 families: on
                # Aspergillus nidulans that took 27 minutes and suppressed 387
                # loci, against ~2 minutes when routed. It is also the wrong
                # experiment -- exhaustive is not the configuration anyone
                # runs, so recall measured that way would not describe real use.
                cmd = [py, "-m", "MATPredict", "detect", "--genome", str(fna),
                       "--out-dir", str(out), "--exclude-records", ",".join(sorted(drop))]
                if a.get("taxid"):
                    cmd += ["--taxid", str(a["taxid"])]
                env = {"PYTHONPATH": str(wt / "src"), "PATH": f"{env_bin}:/usr/bin:/bin"}
                rep.unlink(missing_ok=True)
                r = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=7200)
                if r.returncode != 0:
                    results.append({"record_id": rid, "radius": radius, "status": "error",
                                    "withheld": len(drop), "stderr": r.stderr.strip()[-300:]})
                    print(f"{rid[:28]:30}{radius:8}ERROR {r.stderr.strip()[-90:]}", file=sys.stderr)
                    continue
                stamp.write_text(commit + "\n")
            doc = yaml.safe_load(rep.read_text()) or {}
            if "suppressed_loci" not in doc:
                sys.exit(f"{rep} lacks `suppressed_loci`: it was not written by {wt}/src")
            det = doc.get("detected") or []
            score = score_holdout(
                record_id=rid, spans=t["spans"], detected=det,
                suppressed=doc.get("suppressed_loci") or [], withheld=drop,
                record_families=record_families, record_idiomorphs=record_idiomorphs,
            )
            loc = score.locus or {}
            results.append({
                "record_id": rid, "species": t["species"], "radius": radius,
                "withheld": len(drop), "loci": len(det),
                "status": score.status,
                "expected_idiomorph": sorted(record_idiomorphs.get(rid, ())),
                "family": loc.get("family"),
                "confidence": loc.get("confidence"),
                "idiomorph": loc.get("idiomorph"),
                "families_attempted": len(doc.get("families_attempted") or []),
                "routing_mode": doc.get("routing_mode"),
                "suppressed_unpolished": doc.get("suppressed_unpolished"),
                "source_commit": commit,
            })
            print(f"{rid[:28]:30}{radius:8}withheld={len(drop):3} loci={len(det):3} "
                  f"{score.status} {loc.get('idiomorph') or ''}", file=sys.stderr)
    args.out.write_text(yaml.safe_dump({"results": results}, sort_keys=False))
    print(f"\nwrote {args.out}", file=sys.stderr)
    # Recall = found / (found + suppressed + miss), per radius, never pooled.
    # `no_reference_*` leave the denominator; `wrong_idiomorph` is found but
    # is also printed on its own, because a found locus with the wrong call is
    # not a success.
    for radius in args.radii:
        rows = [r for r in results if r.get("radius") == radius and r["status"] != "no_assembly"]
        c = Counter(r["status"] for r in rows)
        denom = sum(n for s_, n in c.items() if s_ not in NO_REFERENCE and s_ != "error")
        found = sum(c[s_] for s_ in FOUND)
        if denom:
            print(f"  {radius:8} found {found}/{denom}  right idiomorph {c['hit']}  "
                  f"wrong {c['wrong_idiomorph']}  undetermined {c['hit_undetermined']}  "
                  f"bar-withheld {c['suppressed']} (+{c['suppressed_wrong_idiomorph']} wrong)  {dict(c)}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
