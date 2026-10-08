"""The Zygo 23 regression, on BOTH inputs, scored as two separate numbers.

Curator's ruling, 2026-09-26: run both and report two scores.

* `scaffold` -- the funannotate `.gbk` converted to scaffolds (`.fna`) plus
  the annotated proteome (`.faa`), as the 2026-09-21 baseline ran. The truth
  table is in these scaffold coordinates, so calls are scored directly.
* `contig` -- funannotate's `*.contigs.fsa` (the scaffolds split at every N
  gap, written for NCBI submission) with NO proteome, as every run since
  2026-09-23 searched. Harder: assembly gaps split loci, as in many BFD
  genomes. Each call is mapped back to scaffold coordinates through the
  `*.agp` beside the contigs file before scoring.

A genome scores "locus" when a reported call, in scaffold coordinates, lies on
the truth scaffold and overlaps the truth span; "idiomorph" when such a call
also carries the truth idiomorph.

    python scripts/zygo_regression.py run --src <worktree>/src --work $SCRATCH/zygo \\
        --out results/<date>_zygo23 [--inputs scaffold contig] [--jobs 8]
    python scripts/zygo_regression.py score --out results/<date>_zygo23

`run` needs `--src` and `--db` from ONE frozen worktree (a run from a tree
being edited has cost genomes before). Converted FASTA goes to `--work`
(node-local scratch); only reports go to `--out`.
"""
from __future__ import annotations

import argparse
import glob
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import yaml

TRUTH = Path("/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-23_zygo_bar/zygo_truth.tsv")
INPUTS = ("scaffold", "contig")


def read_truth(path: Path) -> list[dict]:
    rows = []
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 6:
            continue
        rows.append(dict(org=f[0], scaffold=f[1], start=int(f[2]), end=int(f[3]),
                         idiomorph=f[4], contigs_fsa=Path(f[5])))
    return rows


def agp_map(contigs_fsa: Path) -> dict[str, tuple[str, int, int, int, str]]:
    """component -> (scaffold, obj_start, obj_end, comp_start, orientation)."""
    m = {}
    for agp in glob.glob(str(contigs_fsa.parent / "*.agp")):
        for line in open(agp):
            f = line.rstrip("\n").split("\t")
            if len(f) >= 9 and f[4] == "W":
                m[f[5]] = (f[0], int(f[1]), int(f[2]), int(f[6]), f[8])
    return m


def to_scaffold(contig: str, start: int, end: int, amap: dict) -> tuple[str, int, int]:
    if contig not in amap:
        return contig, start, end          # already a scaffold name
    scaf, o_start, o_end, c_start, strand = amap[contig]
    if strand == "-":
        a, b = o_end - (end - c_start), o_end - (start - c_start)
    else:
        a, b = o_start + (start - c_start), o_start + (end - c_start)
    return scaf, a, b


def _detect(src: Path, db: Path, genome: Path, proteins: Path | None, out: Path) -> str:
    if (out / "detection_report.yaml").exists() and (out / "detection_report.yaml").stat().st_size:
        return "skip"
    out.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable, "-m", "MATPredict", "detect", "--genome", str(genome),
           "--phylum", "Mucoromycota", "--out-dir", str(out),
           "--evidence-diagnostics", str(out / "evidence_diagnostics.jsonl")]
    if proteins is not None:
        cmd += ["--proteins", str(proteins)]
    # MATPREDICT_HTML=0: no report.html per genome; an env var because `src` may be an older worktree.
    env = dict(os.environ, PYTHONPATH=str(src), MATPREDICT_DB_ROOT=str(db),
               MATPREDICT_HTML=os.environ.get("MATPREDICT_HTML", "0"))
    with open(out / "stdout.log", "w") as so, open(out / "stderr.log", "w") as se:
        rc = subprocess.run(cmd, env=env, stdout=so, stderr=se, timeout=3600).returncode
    return "ok" if rc == 0 else f"rc={rc}"


def cmd_run(args) -> int:
    sys.path.insert(0, str(args.src))
    from MATPredict.detect.annotation_export import convert_genbank

    truth = read_truth(args.truth)
    work = Path(args.work)
    jobs = []
    for row in truth:
        if "scaffold" in args.inputs:
            gbk = next(iter(sorted(row["contigs_fsa"].parent.glob("*.gbk"))), None)
            if gbk is None:
                print(f"FAIL {row['org']}: no .gbk beside {row['contigs_fsa']}", file=sys.stderr)
            else:
                d = work / "scaffold"
                d.mkdir(parents=True, exist_ok=True)
                fna, faa = d / f"{gbk.stem}.fna", d / f"{gbk.stem}.faa"
                if not (fna.exists() and faa.exists()):
                    convert_genbank(gbk, d)
                jobs.append((fna, faa, Path(args.out) / "scaffold" / "runs" / row["org"]))
        if "contig" in args.inputs:
            jobs.append((row["contigs_fsa"], None, Path(args.out) / "contig" / "runs" / row["org"]))
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        for (genome, prot, out), status in zip(jobs, ex.map(
                lambda j: _detect(args.src, args.db, *j), jobs)):
            if status not in ("ok", "skip"):
                print(f"FAIL {out.parent.parent.name}/{out.name}: {status}", file=sys.stderr)
    return cmd_score(args)


def score(runs: Path, truth: list[dict], mapped: bool) -> tuple[int, int, int, list[str]]:
    loc = idio = 0
    misses = []
    for row in truth:
        amap = agp_map(row["contigs_fsa"]) if mapped else {}
        p = runs / row["org"] / "detection_report.yaml"
        det = ((yaml.safe_load(open(p)) or {}).get("detected") or []) if p.exists() else []
        hits = []
        for d in det:
            sc, a, b = to_scaffold(d["contig"], d["start"], d["end"], amap)
            if sc == row["scaffold"] and a <= row["end"] and b >= row["start"]:
                hits.append(d)
        ok_loc = bool(hits)
        ok_id = any(d["idiomorph"] == row["idiomorph"] for d in hits)
        loc += ok_loc
        idio += ok_id
        if not (ok_loc and ok_id):
            misses.append(f"MISS {row['org']} truth={row['idiomorph']} "
                          f"{'no report' if not p.exists() else [(d['contig'], d['start'], d['idiomorph'], d.get('detection_pass')) for d in det]}")
    return loc, idio, len(truth), misses


def cmd_score(args) -> int:
    truth = read_truth(args.truth)
    lines = []
    for name in args.inputs:
        runs = Path(args.out) / name / "runs"
        if not runs.exists():
            continue
        loc, idio, n, misses = score(runs, truth, mapped=(name == "contig"))
        lines.append(f"{name:8s} locus on truth scaffold {loc}/{n}; idiomorph correct {idio}/{n}")
        lines += [f"  {m}" for m in misses]
    text = "\n".join(lines)
    print(text)
    Path(args.out, "score.txt").write_text(text + "\n")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("run", "score"):
        p = sub.add_parser(name)
        p.add_argument("--out", required=True)
        p.add_argument("--truth", type=Path, default=TRUTH)
        p.add_argument("--inputs", nargs="+", choices=INPUTS, default=list(INPUTS))
        if name == "run":
            p.add_argument("--src", type=Path, required=True)
            p.add_argument("--db", type=Path, help="default: <src>/../db")
            p.add_argument("--work", required=True, help="scratch dir for converted FASTA")
            p.add_argument("--jobs", type=int, default=6)
    args = ap.parse_args(argv)
    if args.cmd == "run":
        args.db = args.db or args.src.parent / "db"
        return cmd_run(args)
    return cmd_score(args)


if __name__ == "__main__":
    raise SystemExit(main())
