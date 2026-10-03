"""Submit a detect panel split into genome-size bins (curator ruling 2026-10-01).

Small genomes (< --large-bp, default 500 Mb) go to `short` in chunks sized to
~1-1.5 h of real runtime each. Large genomes go to one or more `epyc` jobs with
GENOME_TIMEOUT=14400 (4 h per genome). Each chunk is one run of
scripts/run_clade_panel.slurm with its own SAMPLE_LIST and OUT_DIR
(OUT/small_000, OUT/large_000, ...), because each job writes its own
rollout_summary.yaml and reports.tar.zst at the end.

Why: in the 2026-09-26 Basidiomycota run (results/2026-09-26_basidiomycota_full/
ANALYSIS.md) genomes < 500 Mb had order medians of 91-330 s and a maximum of
1,737 s; the 58 genomes > 500 Mb (57 Pucciniales) had a median of 2,500 s and
one timed out at the old fixed 3600 s limit.

Genome size comes from BFD tables/asm_stats.parquet (total_length_bp) or from
--sizes, a TSV of ASMID<TAB>total_length_bp. A genome with no size goes to the
small bin and is counted in the plan.

    /usr/bin/python3.12 scripts/submit_binned_panel.py \
        --sample-list LIST.tsv --out /bigdata/.../results/<date>_<name> \
        --worktree /bigdata/.../.claude/worktrees/run-<sha> [--dry-run]

Run from the repo root (logs/ must exist). Paths are absolute; the worktree is
passed explicitly, never resolved from this file's location.
"""
from __future__ import annotations

import argparse
import math
import os
import subprocess
import sys

ASM_STATS = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/tables/asm_stats.parquet"


def read_sample_list(path):
    rows = []
    for line in open(path):
        if not line.strip() or line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if f[0] == "asmid" or f[0] == "ASMID":
            continue
        rows.append((f[0], f[1] if len(f) > 1 else ""))
    return rows


def read_sizes(path):
    if path.endswith(".parquet"):
        try:
            import pyarrow.parquet as pq
        except ImportError:
            sys.exit("pyarrow is needed to read parquet; run with /usr/bin/python3.12 or pass --sizes TSV")
        t = pq.read_table(path, columns=["ASMID", "total_length_bp"]).to_pylist()
        return {r["ASMID"]: int(r["total_length_bp"]) for r in t if r["total_length_bp"] is not None}
    out = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2 and f[1].isdigit():
            out[f[0]] = int(f[1])
    return out


def plan(samples, sizes, *, large_bp, small_cpus, small_seconds, target_hours, large_per_job):
    """Split samples into chunks. Pure; returns a list of dicts.

    Small chunk size = genomes that fill `target_hours` at `small_cpus`-way
    parallelism if each takes `small_seconds` (a measured order median, not a
    guess: default 330 s, the highest Basidiomycota order median).
    """
    small = [s for s in samples if sizes.get(s[0], 0) < large_bp]
    large = [s for s in samples if sizes.get(s[0], 0) >= large_bp]
    per_chunk = max(1, math.floor(target_hours * 3600 * small_cpus / small_seconds))
    chunks = []
    for i in range(0, len(small), per_chunk):
        chunks.append({"bin": "small", "rows": small[i:i + per_chunk]})
    for i in range(0, len(large), large_per_job):
        chunks.append({"bin": "large", "rows": large[i:i + large_per_job]})
    for n, c in enumerate([c for c in chunks if c["bin"] == "small"]):
        c["name"] = f"small_{n:03d}"
    for n, c in enumerate([c for c in chunks if c["bin"] == "large"]):
        c["name"] = f"large_{n:03d}"
    return chunks, per_chunk, sum(1 for s in samples if s[0] not in sizes)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample-list", required=True, help="ASMID<TAB>TAXID per line")
    ap.add_argument("--out", required=True, help="results directory on /bigdata")
    ap.add_argument("--worktree", required=True, help="frozen MATPredict worktree on /bigdata")
    ap.add_argument("--repo", default="/bigdata/stajichlab/jstajich/projects/MATPredict")
    ap.add_argument("--sizes", default=ASM_STATS, help="asm_stats.parquet or ASMID<TAB>bp TSV")
    ap.add_argument("--name", default="binned", help="job-name prefix")
    ap.add_argument("--large-bp", type=int, default=500_000_000)
    ap.add_argument("--small-cpus", type=int, default=16)
    ap.add_argument("--small-seconds", type=float, default=330.0)
    ap.add_argument("--target-hours", type=float, default=1.25)
    ap.add_argument("--large-cpus", type=int, default=16)
    ap.add_argument("--large-per-job", type=int, default=64)
    ap.add_argument("--large-time", default="8:00:00")
    ap.add_argument("--large-timeout", type=int, default=14400)
    ap.add_argument("--mem", default="32gb")
    ap.add_argument("--large-mem", default="64gb")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)

    for p in (a.out, a.worktree):
        if not os.path.abspath(p).startswith("/bigdata/"):
            sys.exit(f"must be on /bigdata (/scratch is node-local): {p}")
    slurm = os.path.join(a.worktree, "scripts", "run_clade_panel.slurm")
    if not os.path.exists(slurm):
        sys.exit(f"no run_clade_panel.slurm in worktree: {a.worktree}")

    samples = read_sample_list(a.sample_list)
    sizes = read_sizes(a.sizes)
    chunks, per_chunk, no_size = plan(samples, sizes, large_bp=a.large_bp, small_cpus=a.small_cpus,
                                      small_seconds=a.small_seconds, target_hours=a.target_hours,
                                      large_per_job=a.large_per_job)
    n_small = sum(len(c["rows"]) for c in chunks if c["bin"] == "small")
    n_large = sum(len(c["rows"]) for c in chunks if c["bin"] == "large")
    print(f"{len(samples)} genomes: {n_small} small (< {a.large_bp / 1e6:.0f} Mb; {no_size} with no size), "
          f"{n_large} large; small chunk = {per_chunk} genomes; {len(chunks)} jobs")

    os.makedirs(os.path.join(a.out, "lists"), exist_ok=True)
    jobs = []
    for c in chunks:
        lst = os.path.join(a.out, "lists", f"{c['name']}.tsv")
        with open(lst, "w") as fh:
            for asmid, taxid in c["rows"]:
                fh.write(f"{asmid}\t{taxid}\n")
        large = c["bin"] == "large"
        cmd = ["sbatch", "--parsable",
               f"--partition={'epyc' if large else 'short'}",
               f"--time={a.large_time if large else '2:00:00'}",
               f"--cpus-per-task={a.large_cpus if large else a.small_cpus}",
               f"--mem={a.large_mem if large else a.mem}",
               f"--job-name={a.name}-{c['name']}",
               f"--output={a.repo}/logs/{a.name}_{c['name']}_%j.log",
               "--export=ALL," + ",".join([
                   f"REPO={a.repo}", f"SRC={a.worktree}/src", f"MATPREDICT_DB_ROOT={a.worktree}/db",
                   f"CLADE={a.name}_{c['name']}", f"SAMPLE_LIST={lst}",
                   f"OUT_DIR={a.out}/{c['name']}",
                   f"GENOME_TIMEOUT={a.large_timeout if large else 3600}"]),
               slurm]
        if a.dry_run:
            print(f"{c['name']}\t{len(c['rows'])}\t" + " ".join(cmd))
            continue
        jid = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout.strip()
        jobs.append(f"{jid}\t{c['name']}\t{len(c['rows'])}")
        print(jobs[-1])
    if jobs:
        with open(os.path.join(a.out, "jobs.tsv"), "a") as fh:
            fh.write("\n".join(jobs) + "\n")


if __name__ == "__main__":
    main()
