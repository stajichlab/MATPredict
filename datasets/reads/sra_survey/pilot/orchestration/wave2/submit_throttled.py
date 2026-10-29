#!/usr/bin/env python3
"""Submit the sbatch lines of a campaign so that at most --max-jobs of them run at once (dependency chains).

Input: the output of `scripts/submit_binned_panel.py --dry-run` (tab-separated lines: chunk name, genome count, the
sbatch command). The lines are dealt round-robin into --max-jobs chains; each job in a chain waits for the previous
one with `--dependency=afterany` (so a failed chunk does not stall its chain). The CPUs in flight are
max-jobs x cpus-per-task and must stay under --cpu-limit (the curator's ruling for the v0.7.0 campaign: under 60).

    submit_throttled.py --lines dry.txt --max-jobs 3 --partition highclock --jobs-out jobs.tsv [--dry-run]
"""
import argparse
import re
import shlex
import subprocess
import sys

ap = argparse.ArgumentParser()
ap.add_argument("--lines", required=True)
ap.add_argument("--max-jobs", type=int, default=3)
ap.add_argument("--cpu-limit", type=int, default=59, help="most CPUs in flight (default 59: under 60)")
ap.add_argument("--partition", help="replace the partition in each sbatch line")
ap.add_argument("--time", help="replace --time in each sbatch line (for example 2:00:00)")
ap.add_argument("--jobs-out", help="write job id, chunk, chain here")
ap.add_argument("--dry-run", action="store_true")
a = ap.parse_args()

cmds = []
for line in open(a.lines):
    parts = line.rstrip("\n").split("\t", 2)
    if len(parts) == 3 and "sbatch" in parts[2]:
        cmds.append(parts)
if not cmds:
    sys.exit(f"no sbatch lines in {a.lines}")
cpus = {int(m.group(1)) for _, _, c in cmds for m in [re.search(r"--cpus-per-task=(\d+)", c)] if m}
if len(cpus) != 1:
    sys.exit(f"expected one cpus-per-task value, found {sorted(cpus)}")
per = cpus.pop()
if a.max_jobs * per > a.cpu_limit:
    sys.exit(f"{a.max_jobs} jobs x {per} CPUs = {a.max_jobs * per} exceeds the limit of {a.cpu_limit}")
print(f"{len(cmds)} chunks, {a.max_jobs} chains x {per} CPUs = {a.max_jobs * per} CPUs in flight (limit {a.cpu_limit})")

chains = [cmds[i::a.max_jobs] for i in range(a.max_jobs)]
rows = []
for k, chain in enumerate(chains):
    prev = None
    for name, n, cmd in chain:
        if a.partition:
            cmd = re.sub(r"--partition=\S+", f"--partition={a.partition}", cmd)
        if a.time:
            cmd = re.sub(r"--time=\S+", f"--time={a.time}", cmd)
        if prev:
            cmd = cmd.replace("sbatch ", f"sbatch --dependency=afterany:{prev} ", 1)
        if a.dry_run:
            print(f"chain {k}: {name} ({n} genomes), waits for {prev or 'nothing'}")
            prev = name
            rows.append((prev, name, k))
            continue
        out = subprocess.run(shlex.split(cmd), capture_output=True, text=True)
        if out.returncode:
            sys.exit(f"sbatch failed for {name}: {out.stderr.strip()}")
        prev = out.stdout.strip().split(";")[0]
        rows.append((prev, name, k))
        print(f"submitted {name} as {prev} (chain {k})")
if a.jobs_out and not a.dry_run:
    with open(a.jobs_out, "w") as fh:
        for jid, name, k in rows:
            fh.write(f"{jid}\t{name}\tchain{k}\n")
