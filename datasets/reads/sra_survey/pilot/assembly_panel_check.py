#!/usr/bin/env python3
"""Type an ASSEMBLY with a species panel, as if its sequence were reads (two copies, so the flank depth is 2): which idiomorph does the assembly itself carry, by the same k-mer rule the reads were typed with?
Used to look at disagreements between a reads call and the `detect` call on the assembly of the same BioSample: if the panel agrees with `detect` the assembly really carries that idiomorph and the reads disagree
(wrong run, mixed sample, divergent allele); if the panel agrees with the reads, the `detect` call is the odd one.

    assembly_panel_check.py --panel PANEL_DIR --prefix PFX --assembly GENOME.fa.gz [--k 31] [--keep]
Prints one row per sub-panel (or the panel itself). Needs the MATPredict checkout on PYTHONPATH (as the worker does).
"""
import argparse, glob, gzip, os, subprocess, sys, tempfile


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0]); ap.add_argument("--panel", required=True); ap.add_argument("--prefix", required=True); ap.add_argument("--assembly", required=True)
    ap.add_argument("--k", default="31"); ap.add_argument("--chunk", type=int, default=1000)
    a = ap.parse_args()
    op = gzip.open if a.assembly.endswith(".gz") else open
    seqs, cur = [], []
    with op(a.assembly, "rt") as fh:
        for ln in fh:
            if ln.startswith(">"):
                if cur: seqs.append("".join(cur)); cur = []
            else: cur.append(ln.strip().upper())
    if cur: seqs.append("".join(cur))
    tmp = tempfile.mkdtemp(dir=os.environ.get("SCRATCH", "/tmp")); fq = os.path.join(tmp, "asm.fastq")
    with open(fq, "w") as o:
        for copy in (1, 2):
            n = 0
            for s in seqs:
                for i in range(0, len(s), a.chunk):
                    c = s[i:i + a.chunk]
                    if len(c) >= 40: n += 1; o.write(f"@c{copy}_{n}\n{c}\n+\n{'I' * len(c)}\n")
    open(os.path.join(tmp, "samples.tsv"), "w").write(f"sample\treads\nasm\t{fq}\n")
    subs = sorted(glob.glob(os.path.join(a.panel, "L*/"))) or [a.panel + "/"]
    first = True
    for sub in subs:
        out = os.path.join(tmp, "r.tsv")
        cmd = [sys.executable, "-m", "MATPredict", "reads-type", "--k", a.k] + [x for n in ("MAT1-1", "MAT1-2") for x in ("--idiomorph", f"{n}={sub}{a.prefix}_{n}.fasta")] + ["--samples", os.path.join(tmp, "samples.tsv"), "--out", out]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode: print(os.path.basename(sub.rstrip("/")), "failed:", r.stderr[-200:]); continue
        rows = open(out).read().replace("\r", "").splitlines()
        if first: print("panel\t" + rows[0]); first = False
        print(os.path.basename(sub.rstrip("/")) + "\t" + rows[1])
    subprocess.run(["rm", "-rf", tmp])


if __name__ == "__main__":
    main()
