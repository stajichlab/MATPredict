#!/usr/bin/env python3
"""Assembly-level idiomorph truth: BLASTN of the idiomorph-specific regions against one assembly.

usage: assembly_idiomorph_blast.py SPECIFIC_REGIONS.fasta GENOME.fasta WORKDIR
Prints one TSV row: genome, then for each region: covered fraction (union of alignments at >= 90% identity, e-value 1e-10)
and best identity. Call rule (explicit): present >= 0.80 covered, absent <= 0.20, otherwise partial.
"""
import subprocess
import sys
from pathlib import Path


def region_lengths(fa):
    out, name = {}, None
    for line in Path(fa).read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split("_")[0]; out[name] = 0
        else:
            out[name] += len(line.strip())
    return out


def union_len(iv):
    iv.sort(); tot, end = 0, 0
    for s, e in iv:
        s = max(s, end + 1) if end >= s else s
        if e >= s:
            tot += e - s + 1; end = max(end, e)
    return tot


def call(cov):
    return "present" if cov >= 0.80 else "absent" if cov <= 0.20 else "partial"


def main(regions, genome, work):
    work = Path(work); work.mkdir(parents=True, exist_ok=True)
    db = work / Path(genome).stem
    subprocess.run(["makeblastdb", "-in", genome, "-dbtype", "nucl", "-out", str(db)], check=True, capture_output=True)
    res = subprocess.run(["blastn", "-query", regions, "-db", str(db), "-evalue", "1e-10", "-perc_identity", "90",
                          "-outfmt", "6 qseqid pident qstart qend"], check=True, capture_output=True, text=True).stdout
    lens = region_lengths(regions)
    iv = {k: [] for k in lens}; best = {k: 0.0 for k in lens}
    for line in res.splitlines():
        q, pid, qs, qe = line.split("\t"); k = q.split("_")[0]
        iv[k].append((int(qs), int(qe))); best[k] = max(best[k], float(pid))
    row = [Path(genome).name.replace(".sorted.fasta", "")]
    for k in sorted(lens):
        cov = union_len(iv[k]) / lens[k]
        row += [f"{cov:.3f}", f"{best[k]:.1f}", call(cov)]
    print("\t".join(row))
    for f in work.glob(db.name + ".n*"):
        f.unlink()


if __name__ == "__main__":
    main(*sys.argv[1:4])
