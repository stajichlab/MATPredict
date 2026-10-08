#!/usr/bin/env python3
"""Write the read pairs in which either mate has a blastx hit.

usage: recruit_pairs.py R1.fastq.gz R2.fastq.gz HIT_IDS.txt OUT_R1.fastq.gz OUT_R2.fastq.gz
HIT_IDS.txt: one read id per line (first token of the FASTQ header, trailing /1 or /2 removed).
"""
import gzip
import sys


def rid(header: str) -> str:
    r = header[1:].split()[0]
    return r[:-2] if r.endswith(("/1", "/2")) else r


def main(r1, r2, ids, o1, o2):
    keep = {l.strip() for l in open(ids) if l.strip()}
    n = 0
    with gzip.open(r1, "rt") as a, gzip.open(r2, "rt") as b, gzip.open(o1, "wt") as oa, gzip.open(o2, "wt") as ob:
        while True:
            ha = a.readline()
            if not ha:
                break
            ra = [ha, a.readline(), a.readline(), a.readline()]
            rb = [b.readline() for _ in range(4)]
            if rid(ha) in keep or rid(rb[0]) in keep:
                oa.writelines(ra); ob.writelines(rb); n += 1
    print(f"recruited {n} pairs from {len(keep)} hit ids", file=sys.stderr)


if __name__ == "__main__":
    main(*sys.argv[1:6])
