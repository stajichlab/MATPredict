#!/usr/bin/env python3
"""Annotation-free search for HMG-box and homeodomain regions in assemblies without gene models.

Read-only discovery tool. For genomes that BFD did not annotate (fragmentary assemblies), a gene-order scan is not possible. This script
searches the assembly directly with tblastn and reports, per genome, the regions that match HMG-class and homeodomain-class controls and
the contigs that carry both within --pair-dist bp.

Queries (FASTA headers): the class comes from the header. Two header styles are understood:
  curated '>id|gene|class' (class sexM, sexP -> HMG; matMc, mat1-Mc, matMi... -> HMG; matPi, mat1-Pc, mat1-Pi, matPc -> HD)
  free text containing 'hd1' or 'hd2' -> HD
Other classes (flank genes) are ignored.
Per genome (TSV: name, fasta): tblastn (E <= --evalue, SEG off); hits of one class on one contig closer than --merge bp are merged into
a region; a pair is an HMG region and an HD region on one contig with a gap <= --pair-dist bp.
Output (--out): regions.tsv (genome, class, contig, start, end, identity, evalue, bitscore) and summary.tsv (genome, contigs, bp,
hmg_regions, hd_regions, pairs, contigs_with_pair, longest_pair_contig_bp).
Limit: tblastn with controls from Mucoromycota and Taphrinomycotina finds HMG and homeodomain genes of any family; a pair in a short contig
is weak evidence (contigs are shorter than the pair distance), and no pair is not evidence of absence.
"""
import argparse
import collections
import gzip
import os
import re
import subprocess
import sys

HMG = {"sexM", "sexP", "matMc", "mat1-Mc", "matMi", "mat1-Mi", "matMc|"}
HD = {"matPi", "mat1-Pc", "mat1-Pi", "matPc", "HD1", "HD2"}


def query_class(header):
    """'HMG', 'HD' or None from a FASTA header (without '>')."""
    h = header.strip()
    low = h.lower()
    if "|" in h:
        tail = h.split("|")[-1].split()[0]
        if "name=" in h:
            for tok in h.split("|"):
                if tok.startswith("name="):
                    tail = tok[5:]
        if tail in HMG:
            return "HMG"
        if tail in HD:
            return "HD"
        return None
    if "hd1" in low or "hd2" in low:
        return "HD"
    return None


def merge(hits, gap):
    """hits: (contig, start, end, ident, evalue, bits) -> merged regions sorted by contig and start."""
    by = collections.defaultdict(list)
    for c, s, e, p, ev, b in hits:
        by[c].append((min(s, e), max(s, e), p, ev, b))
    out = []
    for c, v in by.items():
        v.sort()
        cur = None
        for s, e, p, ev, b in v:
            if cur and s - cur[2] <= gap:
                cur[2] = max(cur[2], e)
                cur[3] = max(cur[3], p)
                cur[4] = min(cur[4], ev)
                cur[5] += b
            else:
                if cur:
                    out.append(tuple(cur))
                cur = [c, s, e, p, ev, b]
        if cur:
            out.append(tuple(cur))
    return sorted(out)


def pairs(hmg, hd, dist):
    """(hmg_region, hd_region) with both on one contig and a gap <= dist."""
    out = []
    for a in hmg:
        for b in hd:
            if a[0] == b[0] and max(0, max(a[1], b[1]) - min(a[2], b[2])) <= dist:
                out.append((a, b))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--queries", nargs="+", required=True, help="FASTA files; classes from headers")
    ap.add_argument("--out", required=True)
    ap.add_argument("--evalue", type=float, default=1e-5)
    ap.add_argument("--merge", type=int, default=2000)
    ap.add_argument("--pair-dist", type=int, default=20000)
    ap.add_argument("--tmp", default=os.environ.get("SCRATCH", "/tmp"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    qfa = {"HMG": os.path.join(a.tmp, "q_hmg.faa"), "HD": os.path.join(a.tmp, "q_hd.faa")}
    handles = {k: open(v, "w") for k, v in qfa.items()}
    counts = collections.Counter()
    for path in a.queries:
        cur, buf = None, []
        for line in open(path):
            if line[0] == ">":
                if cur:
                    handles[cur].write("".join(buf))
                    counts[cur] += 1
                cls = query_class(line[1:])
                cur, buf = cls, ([line] if cls else [])
            elif cur:
                buf.append(line)
        if cur:
            handles[cur].write("".join(buf))
            counts[cur] += 1
    for h in handles.values():
        h.close()
    print("queries:", dict(counts))
    fr = open(os.path.join(a.out, "regions.tsv"), "w")
    fs = open(os.path.join(a.out, "summary.tsv"), "w")
    fr.write("genome\tclass\tcontig\tstart\tend\tidentity\tevalue\tbitscore\n")
    fs.write("genome\tcontigs\tbp\thmg_regions\thd_regions\tpairs\tcontigs_with_pair\tlongest_pair_contig_bp\n")
    for line in open(a.genomes):
        name, fasta = line.rstrip("\n").split("\t")[:2]
        safe = re.sub(r"[^A-Za-z0-9_.-]", "_", name)  # the label may hold '|' or spaces; file names must be shell-safe
        genome = os.path.join(a.tmp, safe + ".fa")
        subprocess.run(["bash", "-c", f"zcat -f \"$0\" > \"$1\"", fasta, genome], check=True)
        db = os.path.join(a.tmp, "db_" + safe)
        subprocess.run(["makeblastdb", "-in", genome, "-dbtype", "nucl", "-out", db, "-logfile", "/dev/null"], check=True)
        lens, cur = {}, None
        for l in open(genome):
            if l[0] == ">":
                cur = l[1:].split()[0]
                lens[cur] = 0
            else:
                lens[cur] += len(l.strip())
        regs = {}
        for cls in ("HMG", "HD"):
            out = os.path.join(a.tmp, safe + "." + cls + ".tb")
            subprocess.run(["tblastn", "-query", qfa[cls], "-db", db, "-seg", "no", "-evalue", str(a.evalue), "-out", out,
                            "-outfmt", "6 sseqid sstart send pident evalue bitscore"], check=True)
            hits = []
            for l in open(out):
                c, s, e, p, ev, b = l.split("\t")
                hits.append((c, int(s), int(e), float(p), float(ev), float(b)))
            regs[cls] = merge(hits, a.merge)
            for r in regs[cls]:
                fr.write(f"{name}\t{cls}\t{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\t{r[5]}\n")
            os.remove(out)
        pr = pairs(regs["HMG"], regs["HD"], a.pair_dist)
        pc = {p[0][0] for p in pr}
        longest = max((lens[c] for c in pc), default=0)
        fs.write(f"{name}\t{len(lens)}\t{sum(lens.values())}\t{len(regs['HMG'])}\t{len(regs['HD'])}\t{len(pr)}\t{len(pc)}\t{longest}\n")
        os.remove(genome)
    fr.close()
    fs.close()
    print(open(os.path.join(a.out, "summary.tsv")).read())


if __name__ == "__main__":
    sys.exit(main())
