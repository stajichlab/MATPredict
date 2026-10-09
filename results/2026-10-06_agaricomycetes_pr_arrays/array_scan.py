#!/usr/bin/env python3
"""Per-genome receptor-array scan for Agaricomycetes (reuses scan_genome.py of PR #33).

Outputs in OUTDIR:
  <asm>.loci.tsv    STE3-like loci (miniprot of ste3_all.faa, merged as in PR #33) with contig length,
                    strict-CAAX (T) count within 10 kb and nearest T distance, tblastn precursor-homology
                    (Hx, E<=1) count within 10 kb and nearest distance
  <asm>.region.tsv  miniprot hits of region queries (HD, STE20, MIPBF): class, contig, start, end, query, qcov, score, ident
  <asm>.cand.tsv    strict-CAAX ORF positions (T) and Hx positions (for precursor clustering)
  <asm>.chance.tsv  p_random: share of 500 random locus-sized windows with a T ORF within 10 kb
No label or detection result is used. Usage: array_scan.py ASMID OUTDIR
"""
import bisect, importlib.util, os, random, subprocess, sys, tempfile
from collections import defaultdict
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("scan_genome", os.path.join(HERE, "scan_genome.py"))
sg = importlib.util.module_from_spec(spec); spec.loader.exec_module(sg)
sg.STE3 = os.path.join(HERE, "ste3_all.faa")
sg.PHERO = os.path.join(HERE, "pheromones_curated.faa")
REGION = os.path.join(HERE, "region_queries.faa")
THREADS = os.environ.get("SLURM_CPUS_PER_TASK", "8")
BIG = 10_000_000


def nearest(arr, a, b):
    if not arr:
        return BIG
    i = bisect.bisect_left(arr, a)
    best = BIG
    for j in (i - 1, i):
        if 0 <= j < len(arr):
            p = arr[j]
            best = min(best, 0 if a <= p <= b else min(abs(p - a), abs(p - b)))
    return best


def region_hits(gz):
    p = subprocess.run([f"{sg.BIN}/miniprot", "-t", THREADS, "-I", "--outn=100", gz, REGION],
                       stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, check=True, text=True).stdout
    out = []
    for line in p.splitlines():
        f = line.split("\t")
        if len(f) < 12:
            continue
        qlen, qs, qe, mat, blen = int(f[1]), int(f[2]), int(f[3]), int(f[9]), int(f[10])
        tags = dict(x.split(":", 2)[0::2] for x in f[12:] if x.count(":") >= 2)
        qcov, sc, ts, te = (qe - qs) / qlen, int(tags.get("AS", 0)), int(f[7]), int(f[8])
        if qcov >= 0.5 and sc >= 150 and te - ts < 60000:
            out.append((f[0].split("|")[1], f[5], ts, te, f[0], qcov, sc, mat / max(1, blen)))
    return out


def main():
    asm, outdir = sys.argv[1], sys.argv[2]
    os.makedirs(outdir, exist_ok=True)
    tmp = tempfile.mkdtemp(prefix="agari_" + asm + "_", dir=os.environ.get("SCRATCH", "/tmp"))
    gz = os.path.join(sg.LIB, asm + ".fa.gz")
    fa = os.path.join(tmp, asm + ".fa")
    seqs = sg.read_fasta(gz)
    with open(fa, "w") as fo:
        for k, v in seqs.items():
            fo.write(f">{k}\n{v}\n")
    clen = {k: len(v) for k, v in seqs.items()}
    loci = sg.miniprot_loci(gz, tmp, threads=int(THREADS))
    cands = sg.orf_candidates(seqs)
    homs = sg.homology_hits(fa, tmp)
    reg = region_hits(gz)
    T, H = defaultdict(list), defaultdict(list)
    for ctg, pos, strand, cls, tail in cands:
        if cls.split("+")[0] == "T":
            T[ctg].append(pos)
    for ctg, pos, q, ev in homs:
        if ev <= 1:
            H[ctg].append(pos)
    for d in (T, H):
        for k in d:
            d[k].sort()

    def cnt(d, c, a, b, W):
        arr = d.get(c, [])
        return bisect.bisect_right(arr, b + W) - bisect.bisect_left(arr, a - W)

    with open(os.path.join(outdir, asm + ".loci.tsv"), "w") as fo:
        fo.write("\t".join(["asm", "contig", "contig_len", "strand", "start", "end", "best_ident", "best_ref_ident", "n_queries",
                            "d_T", "nT_10kb", "d_Hx", "nHx_10kb"]) + "\n")
        for L in sorted(loci, key=lambda x: (x["contig"], x["start"])):
            c, a, b = L["contig"], L["start"], L["end"]
            fo.write("\t".join(map(str, [asm, c, clen.get(c, 0), L["strand"], a, b, f"{L['best_ident']:.3f}",
                                         f"{L['best_ref_ident']:.3f}", L["n_queries"], nearest(T.get(c), a, b),
                                         cnt(T, c, a, b, 10000), nearest(H.get(c), a, b), cnt(H, c, a, b, 10000)])) + "\n")
    with open(os.path.join(outdir, asm + ".region.tsv"), "w") as fo:
        fo.write("class\tcontig\tcontig_len\tstart\tend\tquery\tqcov\tscore\tident\n")
        for cl, c, ts, te, q, qc, sc, idn in reg:
            fo.write(f"{cl}\t{c}\t{clen.get(c, 0)}\t{ts}\t{te}\t{q}\t{qc:.2f}\t{sc}\t{idn:.3f}\n")
    with open(os.path.join(outdir, asm + ".cand.tsv"), "w") as fo:
        fo.write("kind\tcontig\tpos\n")
        for c in T:
            for p in T[c]:
                fo.write(f"T\t{c}\t{p}\n")
        for c in H:
            for p in H[c]:
                fo.write(f"Hx\t{c}\t{p}\n")
    # chance: random windows of median locus length, away from STE3 loci
    rng = random.Random(1)
    ll = int(np.median([L["end"] - L["start"] for L in loci])) if loci else 2000
    ctgs = [(k, l) for k, l in clen.items() if l > ll + 21000]
    tot = sum(l for _, l in ctgs)
    avoid = defaultdict(list)
    for L in loci:
        avoid[L["contig"]].append((L["start"] - 10000, L["end"] + 10000))
    n = hit = tries = 0
    while ctgs and n < 500 and tries < 10000:
        tries += 1
        r = rng.randrange(tot)
        for k, l in ctgs:
            if r < l:
                break
            r -= l
        s = max(10000, min(r, l - ll - 10000)); e = s + ll
        if any(a < e + 10000 and b > s - 10000 for a, b in avoid.get(k, [])):
            continue
        n += 1
        hit += cnt(T, k, s, e, 10000) > 0
    with open(os.path.join(outdir, asm + ".chance.tsv"), "w") as fo:
        fo.write(f"asm\tn_windows\tn_with_T\tp_random\n{asm}\t{n}\t{hit}\t{hit / max(1, n):.4f}\n")
    print(asm, len(loci), "loci", len(reg), "region hits", len(cands), "cands", flush=True)
    os.remove(fa)
    subprocess.run(["rm", "-r", tmp])


main()
