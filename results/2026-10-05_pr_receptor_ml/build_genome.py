#!/usr/bin/env python3
"""Per genome: STE3-like loci (miniprot of ste3_all.faa), the receptor protein and
intron count of each locus, strict-CAAX / precursor-homology counts and nearest
distances, nearest HD-class and flank MAT gene hit. No annotation used.

Usage: build_genome.py ASMID OUTDIR
Writes OUTDIR/<ASMID>.loci.tsv and OUTDIR/<ASMID>.prot.faa
Reuses the functions of scan_genome.py (PR #33).
"""
import bisect
import importlib.util
import os
import subprocess
import sys
import tempfile
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("scan_genome", os.path.join(HERE, "scan_genome.py"))
sg = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sg)
sg.STE3 = os.path.join(HERE, "ste3_all.faa")
sg.PHERO = os.path.join(HERE, "pheromones_curated.faa")
BIN, LIB = sg.BIN, sg.LIB
MAT = os.path.join(HERE, "matgenes.faa")
THREADS = os.environ.get("SLURM_CPUS_PER_TASK", "8")
BIG = 10_000_000
# Leave-own-species-out: curated precursor and MAT-gene queries from the genome's own species
# (record-id taxid prefix) are not used, so a locus does not "find" the gene its own record was
# built from. Record ids start with the species taxid.
OWN_TAX = {
    "GCA_016772295.1": "5346", "GCF_000143185.2": "5334", "GCF_000271585.1": "5325", "GCA_001683735.1": "5627",
    "GCA_984573805.1": "2830151", "GCF_000320585.1": "984962", "GCF_000328475.2": "5270", "GCF_000091045.1": "40410",
    "GCA_056621545.1": "40410", "GCA_000988875.2": "5286", "GCA_921037615.3": "5286", "GCA_026119225.1": "29898",
    "GCA_920103745.3": "5535", "GCA_024748845.1": "5537", "GCA_023212685.1": "1652704", "GCA_023212835.1": "203535",
    "GCA_023212605.1": "203536", "GCA_023212725.1": "349360", "GCA_023212615.1": "49012", "GCA_023212695.1": "63387",
    "GCA_023212635.2": "84751", "GCA_056320075.1": "1708542", "GCF_000263375.1": "671144",
}


def parse_gff(path):
    out, cur = [], None
    for line in open(path):
        if line.startswith("##STA"):
            if cur is not None:
                cur["prot"] = line.split("\t", 1)[1].strip().replace("*", "").replace("-", "")
            continue
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        if f[2] == "mRNA":
            cur = {"contig": f[0], "strand": f[6], "start": int(f[3]), "end": int(f[4]),
                   "score": int(f[5]) if f[5].isdigit() else 0, "prot": "", "ncds": 0}
            out.append(cur)
        elif f[2] == "CDS" and cur is not None:
            cur["ncds"] += 1
    return out


def paf_hits(genome_gz, query, threads, outn=100):
    p = subprocess.run([f"{BIN}/miniprot", "-t", str(threads), "-I", f"--outn={outn}", genome_gz, query],
                       stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, check=True, text=True).stdout
    hits = []
    for line in p.splitlines():
        f = line.split("\t")
        if len(f) < 12:
            continue
        qlen, qs, qe = int(f[1]), int(f[2]), int(f[3])
        tags = dict(x.split(":", 2)[0::2] for x in f[12:] if x.count(":") >= 2)
        hits.append((f[0], (qe - qs) / qlen, f[5], int(f[7]), int(f[8]), int(tags.get("AS", 0))))
    return hits


def nearest(arr_by_ctg, ctg, a, b):
    arr = arr_by_ctg.get(ctg)
    if not arr:
        return BIG
    i = bisect.bisect_left(arr, a)
    best = BIG
    for j in (i - 1, i):
        if 0 <= j < len(arr):
            p = arr[j]
            d = 0 if a <= p <= b else min(abs(p - a), abs(p - b))
            best = min(best, d)
    return best


def main():
    asm, outdir = sys.argv[1], sys.argv[2]
    os.makedirs(outdir, exist_ok=True)
    scratch = os.environ.get("SCRATCH", "/bigdata/stajichlab/jstajich/prml_work/tmp")
    os.makedirs(scratch, exist_ok=True)
    tmp = tempfile.mkdtemp(prefix="prml_" + asm + "_", dir=scratch)
    gz = os.path.join(LIB, asm + ".fa.gz")
    fa = os.path.join(tmp, asm + ".fa")
    seqs = sg.read_fasta(gz)
    with open(fa, "w") as fo:
        for k, v in seqs.items():
            fo.write(f">{k}\n{v}\n")
    loci = sg.miniprot_loci(gz, tmp, threads=int(THREADS))
    cands = sg.orf_candidates(seqs)
    own = OWN_TAX.get(asm[:15], "NONE") + "_"
    homs = [h for h in sg.homology_hits(fa, tmp) if not h[2].startswith(own)]
    gff = os.path.join(tmp, "ste3.gff")
    with open(gff, "w") as fo:
        subprocess.run([f"{BIN}/miniprot", "-t", THREADS, "-I", "--gff", "--trans", "--outn=100", gz, os.path.join(HERE, "ste3_nonref.faa")],
                       stdout=fo, stderr=subprocess.DEVNULL, check=True)
    mr = parse_gff(gff)
    mh = defaultdict(lambda: defaultdict(list))
    for q, qcov, t, ts, te, sc in paf_hits(gz, MAT, THREADS):
        if qcov >= 0.5 and sc >= 150 and te - ts < 60000 and not q.startswith(own):
            mh[q.split("|")[1]][t].append((ts + te) // 2)
    for c in mh:
        for t in mh[c]:
            mh[c][t].sort()
    T, T2, Lx, H = defaultdict(list), defaultdict(list), defaultdict(list), defaultdict(list)
    for ctg, pos, strand, cls, tail in cands:
        base = cls.split("+")[0]
        if base == "T":
            T[ctg].append(pos); T2[ctg].append(pos); Lx[ctg].append(pos)
        elif "T2" in cls:
            T2[ctg].append(pos)
        if base == "L":
            Lx[ctg].append(pos)
    for ctg, pos, q, ev in homs:
        if ev <= 1:
            H[ctg].append(pos)
    for d in (T, T2, Lx, H):
        for k in d:
            d[k].sort()
    fo = open(os.path.join(outdir, asm + ".loci.tsv"), "w")
    fp = open(os.path.join(outdir, asm + ".prot.faa"), "w")
    cols = ["asm", "contig", "strand", "start", "end", "best_query", "best_qcov", "best_ident", "best_score",
            "best_ref", "best_ref_ident", "n_queries", "prot_len", "n_cds", "aln_score",
            "d_T", "d_T2", "d_L", "d_Hx", "d_HD", "d_FLANK", "nT_10kb", "nT_20kb", "nT_50kb",
            "nHx_10kb", "nT2_10kb", "nL_10kb", "n_HD_20kb", "n_FLANK_50kb"]
    fo.write("\t".join(cols) + "\n")
    for Lc in sorted(loci, key=lambda x: (x["contig"], x["start"])):
        c, a, b = Lc["contig"], Lc["start"], Lc["end"]
        ov = [m for m in mr if m["contig"] == c and m["start"] <= b and m["end"] >= a and m["prot"]]
        m = max(ov, key=lambda m: m["score"]) if ov else None

        def cnt(d, W):
            arr = d.get(c, [])
            return bisect.bisect_right(arr, b + W) - bisect.bisect_left(arr, a - W)

        def cnt_m(cls, W):
            arr = mh[cls].get(c, [])
            return bisect.bisect_right(arr, b + W) - bisect.bisect_left(arr, a - W)

        row = [asm, c, Lc["strand"], a, b, Lc["best_query"], f"{Lc['best_qcov']:.2f}", f"{Lc['best_ident']:.3f}",
               Lc["best_score"], Lc["best_ref"], f"{Lc['best_ref_ident']:.3f}", Lc["n_queries"],
               len(m["prot"]) if m else 0, m["ncds"] if m else 0, m["score"] if m else 0,
               nearest(T, c, a, b), nearest(T2, c, a, b), nearest(Lx, c, a, b), nearest(H, c, a, b),
               nearest(mh["HD"], c, a, b), nearest(mh["FLANK"], c, a, b),
               cnt(T, 10000), cnt(T, 20000), cnt(T, 50000), cnt(H, 10000), cnt(T2, 10000), cnt(Lx, 10000),
               cnt_m("HD", 20000), cnt_m("FLANK", 50000)]
        fo.write("\t".join(map(str, row)) + "\n")
        if m:
            fp.write(f">{asm}|{c}:{a}-{b}\n{m['prot']}\n")
    fo.close()
    fp.close()
    print(asm, len(loci), "loci", len(cands), "cands")
    os.remove(fa)


main()
