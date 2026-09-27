#!/usr/bin/env python3
"""Receptor queue step 2: does a pheromone precursor next to an STE3-like
receptor mark the mating-type receptor?

Per genome:
 1. STE3-like loci from miniprot of the step-1 STE3 set (1,060 BFD hits +
    34 curated mating receptors) against the genome DNA, merged into loci.
 2. Pheromone-precursor candidates from a stop-anchored six-frame ORF scan
    of the genome DNA (no annotation used):
      a candidate is a stop codon whose upstream in-frame segment
      (a) ends in a CAAX box just before the stop, and
      (b) has an in-frame Met 20-130 codons upstream of the stop.
    Two CAAX alphabets:
      T  tight   C[VI][IV][AVMG]   (tightened on curated Basidiomycota
                                    precursors, docs/notes/2026-09-20)
      L  textbook C[AVLIM][AVLIM]X
    plus homology: tblastn of 31 curated precursors (-seg no, word 2), E<=1.
 3. Counts of each candidate class within +-W of each locus (W=10, 20 kb),
    and the same counts in random windows of equal size (chance).

Usage: scan_genome.py ASMID OUTDIR [--labels labels.tsv]
Writes OUTDIR/<ASMID>.loci.tsv, <ASMID>.random.tsv, <ASMID>.cand.tsv.
"""
import bisect
import gzip
import os
import random
import re
import subprocess
import sys
from collections import defaultdict

import numpy as np

LIB = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes"
HERE = os.path.dirname(os.path.abspath(sys.argv[0])) if os.path.dirname(sys.argv[0]) else os.getcwd()
BIN = "/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin"
STE3 = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-27_receptor_explore/ste3_all.faa"
PHERO = os.path.join(HERE, "pheromones_curated.faa")
WINDOWS = (10000, 20000, 50000)
N_RANDOM = 1000
MOTIF_T = re.compile(r"C[VI][IV][AVMG]$")
MOTIF_L = re.compile(r"C[AVLIM][AVLIM].$")
MOTIF_T2 = re.compile(r"C[VITE][IVT][AVMG]$")   # tight + T/E at -3, T at -2 (bbp2_a CTIA, bbp2-6 CEVM, Rhodotorula CTIA/CTVA)

CODON = {}
_b = "TCAG"
_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, a in enumerate(_b):
    for j, b in enumerate(_b):
        for k, c in enumerate(_b):
            CODON[a + b + c] = _aa[16 * i + 4 * j + k]
LUT = np.full(125, ord("X"), dtype=np.uint8)
_code = {"T": 0, "C": 1, "A": 2, "G": 3}
for cod, aa in CODON.items():
    LUT[_code[cod[0]] * 25 + _code[cod[1]] * 5 + _code[cod[2]]] = ord(aa)
ENC = np.full(256, 4, dtype=np.uint8)
for ch, v in _code.items():
    ENC[ord(ch)] = v
    ENC[ord(ch.lower())] = v
COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def read_fasta(path):
    seqs, name, buf = {}, None, []
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def translate(seq):
    a = ENC[np.frombuffer(seq.encode(), dtype=np.uint8)]
    n = (len(a) // 3) * 3
    a = a[:n].reshape(-1, 3).astype(np.int32)
    idx = a[:, 0] * 25 + a[:, 1] * 5 + a[:, 2]
    bad = (a >= 4).any(axis=1)
    idx[bad] = 0
    prot = LUT[idx]
    prot[bad] = ord("X")
    return prot.tobytes().decode()


def orf_candidates(seqs):
    """Yield (contig, stop_nt_pos, strand, class, protein_tail) for CAAX ORFs."""
    out = []
    for ctg, s in seqs.items():
        L = len(s)
        rc = s.translate(COMP)[::-1]
        for strand, seq in (("+", s), ("-", rc)):
            for f in range(3):
                p = translate(seq[f:])
                start = 0
                for m in re.finditer(r"\*", p):
                    seg = p[start:m.start()]
                    start = m.end()
                    if len(seg) < 20:
                        continue
                    if seg[-4] != "C":
                        continue
                    if "M" not in seg[max(0, len(seg) - 130):len(seg) - 19]:
                        continue
                    cls = "T" if MOTIF_T.search(seg) else ("L" if MOTIF_L.search(seg) else "C4")
                    if cls != "T" and MOTIF_T2.search(seg):
                        cls = cls + "+T2"
                    aa_stop = m.start()
                    nt = f + 3 * aa_stop  # stop codon start in `seq` coords
                    pos = nt if strand == "+" else L - nt - 3
                    out.append((ctg, pos, strand, cls, seg[-12:]))
    return out


def run(cmd, **kw):
    return subprocess.run(cmd, check=True, **kw)


def miniprot_loci(genome_gz, tmp, threads=16):
    paf = os.path.join(tmp, "ste3.paf")
    with open(paf, "w") as fo:
        run([f"{BIN}/miniprot", "-t", str(threads), "-I", "--outn=100", genome_gz, STE3],
            stdout=fo, stderr=subprocess.DEVNULL)
    alns = []
    for line in open(paf):
        f = line.rstrip("\n").split("\t")
        if len(f) < 12:
            continue
        q, qlen, qs, qe, strand, t, tlen, ts, te, mat, blen = f[0], int(f[1]), int(f[2]), int(f[3]), f[4], f[5], int(f[6]), int(f[7]), int(f[8]), int(f[9]), int(f[10])
        tags = dict(x.split(":", 2)[0::2] for x in f[12:] if x.count(":") >= 2)
        score = int(tags.get("AS", 0))
        qcov = (qe - qs) / qlen
        ident = mat / max(1, blen)
        if qcov < 0.5 or (te - ts) > 8000:
            continue
        alns.append((t, strand, ts, te, q, qcov, ident, score))
    loci = []
    by = defaultdict(list)
    for a in alns:
        by[(a[0], a[1])].append(a)
    for (t, strand), lst in by.items():
        lst.sort(key=lambda x: x[2])
        cur = None
        for a in lst:
            if cur and a[2] <= cur["end"] - 100:
                cur["end"] = max(cur["end"], a[3])
                cur["alns"].append(a)
            else:
                if cur:
                    loci.append(cur)
                cur = {"contig": t, "strand": strand, "start": a[2], "end": a[3], "alns": [a]}
        if cur:
            loci.append(cur)
    for L in loci:
        best = max(L["alns"], key=lambda a: a[7])
        L["best_query"], L["best_qcov"], L["best_ident"], L["best_score"] = best[4], best[5], best[6], best[7]
        # best curated-reference hit
        refs = [a for a in L["alns"] if a[4].startswith("REF|")]
        if refs:
            r = max(refs, key=lambda a: a[7])
            L["best_ref"], L["best_ref_ident"] = r[4], r[6]
        else:
            L["best_ref"], L["best_ref_ident"] = "", 0.0
        L["n_queries"] = len({a[4] for a in L["alns"]})
    return loci


def homology_hits(genome_fa, tmp):
    db = os.path.join(tmp, "gdb")
    run([f"{BIN}/makeblastdb", "-in", genome_fa, "-dbtype", "nucl", "-out", db],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    out = os.path.join(tmp, "phero.tsv")
    run([f"{BIN}/tblastn", "-query", PHERO, "-db", db, "-seg", "no", "-word_size", "2",
         "-evalue", "1", "-outfmt", "6 qseqid sseqid pident length evalue bitscore sstart send",
         "-num_threads", "16", "-max_target_seqs", "500", "-out", out])
    hits = []
    for line in open(out):
        q, s, pid, ln, ev, bs, ss, se = line.split("\t")
        hits.append((s, min(int(ss), int(se)), q.split("|")[0] + "|" + q.split("name=")[-1].split("|")[0], float(ev)))
    return hits


def main():
    asm, outdir = sys.argv[1], sys.argv[2]
    # record-id prefixes (e.g. "5346_,5334_") whose pheromones are excluded from
    # the held-out homology count Hx (same genus as this genome)
    excl = tuple(x for x in (sys.argv[3].split(",") if len(sys.argv) > 3 else []) if x)
    os.makedirs(outdir, exist_ok=True)
    tmp = os.path.join(os.environ.get("SCRATCH", "/tmp"), "phero_" + asm)
    os.makedirs(tmp, exist_ok=True)
    gz = os.path.join(LIB, asm + ".fa.gz")
    fa = os.path.join(tmp, asm + ".fa")
    seqs = read_fasta(gz)
    with open(fa, "w") as fo:
        for k, v in seqs.items():
            fo.write(f">{k}\n{v}\n")
    loci = miniprot_loci(gz, tmp)
    cands = orf_candidates(seqs)
    homs = homology_hits(fa, tmp)
    # index
    idx = defaultdict(lambda: defaultdict(list))
    tails = defaultdict(list)   # contig -> sorted (pos, tail8) of all C-4 ORFs
    for ctg, pos, strand, cls, tail in cands:
        base = cls.split("+")[0]
        if base == "T":
            idx["T"][ctg].append(pos); idx["T2"][ctg].append(pos)
        if "T2" in cls:
            idx["T2"][ctg].append(pos)
        if base in ("T", "L"):
            idx["L"][ctg].append(pos)
        tails[ctg].append((pos, tail[-8:]))
    for ctg in tails:
        tails[ctg].sort()
    for ctg, pos, q, ev in homs:
        idx["H1"][ctg].append(pos)
        if ev <= 0.01:
            idx["H01"][ctg].append(pos)
        if not (excl and q.startswith(excl)):
            idx["Hx"][ctg].append(pos)
    for cls in idx:
        for ctg in idx[cls]:
            idx[cls][ctg].sort()

    def count(cls, ctg, a, b):
        if cls == "R":   # max multiplicity of an identical C-4 ORF tail (8 aa) in the window
            arr = tails.get(ctg, [])
            i = bisect.bisect_left(arr, (a, "")); j = bisect.bisect_right(arr, (b, "~"))
            from collections import Counter
            c = Counter(t for _, t in arr[i:j])
            return max(c.values()) if c else 0
        arr = idx[cls].get(ctg, [])
        return bisect.bisect_right(arr, b) - bisect.bisect_left(arr, a)

    classes = ["T", "T2", "L", "R", "H1", "H01", "Hx"]
    with open(os.path.join(outdir, asm + ".loci.tsv"), "w") as fo:
        hdr = ["asm", "contig", "strand", "start", "end", "best_query", "best_qcov", "best_ident",
               "best_score", "best_ref", "best_ref_ident", "n_queries"]
        for W in WINDOWS:
            hdr += [f"{c}_{W//1000}kb" for c in classes]
        fo.write("\t".join(hdr) + "\n")
        for L in sorted(loci, key=lambda x: (x["contig"], x["start"])):
            row = [asm, L["contig"], L["strand"], L["start"], L["end"], L["best_query"],
                   f"{L['best_qcov']:.2f}", f"{L['best_ident']:.3f}", L["best_score"],
                   L["best_ref"], f"{L['best_ref_ident']:.3f}", L["n_queries"]]
            for W in WINDOWS:
                row += [count(c, L["contig"], L["start"] - W, L["end"] + W) for c in classes]
            fo.write("\t".join(map(str, row)) + "\n")
    # random windows of median locus length + 2W, avoiding STE3 loci
    rng = random.Random(1)
    loc_len = int(np.median([L["end"] - L["start"] for L in loci])) if loci else 2000
    ctgs = [(k, len(v)) for k, v in seqs.items() if len(v) > loc_len + 2 * max(WINDOWS) + 1000]
    tot = sum(l for _, l in ctgs)
    avoid = defaultdict(list)
    for L in loci:
        avoid[L["contig"]].append((L["start"] - max(WINDOWS), L["end"] + max(WINDOWS)))
    with open(os.path.join(outdir, asm + ".random.tsv"), "w") as fo:
        hdr = ["asm", "contig", "start", "end"]
        for W in WINDOWS:
            hdr += [f"{c}_{W//1000}kb" for c in classes]
        fo.write("\t".join(hdr) + "\n")
        n = tries = 0
        while n < N_RANDOM and tries < N_RANDOM * 20 and ctgs:
            tries += 1
            r = rng.randrange(tot)
            for k, l in ctgs:
                if r < l:
                    break
                r -= l
            s = max(max(WINDOWS), min(r, l - loc_len - max(WINDOWS)))
            e = s + loc_len
            if any(a < e + max(WINDOWS) and b > s - max(WINDOWS) for a, b in avoid.get(k, [])):
                continue
            row = [asm, k, s, e]
            for W in WINDOWS:
                row += [count(c, k, s - W, e + W) for c in classes]
            fo.write("\t".join(map(str, row)) + "\n")
            n += 1
    with open(os.path.join(outdir, asm + ".cand.tsv"), "w") as fo:
        fo.write("contig\tpos\tstrand\tclass\ttail\n")
        for c in cands:
            if not c[3].startswith("C4") or "T2" in c[3]:
                fo.write("\t".join(map(str, c)) + "\n")
        for ctg, pos, q, ev in homs:
            fo.write(f"{ctg}\t{pos}\t.\tH:{q}:{ev:g}\t\n")
    gsize = sum(len(v) for v in seqs.values())
    print(f"{asm}\t{gsize}\tloci={len(loci)}\tC4={len(cands)}\tT={sum(1 for c in cands if c[3]=='T')}\tH1={len(homs)}")
    for p in (fa,):
        os.remove(p)


if __name__ == "__main__":
    main()
