#!/usr/bin/env python3
"""Pipeline-independent scan for the genes of a mating locus in one genome.

Looks for three gene classes without using the detect pipeline:
  STE3  pheromone receptor: miniprot of 1,095 STE3 queries (ste3_all.faa,
        results/2026-10-05_caax_receptor_test) and pyhmmer on a six-frame
        stop-to-stop translation with Pfam PF02076.
  HD    homeodomain: miniprot of hd_queries.faa (curated db HD1/HD2/SXI/bE/bW
        records plus Phaffia rhodozyma HD1/HD2) and pyhmmer on the six-frame
        translation with Pfam PF00046 (Homeobox) and PF05920 (Homeobox_KN).
  CAAX  strict-CAAX short ORFs (scan_genome.orf_candidates, motif
        C[VI][IV][AVMG]) so a receptor can be tied to a precursor.

Usage: scan_ste3_hd.py ASMID OUTDIR
Writes OUTDIR/<ASMID>.genes.tsv: kind, contig, strand, start, end, source, score, detail.
Run with the pixi python (pyhmmer, miniprot, numpy on PATH).
"""
import gzip
import os
import re
import subprocess
import sys
from collections import defaultdict

import numpy as np
import pyhmmer

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import scan_genome as sg  # copy of results/2026-10-05_caax_receptor_test/scan_genome.py

LIB = sg.LIB
BIN = sg.BIN
STE3Q = os.path.join(HERE, "ste3_all.faa")
HDQ = os.path.join(HERE, "hd_queries.faa")
PFAM = "/bigdata/stajichlab/shared/lib/funannotate_db/Pfam-A.hmm"
WANT = {"PF02076": "STE3", "PF00046": "HD", "PF05920": "HD"}
MIN_ORF = 30
DOM_E = 1e-4


def _s(x):
    if x is None:
        return ""
    return x.decode() if isinstance(x, bytes) else str(x)


def load_hmms():
    """The three Pfam models; cached in pfam3.hmm next to this script after the
    first run (the full Pfam-A file takes about 30 s to scan)."""
    cache = os.path.join(HERE, "pfam3.hmm")
    src = cache if os.path.exists(cache) else PFAM
    out = []
    with pyhmmer.plan7.HMMFile(src) as hf:
        for hmm in hf:
            acc = _s(hmm.accession).split(".")[0]
            if acc in WANT:
                out.append((acc, hmm))
    if src == PFAM and len(out) == len(WANT):
        tmpf = cache + ".%d" % os.getpid()
        with open(tmpf, "wb") as fo:
            for _, hmm in out:
                hmm.write(fo)
        os.replace(tmpf, cache)
    return out


def six_frame_segments(seqs):
    """Yield (name, protein, contig, strand, frame, aa_start) for stop-to-stop ORF segments."""
    for ctg, s in seqs.items():
        L = len(s)
        rc = s.translate(sg.COMP)[::-1]
        for strand, seq in (("+", s), ("-", rc)):
            for f in range(3):
                p = sg.translate(seq[f:])
                start = 0
                for m in re.finditer(r"\*", p):
                    seg = p[start:m.start()]
                    if len(seg) >= MIN_ORF:
                        yield ctg, strand, f, start, seg, L
                    start = m.end()
                seg = p[start:]
                if len(seg) >= MIN_ORF:
                    yield ctg, strand, f, start, seg, L


def hmm_hits(seqs):
    hmms = load_hmms()
    segs = []
    meta = []
    for i, (ctg, strand, f, aa0, seg, L) in enumerate(six_frame_segments(seqs)):
        name = f"s{i}".encode()
        segs.append(pyhmmer.easel.TextSequence(name=name, sequence=seg).digitize(pyhmmer.easel.Alphabet.amino()))
        meta.append((ctg, strand, f, aa0, L))
    hits = []
    if not segs:
        return hits
    sdb = pyhmmer.easel.DigitalSequenceBlock(pyhmmer.easel.Alphabet.amino(), segs)
    for acc, hmm in hmms:
        for top in pyhmmer.hmmsearch(hmm, sdb, E=1e-2, domE=DOM_E, cpus=1):
            for hit in top:
                i = int(_s(hit.name)[1:])
                ctg, strand, f, aa0, L = meta[i]
                for dom in hit.domains:
                    if dom.i_evalue > DOM_E:
                        continue
                    a = aa0 + dom.env_from - 1          # 0-based aa start in frame
                    b = aa0 + dom.env_to                # exclusive
                    nt0 = f + 3 * a
                    nt1 = f + 3 * b
                    if strand == "+":
                        s, e = nt0 + 1, nt1
                    else:
                        s, e = L - nt1 + 1, L - nt0
                    hits.append((WANT[acc], ctg, strand, s, e, "hmm:" + acc, dom.score, f"{dom.i_evalue:.1e}"))
    return hits


def merge_loci(hits, gap):
    by = defaultdict(list)
    for h in hits:
        by[(h[0], h[1], h[2])].append(h)
    out = []
    for (kind, ctg, strand), lst in by.items():
        lst.sort(key=lambda h: h[3])
        cur = None
        for h in lst:
            if cur and h[3] <= cur[4] + gap:
                cur[4] = max(cur[4], h[4])
                if h[6] > cur[6]:
                    cur[6], cur[5], cur[7] = h[6], h[5], h[7]
            else:
                if cur:
                    out.append(tuple(cur))
                cur = list(h)
        if cur:
            out.append(tuple(cur))
    return out


def miniprot_paf(genome_gz, queries, tmp, tag, threads=2):
    paf = os.path.join(tmp, tag + ".paf")
    with open(paf, "w") as fo:
        sg.run([f"{BIN}/miniprot", "-t", str(threads), "-I", "--outn=50", genome_gz, queries],
               stdout=fo, stderr=subprocess.DEVNULL)
    res = []
    for line in open(paf):
        f = line.rstrip("\n").split("\t")
        if len(f) < 12:
            continue
        q, qlen, qs, qe, strand, t, tlen, ts, te, mat, blen = (f[0], int(f[1]), int(f[2]), int(f[3]), f[4], f[5],
                                                                 int(f[6]), int(f[7]), int(f[8]), int(f[9]), int(f[10]))
        tags = dict(x.split(":", 2)[0::2] for x in f[12:] if x.count(":") >= 2)
        res.append((t, strand, ts + 1, te, q, (qe - qs) / qlen, mat / max(1, blen), int(tags.get("AS", 0)), qe - qs))
    return res


def main():
    asm, outdir = sys.argv[1], sys.argv[2]
    os.makedirs(outdir, exist_ok=True)
    tmp = os.path.join(os.environ.get("SCRATCH", "/tmp"), "scanhd_" + asm)
    os.makedirs(tmp, exist_ok=True)
    gz = os.path.join(LIB, asm + ".fa.gz")
    seqs = sg.read_fasta(gz)
    rows = []

    # miniprot
    for kind, q, minq, tag in (("STE3", STE3Q, 0.5, "ste3"), ("HD", HDQ, 0.0, "hd")):
        hits = []
        for t, strand, s, e, qn, qcov, ident, score, aalen in miniprot_paf(gz, q, tmp, tag):
            if kind == "STE3" and (qcov < minq or (e - s) > 8000):
                continue
            if kind == "HD" and (aalen < 40 or (e - s) > 8000):
                continue
            hits.append((kind, t, strand, s, e, "miniprot", score, f"{qn.split('|')[-1]}:{ident:.2f}:{qcov:.2f}"))
        rows += merge_loci(hits, 100)

    # hmm six frame
    rows += merge_loci(hmm_hits(seqs), 3000)

    # strict CAAX short ORFs
    for ctg, pos, strand, cls, tail in sg.orf_candidates(seqs):
        if cls.split("+")[0] == "T":
            rows.append(("CAAX", ctg, strand, pos + 1, pos + 3, "orf", 0, tail))

    with open(os.path.join(outdir, asm + ".genes.tsv"), "w") as fo:
        fo.write("kind\tcontig\tstrand\tstart\tend\tsource\tscore\tdetail\n")
        for r in sorted(rows, key=lambda r: (r[1], r[3])):
            fo.write("\t".join(map(str, r)) + "\n")
    with open(os.path.join(outdir, asm + ".contigs.tsv"), "w") as fo:
        for k, v in seqs.items():
            fo.write(f"{k}\t{len(v)}\n")
    print(asm, sum(len(v) for v in seqs.values()), len(rows))


if __name__ == "__main__":
    main()
