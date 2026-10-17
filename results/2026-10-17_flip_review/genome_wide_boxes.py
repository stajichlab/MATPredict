"""Genome-wide search for alpha-box and HMG-box proteins: six-frame ORFs (>= 70 aa between stops) scored on the
MAT_core profiles; prints the best hits (score >= MIN) per genome. Tells whether a genome has an alpha box (or HMG box)
anywhere beyond the locus examined.
    genome_wide_boxes.py JOBS_TSV PROFILES_DIR OUT.tsv GENOME_ID...
"""
import csv, gzip, sys, re
from pathlib import Path
import pyhmmer
from Bio import SeqIO
from Bio.Seq import Seq
sys.path.insert(0, str(Path(__file__).resolve().parent))
jobs, prof, out, *ids = sys.argv[1:]
fa = {r[0]: r[4] for r in csv.reader(open(jobs), delimiter="\t")}
hm = {}
for gname in ("MAT1-1-1", "MAT1-2-1"):
    with pyhmmer.plan7.HMMFile(f"{prof}/{gname}.hmm") as fh:
        hm[gname] = fh.read()
ab = pyhmmer.easel.Alphabet.amino()
w = csv.writer(open(out, "w"), delimiter="\t", lineterminator="\n")
w.writerow("genome profile score contig frame_start orf_aa".split())
for gid in ids:
    op = gzip.open if fa[gid].endswith(".gz") else open
    seqs, meta = [], {}
    with op(fa[gid], "rt") as fh:
        for r in SeqIO.parse(fh, "fasta"):
            s = r.seq
            for strand, ss in (("+", s), ("-", s.reverse_complement())):
                for f in range(3):
                    prot = str(ss[f:][: (len(ss) - f) // 3 * 3].translate())
                    pos = 0
                    for orf in prot.split("*"):
                        if len(orf) >= 70:
                            k = f"{r.id}|{strand}|{f}|{pos}"
                            seqs.append(pyhmmer.easel.TextSequence(name=k.encode(), sequence=orf.replace("X", "A")).digitize(ab))
                            meta[k] = len(orf)
                        pos += len(orf) + 1
    block = pyhmmer.easel.DigitalSequenceBlock(ab, seqs)
    for gname, h in hm.items():
        pl = pyhmmer.plan7.Pipeline(ab, Z=1, E=1e9, domE=1e9, bias_filter=False)
        hits = sorted(pl.search_hmm(h, block), key=lambda x: -x.score)[:4]
        for x in hits:
            n = x.name.decode() if isinstance(x.name, bytes) else x.name
            if x.score >= 20:
                w.writerow([gid, gname, round(x.score, 1), n, "", meta[n]])
    print(gid, file=sys.stderr)
