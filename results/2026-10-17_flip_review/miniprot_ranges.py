"""Where the MAT1-1-1 and MAT1-2-1 references align on the locus (aligned reference residue range), by miniprot.
    miniprot_ranges.py LOCI_TSV JOBS_TSV DB_HELOTIALES REF_RECORD... > out.tsv   (run in the pixi env, miniprot on PATH)
"""
import csv, gzip, subprocess, sys, tempfile, os
from Bio import SeqIO
loci, jobs, db, *refs = sys.argv[1:]
fa = {r[0]: r[4] for r in csv.reader(open(jobs), delimiter="\t")}
print("genome\tref\tquery\tqlen\tqstart\tqend\tident\tcontig_start\tcontig_end\tstrand")
for l in csv.DictReader(open(loci), delimiter="\t"):
    g = l["genome"]
    op = gzip.open if fa[g].endswith(".gz") else open
    seq = None
    with op(fa[g], "rt") as fh:
        for r in SeqIO.parse(fh, "fasta"):
            if r.id == l["contig"]:
                seq = str(r.seq); break
    a, b = max(0, int(l["locus_start"]) - 1 - 2000), int(l["locus_end"]) + 2000
    with tempfile.TemporaryDirectory() as t:
        open(f"{t}/r.fa", "w").write(f">r\n{seq[a:b]}\n")
        for ref in refs:
            p = subprocess.run(["miniprot", "-u", "-N", "3", "--outs=0.9", "--paf", f"{t}/r.fa", f"{db}/{ref}/proteins.faa"] if False else
                               ["miniprot", "-u", "-N", "3", "--outs=0.9", f"{t}/r.fa", f"{db}/{ref}/proteins.faa"], capture_output=True, text=True)
            # default miniprot output is PAF
            for ln in p.stdout.splitlines():
                f = ln.split("\t")
                if len(f) < 12: continue
                q, ql, qs, qe, st, _, _, ts, te = f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7], f[8]
                ident = f[9] and round(int(f[9]) / max(1, int(f[10])), 3)
                print("\t".join(map(str, [g, ref, q, ql, qs, qe, ident, a + int(ts) + 1, a + int(te), st])))
