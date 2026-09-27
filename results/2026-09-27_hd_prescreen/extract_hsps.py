"""Step 3a: translated tblastn HSP segments for each candidate cluster -- the
evidence available BEFORE polishing. Admitted clusters (all) plus a random
sample of non-admitted negatives. Queries = the 10 redHD proteins; same
tblastn settings as detect (-seg no, evalue 10). Region = cluster +-500 bp."""
import csv, gzip, os, random, subprocess, sys, tempfile, collections
from concurrent.futures import ThreadPoolExecutor
E = "/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin"
GEN = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes"
SCR = sys.argv[1]
random.seed(1)
rows = list(csv.DictReader(open("clusters.tsv"), delimiter="\t"))
keep = [r for r in rows if r["admitted"] == "True"]
nonadm = [r for r in rows if r["admitted"] != "True"]
keep += random.sample(nonadm, 1500)
by_g = collections.defaultdict(list)
for i, r in enumerate(keep):
    r["cid"] = f"C{i:05d}"; by_g[r["genome"]].append(r)
refs = open("refs_hd.faa").read()
q = "".join(b for b in (">" + x for x in refs.split(">")[1:]) if "redHD" in b.split("\n")[0])
open(f"{SCR}/q.faa", "w").write(q)

def fasta(path):
    name, seq = None, []
    with gzip.open(path, "rt") as fh:
        for l in fh:
            if l.startswith(">"):
                if name: yield name, "".join(seq)
                name, seq = l[1:].split()[0], []
            else: seq.append(l.strip())
    if name: yield name, "".join(seq)

def one(g):
    cl = by_g[g]; need = {r["contig"] for r in cl}
    seqs = {n: s for n, s in fasta(f"{GEN}/{g}.fa.gz") if n in need}
    fa = f"{SCR}/{g}.regions.fa"
    with open(fa, "w") as fo:
        for r in cl:
            s = seqs[r["contig"]]; a = max(0, int(r["start"]) - 500); b = min(len(s), int(r["end"]) + 500)
            fo.write(f">{r['cid']}|{a}\n{s[a:b]}\n")
    out = subprocess.run([f"{E}/tblastn", "-query", f"{SCR}/q.faa", "-subject", fa, "-seg", "no", "-evalue", "10",
                          "-outfmt", "6 qseqid sseqid sstart send evalue bitscore pident sseq"],
                         capture_output=True, text=True, check=True).stdout
    os.remove(fa)
    return out

with ThreadPoolExecutor(4) as ex, open("hsps.tsv", "w") as fo:
    for out in ex.map(one, sorted(by_g)):
        fo.write(out)
with open("clusters_sampled.tsv", "w", newline="") as fo:
    w = csv.DictWriter(fo, fieldnames=list(keep[0]), delimiter="\t"); w.writeheader(); w.writerows(keep)
print(len(keep), "clusters;", sum(1 for _ in open("hsps.tsv")), "HSP rows")
