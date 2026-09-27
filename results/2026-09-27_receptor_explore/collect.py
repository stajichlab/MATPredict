"""Collect PF02076 hits from sampled proteomes + curated mating receptors."""
import csv, glob, os, collections, re
samp = {r["ASMID"]: r for r in csv.DictReader(open("sample.tsv"), delimiter="\t")}
hits = collections.defaultdict(dict)  # asm -> prot -> (score, envfrom, envto)
for f in glob.glob("hmm/*.domtbl"):
    asm = os.path.basename(f)[:-7]
    for l in open(f):
        if l.startswith("#"): continue
        x = l.split()
        pid, score, ef, et = x[0], float(x[7]), int(x[19]), int(x[20])
        if pid not in hits[asm] or score > hits[asm][pid][0]:
            hits[asm][pid] = (score, ef, et)
def fasta(path):
    name = None; seq = []
    for l in open(path):
        if l.startswith(">"):
            if name: yield name, "".join(seq)
            name = l[1:].split()[0]; seq = []
        else: seq.append(l.strip())
    if name: yield name, "".join(seq)
out = open("ste3_all.faa", "w"); tab = open("ste3_hits.tsv", "w")
tab.write("seq_id\tASMID\tspecies\tstrain\tclass\torder\tprotein\tlength\tpfam_score\tkind\n")
n = 0
for asm, prots in hits.items():
    r = samp[asm]
    for pid, s in fasta(r["proteome"]):
        if pid in prots:
            sid = f"{asm}|{pid}"
            out.write(f">{sid}\n{s.rstrip('*')}\n"); n += 1
            tab.write("\t".join([sid, asm, r["SPECIES"], r["STRAIN"], r["CLASS"], r["ORDER"], pid, str(len(s)), str(prots[pid][0]), "query"]) + "\n")
known = 0
for name, s in fasta("db_basidio_proteins.faa"):
    pass
cur = None
for l in open("db_basidio_proteins.faa"):
    pass
seqs = list(fasta("db_basidio_proteins.faa"))
# fasta() drops header fields after the first space; re-read headers to get gene names
hdr = [l[1:].strip() for l in open("db_basidio_proteins.faa") if l.startswith(">")]
for h, (name, s) in zip(hdr, seqs):
    g = re.search(r"name=([^|]+)", h).group(1)
    if g in ("STE3", "STE3a1", "STE3a2", "STE3v2", "pra1", "pra2", "bar3", "bbr2", "pheromone_receptor"):
        rec = h.split("|")[0]; gi = re.search(r"gene_index=(\d+)", h).group(1)
        sid = f"REF|{rec}|{g}|{gi}"
        out.write(f">{sid}\n{s.rstrip('*')}\n"); known += 1
        tab.write("\t".join([sid, rec, rec, "", "", "", g, str(len(s)), "", "known_mating"]) + "\n")
print("query STE3", n, "known mating receptors", known)
