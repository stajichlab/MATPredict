"""Genome-wide homeodomain scan (Pfam Homeodomain PF00046 and Homeobox_KN PF05920 via pyhmmer) on six-frame ORFs >= 50 aa,
then divergent HD-domain pairs within 4 kb. usage: hd_scan.py ACC FASTA.gz OUT.json"""
import gzip, json, os, sys, shutil, tempfile
import pyhmmer
from Bio.Seq import Seq

acc, fa, out = sys.argv[1:4]
abc = pyhmmer.easel.Alphabet.amino()
nm = lambda x: x.decode() if isinstance(x, bytes) else x
hmms = []
with pyhmmer.plan7.HMMFile("/bigdata/stajichlab/shared/lib/funannotate_db/Pfam-A.hmm") as hf:
    for h in hf:
        if nm(h.name) in ("Homeodomain", "Homeobox_KN"): hmms.append(h)
        if len(hmms) == 2: break


def mkseq(name, o):
    try: return pyhmmer.easel.TextSequence(name=name.encode(), sequence=o).digitize(abc)
    except TypeError: return pyhmmer.easel.TextSequence(name=name, sequence=o).digitize(abc)


hits = []  # (contig, strand, gstart, gend, domain, score)
meta = {}
batch = []


def flush():
    global batch
    if not batch: return
    for th in pyhmmer.hmmsearch(hmms, batch, E=1e-4, cpus=2):
        dom = nm(th.query.name)
        for hit in th:
            d = hit.best_domain
            name = nm(hit.name)
            contig, strand, gs, ge, orfaa0 = meta[name]
            # domain aa coords inside the ORF -> genomic nt
            a0, a1 = d.env_from - 1, d.env_to
            if strand == "+": s, e = gs + 3 * a0, gs + 3 * a1 - 1
            else: s, e = ge - 3 * a1 + 1, ge - 3 * a0
            hits.append((contig, strand, s, e, dom, round(hit.score, 1)))
    batch = []


def contigs(path):
    cur = None; buf = []
    for line in gzip.open(path, "rt"):
        if line.startswith(">"):
            if cur is not None: yield cur, "".join(buf).upper()
            cur = line[1:].split()[0]; buf = []
        else: buf.append(line.strip())
    if cur is not None: yield cur, "".join(buf).upper()


n = 0
for contig, seq in contigs(fa):
    L = len(seq)
    if L < 200: continue
    for strand in "+-":
        s_ = seq if strand == "+" else str(Seq(seq).reverse_complement())
        for f in range(3):
            prot = str(Seq(s_[f: L - ((L - f) % 3)]).translate())
            pos = 0
            for orf in prot.split("*"):
                if len(orf) >= 50:
                    # ORF aa offset pos -> nt start in s_ coordinates f + 3*pos
                    ns = f + 3 * pos
                    if strand == "+": gs, ge = ns + 1, ns + 3 * len(orf)
                    else: gs, ge = L - (ns + 3 * len(orf)) + 1, L - ns
                    name = f"o{n}"; n += 1
                    meta[name] = (contig, strand, gs, ge, pos)
                    batch.append(mkseq(name, orf))
                    if len(batch) >= 40000: flush(); meta.clear()
                pos += len(orf) + 1
flush()
# dedupe hits per domain locus, then pairs
hits.sort()
pairs = []
byc = {}
for h in hits: byc.setdefault(h[0], []).append(h)
for c, hs in byc.items():
    for i, a in enumerate(hs):
        for b in hs[i + 1:]:
            if b[2] - a[3] > 4000: break
            if a[1] == "-" and b[1] == "+" and a[4] != "" and not (a[2] == b[2]):
                pairs.append(dict(contig=c, left=[a[2], a[3], a[4], a[5]], right=[b[2], b[3], b[4], b[5]], gap=max(0, b[2] - a[3])))
json.dump(dict(acc=acc, orfs=n, hd_hits=len(hits), hits=hits[:2000], divergent_pairs=pairs), open(out, "w"))
