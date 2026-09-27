"""Steps 2-4: leave-one-species-out HD1/HD2 HMMs from the redHD records,
scored against each candidate cluster's tblastn HSP translations (pre-polish
evidence). Record-strain genomes are excluded from evaluation."""
import csv, collections, subprocess, io, time
import pyhmmer
from pyhmmer.easel import Alphabet, TextSequence, TextMSA, DigitalSequenceBlock
AB = Alphabet.amino()
SAMP = {r["ASMID"]: (r["SPECIES"] or r["SPECIES_IN"]) for r in csv.DictReader(open("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"))}
RECORD_CONTIG = ("LCTV02", "CAKKSX03", "JANBVD01", "CAKLCE03", "PQDI01")
# references
refs = []
name = None
for l in open("refs_hd.faa"):
    if l.startswith(">"):
        f = l[1:].strip().split("|"); name = f; refs.append([f, ""])
    else: refs[-1][1] += l.strip()
redhd = [(f[3].split("=")[1], f[-1], s) for f, s in refs if "redHD" in f[0]]
def build(gene, exclude_species):
    seqs = [(f"{sp}_{i}", s) for i, (g, sp, s) in enumerate(redhd) if g == gene and sp != exclude_species]
    fa = "".join(f">{n.replace(' ','_')}\n{s}\n" for n, s in seqs)
    aln = subprocess.run(["mafft", "--quiet", "--localpair", "--maxiterate", "1000", "-"], input=fa, capture_output=True, text=True, check=True).stdout
    names, cur = [], []
    recs = {}
    for l in aln.splitlines():
        if l.startswith(">"): n = l[1:]; recs[n] = ""; names.append(n)
        else: recs[n] += l.strip()
    msa = TextMSA(name=f"{gene}".encode(), sequences=[TextSequence(name=n.encode(), sequence=recs[n].upper()) for n in names]).digitize(AB)
    hmm, _, _ = pyhmmer.plan7.Builder(AB).build_msa(msa, pyhmmer.plan7.Background(AB))
    return hmm, len(seqs)
species_set = sorted({sp for _, sp, _ in redhd})
HMMS = {}
for ex in [None] + species_set:
    HMMS[ex] = {g: build(g, ex) for g in ("HD1", "HD2")}
# clusters + hsps
cl = {r["cid"]: r for r in csv.DictReader(open("clusters_sampled.tsv"), delimiter="\t")}
segs = collections.defaultdict(list)
for l in open("hsps.tsv"):
    q, s, ss, se, ev, bs, pid, sseq = l.rstrip("\n").split("\t")
    cid = s.split("|")[0]
    seq = sseq.replace("-", "").replace("*", "X")
    if len(seq) >= 10: segs[cid].append((q.split("|")[3].split("=")[1], seq, float(bs)))
t0 = time.time(); out = []
for cid, r in cl.items():
    g = r["genome"]; sp = SAMP.get(g, "?")
    is_record = r["contig"].startswith(RECORD_CONTIG)
    ex = sp if sp in species_set else None
    best = {"HD1": -1e9, "HD2": -1e9}
    uniq = {seq for _, seq, _ in segs.get(cid, [])}
    if uniq:
        block = DigitalSequenceBlock(AB, [TextSequence(name=f"s{i}".encode(), sequence=s).digitize(AB) for i, s in enumerate(uniq)])
        for gene in ("HD1", "HD2"):
            hmm = HMMS[ex][gene][0]
            for hits in pyhmmer.hmmsearch([hmm], block, E=1e9, T=None, Z=1):
                for h in hits:
                    best[gene] = max(best[gene], h.score)
    out.append(dict(cid=cid, genome=g, species=sp, label=r["label"], admitted=r["admitted"], polish_capped=r["polish_capped"],
                    record_genome=is_record, n_hsp=len(segs.get(cid, [])), max_blast_bits=max([b for *_, b in segs.get(cid, [])] or [0]),
                    hd1=round(best["HD1"], 1), hd2=round(best["HD2"], 1), score=round(max(best.values()), 1), loso_excluded=ex or ""))
with open("scores.tsv", "w", newline="") as fo:
    w = csv.DictWriter(fo, fieldnames=list(out[0]), delimiter="\t"); w.writeheader(); w.writerows(out)
print("scored", len(out), "clusters in", round(time.time() - t0, 1), "s; HMM training sizes", {k: {g: v[1] for g, v in d.items()} for k, d in HMMS.items()})
