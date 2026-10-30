"""Parse AUGUSTUS GFF3 predictions on the MAT windows (windows.fa), convert to genome coordinates, and classify each predicted protein
(Pfam MAT domains, phmmer best reference, Saccharomycotina MATa2 profile). Usage: classify_aug.py data/aug_*.gff3 > table"""
import re, sys, pyhmmer
from pyhmmer.easel import SequenceFile, Alphabet, TextSequence
def _s(x): return x.decode() if isinstance(x, bytes) else x
win = {}
for l in open("data/windows.fa"):
    if l.startswith(">"):
        n, c, a, b = l[1:].split(); win[n] = (c, int(a), int(b))
genes = []  # (species, window, id, strand, [cds], protein)
for fn in sys.argv[1:]:
    sp = fn.split("aug_")[1].rsplit(".gff3", 1)[0]
    cur = None; prot = None
    for l in open(fn):
        if l.startswith("#"):
            if "protein sequence = [" in l:
                prot = l.split("[", 1)[1].strip()
                if prot.endswith("]"): prot = prot[:-1]; cur_p = prot; genes[-1][-1] = cur_p; prot = None
                else: prot = prot
            elif prot is not None:
                t = l[1:].strip()
                if t.endswith("]"): prot += t[:-1]; genes[-1][-1] = prot; prot = None
                else: prot += t
            continue
        f = l.rstrip("\n").split("\t")
        if len(f) < 9: continue
        if f[2] == "transcript":
            genes.append([sp, f[0], re.search(r"ID=([^;]+)", f[8]).group(1), f[6], [], ""])
        elif f[2] == "CDS":
            genes[-1][4].append((int(f[3]), int(f[4])))
alpha = Alphabet.amino()
ts = [TextSequence(name=f"{g[0]}|{g[1]}|{g[2]}".encode(), sequence=g[5]).digitize(alpha) for g in genes if g[5]]
with SequenceFile("data/all_refs.faa", digital=True, alphabet=alpha) as sf: db = sf.read_block()
best = {}
for hits in pyhmmer.phmmer(ts, db, E=1e-3):
    if len(hits): best[_s(hits.query.name)] = (_s(hits[0].name), hits[0].score)
want = {"HMG_box", "HMG_box_2", "Homeodomain", "Homeobox_KN", "MATalpha_HMGbox", "I_LWEQ", "ANTH", "ENTH"}
dom = {}
with pyhmmer.plan7.HMMFile("/Users/jstajich/db/funannotate/Pfam-A.hmm") as hf:
    for hits in pyhmmer.hmmscan(ts, hf, E=1e-1):
        dom[_s(hits.query.name)] = [(_s(h.name), round(h.score, 1)) for h in hits if _s(h.name) in want]
with pyhmmer.plan7.HMMFile("/Users/jstajich/projects/MATPredict/db/Ascomycota/classifiers/MATa2_saccharomycotina/MATa2_saccharomycotina.hmm") as hf:
    a2 = {}
    for hits in pyhmmer.hmmsearch(list(hf), ts, E=1e-1):
        for h in hits: a2[_s(h.name)] = round(h.score, 1)
print("species\twindow\tgene\tgstart\tgend\tstrand\tcds\taa\tbest_ref\tbits\tpfam\tMATa2_profile")
for g in genes:
    k = f"{g[0]}|{g[1]}|{g[2]}"
    c, a, b = win[g[1]]
    cds = sorted((a + x - 1, a + y - 1) for x, y in g[4])
    br = best.get(k); nm = ""
    if br: nm = (br[0].split("name=")[1].split("|")[0] if "name=" in br[0] else br[0])[-40:]
    print("\t".join(map(str, [g[0], g[1], g[2], cds[0][0], cds[-1][1], g[3], ";".join(f"{x}-{y}" for x, y in cds), len(g[5]), nm, f"{br[1]:.0f}" if br else "", dom.get(k, ""), a2.get(k, "")])))
