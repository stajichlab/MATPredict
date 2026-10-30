"""Leave-one-clade-out and leave-one-species-out evaluation of a1 vs alpha2 (and a2, alpha1) profiles; builds the final profiles.
Run on the HPCC (mafft module): python loco.py"""
import subprocess, sys, tempfile, collections, io, json
import pyhmmer
from pyhmmer.easel import Alphabet, TextSequence, TextMSA, DigitalSequenceBlock
from pyhmmer.plan7 import Builder, Background
alpha = Alphabet.amino()
meta = [l.rstrip("\n").split("\t") for l in open("data/set_meta.tsv")][1:]
seqs = {}   # id -> (class, species, clade, seq)
fa = {}
for c in ("a1", "alpha2", "a2", "alpha1"):
    name = None
    for l in open(f"data/set_{c}.faa"):
        if l.startswith(">"): name = l[1:].strip(); fa[name] = ""
        elif name: fa[name] += l.strip()
for c, sp, cl, src, ident, ln in meta:
    seqs[ident.replace(" ", "_")] = (c, sp, cl, fa[ident.replace(" ", "_")])
def mafft(items):
    with tempfile.NamedTemporaryFile("w", suffix=".faa", delete=False) as fh:
        for i, s in items: fh.write(f">{i}\n{s}\n")
        fh.flush()
    out = subprocess.run(["mafft", "--localpair", "--maxiterate", "1000", "--quiet", fh.name], capture_output=True, text=True)
    if out.returncode: sys.exit("mafft failed: " + out.stderr[-500:])
    out = out.stdout
    d = {}; n = None
    for l in out.split("\n"):
        if l.startswith(">"): n = l[1:]; d[n] = ""
        elif n: d[n] += l.strip()
    return d
def build(name, items):
    aln = mafft(items)
    msa = TextMSA(name=name.encode(), sequences=[TextSequence(name=i.encode(), sequence=s.replace(".", "-")) for i, s in aln.items()]).digitize(alpha)
    hmm, _, _ = Builder(alpha, seed=42).build_msa(msa, Background(alpha))
    return hmm
def score(hmm, seq):
    q = TextSequence(name=b"q", sequence=seq).digitize(alpha)
    best = 0.0
    for hits in pyhmmer.hmmsearch([hmm], [q], E=1e6, Z=1000):
        for h in hits: best = max(best, h.score)
    return best
CL = ["Saccharomycetaceae", "Serinales", "Pichiaceae", "Dipodascales", "other"]
def cls_items(c, exclude=lambda v: False):
    return [(i, v[3]) for i, v in seqs.items() if v[0] == c and not exclude(v)]
res = []
for scheme in ("clade", "species"):
    groups = sorted({(v[2] if scheme == "clade" else v[1]) for v in seqs.values() if v[0] in ("a1", "alpha2")})
    for g in groups:
        ex = (lambda v, g=g: (v[2] if scheme == "clade" else v[1]) == g)
        ha1 = build("a1", cls_items("a1", ex)); hal2 = build("alpha2", cls_items("alpha2", ex))
        n1 = len(cls_items("a1", ex)); n2 = len(cls_items("alpha2", ex))
        for i, v in seqs.items():
            if v[0] in ("a1", "alpha2") and ex(v):
                s1, s2 = score(ha1, v[3]), score(hal2, v[3])
                res.append((scheme, g, i, v[0], v[1], round(s1, 1), round(s2, 1), "a1" if s1 > s2 else "alpha2", n1, n2))
with open("data/loco.tsv", "w") as fh:
    fh.write("scheme\theldout\tid\ttrue\tspecies\ta1_score\talpha2_score\tcall\tn_a1_train\tn_alpha2_train\n")
    for r in res: fh.write("\t".join(map(str, r)) + "\n")
import os; os.makedirs("hmm", exist_ok=True)
for c in ("a1", "alpha2", "a2", "alpha1"):
    h = build(c, cls_items(c))
    with open(f"hmm/{c}.hmm", "wb") as fh: h.write(fh)
    print(c, len(cls_items(c)), "sequences, model length", h.M)
