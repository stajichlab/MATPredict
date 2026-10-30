"""Six-frame scan of the whole S. schoenii genome (GCA_054953955.1) for MATa2-like (HMG-like) and alpha2-like / a1-like ORFs (stop-to-stop segments >= 60 aa,
default HMMER filters, E <= 1e-2). Lists hits within 15 kb of a MAT homeodomain/alpha-box hit as 'cassette neighbours'."""
import itertools, pyhmmer, sys, collections
from pyhmmer.easel import Alphabet, TextSequence
CODON = {a + b + c: aa for (a, b, c), aa in zip(itertools.product("TCAG", repeat=3), "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}
def rc(s): return s[::-1].translate(str.maketrans("ACGTacgt", "TGCAtgca"))
def tr(s): return "".join(CODON.get(s[i:i + 3].upper(), "X") for i in range(0, len(s) - 2, 3))
g = {}; n = None
for l in open("work/ss_genome.fa"):
    if l.startswith(">"): n = l[1:].split()[0]; g[n] = []
    else: g[n].append(l.strip())
g = {k: "".join(v) for k, v in g.items()}
alpha = Alphabet.amino(); H = {}
for n in ("a2", "a2_no_Pichiaceae", "a2_no_Saccharomycetaceae", "sacch_a2", "a1", "alpha2", "alpha1"):
    p = "../../db/Ascomycota/classifiers/MATa2_saccharomycotina/MATa2_saccharomycotina.hmm" if n == "sacch_a2" else f"hmm/{n}.hmm"
    with pyhmmer.plan7.HMMFile(p) as hf: H[n] = hf.read()
segs = []
for c, s in g.items():
    if len(s) < 200000: continue      # chromosomes only
    for strand, t in (("+", s), ("-", rc(s))):
        for f in range(3):
            p = tr(t[f:]); pos = 0
            for seg in p.split("*"):
                if len(seg) >= 60:
                    n0 = f + 3 * pos; n1 = n0 + 3 * len(seg) - 1
                    a, b = (n0 + 1, n1 + 1) if strand == "+" else (len(s) - n1, len(s) - n0)
                    segs.append((f"{c}:{a}-{b}({strand})", seg))
                pos += len(seg) + 1
print("segments", len(segs), file=sys.stderr)
q = [TextSequence(name=k.encode(), sequence=v).digitize(alpha) for k, v in segs]
hits = collections.defaultdict(dict)
for n, h in H.items():
    for th in pyhmmer.hmmsearch([h], q, E=1e-2, cpus=4):
        for x in th: hits[str(x.name)][n] = round(x.score, 1)
rows = sorted(hits.items(), key=lambda kv: -max(kv[1].get("a2", 0), kv[1].get("sacch_a2", 0)))
with open("data/ss_genome_scan.tsv", "w") as fh:
    fh.write("segment\t" + "\t".join(H) + "\n")
    for k, v in rows: fh.write(k + "\t" + "\t".join(str(v.get(n, "")) for n in H) + "\n")
print("top a2 / sacch_a2:")
for k, v in rows[:12]: print(k, v)
