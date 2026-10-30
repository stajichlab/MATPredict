"""Six-frame segments (stop to stop) of B. californica PB4207 JAMAIE010000005.1 618,500-623,500 and of the other three californica assemblies'
homologous region scored with the alpha2, a1 and alpha1 profiles; coordinates 1-based."""
import itertools, pyhmmer, sys
from pyhmmer.easel import Alphabet, TextSequence
CODON = {a + b + c: aa for (a, b, c), aa in zip(itertools.product("TCAG", repeat=3), "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}
def rc(s): return s[::-1].translate(str.maketrans("ACGTacgt", "TGCAtgca"))
def tr(s): return "".join(CODON.get(s[i:i + 3].upper(), "X") for i in range(0, len(s) - 2, 3))
def fasta(p):
    d = {}; n = None
    for l in open(p):
        if l.startswith(">"): n = l[1:].split()[0]; d[n] = ""
        elif n: d[n] += l.strip()
    return d
alpha = Alphabet.amino(); H = {}
for n in ("alpha2", "a1", "alpha1"):
    with pyhmmer.plan7.HMMFile(f"hmm/{n}.hmm") as hf: H[n] = hf.read()
def segs(seq, a, b, off):
    sub = seq[a - off:b - off + 1]; out = []
    for strand, s in (("+", sub), ("-", rc(sub))):
        for f in range(3):
            p = tr(s[f:]); pos = 0
            for seg in p.split("*"):
                if len(seg) >= 40:
                    n0 = f + 3 * pos; n1 = n0 + 3 * len(seg) - 1
                    g0, g1 = (a + n0, a + n1) if strand == "+" else (b - n1, b - n0)
                    out.append((f"{strand}{f}:{g0}-{g1}", seg))
                pos += len(seg) + 1
    return out
def score(seg):
    q = [TextSequence(name=k.encode(), sequence=v).digitize(alpha) for k, v in seg]
    res = {k: {} for k, _ in seg}
    for n, h in H.items():
        for hits in pyhmmer.hmmsearch([h], q, E=1e9, domE=1e9, F1=1.0, F2=1.0, F3=1.0, bias_filter=False):
            for x in hits: res[str(x.name)][n] = x.score
    return res
D = "../2026-10-30_barnettozyma/data/"
for fa, off, a, b in ((D + "GCA_023627795.1_JAMAIE010000005.1_589225-655853.fa", 589225, 618500, 623500),):
    seq = next(iter(fasta(fa).values())); s = segs(seq, a, b, off); r = score(s)
    top = sorted(r.items(), key=lambda kv: -max(kv[1].values() or [0]))[:6]
    for k, v in top: print(k, len(dict(s)[k]), {n: round(x, 1) for n, x in v.items()})
print("--- alignment span of the best alpha2 segment, and the same search in the other californica windows")
seqs_ = dict(s)
best = max(r.items(), key=lambda kv: kv[1].get("alpha2", -9))[0]
q = [TextSequence(name=best.encode(), sequence=seqs_[best]).digitize(alpha)]
for hits in pyhmmer.hmmsearch([H["alpha2"]], q, E=1e9, domE=1e9, F1=1.0, F2=1.0, F3=1.0, bias_filter=False):
    for x in hits:
        for d in x.domains:
            a = d.alignment
            print(best, "hmm", a.hmm_from, a.hmm_to, "of", H["alpha2"].M, "target", a.target_from, a.target_to, "of", len(seqs_[best]), round(d.score, 1))
            print(seqs_[best][:80], "... starts with", seqs_[best][:12])
for fa in sorted(__import__("glob").glob(D + "GCA_*_*.fa")):
    seq = next(iter(fasta(fa).values())); sg = segs(seq, 1, len(seq), 1)
    # keep segments >= 90 aa to limit work
    sg = [(k, v) for k, v in sg if len(v) >= 90]
    rr = score(sg)
    topk = sorted(rr.items(), key=lambda kv: -kv[1].get("alpha2", -9))[:2]
    print(fa.split("/")[-1][:30], [(k, len(dict(sg)[k]), round(v.get("alpha2", 0), 1), round(v.get("a1", 0), 1)) for k, v in topk])
