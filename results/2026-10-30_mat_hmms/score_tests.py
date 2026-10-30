"""Score the test proteins with the new a1 / alpha2 / a2 / alpha1 profiles and the Saccharomycotina MATa2 profile, against a background
of the S. schoenii proteome (TransDecoder ORFs, a non-MAT-biased Saccharomycotina set)."""
import sys, pyhmmer
from pyhmmer.easel import Alphabet, TextSequence
alpha = Alphabet.amino()
HM = {n: f"hmm/{n}.hmm" for n in ("a1", "alpha2", "a2", "alpha1", "a2_no_Pichiaceae", "a2_no_Saccharomycetaceae", "a2_no_Serinales")}
HM["sacch_a2(existing)"] = "../../db/Ascomycota/classifiers/MATa2_saccharomycotina/MATa2_saccharomycotina.hmm"
import os
hmms = {}
for n, p in HM.items():
    if os.path.exists(p):
        with pyhmmer.plan7.HMMFile(p) as hf: hmms[n] = hf.read()
def fa(path):
    d = {}; n = None
    for l in open(path):
        if l.startswith(">"): n = l[1:].split()[0]; d[n] = ""
        elif n: d[n] += l.strip().rstrip("*")
    return d
test = fa("data/test_proteins.faa")
bg = {k: v for k, v in fa("/Users/jstajich/projects/MATPredict/results/2026-10-29_ss_records/data/transcripts.fa.transdecoder.pep").items() if len(v) >= 60}
print("background proteins", len(bg), file=sys.stderr)
def block(d): return [TextSequence(name=k.encode(), sequence=v).digitize(alpha) for k, v in d.items()]
def run(h, d):
    out = {}
    for hits in pyhmmer.hmmsearch([h], block(d), E=1e9, domE=1e9, incE=1e9, incdomE=1e9, Z=3067, F1=1.0, F2=1.0, F3=1.0, bias_filter=False):
        for x in hits: out[str(x.name)] = max(out.get(str(x.name), -9), x.score)
    return out
BG = {n: run(h, bg) for n, h in hmms.items()}
T = {n: run(h, test) for n, h in hmms.items()}
print({n: len(BG[n]) for n in hmms}, file=sys.stderr)
cols = list(hmms)
with open("data/test_scores.tsv", "w") as fh:
    fh.write("protein\tlen\t" + "\t".join(cols) + "\n")
    for k, v in test.items():
        fh.write(f"{k}\t{len(v)}\t" + "\t".join(f"{T[n].get(k, 0):.1f}" for n in cols) + "\n")
    fh.write("BACKGROUND_max\t\t" + "\t".join(f"{max(BG[n].values()):.1f}" for n in cols) + "\n")
    fh.write("BACKGROUND_p99\t\t" + "\t".join(f"{sorted(BG[n].values())[int(0.99*len(BG[n]))]:.1f}" for n in cols) + "\n")
    # rank of each test protein among the background for the a2 profiles
    for n in cols:
        pass
rank = {}
with open("data/test_rank.tsv", "w") as fh:
    fh.write("protein\t" + "\t".join(f"{n}:bits(bg proteins >=)" for n in cols) + "\n")
    for k in test:
        fh.write(k + "\t" + "\t".join(f"{T[n].get(k,0):.0f}({sum(1 for s in BG[n].values() if s >= T[n].get(k,0))})" for n in cols) + "\n")
print(open("data/test_scores.tsv").read()); print(open("data/test_rank.tsv").read())
