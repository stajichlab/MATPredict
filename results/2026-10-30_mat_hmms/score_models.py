"""Score S. schoenii candidate protein models (stdin-free; edit MODELS) against the a1, alpha2, a2 profiles."""
import pyhmmer
from pyhmmer.easel import Alphabet, TextSequence
seq = "".join(l.strip() for l in open("../2026-10-30_saccharomycopsis/chrIV_463k-542k.fa") if not l.startswith(">")).upper()
W = 463000
from itertools import product
bases = "TCAG"; aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
codon = {a + b + c: aa[i] for i, (a, b, c) in enumerate(product(bases, repeat=3))}
def sub(a, b): return seq[a - W:b - W + 1]
def rc(s): return s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))
def tr(s): return "".join(codon.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))
def prot(exons, strand):
    s = "".join(sub(a, b) for a, b in sorted(exons)); return tr(s if strand == "+" else rc(s))
MODELS = {
 "left_a1_1exon_ATG471961-472464": ([(471961, 472464)], "+"),
 "left_a1_p5_2exon_103aa": ([(471961, 472227), (472445, 472489)], "+"),
 "left_a1_agent_471889-472461": ([(471889, 472461)], "+"),
 "right_a1_536074-536634_minus": ([(536074, 536634)], "-"),
 "right_a1_3exon_211aa": ([(536046, 536090), (536296, 536634), (536693, 536944)], "-"),
 "left_alpha2_470349-470951": ([(470349, 470951)], "+"),
 "right_a2_537367-538167": ([(537367, 538167)], "+"),
}
alpha = Alphabet.amino()
hmms = {n: next(pyhmmer.plan7.HMMFile(f"hmm/{n}.hmm")) for n in ("a1", "alpha2", "a2", "alpha1")}
for name, (ex, st) in MODELS.items():
    p = prot(ex, st); q = TextSequence(name=name.encode(), sequence=p.replace("*", "")).digitize(alpha)
    out = []
    for hn, h in hmms.items():
        hits = list(pyhmmer.hmmsearch(h, [q], E=10))[0]
        out.append(f"{hn} {hits[0].score:.0f}" if len(hits) else f"{hn} -")
    print(name, len(p.rstrip('*')), "aa", "stop" if p.endswith("*") else "nostop", "start", p[0], " | ", "  ".join(out))
