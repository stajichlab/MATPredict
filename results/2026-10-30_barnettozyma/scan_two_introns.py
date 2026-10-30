"""Exhaustive search for a gene with up to two introns in a window region, for the divergent homeodomain MATalpha2 gene.
Donor GTA[ACT][GT] (yeast consensus GTATGT, relaxed), acceptor [CTA]AG, intron 35-200 bp. ORF must start with M, end at a stop,
be >= minaa, and contain every junction. Scored by a Pfam profile (default Homeodomain) and by phmmer against the curated alpha2 proteins.
Usage: scan_two_introns.py window.fa strand lo hi [minaa] [pfam name]"""
import sys, re, bisect, pyhmmer
from Bio import SeqIO
from Bio.Seq import Seq
from pyhmmer.easel import Alphabet, TextSequence, SequenceFile
def nm(x): return x.decode() if isinstance(x, bytes) else x
fa, strand, lo, hi = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
minaa = int(sys.argv[5]) if len(sys.argv) > 5 else 120
pf = sys.argv[6] if len(sys.argv) > 6 else "Homeodomain"
S = str(next(SeqIO.parse(fa, "fasta")).seq).upper()
R = S[lo - 1:hi]; T = R if strand == "+" else str(Seq(R).reverse_complement())
gpos = (lambda i: lo + i) if strand == "+" else (lambda i: hi - i)
don = [m.start() for m in re.finditer(r"(?=GTA[ACT][GT])", T)]
acc = [m.end() for m in re.finditer(r"(?=[CTA]AG)", T)]           # index just after AG = m.start()+3 ; use end of lookahead offset
acc = [m.start() + 3 for m in re.finditer(r"(?=[CTA]AG)", T)]     # exclusive end of intron
introns = [(d, a) for d in don for a in acc[bisect.bisect_left(acc, d + 35):bisect.bisect_right(acc, d + 200)]]
alpha = Alphabet.amino()
def translate_orfs(chain, juncs):
    out = []
    for fr in range(3):
        L = len(chain) - ((len(chain) - fr) % 3)
        aa = str(Seq(chain[fr:L]).translate())
        pos = 0
        for part in aa.split("*"):
            k = part.find("M")
            if k >= 0 and len(part) - k >= minaa:
                a0 = fr + (pos + k) * 3; a1 = fr + (pos + len(part)) * 3
                if all(a0 < j < a1 for j in juncs): out.append(part[k:])
            pos += len(part) + 1
    return out
cands = {}
def add(prot, spec): cands.setdefault(prot, spec)
for (d, a) in introns:                                                    # one intron
    ch = T[:d] + T[a:]
    for p in translate_orfs(ch, [d]): add(p, ((d, a),))
# two introns: restrict to pairs of introns that do not overlap, spacing >= 30
import itertools
sel = introns
for i, (d1, a1) in enumerate(sel):
    for (d2, a2) in sel:
        if d2 >= a1 + 30 and d2 - a1 < 800:
            ch = T[:d1] + T[a1:d2] + T[a2:]
            for p in translate_orfs(ch, [d1, d1 + (d2 - a1)]): add(p, ((d1, a1), (d2, a2)))
keys = list(cands)
seqs = [TextSequence(name=str(k).encode(), sequence=p).digitize(alpha) for k, p in enumerate(keys)]
with pyhmmer.plan7.HMMFile("/Users/jstajich/db/funannotate/Pfam-A.hmm") as hf:
    hmm = next(h for h in hf if nm(h.name) == pf)
res = []
for hits in pyhmmer.hmmsearch(hmm, seqs, E=10):
    for h in hits: res.append((h.score, int(nm(h.name))))
res.sort(reverse=True)
for sc, k in res[:int(sys.argv[7]) if len(sys.argv) > 7 else 10]:
    sp = cands[keys[k]]
    print(f"{sc:.1f} {pf} bits, {len(keys[k])} aa, introns (window): " + "; ".join(f"{min(gpos(d),gpos(a-1))}-{max(gpos(d),gpos(a-1))} ({a-d} bp {T[d:d+6]}..{T[a-3:a]})" for d, a in sp))
    print("   ", keys[k])
