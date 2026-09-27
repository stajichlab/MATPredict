"""Step 5: build two HMG-box alignments and keep the better one.

A: hmmalign --trim to PF00505; keep match columns only (no insert columns).
B: MAFFT (L-INS-i if n < 1000, else --auto), then ClipKIT smart-gap.
Metrics per alignment: columns, parsimony-informative columns (>= 2 residues each
seen in >= 2 sequences, gaps and X ignored), gap fraction.
Also reports how much of the HMG box each MAT1-2-1 outgroup covers in A.
"""
import collections, subprocess, sys
def read_afa(f):
    d, n = {}, None
    for l in open(f):
        l = l.rstrip()
        if l.startswith(">"): n = l[1:].split()[0]; d[n] = []
        elif n: d[n].append(l)
    return {k: "".join(v) for k, v in d.items()}
def write(d, f):
    with open(f, "w") as fo:
        for k, v in d.items(): fo.write(f">{k}\n{v}\n")
def metrics(d):
    seqs = list(d.values()); L = len(seqs[0]); pi = 0; gaps = 0
    for i in range(L):
        col = [s[i] for s in seqs]
        gaps += sum(c in "-." for c in col)
        cnt = collections.Counter(c for c in col if c not in "-.X")
        if sum(1 for v in cnt.values() if v >= 2) >= 2: pi += 1
    return L, pi, round(gaps / (L * len(seqs)), 3)
n = sum(1 for l in open("domains.faa") if l.startswith(">"))
subprocess.run("hmmalign --trim --outformat afa PF00505.hmm domains.faa > aln_hmm_raw.afa", shell=True, check=True)
raw = read_afa("aln_hmm_raw.afa")
L = len(next(iter(raw.values())))
keep = [i for i in range(L) if all(not (s[i].islower() or s[i] == ".") for s in raw.values())]
A = {k: "".join(v[i] for i in keep).upper() for k, v in raw.items()}
write(A, "aln_hmm.afa")
mode = "--localpair --maxiterate 1000" if n < 1000 else "--auto"
subprocess.run(f"mafft {mode} --thread 4 --quiet domains.faa > aln_mafft_raw.afa", shell=True, check=True)
subprocess.run("clipkit aln_mafft_raw.afa -m smart-gap -o aln_mafft.clipkit.afa > /dev/null", shell=True, check=True)
B = {k: v.upper() for k, v in read_afa("aln_mafft.clipkit.afa").items()}
mA, mB = metrics(A), metrics(B)
with open("alignment_comparison.txt", "w") as fo:
    fo.write(f"sequences {n}\n")
    fo.write(f"A hmmalign PF00505 match cols: columns {mA[0]} parsimony-informative {mA[1]} gap_fraction {mA[2]}\n")
    fo.write(f"B mafft {mode} + clipkit smart-gap: columns {mB[0]} parsimony-informative {mB[1]} gap_fraction {mB[2]}\n")
    for k, v in A.items():
        if k.startswith("OUT|"):
            fo.write(f"MAT1-2-1 coverage of HMG box (A): {k} {sum(c not in '-.' for c in v)}/{len(v)}\n")
print(open("alignment_comparison.txt").read())
