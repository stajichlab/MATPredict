"""Candidate proteins to score, and six-frame ORFs of the B. californica region 620,300-620,800 (unmodelled homeodomain)."""
import re, itertools
CODON = {a + b + c: aa for (a, b, c), aa in zip(itertools.product("TCAG", repeat=3), "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}
def rc(s): return s[::-1].translate(str.maketrans("ACGTacgt", "TGCAtgca"))
def tr(s): return "".join(CODON.get(s[i:i + 3].upper(), "X") for i in range(0, len(s) - 2, 3))
def fasta(path):
    d = {}; n = None
    for l in open(path):
        if l.startswith(">"): n = l[1:].split()[0]; d[n] = ""
        elif n: d[n] += l.strip()
    return d
out = []
cur = fasta("data/curated_mat_proteins.faa")
for k, v in cur.items():
    if "4909_cbs-573_MATpk_MATa|MATa2" in k: out.append(("Pkud_MATa2_control", v))
pf = fasta("../2026-10-30_petasospora/data/models_proteins.faa")
for k, v in pf.items():
    if k.startswith("MATa_JOY008|MATa") or k.startswith("MATalpha_TBRC4498|MATalpha1"): out.append(("PF_" + k.split("|")[0].split("_")[1] + "_" + k.split("|")[1], v))
bc = fasta("../2026-10-30_barnettozyma/models_proteins.faa")
for k, v in bc.items():
    if k in ("BC_MATa2", "BC_MATa1", "BC_MATalpha1"): out.append((k, v))
ss = fasta("/Users/jstajich/projects/MATPredict/results/2026-10-29_ss_records/data/longest_orfs.pep")
for k in ("SS.1982.2.p3", "SS.1982.1.p3", "SS.1982.2.p4", "SS.1982.2.p5", "SS.2002.2.p2", "SS.2001.1.p4", "SS.2002.1.p2"):
    out.append((k, ss[k].rstrip("*")))
# B. californica unmodelled region
g = fasta("../2026-10-30_barnettozyma/data/GCA_023627795.1_JAMAIE010000005.1_589225-655853.fa")
seq = next(iter(g.values())); OFF = 589225
for a, b in ((620000, 621100),):
    sub = seq[a - OFF:b - OFF + 1]
    for strand, s in (("+", sub), ("-", rc(sub))):
        for f in range(3):
            p = tr(s[f:])
            pos = 0
            for seg in p.split("*"):
                if len(seg) >= 50: out.append((f"BCreg_{strand}{f}_{pos}", seg))
                pos += len(seg) + 1
with open("data/test_proteins.faa", "w") as fh:
    for n, s in out: fh.write(f">{n}\n{s}\n")
print(len(out))
