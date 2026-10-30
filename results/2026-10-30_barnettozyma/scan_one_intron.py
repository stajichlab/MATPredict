"""Exhaustive one-intron (and zero-intron) gene-model search in a window region, scored with the curated Saccharomycotina
MATa2 profile and with phmmer against the curated MATa2 proteins. Usage: scan_one_intron.py window.fa strand lo hi [minaa] [model: HMM path or pfam:NAME]
Introns: GT..AG, 40-250 bp. Only ORFs that contain the splice junction and have >= minaa aa are scored. Prints the top 12."""
import sys, pyhmmer
def nm(x): return x.decode() if isinstance(x, bytes) else x
from Bio import SeqIO
from Bio.Seq import Seq
from pyhmmer.easel import Alphabet, TextSequence
fa, strand, lo, hi = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
minaa = int(sys.argv[5]) if len(sys.argv) > 5 else 70
MODEL = sys.argv[6] if len(sys.argv) > 6 else "/Users/jstajich/projects/MATPredict/db/Ascomycota/classifiers/MATa2_saccharomycotina/MATa2_saccharomycotina.hmm"   # path, or pfam:NAME
S = str(next(SeqIO.parse(fa, "fasta")).seq).upper()
R = S[lo - 1:hi]
T = R if strand == "+" else str(Seq(R).reverse_complement())
def gpos(i):  # index on T -> window coordinate
    return lo + i if strand == "+" else hi - i
alpha = Alphabet.amino()
def load_model():
    if MODEL.startswith("pfam:"):
        with pyhmmer.plan7.HMMFile("/Users/jstajich/db/funannotate/Pfam-A.hmm") as hf:
            for h in hf:
                if nm(h.name) == MODEL[5:]: return h
        raise SystemExit("no such Pfam profile")
    with pyhmmer.plan7.HMMFile(MODEL) as hf: return hf.read()
if MODEL.startswith("refs:"):
    import re as _re
    _, rfa, rx = MODEL.split(":", 2)
    from pyhmmer.easel import SequenceFile
    with SequenceFile(rfa, digital=True, alphabet=alpha) as sf: _db = [x for x in sf.read_block() if _re.search(rx, nm(x.name))]
    from pyhmmer.easel import DigitalSequenceBlock
    _db = DigitalSequenceBlock(alpha, _db)
    hmm = None
else:
    hmm = load_model()
cands = []
CURJ = [None]
def orfs(spliced, junction):
    for fr in range(3):
        aa = str(Seq(spliced[fr:len(spliced) - ((len(spliced) - fr) % 3)]).translate())
        pos = 0
        for part in aa.split("*"):
            a0 = fr + pos * 3; a1 = a0 + len(part) * 3
            if len(part) >= minaa and (junction is None or a0 < junction < a1):
                cands.append((spliced, part, a0, junction, CURJ[0]))
            pos += len(part) + 1
orfs(T, None)
n = len(T)
gts = [i for i in range(n - 1) if T[i:i + 2] == "GT"]
ags = [j for j in range(1, n) if T[j - 1:j + 1] == "AG"]
import bisect
for i in gts:
    for j in ags[bisect.bisect_left(ags, i + 40):bisect.bisect_right(ags, i + 250)]:
        CURJ[0] = j
        orfs(T[:i] + T[j + 1:], i)
        # record intron with the candidate
        cands[-1:] = cands[-1:] if cands and cands[-1][3] == i else cands[-1:]
# scoring: dedupe by protein
seen = {}
for sp, part, a0, junc, jend in cands:
    seen.setdefault(part, (sp, a0, junc, jend))
seqs = [TextSequence(name=str(k).encode(), sequence=p).digitize(alpha) for k, p in enumerate(seen)]
res = []
keys = list(seen)
if hmm is None:
    for hits in pyhmmer.phmmer(seqs, _db, E=1e-1):
        if len(hits): res.append((hits[0].score, int(nm(hits.query.name))))
else:
  for hits in pyhmmer.hmmsearch(hmm, seqs, E=10):
    for h in hits:
        k = int(h.name if isinstance(h.name, str) else h.name.decode()); res.append((h.score, k))
res.sort(reverse=True)
for sc, k in res[:12]:
    p = keys[k]; sp, a0, junc, jend = seen[p]
    print(f"{sc:.1f} bits len={len(p)} junction_idx_on_T={junc} cds_start_idx={a0}", ("intron window %d-%d len %d seq %s...%s" % (min(gpos(junc), gpos(jend)), max(gpos(junc), gpos(jend)), jend - junc + 1, T[junc:junc+8], T[jend-5:jend+1])) if junc is not None else "no intron")
    print("   ", p)
