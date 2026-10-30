"""ORFs of the StringTie transcripts in a window of S. schoenii chromosome IV with genomic CDS coordinates, Pfam domains
(HMG_box, Homeodomain, MATalpha_HMGbox, ANTH/ENTH, I/LWEQ) and the best curated MAT reference hit.
Usage: chr4_orfs.py data refs.faa Pfam-A.hmm A B"""
import re, sys, pyhmmer
from pyhmmer.easel import SequenceFile, Alphabet, TextSequence
d, refs, pfam, A, B = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5])
tx = {}
for l in open(f"{d}/stringtie.gtf"):
    f = l.split("\t")
    if f[0] != "CM143129.1" or f[2] != "exon": continue
    tid = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
    tx.setdefault(tid, {"strand": f[6], "ex": []})["ex"].append((int(f[3]), int(f[4])))
tx = {k: v for k, v in tx.items() if A <= min(a for a, b in v["ex"]) and max(b for a, b in v["ex"]) <= B}
def to_genome(tid, s, e):  # 1-based transcript coords (transcript orientation = genome + order of exons, '+' transcripts; '-' reversed)
    ex = sorted(tx[tid]["ex"]); strand = tx[tid]["strand"]
    blocks, off = [], 0
    seq = ex if strand != "-" else ex[::-1]
    for a, b in seq:
        L = b - a + 1
        lo, hi = max(s, off + 1), min(e, off + L)
        if lo <= hi:
            blocks.append((a + (lo - off - 1), a + (hi - off - 1)) if strand != "-" else (b - (hi - off - 1), b - (lo - off - 1)))
        off += L
    return sorted(blocks)
orfs = []; name = None
for l in open(f"{d}/longest_orfs.pep"):
    if l.startswith(">"):
        m = re.match(r">(\S+) type:(\S+) len:(\d+) \S+ (\S+):(\d+)-(\d+)\(([+-])\)", l)
        name = m.group(1); tid = name.rsplit(".p", 1)[0]
        orfs.append([name, tid, m, ""])
    else: orfs[-1][3] += l.strip()
orfs = [o for o in orfs if o[1] in tx and len(o[3].rstrip("*")) >= 60]
alpha = Alphabet.amino()
seqs = {o[0]: TextSequence(name=o[0].encode(), sequence=o[3].rstrip("*")).digitize(alpha) for o in orfs}
with SequenceFile(refs, digital=True, alphabet=alpha) as sf: db = sf.read_block()
best = {}
for hits in pyhmmer.phmmer(list(seqs.values()), db, E=1e-2):
    if len(hits):
        h = hits[0]; best[hits.query.name.decode() if isinstance(hits.query.name, bytes) else hits.query.name] = (h.name.decode() if isinstance(h.name, bytes) else h.name, h.score)
dom = {}
want = {"HMG_box", "Homeodomain", "MATalpha_HMGbox", "ANTH", "ENTH", "I_LWEQ", "HMG_box_2", "Homeobox_KN"}
with pyhmmer.plan7.HMMFile(pfam) as hf:
    for hits in pyhmmer.hmmscan(list(seqs.values()), hf, E=1e-3):
        q = hits.query.name.decode() if isinstance(hits.query.name, bytes) else hits.query.name
        dom[q] = [((h.name.decode() if isinstance(h.name, bytes) else h.name), round(h.score, 1)) for h in hits if (h.name.decode() if isinstance(h.name, bytes) else h.name) in want]
print("orf\ttranscript\tstrand\tlen\tcds_blocks\tbest_ref\tbits\tpfam")
for name, tid, m, s in orfs:
    ts, te, st = min(int(m.group(5)), int(m.group(6))), max(int(m.group(5)), int(m.group(6))), m.group(7)
    blocks = to_genome(tid, ts, te) if True else []
    br = best.get(name)
    ref = ""
    if br:
        ref = br[0].split("|")[0][-40:] + ("|" + br[0].split("name=")[1].split("|")[0] if "name=" in br[0] else "|" + br[0].split("|")[-1])
    gs = {("+","+"):"+",("+","-"):"-",("-","+"):"-",("-","-"):"+"}.get((tx[tid]["strand"], st), "?")
    print(name, tid, gs, len(s.rstrip('*')), ";".join(f"{a}-{b}" for a, b in blocks), ref, f"{br[1]:.0f}" if br else "", dom.get(name, ""), sep="\t")
