"""Six-frame stop-to-stop segments of a locus window, searched with the Pfam MAT-relevant profiles
(HMG_box, HMG_box_2, Homeodomain, Homeobox_KN, MATalpha_HMGbox, ANTH, ENTH, I_LWEQ) and with phmmer against the curated MAT proteins.
Usage: orf_scan.py window.fa refs.faa Pfam-A.hmm lo hi   (lo/hi: 1-based window coordinates to scan)
Prints a TSV: seg_id, strand, window_start, window_end, aa, pfam hits (name:bits:dom_start-dom_end in aa), best ref (name, bits)."""
import sys, pyhmmer
def nm(x): return x.decode() if isinstance(x, bytes) else x
from Bio import SeqIO
from Bio.Seq import Seq
from pyhmmer.easel import SequenceFile, Alphabet, TextSequence
fa, refs, pfam, lo, hi = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5])
MINAA = int(sys.argv[6]) if len(sys.argv) > 6 else 40
rec = next(SeqIO.parse(fa, "fasta")); S = str(rec.seq).upper()
alpha = Alphabet.amino()
segs = []   # (id, strand, wstart, wend, aa)
for strand in "+-":
    s = S if strand == "+" else str(Seq(S).reverse_complement())
    L = len(S)
    for fr in range(3):
        aa = str(Seq(s[fr:len(s) - ((len(s) - fr) % 3)]).translate())
        pos = 0
        for part in aa.split("*"):
            n = len(part)
            if n >= MINAA:
                a0 = fr + pos * 3; a1 = a0 + n * 3          # 0-based [a0, a1) on the oriented strand
                ws, we = (a0 + 1, a1) if strand == "+" else (L - a1 + 1, L - a0)
                if we >= lo and ws <= hi:
                    segs.append((f"s{len(segs)+1}", strand, ws, we, part))
            pos += n + 1
seqs = [TextSequence(name=i.encode(), sequence=p).digitize(alpha) for i, st, ws, we, p in segs]
want = {"HMG_box", "HMG_box_2", "Homeodomain", "Homeobox_KN", "MATalpha_HMGbox", "ANTH", "ENTH", "I_LWEQ", "HMG_box_5", "HMG_box_4"}
profs = []
with pyhmmer.plan7.HMMFile(pfam) as hf:
    for h in hf:
        if nm(h.name) in want: profs.append(h)
        if len(profs) == len(want): break
dom = {}
for hits in pyhmmer.hmmsearch(profs, seqs, E=1e-2, bit_cutoffs=None):
    pn = nm(hits.query.name)
    for h in hits:
        for d in h.domains.included if len(h.domains.included) else []:
            dom.setdefault(nm(h.name), []).append(f"{pn}:{d.score:.0f}:{d.alignment.target_from}-{d.alignment.target_to}")
        if not len(h.domains.included) and h.score > 12:
            dom.setdefault(nm(h.name), []).append(f"{pn}~:{h.score:.0f}")
with SequenceFile(refs, digital=True, alphabet=alpha) as sf: db = sf.read_block()
best = {}
for hits in pyhmmer.phmmer(seqs, db, E=1e-2):
    if len(hits):
        h = hits[0]; best[nm(hits.query.name)] = (nm(h.name), h.score)
print("seg\tstrand\twstart\twend\taa\tpfam\tbest_ref")
for i, st, ws, we, p in segs:
    br = best.get(i); d = dom.get(i, [])
    if br or d:
        nm = br[0].split("name=")[1].split("|")[0] if br and "name=" in br[0] else (br[0].split("|")[-1] if br else "")
        print(i, st, ws, we, len(p), ",".join(d), f"{nm}:{br[1]:.0f}" if br else "", sep="\t")
