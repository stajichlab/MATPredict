"""Six-frame stop-to-stop fragments (>= MINLEN aa) of a genomic window, scanned with Pfam-A (MAT-relevant domains) and
phmmer against the MAT reference proteins. Usage: scan_window.py FASTA CONTIG START END OUT.tsv [MINLEN]"""
import sys, subprocess, pyhmmer
from pyhmmer.easel import SequenceFile, Alphabet, TextSequence
fa, contig, a, b, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
MINLEN = int(sys.argv[6]) if len(sys.argv) > 6 else 40
seq = "".join(subprocess.check_output(["samtools", "faidx", fa, f"{contig}:{a}-{b}"], text=True).split("\n")[1:])
CODON = {}
bases = "TCAG"; aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, x in enumerate(bases):
    for j, y in enumerate(bases):
        for k, z in enumerate(bases): CODON[x + y + z] = aas[16 * i + 4 * j + k]
def rc(s): return s[::-1].translate(str.maketrans("ACGTacgtN", "TGCAtgcaN"))
def tr(s): return "".join(CODON.get(s[i:i + 3].upper(), "X") for i in range(0, len(s) - 2, 3))
frags = []
for strand, s in (("+", seq), ("-", rc(seq))):
    for f in range(3):
        p = tr(s[f:]); pos = 0
        for piece in p.split("*"):
            if len(piece) >= MINLEN:
                st = f + 3 * pos; en = st + 3 * len(piece)  # nt offsets in s
                if strand == "+": g0, g1 = a + st, a + en - 1
                else: g0, g1 = b - en + 1, b - st
                frags.append((strand, f, g0, g1, piece))
            pos += len(piece) + 1
def _s(x): return x.decode() if isinstance(x, bytes) else x
alpha = Alphabet.amino()
ts = [TextSequence(name=f"{s}{f}_{g0}-{g1}".encode(), sequence=p).digitize(alpha) for s, f, g0, g1, p in frags]
dom = {}
want = {"HMG_box", "HMG_box_2", "Homeodomain", "Homeobox_KN", "MATalpha_HMGbox", "I_LWEQ", "ANTH", "ENTH"}
with pyhmmer.plan7.HMMFile("/Users/jstajich/db/funannotate/Pfam-A.hmm") as hf:
    for hits in pyhmmer.hmmscan(ts, hf, E=1e-2):
        q = _s(hits.query.name)
        dom[q] = [(_s(h.name), round(h.score, 1)) for h in hits if _s(h.name) in want]
with SequenceFile(sys.argv[0].rsplit("/", 1)[0] + "/../2026-10-29_ss_records/refs.faa" if False else "/Users/jstajich/projects/MATPredict/results/2026-10-30_petasospora/data/all_refs.faa", digital=True, alphabet=alpha) as sf: db = sf.read_block()
best = {}
for hits in pyhmmer.phmmer(ts, db, E=1e-3):
    if len(hits):
        h = hits[0]; best[_s(hits.query.name)] = (_s(h.name), h.score)
with open(out, "w") as o:
    o.write("frag\tstrand\tstart\tend\taa\tpfam\tbest_ref\tbits\tseq\n")
    for (s_, f, g0, g1, p), t in zip(frags, ts):
        n = _s(t.name); d = dom.get(n, []); br = best.get(n)
        if d or br:
            nm = br[0] if br else ""
            nm = (nm.split("name=")[1].split("|")[0] if "name=" in nm else nm)[-45:]
            o.write(f"{n}\t{s_}\t{g0}\t{g1}\t{len(p)}\t{d}\t{nm}\t{br[1]:.0f}\t{p}\n" if br else f"{n}\t{s_}\t{g0}\t{g1}\t{len(p)}\t{d}\t\t\t{p}\n")
