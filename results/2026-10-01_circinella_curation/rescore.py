"""Score the Circinella-group HMG proteins (label-tree rnhA_models.faa) with
the classifier before (polish-scope-cuts a8863f3) and after the rebuild.
Usage: python rescore.py FAA OLD_DIR NEW_DIR > rescore.tsv"""
import sys
from pathlib import Path
import pyhmmer
from Bio import SeqIO
faa, old, new = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
al = pyhmmer.easel.Alphabet.amino()
def load(d):
    return {n: pyhmmer.plan7.HMMFile(d / f"{n}.hmm").read() for n in ("sexM", "sexP")}
H = {"old": load(old), "new": load(new)}
seqs = [pyhmmer.easel.TextSequence(name=r.id.encode(), sequence=str(r.seq).rstrip("*")).digitize(al)
        for r in SeqIO.parse(faa, "fasta")]
sc = {}
for tag, hm in H.items():
    for g, h in hm.items():
        for hits in pyhmmer.hmmsearch([h], seqs, E=1e9, domE=1e9):
            for x in hits:
                sc[(tag, g, x.name if isinstance(x.name, str) else x.name.decode())] = x.score
print("id\told_sexM\told_sexP\tnew_sexM\tnew_sexP")
for s in seqs:
    n = s.name if isinstance(s.name, str) else s.name.decode()
    print(n, *(round(sc.get((t, g, n), 0.0), 1) for t in ("old", "new") for g in ("sexM", "sexP")), sep="\t")
