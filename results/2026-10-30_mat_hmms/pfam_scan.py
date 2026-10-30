"""Pfam HMG-box and homeodomain models against the test proteins (all hits, E <= 10, with the Pfam GA bit cutoff)."""
import pyhmmer
from pyhmmer.easel import SequenceFile, Alphabet
PF = "/Users/jstajich/db/funannotate/Pfam-A.hmm"
WANT = {"HMG_box", "HMG_box_2", "HMG_box_5", "HMGL-like", "Homeodomain", "Homeobox_KN", "MATalpha_HMGbox", "SoxE", "Sox17_18_mid"}
alpha = Alphabet.amino()
with SequenceFile("data/test_proteins.faa", digital=True, alphabet=alpha) as sf: q = sf.read_block()
models = []
with pyhmmer.plan7.HMMFile(PF) as hf:
    for h in hf:
        if str(h.name) in WANT: models.append(h)
print("models", [str(m.name) for m in models])
rows = {}
for hits in pyhmmer.hmmsearch(models, q, E=10):
    m = str(hits.query.name); ga = hits.query.cutoffs.gathering[0] if hits.query.cutoffs.gathering else None
    for h in hits:
        rows.setdefault(str(h.name), []).append(f"{m} {h.score:.1f}/{ga:.1f}" if ga else f"{m} {h.score:.1f}")
with open("data/pfam_hits.tsv", "w") as fh:
    for k in [str(s.name) for s in q]:
        fh.write(k + "\t" + "; ".join(rows.get(k, [])) + "\n")
print(open("data/pfam_hits.tsv").read())
