import sys
import pyhmmer
from pyhmmer.plan7 import HMMFile
from pyhmmer.easel import SequenceFile, Alphabet
hmm_path, seq_path, out = sys.argv[1], sys.argv[2], sys.argv[3]
with HMMFile(hmm_path) as hf:
    hmm = hf.read()
seqs = SequenceFile(seq_path, digital=True, alphabet=Alphabet.amino()).read_block()[:25]
rows = []
for top in pyhmmer.hmmsearch(hmm, seqs, cpus=1):
    for h in top:
        rows.append(f"{h.name}\t{h.score:.3f}")
open(out, 'w').write("\n".join(rows) + "\n")
print(len(rows), "hits written to", out)
