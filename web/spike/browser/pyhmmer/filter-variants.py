import sys
sys.path.insert(0, "/psutil_stub")
from pyhmmer.plan7 import HMMFile, Pipeline
from pyhmmer.easel import SequenceFile, Alphabet
alpha = Alphabet.amino()
with HMMFile("/sexM.hmm") as hf:
    hmm = hf.read()
seqs = SequenceFile("/train.faa", digital=True, alphabet=alpha).read_block()[:3]
variants = {"default": {}, "F1=1": {"F1": 1.0}, "F2=1": {"F2": 1.0}, "F3=1": {"F3": 1.0}, "bias off": {"bias_filter": False}}
for name, kw in variants.items():
    pl = Pipeline(alpha, **kw)
    hits = pl.search_hmm(hmm, seqs)
    print("variant", name, "hits", len(hits), flush=True)
