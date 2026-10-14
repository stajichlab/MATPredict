import sys
sys.path.insert(0, "/psutil_stub")
from pyhmmer.plan7 import HMMFile
with HMMFile(sys.argv[1]) as hf:
    hmm = hf.read()
print("params", list(hmm.evalue_parameters.as_vector()), flush=True)
