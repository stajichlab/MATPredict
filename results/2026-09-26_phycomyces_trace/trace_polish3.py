"""Log creation of every GeneCluster on contig_552 (id + creating line), and
the cluster id each classify() call is keyed under."""
import sys, inspect
import MATPredict.detect.pipeline as P
Orig = P.GeneCluster
class Logged(Orig):
    def __init__(self, *a, **k):
        super().__init__(*a, **k)
        if self.contig == "contig_552":
            fr = inspect.currentframe().f_back
            print("TRACE new", id(self), self.start, self.end,
                  f"{fr.f_code.co_filename.split('/')[-1]}:{fr.f_lineno} {fr.f_code.co_name}", file=sys.stderr)
P.GeneCluster = Logged
import MATPredict.detect.scoring as S, MATPredict.detect.flank_carried as F
for mod in (S, F):
    if hasattr(mod, "GeneCluster"):
        mod.GeneCluster = Logged
ocl = P.classify
def cl(*a, **k):
    r = ocl(*a, **k)
    f = inspect.currentframe().f_back.f_locals
    c = f.get("cluster")
    print("TRACE classify", id(c), c.contig, f.get("gene_name"), r.status, file=sys.stderr)
    return r
P.classify = cl
from MATPredict.__main__ import main
sys.exit(main(sys.argv[1:]))
