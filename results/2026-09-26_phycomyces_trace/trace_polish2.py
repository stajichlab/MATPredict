"""Log, per cluster, the admitted families and the genes sent to polishing.

Wraps pipeline globals that the polish loop calls by module name:
_families_meeting_evidence_floor, _padded_window, classify.
"""
import sys, inspect
import MATPredict.detect.pipeline as P
ofl, opw, ocl = P._families_meeting_evidence_floor, P._padded_window, P.classify
def fl(cluster, families, floor):
    r = ofl(cluster, families, floor)
    if cluster.contig in ("contig_552", "scaffold_145"):
        genes = sorted({(str(h.family_key), h.gene_name) for h in cluster.hits})
        print("TRACE floor", cluster.contig, cluster.start, cluster.end, "floor=", floor,
              "admitted=", [str(f.key) for f in r], "hits=", genes, file=sys.stderr)
    return r
def pw(cluster, padding, lengths):
    r = opw(cluster, padding, lengths)
    g = inspect.currentframe().f_back.f_locals.get("gene_name")
    print("TRACE window", cluster.contig, g, r, file=sys.stderr)
    return r
def cl(*a, **k):
    r = ocl(*a, **k)
    f = inspect.currentframe().f_back.f_locals
    print("TRACE classify", f.get("cluster").contig, f.get("gene_name"), r.status, file=sys.stderr)
    return r
P._families_meeting_evidence_floor, P._padded_window, P.classify = fl, pw, cl
from MATPredict.__main__ import main
sys.exit(main(sys.argv[1:]))
