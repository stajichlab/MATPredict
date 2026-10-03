"""Run `matpredict detect` with the pipeline's polish and classify calls logged.

Usage: trace_polish.py <detect args...>   (logs to stderr lines starting TRACE)
"""
import sys
import MATPredict.detect.pipeline as P

def wrap(name):
    orig = getattr(P, name)
    def f(*a, **k):
        r = orig(*a, **k)
        if name.startswith("polish"):
            print("TRACE", name, k.get("gene_name"), k.get("window"),
                  None if r is None else (r.contig, r.start, r.end), file=sys.stderr)
        else:
            print("TRACE classify ->", getattr(r, "status", r), file=sys.stderr)
        return r
    setattr(P, name, f)

for n in ("polish_with_miniprot", "polish_with_exonerate", "classify"):
    wrap(n)
from MATPredict.__main__ import main
sys.argv = ["matpredict"] + sys.argv[1:]
sys.exit(main())
