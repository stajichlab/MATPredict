import os, pickle, runpy, sys
from MATPredict.detect import pipeline
orig = pipeline.build_receptor_arrays
def wrap(fams, hits):
    pickle.dump([(h.family_key.locus_name,h.gene_name,h.method,h.contig,h.start,h.end,h.strand,h.identity,h.coverage,h.evalue,h.bitscore,h.align_length_aa,h.reference_length_aa,h.reference_record_id,h.superseded_by) for h in hits], open(os.environ["DUMP"], "wb"))
    return orig(fams, hits)
pipeline.build_receptor_arrays = wrap
sys.argv = ["MATPredict", "detect"] + sys.argv[1:]
runpy.run_module("MATPredict", run_name="__main__")
