"""Run `matpredict detect` with exonerate SKIPPED for the PAP1/OBP1/PIK1 flanks.

Experiment only (flank ablation, 2026-09-26): measures what polishing the
flanks with miniprot alone costs in calls, since exonerate on these long
proteins is ~90% of the new roster's runtime. No source file is changed; the
default polisher bound into run_pipeline's signature is swapped at runtime.
"""
import runpy, sys
from MATPredict.detect import pipeline

FLANKS = {"PAP1", "OBP1", "PIK1"}
f = pipeline.run_pipeline
names = f.__code__.co_varnames[: f.__code__.co_argcount]
i = names.index("polish_with_exonerate") - (len(names) - len(f.__defaults__))
defaults = list(f.__defaults__)
original = defaults[i]

def exonerate_except_flanks(*args, gene_name=None, **kwargs):
    if gene_name in FLANKS:
        return None          # classify() then uses miniprot's model alone
    return original(*args, gene_name=gene_name, **kwargs)

defaults[i] = exonerate_except_flanks
f.__defaults__ = tuple(defaults)
sys.argv = ["MATPredict"] + sys.argv[1:]
runpy.run_module("MATPredict", run_name="__main__")
