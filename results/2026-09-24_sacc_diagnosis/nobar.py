# Run `detect` with the >=2 modelled-gene bar disabled (min_polished_genes=0), no source edits.
import sys, functools
import MATPredict.detect.cli as cli
import MATPredict.detect.pipeline as P
orig=P.run_pipeline
cli.run_pipeline=functools.partial(orig, min_polished_genes=0)
from MATPredict.__main__ import main
sys.exit(main(sys.argv[1:]))
