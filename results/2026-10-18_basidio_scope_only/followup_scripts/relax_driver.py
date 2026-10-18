"""Analysis only: `matpredict detect` with the modelled-gene bar changed through run_pipeline's own parameter
(min_polished_genes), no change to src. usage: relax_driver.py BAR genome.fa taxid out_dir"""
import argparse, sys
from MATPredict.detect import cli

BAR = int(sys.argv[1])
REAL = cli.run_pipeline


def wrapped(*a, **k):
    k["min_polished_genes"] = BAR
    return REAL(*a, **k)


cli.run_pipeline = wrapped
p = argparse.ArgumentParser()
cli.register_subcommands(p.add_subparsers(dest="cmd"))
args = p.parse_args(["detect", "--genome", sys.argv[2], "--taxid", sys.argv[3], "--out-dir", sys.argv[4]])
cli._cmd_detect(args)
