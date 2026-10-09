#!/bin/bash
# Pull the HPCC scan output into scan_out/ (per-genome tsv) and pack it. Run from this directory.
mkdir -p scan_out
ssh hpcc2 'cd /bigdata/stajichlab/jstajich/agari_work/out && tar cf - *.loci.tsv *.region.tsv *.cand.tsv *.chance.tsv' 2>/dev/null | tar xf - -C scan_out
ls scan_out/*.chance.tsv | wc -l
