Stage 0b rescue test (design note section 3.4b / 9). Run 2026-10-21, HPCC job 30234460, exonerate 2.4.0.
truth.py -> truth.tsv (34 curated MAT1-1-1 / MAT1-2-1 genes with deposited exons); fetch.py fetched +-3 kb windows from NCBI (the windows FASTA, ~0.5 MB, is not committed: rerun fetch.py);
build_refs.py -> refs.faa (55 curated proteins); rescue_test.py (withhold record or genus, truncate references to 30/50/70/100%, production exonerate settings, extend, compare with the deposited end); analyze*.py -> tables.
s0b_rescue.tsv: all 272 runs. truncated_rows2.tsv: the 189 truncated alignments that reached a stop, with the acceptance rule applied.
