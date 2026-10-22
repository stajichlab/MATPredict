Stage 0c held-out set (design note section 13), built 2026-10-22 from NCBI protein records (GenPept text, not committed: rerun search.py + fetch_gp.py).
search.py -> ids; fetch_gp.py -> GenPept; parse_gp.py / filter2.py / choose_set.py -> candidates and the 106 selected genes (selected.tsv; 333 direct-deposit candidates in candidates_direct.tsv, of which 207 have >= 1 kb downstream and 101 are NOT used and are reserved for any future validation);
make_truth.py -> truth_heldout.tsv; windows fetched with ../stage0b_rescue/fetch.py (heldout_windows.fa not committed); rescue_test3.py (HPCC job 30258079) -> heldout_rescue.tsv; analyze_ho.py -> heldout_truncated_rows.tsv.
Truth = the depositors' annotated CDS structures (not RNA-supported; some may be predictions).
