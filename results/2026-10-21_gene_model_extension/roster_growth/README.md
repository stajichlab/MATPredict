Roster-growth experiment (design note section 14b), 2026-10-22, HPCC array 30261437.
Augmented references = the 55 curated proteins + the 106 held-out deposited proteins (refs_aug.faa is rebuilt by build_aug.py from GenPept; not committed). Leave-genus-out: for each of the 140 genes (34 curated + 106 deposits) every reference of the gene's own genus is withheld.
rescue_test4.py (radius = genus only), analyze_rg.py compares with the baseline (55 references: s0b_rescue2 genus rows for the 34 genes, heldout_rescue for the 106). augmented_rows.tsv / baseline_rows.tsv: the compared truncated rows.
The 106 genes are development data; this experiment changes only the reference set, but is not a validation.
