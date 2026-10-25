# Intron cap and the Basidiomycota b-locus (2026-10-25)

Evidence for the per-family `max_intron_bp` in `db/Basidiomycota/order.yml` (bLocus, 10 kb). Frozen clean-run
code and database (`cleanrun-2026-10b`) unless stated; clean-run settings.

PR #92 capped polish alignments at a 5 kb intron. Between v0.6.0 and the clean run, 7 Basidiomycota genomes lost a
call to that merge (bisect: first bad merge 295df47). The 7 genomes were re-run at several caps.

| genome | what v0.6.0 called | returns at |
|---|---|---|
| GCA_001444555.1, GCA_022813745.1 | bE/bW, bW 3.6-3.8 kb, largest intron 2.7 kb | 8 kb |
| GCA_015039315.1 | bE/bW, bW 2.5 kb (110 aa of a 647 aa reference), largest intron 2.2 kb | 10 kb |
| GCA_049862275.1, GCA_002897395.1, GCA_001600235.1 | bE/bW on 26-78 kb models with 23-67 kb "introns" (chimeric) | 50-200 kb only |
| GCA_900009495.1 (pheromone receptor) | PR, 2.7 kb model at 20 kb, 10.9 kb model with an 8 kb intron at 50 kb | 20 kb |

So 10 kb for the b-locus restores the first three genomes, the three chimeric ones stay lost, and the
pheromone-receptor genome stays lost (PR #92 found long-intron pheromone-receptor models to be chimeric; PR and HD keep 5 kb).

Files:
- `out_cap_7genomes.txt`: the 7 genomes at 5, 20, 50 and 200 kb (models and largest introns).
- `out_cap_4genomes_mid.txt`: the four ordinary-model genomes at 5, 6, 8, 10, 15, 20 and 50 kb.
- `out_panels_family_cap_vs_clean_run.txt`: the real family cap (new main code) on 216 genomes (all Ustilaginales and
  Malasseziales plus the 7) and on a seeded random sample of 500 other Basidiomycota genomes, against the clean run:
  3 b-locus calls restored and no new model with an intron over 5 kb in the first; in the second 490 identical, 6 changed, 1 gained.
  (The script's own "intron > 5 kb" count ignores changed loci; read the model lists, where one bW model is 10.5 kb with a 7.7 kb intron.)
- `out_control_5kb_vs_family_cap.txt`: the 7 changed or gained genomes re-run with the cap forced back to 5 kb. With
  the new code at 5 kb their loci and model sizes equal the clean run, so every change is the cap's: one bW model
  770 bp -> 10,476 bp (7.7 kb intron, call unchanged), two small bW boundary shifts, one new low-confidence HD locus.
  The three `MAT`-family changes flagged in the panel are `MFa` status labels (`polished_disagree` -> `polished_agree`),
  not loci; confidence is unchanged.
