# `nearest_reference` on the v0.9.0 campaign: real values (#200)

Date 2026-11-21. Job 32359463; `label_loci_table.py` at `main` f901c929 on `results/2026-11-09_v090_ascomycota/loci.tsv` (20,844 loci). Counts in `results/2026-11-21_nearest_reference_200/summary.tsv`; the full table stays on the HPCC. This is the step PR #311 left open ("tables not regenerated, curator has not seen real values").

## Tests from the issue, on real data
- Every locus has a value: 0 blank, 1 explicit `none` (an undetermined-allele locus), 20,844 of 20,844 labelled.
- Same input, same output was tested in #311 on fixtures; not re-run here.
- 63 `NoRef` loci: `Peziz__NoRef_MAT1-1` 31, `_MAT1-2` 31, `_undet` 1. The v0.7.0 count was also 63 (31 / 32); I compared the counts only, not the loci one by one.

## What the column says
- **NoRef loci** sit nearest the frozen references XanPar (16 loci have it first in MAT1-1-1 or MAT1-2-1; I did not total the other groups), NeuCra, AspFum, BipMay and BotCin. Examples: `MAT1-2-1=NeuCra` 9, `MAT1-2-1=XanPar` 9, `MAT1-1-1=XanPar` 7, `MAT1-2-1=AspFum` 5, `MAT1-2-1=BipMay` 4. Multi-gene loci list each core gene, for example `MAT1-1-1=BipMay|MAT1-1-4=AspFum`, and genes of one locus can point to different references.
- **Undetermined-allele loci (976):** mostly NeuCra (460) and BipMay (133) by MAT1-1-1, then the yeast records.
- **Commonest values overall:** `MAT1-2-1=NeuCra` 1,829; the three-gene NeuCra locus 1,485; the S. cerevisiae MAT locus 1,101 and 987; the A. fumigatus locus 916; the C. auris MTL-A locus 663.

## Question for the curator
Yeast loci (Saccharomycotina) show the matched record id (for example `MATA2=4932_s288c_MATsc_HMRa`, `MTLA2=4959_cbs767_MTL_A`), not a system code. That is the fallback the definition gives for a record whose class has no registry code, because `pezizomycotina_class_systems` only covers Pezizomycotina classes. Should yeast records map to their system code (for example MATsc to SacCer, MTL to CanAlb), so one column reads the same across the tree? Without the change, the column is unambiguous but mixes codes and record ids.

## Caveat
Record classes are read from the current `db/`, which has records added after the campaign's frozen copy (MTLA2, the two Exophiala records). I did not check whether these records change any locus's matched reference, and I did not diff the column against the frozen db.
