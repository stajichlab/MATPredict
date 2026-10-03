# Second Wallemia record (v2) and the 51-genome re-run

2026-09-27, branch `curation-puccinio` (local), commits `3d8a755` (record +
roster), `4717ffd` (publication-highlights). Frozen run worktree `run-3d8a755`,
SLURM job 29125916. Build script: `build_wallemia_v2.py`; comparison:
`compare.py`, `per_genome.tsv`; run output `run_3d8a755/`.

## Genome choice

The 18 BFD Wallemiales genomes carrying the other version (`versions.tsv` of the
follow-up): W. mellicola 13, W. ichthyophaga 4, W. canadensis 1. Chosen:
*W. canadensis* EXF-10342 (GCA_056320075.1). It is the best-assembled (SPAdes,
159 contigs, contig N50 186 kb; the others are IDBA-Hybrid, 202-298 contigs, N50
112-199 kb) and the only one whose assembly annotates the locus receptor (no
W. mellicola or W. ichthyophaga other-version genome has a CDS at its STE3
hit). Its receptor is 80% (tblastn, 136 aa) to the W. mellicola version and
65% to the W. ichthyophaga version.

## Record `1708542_exf-10342_wallMAT_v2` (tier 2, putative, pending sign-off)

JBHFOR010000004.1:328,288-351,767. STE3v2 = KAN3001147.1, CAF1 = KAN3001148.1
(99.6% to v1), BAP31 = KAN3001158.1 (98.4% to v1): all from the assembly's
CDS features, 100% match on validation. HMG present: false (only a 43-aa
HMG-box fragment, 60.5% to the v1 HMG; Gostincar et al. 2019 report a truncated
HMG in the inverted version). SXI1 present: false (absent from that version per
Gostincar 2019; no homeodomain hit). The annotated STE3v2 (192 aa) is probably
3'-truncated (B5.5); the record keeps the annotated protein.

## Roster change (`wallMAT`)

- Vocabulary pattern -> enum `["v1", "v2"]`.
- Idiomorph labels: `v1` = W. mellicola CBS 633.66 version (the first record's
  label, unchanged), `v2` = this version. Placeholders; which is which mating
  type is not known.
- The receptor is split by version: `STE3` (present_in v1, the signed-off
  record's gene name kept unchanged) and `STE3v2` (present_in v2). SXI1 and HMG
  are v1-only. BAP31/CAF1 shared.
- max_cluster_gap_bp 20,000 -> 25,000 (v2 CAF1-BAP31 gap 20,489 bp).

## 51-genome re-run (vs the one-record run `run-a6770ca`)

| | before | after |
|---|---|---|
| v1-type genomes (33) | undetermined, high | **v1, high** (33/33; 0 changes of confidence or class) |
| other-version genomes (18) | undetermined, medium | **v2, medium** (18/18) |

- v2 vote: v2 178-234 vs v1 26-74 (best-bitscore). STE3v2 identity: W. mellicola
  93.0-94.6%, W. ichthyophaga 76.5%, W. canadensis 100%.
- No genome gets both v1 and v2.
- Per species: W. mellicola 14 v1 / 13 v2; W. ichthyophaga 18 / 4 (the 4 =
  EXF-759, EXF-3555, EXF-8622, EXF-8623, the published inverted strains);
  W. hederae 1 / 0; W. canadensis 0 / 1.

## Why v2 stays medium (0 of 18 reach high)

Not the record: `tiering.assign_tier` returns medium whenever
`any_gene_unpolished`. In every v2 genome the cluster also holds a weak,
unmodelled hit to the **v1-only** HMG gene (33-37%; the truncated HMG-box
fragment), so the rule fires even though the gene is not expected for the
called version. STE3v2, CAF1 and BAP31 are all modelled. A code change would
be needed: ignore unpolished genes that are not expected for the called
idiomorph (`expected_genes_for_idiomorph`) when setting the tier. Not changed
here (out of scope).
