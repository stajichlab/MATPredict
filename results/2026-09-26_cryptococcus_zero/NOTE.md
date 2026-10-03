# Why C. neoformans got 0 calls in the Basidiomycota pilot

2026-09-26. Code checked: `f81dad1` (frozen worktree `run-f81dad1`, the pilot's
baseline). No code changed.

## Short answer

This is not a regression, and not a routing, cap, flank-rule or anchor effect.

- The pilot held ONE C. neoformans genome: GCA_002221985.1 (A1-35-8).
- That genome is one of 4 C. neoformans genomes uncalled since the
  2026-09-23 Tremellales bar run (37 of 43 Cryptococcus genomes called there).
  The other three: GCA_002217545.1 (cng10), GCA_002220035.1 (MW_RSA852),
  GCA_002222025.1 (cng4).
- On 2026-09-21 (before the modelled-gene bar) the same 4 genomes got only
  scattered single-gene `idiomorph_gene_only` calls that named both a and
  alpha. They were never called correctly.

## Cause

All 4 genomes are MATalpha, and their MAT locus is split over 5-6 contigs.

Direct tblastn with the 13 JEC20/JEC21 reference proteins, A1-35-8:

| gene | contig | identity |
|---|---|---|
| SXI1 | APKI01000103.1 (93.5 kb) | 73.9% |
| FAO1 | APKI01000103.1 | 75.3% |
| STE3 (alpha) | APKI01000104.1 | 83.7% |
| MFalpha | APKI01000105.1 / 107.1 | 86.8% |
| RPL39 | APKI01000108.1 | 94.1% |
| PAN6 | APKI01000109.1 | 79.1% |

The other three genomes show the same split: SXI1 and FAO1 together on one
contig, STE3, PAN6, MFalpha and RPL39 each on other contigs. In the called
control cng8 (GCA_002222015.1), SXI1, FAO1, STE3 and MFalpha are on one
contig.

Why the real locus is lost (in `run_pipeline`, `src/MATPredict/detect/pipeline.py`):

1. `fraction_found` counts only the family's non-optional, searchable genes.
   In `db/Basidiomycota/order.yml` (MAT locus, scope 5234) SXI1 and SXI2 are
   `optional: true`; MFalpha, MFa and RPL39 are too short to search. So the
   denominator is 3: FAO1, PAN6, STE3.
2. The true-locus cluster on APKI01000103.1 has SXI1 + FAO1. SXI1 does not
   count, so it scores 1/3 = 0.33, below the 0.5 ambiguity floor. It is
   dropped before it is ever built (`score.fraction_found < ambiguity_floor`
   in the per-cluster loop). It never appears in `suppressed_loci` either.
   Its evidence row shows `admitted: true`, best identity 89.9%, 7 polish
   pairs attempted.
3. A paralog cluster on APKI01000022.1 (PAN6 + STE3 hits, best identity
   40.7%) scores 2/3 and passes the floor, but no gene models, so the
   modelled-gene bar withholds it. That is the only entry the report shows.
4. Cross-contig merging (`allow_cross_contig_fragments`) is off by default,
   by design (2026-09-2x ruling), so the split locus cannot be rejoined.

## Test of the fix (scratch DB copy, no code change)

With SXI1 counted as a required gene (DB copy only; SXI2 left optional,
which is valid for MATalpha genomes only):

| genome | before | after |
|---|---|---|
| A1-35-8 | none | alpha, medium, partial_locus, 2 modelled (SXI1, FAO1) |
| cng10 | none | same |
| MW_RSA852 | none | same |
| cng4 | none | same |
| cng8 (control) | alpha high, 4 modelled | alpha high mat_locus, 4 modelled (unchanged) |

## Recommended fix

Count an idiomorph-restricted HD pair as ONE required slot, satisfied by
either member: SXI1 (alpha) or SXI2 (a). This is the gene that defines the
idiomorph, so a cluster holding it should count it.

- Where: `score_cluster` in `src/MATPredict/detect/scoring.py` (the
  `optional` handling near lines 84-125), plus a roster field in
  `db/Basidiomycota/order.yml`, e.g. `idiomorph_alternatives: [SXI1, SXI2]`
  or a shared `slot: SXI` on both entries. Validate the schema in
  `db/_schema/order.schema.yaml`.
- Do not simply drop `optional: true` from both: that makes the denominator 5
  and the true cluster still scores 2/5 = 0.4.
- Check other families with idiomorph-restricted optional pairs (e.g.
  sexM/sexP, MAT1-1-1/MAT1-2-1) before applying the slot generically; the
  change raises `fraction_found` wherever such a gene is found.
- Re-run the 43 Cryptococcus genomes of the 2026-09-23 bar run and the
  23-genome Zygo regression after the change.

## Also seen (not investigated)

C. gattii LA55 (FAO1, MFalpha, SXI1; 0 modelled) and C. decagattii WM1802
(STE3, SXI2, RPL39; 0 modelled) are also uncalled on 2026-09-23. They may
share this cause; not checked.
