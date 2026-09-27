# Replay: ignore allele-absent genes when assigning confidence

2026-09-27, read-only replay of existing reports. Scripts: `replay.py`
(reimplements `tiering.assign_tier` and the post-tier caps in
`pipeline._build` / `flank_carried`), `summarise.py`. Tables: `calls_*.tsv`,
`risers.tsv`, `summary.txt`, `closeness.txt`. Each run is replayed with the
code and roster of the tree that produced it (frozen worktrees; 882aa01,
f7b9773, 634dda4 re-extracted with `git archive`).

## Reproduction

100% of calls reproduced in 7 of 8 runs; Mucoromycota 1a00b0a 241/245
(98.4%). The 4 misses are excluded from all counts. Getting there needed the
existing `model_idiomorph_alternatives` rule: for Mucoromycota the current
code ALREADY ignores the losing (other-allele) half of sexM/sexP in the
unpolished check. The proposed rule generalises that.

## Variants

- **A** (as proposed): ignore UNPOLISHED genes whose `present_in_idiomorphs`
  excludes the called idiomorph when deciding `any_gene_unpolished`.
- **B**: A + drop allele-absent genes (any status) from the found set used for
  the expected core (`expected_genes_for_idiomorph`).
- **B'**: A + drop only allele-absent genes below 50% identity from that set.

Why B matters: the 18 Wallemia v2 calls have NO unpolished gene. They are
medium because weak v1-gene cross-hits (HMG, STE3) widen the expected core to
the full roster, and SXI1 is then "missing" (`core_not_found`). Rule A cannot
reach them.

## Calls rising medium -> high

| Run | Family | medium | A | B | B' |
|---|---|---:|---:|---:|---:|
| basidio_ad1f865 | Basidiomycota:MAT | 209 | 3 | 3 | 3 |
| cap6_f7b9773 | Ascomycota:MAT | 156 | 14 | 25 | 24 |
| cap6_f7b9773 | Ascomycota:MATsc | 177 | 2 | 78 | 7 |
| cap6_f7b9773 | Ascomycota:MTL | 25 | 4 | 6 | 4 |
| dothideo_d18ec5a | Ascomycota:MAT | 43 | 4 | 12 | 12 |
| dothideo_52cd292 | Ascomycota:MAT | 46 | 9 | 17 | 17 |
| early_634dda4 | Mucoromycota:MAT | 117 | 0 | 1 | 1 |
| mucoro_1a00b0a | Mucoromycota:MAT | 85 | 0 | 29 | 0 |
| serinales_882aa01 | Ascomycota:MTL | 625 | 16 | 94 | 82 |
| wallemia_3d8a755 | Basidiomycota:wallMAT | 18 | 0 | 18 | 17 |
| **total** | | | **52** | **283** | **167** |

No other tier changes under A. B lowers 29 MATsc and 1 PM call (the found set
changes the roster), B' none recorded as lowered in counted runs.
All other Basidiomycota families (HD 648 medium, Bbeta, Balpha, PR, aLocus):
0 risers under every variant.

## Strong ignored genes (>60% identity): 108 risers, all under B only

- MATsc S. cerevisiae (71): MATA2 at ~100% ignored in MATalpha calls — the
  a2/alpha2 X-region homology shared by all three cassettes.
- Mucoromycota Rhizopus (29): btbA at 66-99% ignored in Minus calls. btbA is
  rostered `present_in_idiomorphs: [Plus]` but is clearly present in Minus
  loci — a roster error, not a reason to promote.
- Serinales MTL (7): MTLA1 at 60-100% in alpha calls (likely collapsed a/alpha,
  cf. the C. albicans read study).
- Wallemia (1): HMG 60.5%.
B would promote all of these. B' promotes none of them by construction.

## Spot checks (~25 risers)

- A, Ascomycota MAT: called MAT1-2-1 modelled at 40-57%, ignored MAT1-1-x
  cross-hits 24-33% unpolished, SLA2/APN2/COX13 flanks present. Plausible.
- A, Basidiomycota MAT (Cryptococcus): SXI1 83%, STE3 94% modelled; ignored
  SXI2 26-32% unpolished. Plausible.
- A, MTL: mixed. Some alpha calls have no flanks and all core genes 22-50%;
  one ignores MTLA2 52% + MTLA1 50% against a called MTLalpha1 46%. Not safe.
- B', Wallemia v2: STE3v2 76-95% modelled, BAP31/CAF1 flanks, ignored HMG/STE3
  33-44%. Plausible.
- B', MATsc (7): every core gene 29-42%, no flanks. Not plausible as high.
- B', Serinales MTL: alpha genes 46-61% vs ignored a-genes 30-48%. Separation
  often thin.

## Closeness guard

Many risers ignore a cross-hit almost as strong as the called gene
(`closeness.txt`): gap between the called allele's best MODELLED core identity
and the strongest ignored gene < 10 points in 29/52 A risers and 78/167 B'
risers; in some the ignored gene is STRONGER than the called one. Requiring a
gap >= 10 points leaves A 23 risers, B' 89 (Serinales 39, Dothideomycetes 20,
cap-panel MAT 8, Wallemia 17, Cryptococcus 3, MTL 2).

## Recommendation

Adopt B' with the closeness guard, not A alone and never B: ignore an
allele-absent gene (for the unpolished check AND the expected-core set) only
when it is below 50% identity AND at least 10 points below the called allele's
best modelled core gene. That reaches the Wallemia v2 calls (17/18), keeps the
S. cerevisiae cassette, Rhizopus btbA and collapsed a/alpha cases out, and
blocks promotions where the "absent" allele is as strong as the called one.
Separately: btbA's `present_in_idiomorphs: [Plus]` is contradicted by 29
strong Minus-locus hits and should be reviewed.
