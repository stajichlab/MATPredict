# Mucoromycota re-score after the 2026-09-26 idiomorph fixes

Code `4a24ffb` (frozen `run-4a24ffb`) vs baseline `f7796ae` (same code before
items 1-3), both on the 293 Mucoromycota genomes of the early-diverging scan,
SLURM short, 293/293 finished each. Changes in `4a24ffb`:
btbA non-informative (08a2616); model both sexM and sexP (6c8050e, 4a24ffb);
relaxed pass counts modelled genes (047b5f2). Full tables:
`compare_f7796ae_vs_4a24ffb.txt` (same-code baseline) and
`compare_634dda4_vs_4a24ffb.txt` (the original scan).

## Label changes vs f7796ae

| change | calls | where |
|---|---:|---|
| unchanged | 185 | |
| Plus -> Minus | 34 | all Rhizopodaceae; includes all 33 sexM-only Rhizopus calls |
| Minus -> Plus | 13 | Backusella 6, Apophysomyces 4 (+1 Mucoraceae), Blakeslea 1, Syncephalastrum 1 |
| new | 13 | 8 relaxed-pass, 5 strict (below) |
| lost | 0 | |

* **33/33 Rhizopus sexM-only calls are now Minus** (btbA no longer votes).
* **11 of the 13 sexP-clade Minus calls are now Plus**, each decided on the two
  exonerate models: sexP scored 264-449 vs sexM 110-147. This agrees with the
  HMG-box tree (sexP clade, UFBoot 99). The 2 Umbelopsis calls stay Minus on
  the first-pass verdict (no model pair to compare).
* **Lichtheimiaceae 1 -> 3 called** (2 Plus, 1 Minus; the 2 new are Rhizomucor
  pusillus relaxed calls). Syncephalastraceae 4/9 unchanged in count, labels now
  2 Plus / 2 Minus.

## New strict calls (5)

* Mucor irregularis B50 and B7584: Plus, 5 modelled genes. Before, withheld as
  flank-carried (core outside flank span). These are the two real loci the
  Ascomycota flank-rule audit said were wrongly withheld.
* Benjaminiella poitrasii (Minus), Cokeromyces recurvatus (Plus), Radiomyces
  spectabilis (Plus): before, withheld with 1 modelled gene; now the second
  core model counts.

## Relaxed pass (the medium-cap exploration)

91 relaxed candidates in 17 genomes (the pass runs only where the strict pass
found nothing). 8 reached 2 modelled genes and were reported, all at medium.
The strict tier rule would give 10 of the 91 `high` and 81 `medium`; of the 8
reported, 1 (GCA_019828915.1, tptA+sexM) would be high. Zygo 23 contig input:
NRRL_1554 comes back Plus through this pass (uncapped: high).

Decision on the cap: **kept at medium.** Uncapping changes 1 of 8 reported
calls here. Every relaxed call is by definition below the fraction floor and
most are 2-gene clusters on short contigs; the data give no reason to let
them reach high.

## Zygo 23 (`zygo23_4a24ffb/score.txt`)

scaffold+proteome 23/23 locus, 23/23 idiomorph; contigs genome-only 23/23,
23/23 (was 22/23 since 6bc985d).

## Limits

No truth exists for most of these genomes; label changes are scored against
the HMG-box tree and the resolution scores, not against known mating types.
The Rhizopus Minus labels are consistent with sexM at ~96% identity and no
sexP, not independently confirmed.
