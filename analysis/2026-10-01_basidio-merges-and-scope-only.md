# Basidiomycota branch merges and scope-only families
Status: decided

## Question
Can the signed-off records on curation-puccinio (Pucciniales, Sporidiobolales,
Wallemia, Rhodotorula; Agaricomycete HD anchors) and basidio-anchors
(Polyporales/Russulales/Russulaceae PR) join PR #9 without harming existing
calls?

## Data and code version
- Baseline: PR #9 c26669c (frozen worktree run-c26669c)
- Step 1: + curation-puccinio, cf57bad (run-cand-puccinio)
- Step 2: + basidio-anchors, 7a9c987 (run-cand-basidio)
- Step 3: + scope-only families and Boletales in PR scope, 5885a7e (cand-fix)
- Inputs: testset/regression_panel.tsv (166 genomes) and Zygo 23 on both inputs

## Method
1. scripts/run_regression_panel.sh for each step against the one before.
2. scripts/regression_check.py diff of baseline vs step 3 (net).
3. Every call-touching change was listed with species, order and routing mode.

## Results
Step 1 (results/2026-10-01_regression_puccinio/diff/): Ascomycota,
Mucoromycota, Zygo 23 unchanged. Basidiomycota: the new redPR and wallMAT
families made calls in seven phylum-fallback genomes outside their lineage
(Rhizoctonia, Auricularia, Fomitiporia, Dacryopinax, Exobasidium,
Microbotryum, Leucosporidium), redPR mostly "A1/A2 high". Cause: fallback
searches every family, and the STE3 references match any basidiomycete
receptor. Serpula lacrymans moved from fallback to lineage routing (new
Boletales HD record) and lost 3 PR calls; PR was not in the Boletales scope.

Net, baseline vs step 3 (results/2026-10-01_regression_net_final/):

| panel | call-touching changes |
|---|---|
| Ascomycota (30), Mucoromycota (80), Zygo 23 (46) | 0 |
| record sources (23) | 5 model/span, no call or label change |
| Basidiomycota (33) | 14 gained, 2 lost, 0 label changes, 19 model/span |

- Gained, all in-lineage: redPR + redHD (R. toruloides), redHD
  (S. johnsonii), wallMAT (W. ichthyophaga), rustHD (M. americana,
  P. graminis), HD x7 from the anchors (Ganoderma 2, Gelatoporia, Trametes,
  Serpula, Stereum), PR (Laccaria).
- Lost: Serpula 1 of 3 and Stereum 1 PR call, both CAAX-only "unverified",
  now withheld below the fraction floor.
- Three PR calls (Ganoderma, Trametes, Heterobasidion) lost the "unverified"
  label: a curated pheromone gene now supports them.
- The seven wrong-lineage calls are gone.

## What changed in detection
- 5885a7e: `fallback_searchable: false` (roster flag) on redPR, redHD, rustHD,
  wallMAT: never searched in phylum fallback or explicit-phylum runs;
  `--exhaustive` still searches them. PR scope adds Boletales (68889).

## Limits
- The panel has 33 Basidiomycota genomes; fallback orders are sparsely
  covered. The full 3,270-genome run has not been repeated.
- The two lost PR calls and the Serpula calls rest on the CAAX scan, whose
  false-positive rate is not bounded.

## Curator decisions
- 2026-10-01: scope-only families (option a); add Boletales to PR scope;
  signed off the net result; merged into PR #9 at 5885a7e.
- 2026-10-01: no new cap-off test; the 2026-09-27 test (2/50 gained) stands.

## Files
results/2026-10-01_regression_{puccinio,basidio_anchors,scope_only,net_both,net_final}/
