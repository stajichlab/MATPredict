# Polish cap: a core pair is not capped out (2026-10-17)

## Question
Seven Komagataella genomes with a real two-gene MTL cluster were silent on main. Is the per-family polish cap the cause, and can the cap be made independent of how many other clusters a genome admits?

## Step 1: cause confirmed
Seven genomes (GCA_023013455.1, GCA_023013825.1, GCA_023013855.1, GCA_023014715.1, GCA_023013895.1, GCA_023014845.1, GCA_030557195.1), own taxid, main (src at origin/main 5b9a9e0), capped vs `--max-polished-clusters-per-family 0`:

- Capped: no call in all 7. MTL admitted 18-24 clusters, 12-18 `polish_capped: true`; the true cluster (two core genes, only `core_MAT` role, best identity 43.4%, 7 hits) ranked 8th by the rank, so it was capped and withheld as `modelled_gene_bar`. The strong-identity tier (50%) did not catch it.
- Uncapped: `MTL A` appears in all 7 (6 high, 1 medium).
- Cost of uncapped: wall 104-160 s -> 178-297 s (about 1.8x).

## Step 2: design
`select_polish_clusters(..., core_genes)` / `has_core_pair`: a cluster with at least two distinct live (not superseded) REQUIRED `core_MAT` genes of one single idiomorph (families: MTL A, MATsc, MATyl, PM, mat1, wallMAT v1) is polished ON TOP of the existing selection. Property of the cluster alone, so the number of clusters or families does not matter. Default on; `--no-polish-core-pair` / `polish_core_pair=False` restores main.

Two earlier versions were measured and rejected:
1. Any two core genes of one idiomorph: counted the optional accessory genes of the Pezizomycotina `MAT` family (MAT1-1-1 + MAT1-1-3 ...). On 576 genomes it polished +30% clusters (3,289 -> 4,284) for 0 calls. Fixed by requiring non-optional genes.
2. Pair takes a slot from the weakest ranked cluster (cost neutral): on 1,505 genomes it displaced two real MTL alpha calls (GCA_030567475.1, GCA_047651925.1). Fixed by adding the pair on top of the cap, never in place of a cluster.

Flag: the design is not risky, so it is on by default, but the switch is kept because it makes the A/B and any rollback one flag.

## Step 4: campaign sample (final code, commit 03d9d61)
Sample: 1,505 Ascomycota genomes from the v0.6.0 campaign, seed 20261017: 110 each from the ten largest orders by non-fallback genomes (Hypocreales, Eurotiales, Serinales, Saccharomycetales, Pleosporales, Botryosphaeriales, Glomerellales, Magnaporthales, Onygenales, Sordariales; Sordariales has 61 (all it had after filtering), Botryosphaeriales 110) plus 405 phylum_fallback genomes (7 Komagataella included, 398 stratified over orders, at least 5 each). Own taxids from assembly_summary_GCA/GCF. Main vs fix, plus the regression panel (166 genomes) and Zygo 23.

| measure | main | fix |
|---|---|---|
| calls (1,505 genomes) | 1,601 | 1,660 |
| genomes with a call | 1,263 | 1,319 |
| calls gained / lost / label or confidence changed | | +59 / 0 / 0 |
| polished (cluster, family) pairs | 14,926 | 15,250 (+324, +2.2%) |
| paired wall time, 230 affected genomes (same task, back to back) | 19,889 s | 19,763 s (-0.6%, median delta 0 s) |
| regression panel (166 genomes) | 156 calls | 156 calls, 0 touch a call, 4 withheld-only changes |
| Zygo 23 (scaffold and contig) | 23/23 | 23/23, 0 changed loci |

Gained calls: 58 `MTL A` (Komagataella: 40 high, 18 medium, including CBS 7435 GCA_900235035.2 and other K. phaffii assemblies; all at 0-2.3 kb contig edges or one 8 kb span, `idiomorph_gene_only`) and 1 `MATsc MATalpha` medium (GCA_054953955.1). 59 of 1,505 genomes gain a call (one each). Gains appear only in MTL and MATsc (yeast families); in the first 1,051 genomes (the lineage-routed orders, interim run) no call changed. Which routing mode each of the 59 had was not tabulated (K. phaffii is phylum_fallback).

`core_model_changed` (2, GCA_030567475.1, MTL/MATsc alpha loci): coordinates and identity are identical in the two outputs; reported by the tool on this genome, not examined further.

## Campaign-wide importance
Small in count (59 of 1,505 genomes, 3.9%) but concentrated: it is the Komagataella (K. phaffii) problem plus about a dozen other Saccharomycotina. The v0.6.0 campaign had the same artefact, so a full Ascomycota rerun would add roughly the same fraction among fallback-routed yeast genomes. Cost is negligible (+2.2% polished pairs, wall time unchanged within noise). Proposal (not started): rerun only the 1,472 phylum_fallback genomes plus Saccharomycotina lineage-routed ones (Saccharomycetales, Serinales), about 5,000 genomes; no need for the whole 19,415.

## Caveats and unverified
- The gained Komagataella calls are `idiomorph_gene_only` at contig edges with no flank gene. The earlier agent's analysis (K. phaffii deposit proteins MATa1/MATa2 100% on that contig, placed at CBS 7435 2.7-5.0 kb) supports them as real; I did not re-verify, and the dipodascales note treated similar calls as suspect. Curator should confirm.
- GCA_054953955.1 `MATsc` gain and the 8 kb-span gains (GCA_023011165.1, GCA_023015205.1, GCA_030557195.1 span 1-8417) were not examined.
- Fix results for the 1,275 genomes where the selection did not change are taken from the intermediate run (same selection as main; polished counts equal), not rerun with the final commit. The final commit was run on the 230 genomes whose selection differed in the intermediate version, plus 3 panel genomes. The intermediate version differs from the final only in dropping displacement, so genomes with no added cluster are identical by construction.
- Unpaired wall times of the whole sample are confounded by cluster load (-7%); only the paired subset is a valid timing.
- Zygo 23 was run on the intermediate (slot-taking) version only; Mucoromycota has no family with two required core genes of an idiomorph, so the final version cannot differ; not re-run.
- Panel runs used own-taxid routing, not `--phylum Mucoromycota`.
- Composition with PR #97 (`--strong-single-gene-policy`): merged in a throwaway worktree; the only conflict is the `select_polish_clusters` call line, resolved by keeping both (`polish_allowed |= ssg_only` after my call); 167 polish/cap/single tests pass there. The merge itself is not committed here.
- The fix runs sit in `/bigdata/stajichlab/jstajich/projects/cap` (work_sample, work_panel, final_sample, final_panel, zygo_main, zygo_fix2). No cap jobs are queued.

## Tests
`tests/detect` and `tests/db` on the HPCC test environment, final code: 1,128 passed, 4 skipped. New: `tests/detect/test_polish_core_pair.py` (12 tests: capped-out pair now polished, unrelated clusters unchanged, many extra families do not change the result, optional genes do not make a pair, flag on/off through the pipeline).
