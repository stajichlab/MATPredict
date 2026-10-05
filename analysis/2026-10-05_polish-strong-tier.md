# Polish cap identity tier (2026-10-05)

Status: implemented on branch `tiered-polish-rank` (PR open); full Dothideomycete re-run pending when this was written.

## Question
After the Dothideomycete records were merged, 15 genomes lost their call to the per-family polish cap, including the genome one record came from (`analysis/2026-10-05_dothideomycetes-full-run.md`). Can the pre-cap ranking be fixed without paying for
a higher cap?

## The change
`select_polish_clusters` now takes an identity tier. With the default cap of 6 per family, every admitted cluster whose best identity (any of the family's hits, superseded cross-hits included, as `_polish_rank` counts them) is at or above 50% is polished,
even past the cap, and the remaining slots up to the cap are filled from the rest in the existing rank (distinct genes, identity, hits, with the V3 strong-fragment preference kept). `run_pipeline(polish_strong_identity=50.0)`; CLI `--polish-strong-identity PCT`, default 50, 0 restores the plain cap.
The cap's "never more than N per family" no longer holds: a genome with more than N strong clusters polishes all of them (2 of 2,722 Dothideomycete genomes in the replay; at most 8 strong clusters).

## Evidence
- **Replay** over 2,722 Dothideomycete genomes and 2,582 true loci (`results/2026-10-05_dothideomycetes_full/cap_rank_replay_summary.txt`, `cap_tiered_replay_summary.txt`): the plain cap loses 17 true loci; this rule loses 0 at 16,284 clusters polished against 16,281 today (1.00x);
  a cap of 15 would need 2.47x the work, no cap 5.59x. The 21 true loci below 50% (mostly *Botryosphaeria dothidea*, 44.5%) are kept by the top-up.
- **Tests:** 14 new (`tests/detect/test_polish_strong_tier.py`): the IPO323 scenario, more strong clusters than the cap, no strong cluster (equals the plain cap), the inclusive 50.0 boundary, per-family independence, V3 priority in the top-up, the pipeline with the tier on and off, the CLI default and zero.
  Two existing cap tests now pin the plain cap (`polish_strong_identity=None`) because they test the cap itself; their fixture has two clusters above 50%. Full suite on HPCC: 1036 passed, 0 failed (the `test` environment needed pyhmmer, installed).
- **The 15 lost genomes, default flags** (`results/2026-10-05_tiered_lost15/`): all 15 recover, each identical to its v0.6.0 call (MAT1-1, medium, partial locus), including IPO323 at `NC_018206.1:616,667-623,138`. 12 minutes for 15 genomes on 8 CPUs, against 43 minutes with the cap off.
- **Regression gate** (`results/2026-10-05_regression_tiered_cap/`; baseline `af3c615`, candidate `0e39730`, same database, the fixed 163-genome panel and Zygo 23):

| panel | genomes | loci changed | touch a call |
|---|---|---|---|
| Ascomycota | 30 | 16 | 0 |
| Basidiomycota | 33 | 0 | 0 |
| Mucoromycota | 80 | 62 | 0 |
| record sources | 23 | 3 | 1 (call gained) |
| Zygo 23 | 46 | 2 | 0 |

  No call is lost and no label, confidence, locus class or gene model changes anywhere. The single call gained is `GCF_000011425.1`, the *A. nidulans* FGSC A4 genome, whose own MAT locus was withheld as `modelled_gene_bar+polish_capped` in the baseline and is now called (MAT1-1, medium): a curated record's own genome
  failing to call was the same cap problem. The other changes are loci withheld on both sides: whether a withheld cluster was capped or polished swaps (Ascomycota 14, Mucoromycota 41, record sources 2, Zygo 1), and clusters the tier now polishes that stay withheld below the fraction floor (Ascomycota 2, Mucoromycota 21 including 3 reversals, Zygo 1).

## Limits and decisions
- The 50% line is calibrated on Dothideomycetes. The regression panels contain few Dothideomycetes (one Ascomycota genome per order), so the full Dothideomycete re-run is the in-domain check; its comparison with the earlier run is added to this note when it finishes.
- Divergent clades with real loci under 50% are unaffected as long as fewer than 6 clusters clear 50%; where more do, the weak ones lose slots. The replay found 2 such genomes in 2,722.
- Rollback: `--polish-strong-identity 0`.
- Curator: confirm 50% as the default threshold, or ask for a per-family value in `order.yml`.
