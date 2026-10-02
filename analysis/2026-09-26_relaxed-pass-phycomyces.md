# Phycomyces NRRL_1554 regression and the relaxed pass
Status: decided

## Question
Why was a correct Plus/high call (2026-09-21) lost by 2026-09-23?

## Results (`results/2026-09-26_phycomyces_trace/NOTE.md`)
- Since 6bc985d, `_relaxed_results` never set `polished_genes`, so every
  relaxed-pass call was withheld. On contig input the locus is split, so only
  the relaxed pass can call it.
- Fix 047b5f2; Zygo 23 back to 23/23 on contigs.
- Later (Fable review F1, fix 4f37cb7): the relaxed gate ran before the
  withhold rules; now it runs on strict results that survive them.

## Files
`results/2026-09-26_phycomyces_trace/`, `scripts/zygo_regression.py`.
