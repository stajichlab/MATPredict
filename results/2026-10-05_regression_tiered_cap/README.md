Regression run-tiered (candidate, commit 0e39730 of branch `tiered-polish-rank`) vs run-af3c615 (baseline, `main` at the branch point). Same database; only `src/MATPredict/detect/pipeline.py` and `cli.py` differ.
Panel: `testset/regression_panel.tsv` (Mucoromycota 80, Basidiomycota 33, Ascomycota 30, record sources 23) and Zygo 23 on both inputs. Jobs: `jobs.tsv`. Summary: `diff/summary.md`; per panel `regression_summary.md` and `regression_withheld_changes.md`.
Result: no call lost and no label changed in any panel; 1 call gained (GCF_000011425.1, the *A. nidulans* FGSC A4 genome a curated record came from). The remaining changes are loci withheld on both sides, or absent on one side and withheld below the fraction floor on the other.
See `analysis/2026-10-05_polish-strong-tier.md`.
