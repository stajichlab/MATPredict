# Why Phycomyces blakesleeanus NRRL_1554 lost its Plus call

2026-09-26. Curator ruling: trace the miss.

## Answer

Two separate changes combine. Only the second one is a defect.

1. **The Zygo harness changed its input on 2026-09-23.** The 2026-09-21 run
   searched scaffolds converted from the annotation `.gbk`, with the annotated
   proteome (`--proteins`). From 2026-09-23 on, the harness searches the
   funannotate `.contigs.fsa`, genome-only. That file splits scaffold_145 at
   its N gaps (AGP: contig_551 | 10 N | contig_552 | 10 N | contig_553). The
   locus then sits on contig_552 (11,567 bp) with only sexP and rnhA; tptA is
   on contig_553. So the cluster scores 2/5 = 0.4, below the 0.5 floor, and
   the strict pass skips it. The relaxed pass is meant for exactly this case.

2. **Defect, commit `6bc985d` (2026-09-22, "a locus must rest on gene models,
   not on bare alignments").** It added the modelled-gene bar (at least 2
   modelled genes) and applies it after the relaxed pass. But
   `_relaxed_results` (`src/MATPredict/detect/pipeline.py`) builds its
   `DetectionResult` without setting `polished_genes`, so the field keeps its
   default of 0. **Every relaxed-pass call is therefore withheld as
   `modelled_gene_bar`, however well its genes are modelled.**

The prior agent's "272 bp from a contig end" is true on the contig input
(contig_552:272-10,473), not on the scaffold the 2026-09-21 run used.

## Evidence

Code x input matrix (`run_matrix.sh`, `summarize.py`, `runs/`):

| code | scaffolds + proteome | scaffolds | contigs.fsa |
|---|---|---|---|
| 382f4b3 (09-21) | Plus/high/mat_locus, 4 genes | Plus/high/mat_locus, 4 genes | Plus/medium/partial_locus (relaxed) |
| eab0a6b (parent of 6bc985d) | | | Plus/medium/partial_locus (relaxed) |
| 6bc985d | | | not detected: "0 modelled gene(s)" |
| 90db4af (PR #9 head) | Plus/high/mat_locus, 4 genes | Plus/high/mat_locus, 4 genes | withheld, `polished_genes: 0` |

- The current code is not worse on the scaffold input. Neither the database
  nor the search changed the result.
- The genes model correctly. Traced inside the 90db4af pipeline
  (`trace_polish2.py`): sexP `polished_agree`, rnhA `polished_disagree`, both
  in window contig_552:1-11,567. Called directly (`probe_polish.py`):
  sexP 9,505-10,476 at 100%, rnhA 2,695-5,289 at 100% (exonerate). The count
  is lost only because the relaxed result never reads `polish_by`.
- Zygo 23 on contigs.fsa (`run_zygo23.sh`, scored with
  `results/2026-09-26_sxi_slot_gateA/score_zygo.py`): eab0a6b 23/23 locus
  and idiomorph, 1 relaxed call; 6bc985d 22/23, 0 relaxed calls. NRRL_1554
  is the only Zygo genome that needed the relaxed pass.
- Current reports hold no relaxed call at all: 0 in 294 early-diverging
  Mucoromycota genomes and 0 in 2,368 Serinales genomes.

## Proposed fix

In `_relaxed_results`, set `polished_genes` from the same rule the strict
path uses (`_modelled_gene_names` in `run_pipeline`: distinct live gene names
with a polished outcome in `polish_by` for this (cluster, family), plus
`diamond_proteome` hits). Move that rule to a module-level helper so both
paths share it. Add a test: a relaxed cluster with 2 polished genes is
reported, and one with 0 is withheld. Then re-run Zygo 23 (expect 23/23 on
contigs.fsa) and check how many relaxed calls return in a Mucoromycota panel.

Separately, decide which input the Zygo regression should use. The truth
table is in scaffold coordinates, and the 2026-09-21 baseline used scaffolds
plus proteome. The contig input is a harder, different test.

## Who else is affected

Any genome where the strict pass finds nothing, and the best cluster is
split by an assembly break but has 2+ well-modelled genes. That is the
fragmented-assembly case the relaxed pass exists for. Since 6bc985d the
relaxed pass has emitted no reportable call anywhere. Not measured: how many
genomes in past panels would regain a call.

No code was changed. Frozen worktrees created: run-382f4b3, run-90db4af,
run-6bc985d, run-eab0a6b.
