# Circinella minor lost its Plus call under the MAT-gene gate
Status: open

## Question
Why did the MAT-gene gate withhold the Plus call in Circinella minor
GCA_016758965.1 on the curation-umbelopsis re-measure, and does the gate
also hit other Circinella genomes?

## Data and code version
- Gated run: code 667252a (curation-umbelopsis rebased on 076afe4, classifier
  rebuilt with the Umbelopsis and S. racemosum records), frozen worktree
  run-667252a. Report:
  results/2026-09-28_umbelopsis_rebased/Mucoromycota_667252a/runs/GCA_016758965.1_ASM1675896v1/
- Last run where it was called WITH the gate: a9fb69c
  (results/2026-09-28_next_fixes/Mucoromycota_a9fb69c/); before the gate:
  f25cf70 (results/2026-09-28_review_fixes/Mucoromycota_f25cf70/).
- Genome: BFD input_clean_genomes/GCA_016758965.1_ASM1675896v1.fa.gz
  (1,332 contigs). BFD family: Lichtheimiaceae.
- Held-out comparison: results/2026-09-28_lcg_holdout/ (code 076afe4).

## Method
1. Read the locus from the three reports and evidence_diagnostics.jsonl.
2. tblastn (-seg no, e<=1e-3) of the run's 51 reference proteins
   (_reference.faa) against the whole genome (work/hits.tsv).
3. Read the 29 LCG Circinella reports (runs2/).

## Results
Locus JAEPRB010000020.1:310,559-318,219 (contig length 325,645):
- Genes at the locus: sexP (single exon 310,595-311,089, 28.5% identity,
  modelled, polished_single) and rnhA (56.3%, 12 exons, modelled). A weak
  sexM hit (33.3%, 27% coverage) overlaps sexP. No other roster gene.
  fraction_found 0.4 (below the 0.5 floor); called by the relaxed pass.
- Classifier on the sexP model: a9fb69c Plus 112.5 vs Minus 30.1 (margin
  82.3); 667252a Plus 85.9 vs Minus 32.3 (margin 53.6).
- Gate: needs absolute >=100 bits OR >=2 roster flanks modelled at >=40%.
  a9fb69c passed on score (112.5). 667252a fails both: 85.9 bits and one
  flank (rnhA).
- Cause of the loss: the classifier REBUILD lowered this protein's absolute
  score by 26.6 bits (the sexP HMM now includes the new, divergent
  Umbelopsis and S. racemosum sexP). The locus and model did not change.
- Flanks elsewhere in the genome (best tblastn hits): tptA 57.7% and algA
  43.1% on JAEPRB010000215.1 (20.2-32.8 kb of a 67,394 bp contig, not at an
  end); glrA 54.4% on JAEPRB010000004.1:313.9 kb. None is within 300 kb
  upstream or 7.4 kb downstream (to the contig end) of sexP-rnhA. So in this
  assembly the sexP-rnhA pair sits apart from tptA/algA/glrA. This matches
  the S. racemosum NRRL 2496 pattern (sexP 69 bp from rnhA, other flanks
  elsewhere). The contig downstream of rnhA ends after 7.4 kb, so glrA could
  lie past the end; tptA upstream cannot be explained by a contig end.
- LCG held-out (076afe4, before the rebuild): 20 of 29 Circinella called.
  18 calls are sex+rnhA partial loci at medium; C. rigida NRRL 2341 and
  C. simplex NRRL 2407 are full loci (tptA/algA/glrA present) at high.
  Of the 9 uncalled, 7 have a sexP+rnhA locus typed Plus from an HSP
  fragment (margins 78.9-83.7, 112-116 bits) withheld by the modelled-gene
  bar and the fraction floor (sexP not modelled). None was withheld by the
  gate at 076afe4.
- Two LCG strains whose names carry a mating-type label are called the
  opposite way: C. angarensis RSA_198_Plus -> Minus (margin 68.8) and
  C. umbellata RSA_505_Plus -> Minus (margin 69.0). Not resolved: label
  convention, strain error, or detection error.

## What changed in detection
None (trace only).

## Limits
- One gated genome traced in depth. The LCG check used 076afe4, not the
  rebuilt classifier, so how many LCG Circinella calls the rebuild would
  push below 100 bits is not measured. Jena has no Circinella.
- Whether sexP-rnhA separated from tptA/algA is biology or an assembly
  artefact is inferred from one assembly (tptA not adjacent despite 310 kb
  of upstream sequence).

## Curator decisions
- Open: gate robustness (see recommendation) and the two name/label
  conflicts. Circinella is in BFD's Lichtheimiaceae, which the curator
  ruled (2026-09-28, option a) to be called only at >=100 bits for now;
  under that ruling the withholding is consistent.

## Recommendation (not implemented)
The 100-bit absolute threshold is not stable across classifier rebuilds
(this protein moved 26.6 bits with no locus change). In
src/MATPredict/detect/mat_gene_gate.py, derive the threshold from the
classifier build (record it in db/<Phylum>/classifiers/<family>/manifest.yaml
from the leave-one-genus-out full-protein score distribution, written by
scripts/build_idiomorph_hmms.py) instead of a fixed roster number, and
re-validate on the F3 set after each rebuild. Separately, the sexP-rnhA-only
arrangement in Circinella/Syncephalastrum argues for counting a strongly
modelled rnhA (>=50%) adjacent to the core as sufficient flank support in
families where tptA/algA are known to be displaced; this belongs with the
Lichtheimiaceae exploration.

## Files
- work/hits.tsv (tblastn of references vs genome)
- LCG table: results/2026-09-28_lcg_holdout/per_genome.tsv
