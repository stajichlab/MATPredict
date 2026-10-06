# Lichtheimiaceae MAT exploration
Status: open

## Question
Can MATPredict call Plus/Minus in Lichtheimiaceae (wide sense, Walther 2019)?
Do hypotheses H1 (rnhA counts as flank support) and H2 (admit a core-only
cluster typed at or above the gate) help, and at what cost?

## Data and code version
- Classifier: sexM/sexP/P1 HMMs from `.claude/worktrees/run-52b3ff9/db/Mucoromycota/classifiers/MAT/`; gate 98.4, margin 25.
- Reports replayed: LCG `2026-10-01_heldout_rerun/lcg_runs` (52b3ff9); BFD `2026-09-30_cap_v3/regression/cand/mucoromycota/runs` (V3); Jena `2026-10-01_heldout_rerun/jena_scaffolds`.
- Inputs: `genomes.tsv` (121 family genomes); `all.tsv` (765 genomes: LCG 621, BFD 80, Jena 64). SLURM job 29324717.
- Check protein: CDS03202.1, *L. ramosa* SexM (Schulz 2016).

## Method
1. tblastn of the 18 sexM/sexP references; merge hits within 3 kb; exonerate protein2genome model per locus; score with the classifier (`worker.py`). Score all annotated proteins too (LCG, Jena).
2. miniprot of tptA, rnhA, glrA, algA, btbA references.
3. Anchor A = best tblastn hit of CDS03202.1. Anchor B = best typed candidate. Neighbours within ±100 kb from gff3; DIAMOND all-vs-all (≥40% id, ≥50% coverage both), single-linkage clusters (`synteny.py`, `synteny_clusters.tsv`).
4. Replay H1 and H2 on the reports (`analyze.py`).

## Results
Per genus (`by_genus.tsv`, `per_genome.tsv`):

| Genus | n | Reported now | H2 gains | H1 gains | Core locus ≤50 kb of rnhA/tptA |
|---|---|---|---|---|---|
| Lichtheimia | 51 | 0 | 50 | 0 | 0 |
| Circinella | 30 | M16 P1 | 2 | 8 | 28 |
| Thamnostylum | 8 | M5 P1 | 0 | 2 | 7 |
| Fennellomyces | 8 | M5 | 0 | 0 | 8 |
| Rhizomucor | 14 | P4 | 10 | 0 | 2 |
| Dichotomocladium | 5 | 0 | 4 | 0 | 0 |

**CDS03202.1 check.**
- The shipped classifier scores CDS03202.1 as sexM 87.0, sexP 72.1, P1 91.5. That is below the gate and untyped.
- Anchor A is a full-length ortholog in 33/33 annotated *Lichtheimia* genomes (blastp to CDS03202.1 is 62–96% identity; to Mucorales sexM, 35–40%).
- At A, the annotated proteins score sexM 50–98. None is typed.
- Anchor B (the typed candidate) is never on A's contig (0/80 genomes with both anchors).
- The H2 models have a mean identity to CDS03202.1 of 36.7% (Minus) and 27.5% (Plus). So the H2 gains in *Lichtheimia* are other HMG genes.

**Synteny.**
- The A gene is in one cluster across all 7 annotated genera (73 genomes). An Hsp90 gene and two hypothetical genes are conserved within 50 kb of it in all genera (*Lichtheimia* 29–31/33, *Circinella* 20–24/28, *Rhizomucor* 7/9).
- No flank gene is on A's contig in any genome (`flank_positions.tsv`).
- In *Circinella* the A ortholog is present in genomes with Minus-typed and Plus-typed rnhA loci. So it is not idiomorph-specific there. This is inferred for *Lichtheimia* too, but not shown.
- In the Circinella group, the sexM/sexP locus sits 2–9 kb from rnhA. tptA sits with algA elsewhere (algA–tptA on the same contig in 21/30 *Circinella*).

**Labels (from filenames).**
- In the Circinella group, 6 of 8 labelled genomes that have a call or an H1 call are inverted: 4 Plus-labelled genomes are called Minus, and 2 Minus-labelled genomes get H1 Plus. 2 agree.

**H1.**
- 11 loci; 9 genomes gain a first call. All gains are Plus with sexP 70–88 (below the gate).
- 2/2 labelled gains contradict the label. 1 gain adds a second idiomorph.
- The 1 gain outside the family, *Rhizopus microsporus* NRRL A-17693, has scores identical to *C. umbellata*. This suggests a misidentified genome (unverified).

**H2.**
- 112 loci. 43 are outside the family; of these, 18 add a second call, and 15 of those are the opposite idiomorph.
- Zygo 23: 0 gains for H1 and for H2 (all 23 analysed).

## Limits
- Exonerate models cover about 60–80 aa (the HMG region). Single-linkage clustering can merge paralogs.
- There are no curator mating-type labels for *Lichtheimia* or *Rhizomucor*. The +/− convention of the RSA strains is not checked against sexP/sexM.

## Curator decisions
Open:
- H1 and H2 are not recommended.
- How to treat CDS03202.1 (the published "SexM").
- The label inversion in the Circinella group.

## Files
`genomes.tsv`, `all.tsv`, `per_genome.tsv`, `candidates.tsv`, `h1_h2_gains.tsv`, `by_genus.tsv`, `summary.txt`, `anchors.tsv`, `flank_positions.tsv`, `synteny_clusters.tsv`, `h2_vs_cds03202.tsv`, `worker.py`, `analyze.py`, `synteny.py`.
