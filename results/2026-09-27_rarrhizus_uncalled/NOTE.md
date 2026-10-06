# Uncalled Rhizopus arrhizus genomes (scan 4457c3a)

12 uncalled: 11 R. arrhizus + 1 R. delemar (GL13), all from the "GL" genome
series (GCA_0117640xx / 011801xxx / 011952xxx; list in uncalled.txt).
Report reason for each: best cluster below the 2-modelled-gene bar.

## Cause: the MAT locus is split across short contigs (assembly), not absent

Direct tblastn of the R. arrhizus records (64495_cbs110-17 Minus,
64495_cbs346-36 Plus) against each genome (best_hits.tsv):
- Every uncalled genome carries sexP at 98.4% (full length) and tptA,
  rnhA, btbA at 98.6-100%. No sexM hit. So all 12 are Plus by sequence.
- In all 12, sexP sits ALONE on its own contig, and almost always within
  ~200 bp of a contig end (e.g. JAANQS010002783.1:2,803-3,741 of 3,932).
  tptA+btbA are on a second contig and rnhA on a third.
- In called GL genomes sexP (or sexM) shares a contig with rnhA
  (e.g. GCA_011951295.2 JAANRC010000468.1), giving 2 modelled genes.
- Contiguity does not separate the groups (asm_stats.tsv): uncalled N50
  3.9-18.3 kb, called N50 6.0-17.1 kb. The break falls between sexP and
  rnhA in these particular assemblies.
- Not a taxon mislabel: tptA/rnhA are 100% identical to the R. arrhizus
  references.
- The "~40% identity" in the 2026-09-26 note came from the best withheld
  cluster being an HMG paralog cluster (e.g. JAANQS010000160.1, 40.5%);
  the true genes were single-gene clusters on separate contigs.

## Classification

All 12: assembly fragmentation at the locus (sexP isolated at a contig end).
Not a detection bug in the narrow sense: cross-contig merging is off by
design, and a single modelled core gene cannot pass the 2-gene bar.

## Implications

- The Rhizopus Plus/Minus tally is missing 12 Plus genomes.
- Same failure mode as the Cryptococcus split loci (SXI slot fix) and the
  C. auris gap cases.
