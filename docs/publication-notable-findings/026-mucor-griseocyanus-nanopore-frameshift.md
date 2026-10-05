# 026. A Nanopore homopolymer frameshift truncates the Mucor griseocyanus sexM and costs its Minus label

- **Category:** assembly-or-annotation artefact; method
- **Status:** artefact-explained (no short reads to confirm the error)
- **Lineage:** Mucorales, Mucoraceae

## Summary
The Mucor griseocyanus CBS 116.08 assembly (Oxford Nanopore MinION only) has a
2-nt frameshift inside its sexM, in an A5 homopolymer. The gene model stops being
translated there (about 112 aa), so the classifier scores a fragment: the call
was Minus by a thin margin (31.1 bits) and becomes undetermined (16.3) after the
2026-10-04 rebuild. The frame-corrected protein (253 aa) scores Minus by +94 to
+116 bits. The locus is otherwise intact.

## Evidence
- Assembly GCA_060309175.1 (CBS116.08v1, Westerdijk Institute; contig-level,
  contig N50 1.73 Mb, 35x); reads: SRR32276708, OXFORD_NANOPORE MinION
  (BioSample SAMN46730503).
- Locus JBLNHG010000007.1:607,370-664,666 (tptA, sexM, rnhA, algA, glrA);
  sexM 659,923-660,683 (-).
- exonerate protein2genome of the curated M. circinelloides R7B sexM
  (36080_r7b_MAT_Minus): 77.7% identity over aa 1-247 of 249, with one 2-nt
  frameshift after "...PKPSR" at the sequence ...CCTT ctAGAAAAAGAAAtctt....
- Classifier on the frame-corrected 253-aa protein: 2026-10-01 build sexM 156.6
  vs sexP 40.9; rebuild 137.9 vs 44.3. On the truncated model: 35.2 vs 4.1, then
  20.1 vs 3.8.
- Results: results/2026-10-04_mucoro_bfd_rebuild/ (compare_vs_campaign.txt);
  analysis/2026-10-04_mycotypha-record-and-classifier-rebuild.md.

## Method that found it
Call-by-call comparison of 288 BFD Mucoromycotina genomes before and after a
classifier rebuild; exonerate alignment of the curated sexM to the locus.

## Limits
Without short reads the frameshift cannot be proven to be an assembly error.
How many other Nanopore-only assemblies lose a label this way is not measured.
