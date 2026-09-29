# Next fixes: MAT-gene gate, F6, CAAX label, record assembly field (2026-09-28)

Branch `next-fixes` from PR #9 3aec88b. Commits: 801c55f (MAT-gene gate),
e41b08b (F6), 2993ee1 (all CAAX-dependent calls unverified), a9fb69c
(`locus.assembly_accession` + backfill), d8f70a5 (self-check reads
assembly-type segments and merged calls). 900 tests pass (882 before).
Frozen worktrees: run-a9fb69c (Mucoromycota, Zygo), run-nf-measure
(merge-measure2's Basidiomycota records rebased onto a9fb69c: Agaricales
panel, self-call panels). Baselines: results/2026-09-28_review_fixes/ (f25cf70
and its Agaricales panel).

## Zygo 23
scaffold 23/23 locus, 23/23 idiomorph; contig 23/23, 23/23 (zygo23_a9fb69c).

## Mucoromycota (293 genomes; compare_mucoro.txt)
Genomes called 258 -> 253. 17 calls withheld by the MAT-gene gate, 3 gained,
12 medium -> high, 0 label changes.
- Gate (item 1), 17 withheld: 16 were `undetermined` (no mating type) and 1
  a Benjaminiella second call (Minus/medium, score 86.9, no flanks; the genome
  keeps its Minus/high call). Scores 34-77 bits, all model-typed, 0-1
  supporting flanks. The HMG-paralog pairs are gone: Syncephalastrum x6
  (score 64-67, glrA alone), Circinella x2.
- 5 genomes lose every call, all of which were `undetermined` only: Rhizomucor
  miehei, R. pusillus FCH_5_7 and Rhipu1 (score 39.5, no flanks >= 40%),
  Circinella umbellata, Phascolomyces articulosus. FLAG: the Lichtheimiaceae
  (Rhizomucor) loci lack Mucorales flanks by biology (flank-synteny test
  failed), so the gate can keep them only at >= 100 bits.
- 3 gained, via the relaxed pass once the strict paralog call is withheld:
  S. racemosum NRRL 2496 MCGN01000004.1:1,755,506 Plus/medium (the curated
  record's own locus, not called before), S. monosporum B8922 Minus/medium,
  S. monosporum PYS2302 Plus/medium.
- F6 (item 2), 12 medium -> high: 11 had only sub-39-bit unmodelled genes; 1
  (GCA_019677225.1) had a sub-floor tptA plus the other-allele sexP already
  ignored by the allele-absent rule.

## Agaricales CAAX panel (128 genomes; compare_agaricales.txt)
Calls 281 -> 281, genomes called 124 -> 124. Unverified calls 29 (21 genomes)
-> 102 (77 genomes): every CAAX-dependent call now labelled. F6: 2
medium -> high (both sub-floor unmodelled genes).

## Record self-call (selfcall_results.tsv)
The 22 backfilled records plus 10 Basidiomycota records whose assembly is in
their segments, run on run-nf-measure's database.
- 27 called.
- Not real misses (all called at the record coordinates, checked by hand):
  C. auris B11221 (RefSeq NW_021640163.1 vs GenBank PGLS01000002.1, identical
  coordinates); Trametes FP-101664 and Heterobasidion TC 32-1 (BFD holds the
  GCF assemblies; PR called on NW_ contigs at the record start). The check
  matches contigs by accession, so GCA/GCF pairs need a map (not built).
- fen-198 (GCA_964255735.1): record has no segment coordinates; both
  MAT1-1 and MAT1-2 are called in the genome.
- 5011_unknown-1_PM_combined: withheld (modelled_gene_bar) -- a real
  withhold, not diagnosed.

## Assembly backfill (backfill_written.tsv)
115 records: 22 resolved (NCBI elink nuccore->assembly + assembly esummary),
92 locus-specific deposits in no assembly (null), 1 without a segment
accession, 0 ambiguous. Every record still validates.
