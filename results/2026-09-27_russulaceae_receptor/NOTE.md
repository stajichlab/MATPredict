# Russulaceae B (PR) record: Russula nobilis gfRusNobi1 (2026-09-27)

Record 2830151_kdtol00553_PR_B1 (branch basidio-anchors, commit 4a0f0eb),
tier 2, pending sign-off.

## Choice of genome
Candidates with a strict-CAAX receptor-precursor cluster in the positional
scan: Russula nobilis (2 strict clusters), R. versicolor (2), Lactarius
sanguifluus (1). R. nobilis GCA_984573805.1: chromosome level, 123 contigs,
contig N50 1.48 Mb; three tandem STE3 receptors (Pfam PF02076 at GA) on
OZ475200.1:1,514,818-1,522,262, best blastp hits the Heterobasidion B
receptors (38.5-53.5%, E <= 1.7e-100). R. versicolor has only one full-length
receptor at its cluster (the other hit is 243 aa); L. sanguifluus is a
978-contig assembly. No Russulaceae mating-type deposit exists (NCBI
nucleotide: WGS scaffolds and mip partial CDS only; PubMed: no paper).
The public assembly has no gene annotation: receptor exons are the BFD
funannotate models FD20D505_007127/8/9, recorded by coordinates. Precursor:
unannotated ORF 1,529,356-1,529,433(-), 25 aa MVGKPTSRKRKRYEIRNLGMCCVVG,
strict CAAX CVVG, 7.1 kb from the third receptor; no homology to curated
pheromones (tblastn E >= 0.1).

## Validation (15 Russulales genomes; before run-29818af, after run-4a0f0eb; db-only difference)
- PR genomes called: 1 -> 1 (Heterobasidion annosum only, unchanged). 0 of 8
  Russulaceae genomes called, INCLUDING R. nobilis itself.
- HD: 14/15 called. Median wall 153 -> 161 s.

## Why the record does not call its own genome
- At OZ475200.1 the three receptors form one cluster (1,510,837-1,522,259,
  126 hits, best identity 97.2%) with gene_count 1: all three receptors carry
  the family's generic name `pheromone_receptor`, which counts as ONE
  distinct gene. The admission floor is >= 2 distinct genes (pipeline.py
  `>=2 distinct genes`), so the cluster is never admitted.
- The 25-aa precursor gets no search hit at all (no evidence row near
  1,529,356), although it is in the reference set: a 25-aa tblastn query
  cannot reach significance genome-wide. So the second distinct gene never
  appears.
- The four withheld PR loci in this genome are elsewhere (0 modelled genes).

## Consequence
Adding Russulaceae references cannot help until detection can admit a
receptor-only cluster or find a tiny precursor. That is the positional
CAAX filter / short-pheromone search already queued (receptor queue b), or a
rule that counts distinct receptor copies. The record itself is sound (all
genes re-derived; tests pass) but has no detection effect today.
