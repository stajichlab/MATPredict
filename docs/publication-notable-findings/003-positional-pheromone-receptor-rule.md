# 003. Mating pheromone receptors identified by position, from genome DNA alone

- **Category:** method
- **Status:** candidate (small ground truth: 9 receptors in 6 genomes)
- **Lineage:** Basidiomycota (works in Agaricomycetes and Ustilago; not in rusts)

## Summary
Agaricomycete genomes carry 3-7 STE3-like receptor copies, and a tree does not
separate the mating receptors from the others. A short ORF ending in a strict
CAAX motif within 10 kb of an STE3 locus flags every known Agaricomycete mating
receptor and none of the other copies. The scan uses genome DNA only, because
annotations miss tiny pheromone genes.

## Evidence
- Tree result: `results/2026-09-27_receptor_explore/NOTE.md`, `receptor_tree.pdf`
  (238 proteomes, 1,060 PF02076 hits + 34 curated receptors; FastTree).
  Sporidiobolales A1 and A2 each form a clade (local support 0.94, 0.88);
  Agaricomycete mating receptors do not.
- Rule and results: `results/2026-09-27_pheromone_positional/NOTE.md`,
  `gt_eval.txt`, `uncurated_loci.tsv`.
  - Ground truth (Coprinopsis, Schizophyllum, C. neoformans JEC21, U. maydis
    521, R. toruloides CBS 14 and NBRC 0880; 9 mating, 25 other STE3 loci):
    strict `C[VI][IV][AVMG]` within 10 kb flags 6/9 mating and 0/25 other;
    random same-size windows 2.5%. Coprinopsis 3/3, Schizophyllum 2/2.
    With homology to curated precursors: 4/9 and 0/25, random ~0.1%.
  - The textbook CAAX alphabet `C[AVLIM][AVLIM]X` flags 42% of random windows.
  - 52 uncurated Agaricomycete genomes, 566 STE3 loci: strict 75 flagged vs
    13.6 expected (5.5x); strict + homology 36 vs 0.3. Polyporales 13/15
    genomes, Russulales 8/12, Boletales 2/15, Hymenochaetales 1/10.
  - Flagged copies sit nearer the curated mating receptors in the tree (median
    patristic distance 1.12-1.19 vs 2.20).
  - Rusts (4 genomes, 16 loci): 1 flagged, 0 with homology.
- Misses: Cryptococcus (pheromone 45 kb away); Rhodotorula, whose candidate
  precursors end CTIA/CTVA (caught only by the relaxed class T2).
  CTVA is the C-terminus of the rhodotorucine A precursor (Akada et al. 1989;
  literature, not re-checked in this project).

## Method that found it
Six-frame stop-anchored ORF scan (in-frame Met 20-130 codons upstream, Cys at
-4); miniprot of 1,094 STE3 proteins to locate receptor loci; chance from 1,000
random windows per genome (`scan_genome.py`, `evaluate_gt.py`,
`analyze_uncurated.py`). Committed in polish-scope-cuts `2a310aa`.

## Verification done / still open
- Done: tested against curated receptors; chance rate measured.
- Open: curation of Polyporales and Russulales receptor-precursor clusters
  (receptor queue step c, running 2026-09-27); a lineage motif for
  Sporidiobolales (e.g. `CT[IV]A`); rust precursors.

## Limits
Ground truth is small; unlabelled copies are assumed non-mating; the homology
rule is weakest where precursors diverge (Boletales, Hymenochaetales, rusts,
Rhodotorula).

## Related
`docs/notes/2026-09-20_short-pheromone-orf-detection.md`; annotation report
B1.6 (U. maydis mfa1 missing from the modern annotation).
