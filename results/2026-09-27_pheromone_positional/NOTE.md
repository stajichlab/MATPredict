# Receptor queue, step 2: does a nearby pheromone precursor mark the mating receptor?

2026-09-27. Read-only exploration. Scripts: `scan_genome.py`, `run_gt.sh`,
`evaluate_gt.py`, `run_uncur.sh`, `analyze_uncurated.py`. Outputs: `out_gt/`,
`out_uncur/`, `out_rust/`, `gt_eval.txt`, `gt_loci_labelled.tsv`,
`uncurated_loci.tsv`, `uncurated_per_genome.tsv`, `uncurated_strict_summary.txt`.

## Method (genome DNA, no annotation used)

- STE3-like loci: miniprot of the step-1 set (1,060 BFD STE3 hits + 34 curated
  mating receptors) against the genome; alignments with >=50% query cover and
  genomic span <=8 kb, merged when they overlap.
- Precursor candidates: stop-anchored six-frame scan. A candidate is a stop
  codon whose upstream in-frame segment has an in-frame Met 20-130 codons
  before the stop and a cysteine at -4. Classes:
  - T strict: `C[VI][IV][AVMG]` (tightened on curated Basidiomycota
    precursors, docs/notes/2026-09-20_short-pheromone-orf-detection.md)
  - T2: `C[VITE][IVT][AVMG]` (adds T/E at -3, T at -2)
  - L textbook: `C[AVLIM][AVLIM]X`
  - R: >=2 short ORFs with an identical 8-aa C-terminal tail (tandem
    precursor copies)
  - Hx: tblastn hit of the 31 curated precursors (-seg no, word 2, E<=1),
    excluding precursors of the genome's own genus (held out)
- Counts within +-10, 20, 50 kb of each locus; chance = the same counts in
  1,000 random windows of equal size per genome, avoiding STE3 loci.

## Ground truth (6 genomes: Coprinopsis, Schizophyllum, C. neoformans JEC21,
## U. maydis 521, R. toruloides CBS 14 and NBRC 0880)

34 STE3-like loci; 9 labelled mating (record coordinates, or >=95% to the
strain's own record receptor); 25 other.

| rule | mating | other | random window |
|---|---|---|---|
| T strict, +-10 kb | 6/9 | 0/25 | 2.5% |
| T strict + Hx, +-10/20 kb | 4/9 | 0/25 | ~0.1% |
| T2 or R, +-20 kb | 9/9 | 0/25 | 13.6% |
| L textbook, +-10 kb | 6/9 | 9/25 | 42% |

- The textbook CAAX alphabet is useless (random windows 42%).
- Strict T at 10 kb: all Agaricomycete mating receptors flagged
  (Coprinopsis 3/3, Schizophyllum 2/2) and Ustilago; misses Cryptococcus
  (MFalpha is 45 kb away, outside the window) and both Rhodotorula (their
  precursors end CTIA/CTVA, caught only by T2 or R).

## Uncurated Agaricomycete orders (52 genomes, 566 STE3-like loci)

| rule | observed | expected by chance | ratio |
|---|---|---|---|
| T strict, +-10 kb | 75 | 13.6 | 5.5 |
| T strict + Hx | 36 | 0.3 | ~130 |
| T2 or R, +-20 kb | 141 | 101 | 1.4 |

Strict rule per order (genomes with >=1 flagged copy / flagged copies;
T+Hx clusters within 50 kb):
- Polyporales 13/15 genomes, 45 copies; T+Hx in 11/15 genomes, 13 clusters.
- Russulales 8/12, 20 copies; T+Hx in 4/12, 4 clusters.
- Boletales 2/15, 3 copies; Hymenochaetales 1/10, 7 copies; T+Hx 0.

Tree (step-1 FastTree): copies flagged by the strict rule sit closer to the
curated Coprinopsis/Schizophyllum mating receptors than unflagged copies
(median patristic distance 1.12-1.19 vs 2.20).

## Rusts (4 genomes: P. graminis, M. larici-populina, P. triticina Pt76,
## P. striiformis)

16 STE3-like loci; 1 flagged by strict T, none by T+Hx. The rule finds
essentially nothing in rusts.

## Limits

- Ground truth is 9 mating receptors in 6 genomes; the other-copy labels
  assume unlabelled copies are non-mating.
- The best-specificity rule (T+Hx) uses homology to curated precursors, which
  are mostly Agaricales/Tremellales/Ustilago; it is least sensitive where
  precursors diverge (Boletales, Hymenochaetales, rusts, Rhodotorula).
- Loci merge nearby receptors; multi-copy B loci are normal.
- Assembly gaps at loci and genome-specific CAAX rates vary (random-window
  rate 0-6% for strict T).
