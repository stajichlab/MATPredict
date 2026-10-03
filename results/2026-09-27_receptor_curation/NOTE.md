# Polyporales and Russulales B (PR) locus records (2026-09-27)

Receptor queue step (c), curator ruling 2026-09-27. Targets: the
receptor + strict-CAAX precursor clusters found by the positional scan
(results/2026-09-27_pheromone_positional/). Commit `29818af`, branch
basidio-anchors. All records tier 2 (genome-derived), pending sign-off.

## Records

| record | assembly / contig : span | receptors | precursors (strict CAAX) | source |
|---|---|---|---|---|
| 984962_tc-32-1_PR_B1 (Heterobasidion irregulare, Russulales) | GCA_000320585.2 KI925460.1:825,882-849,398 | 5 (ETW79702/03/04/06/09) | 3 (ETW79705 CILF, ETW79710 CTIA, ETW79712 CVVA; 1 strict) | JGI annotation, Olson et al. 2012 |
| 5325_fp-101664-ss1_PR_B1 (Trametes versicolor, Polyporales) | GCA_000271585.1 JH711790.1:1,556,459-1,581,006 | 3 (EIW56739/40/41) | 9 (7 strict: CIIA/CVIA) | JGI annotation, Floudas et al. 2012 |
| 5627_9006-11_PR_B1 (Grifola frondosa, Polyporales) | GCA_001683735.1 LUGG01000005.1:512,153-547,096 | 5 (OBZ74837/42, OBZ74475, OBZ74569, OBZ74474) | 3 unannotated ORFs by coordinates (CVIA, CIIG, CIIA; all strict); the last is the ortholog of WM1-25 ph1 LC706365.1 (tblastn 55%, E=2e-11) | locus existence tier 1: tetrapolar crosses, Zhang et al. 2023 (PMID 37888215) |

- Every annotated gene's translation from the recorded coordinates equals
  its annotated protein; unannotated precursors start with Met and have no
  internal stop. `curate-db build-gff` derives the same proteins.
- Segment sources are the INSDC contigs, so coordinate-only precursors are
  derivable. Idiomorph `B1` is a placeholder (alleles counted, curator
  convention 2026-09-26).

## Family choice

Kept in `Basidiomycota:PR` with scope widened to Polyporales (5303) and
Russulales (452342). The family's roster is generic (`pheromone_receptor`,
`fungal_mating_type_pheromone`), so no new family or gene names are needed,
and the same family then covers every Agaricomycete order with a curated B
locus.

## Validation (35 BFD genomes; before run-f26c397, after run-29818af)

| | Polyporales (20) | Russulales (15) |
|---|---|---|
| PR called, before -> after | 0 -> 9 genomes (11 calls) | 0 -> 1 (1 call) |
| PR calls at a strict-CAAX flagged STE3 locus | 7 of 11 | 1 of 1 |
| PR call on an HD contig | 1 of 11 | 0 |
| HD calls | 20 -> 20 | 14 -> 14 |
| median wall per genome | 89 -> 364 s (4.1x) | 69 -> 153 s (2.2x) |

- All 12 PR calls are medium, `idiomorph_gene_only`, each with a modelled
  receptor (58-84%) and a modelled precursor (64-71% in the calls checked).
- 176 PR clusters are withheld (37 with one modelled gene, 139 with none).
- PR and HD are on different contigs in 11 of 12 calls (tetrapolar layout);
  the one exception (Trametes meyenii) is not checked.
- Russulales gains little: one record from Heterobasidion does not reach
  Lactarius, Russula, Hericium, Auriscalpium.
- Multi-copy counts are not measurable from the reports: gene_evidence keeps
  one entry per gene name, so every call shows 1 receptor + 1 precursor.

## Not usable

- Coriolopsis trogii b1-1/b1-2/b2-1/b2-2 (MF990238-41): receptor genes only,
  a direct submission with no paper and no precursors.
- Phanerochaete STE3.1-3.5 (HQ188385-421), Wolfiporia STE3 (KF044372-81,
  LC672030-2), Amylostereum RAB1 (EU380312-3), Lenzites (AY226021): partial
  receptor CDS only; the B locus cannot be delimited from them.
- Ganoderma boninense (LR7255xx, ON0599xx): receptor fragments and precursor
  mRNAs without a genomic locus.
- Annotation errors logged: B6.1 (Grifola precursors unannotated), B6.2
  (Trametes EIW57105.1 ends CTIAW, not CAAX).

## Curator decisions

1. Sign off the three records.
2. Runtime: PR search raises Polyporales wall time 4.1x; accept or cap.
3. Russulales needs a second record from Russulaceae (Lactarius/Russula) to
   be covered.
