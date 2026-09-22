# Tremellales at scale, and what the figures say about classification (2026-09-21)

334 Tremellales genomes from BFD, genome-only, `--taxid` routing, cross-contig
fragment merging off. 334/334 completed, 0 failures, ~11 s/genome. Plus the
*C. deneoformans* JEC21/JEC20 pair run in both merging configurations.

Reports: `results/2026-09-21_tremellales334/` (288 KB + 3 KB, zstd).
Figures are rendered from those reports; none is constructed.

## 1. The 120 kb Tremellales gap works

The curator set `max_cluster_gap_bp: 120000` for Basidiomycota `MAT` on
2026-09-21 from the two curated *C. deneoformans* records. It was untested. It
works, and the locus lands within 2 bp of the curated span:

| | curated record | detected | idiomorph | conf |
|---|---|---|---|---|
| JEC21 | 103,914 bp | **103,912 bp** | alpha ✓ | high |
| JEC20 | 96,538 bp | **96,539 bp** | a ✓ | high |

Each found exactly its own idiomorph's complete gene set. Under the previous
25 kb default this locus could not be clustered at all — its largest internal
gene gap is 74,899 bp. Routing resolved `lineage` → 1 family, so the
order-level scope set the same day also works.

## 2. Cross-contig merging: OFF is correct, and the two changes do different jobs

| run | loci | idiomorph | confidence |
|---|---|---|---|
| JEC21 merge on | 11 | alpha ✓ | high |
| JEC21 merge off | 11 | alpha ✓ | high |
| JEC20 merge **on** | 10 | **undetermined ✗** | **low** |
| JEC20 merge off | 11 | a ✓ | high |

Figures D and E are the same genome and differ only in that flag. **The gap
finds the locus; the merging default calls it correctly.** PR #4's choice of
`off` by default is supported by this pair.

Figure D also shows the fragmented-reporting defect the Basidiomycota note
predicted: the entry's top-level `contig/start/end` claims 96,539 bp on
CM152755.1 while its evidence actually spans 1,262,069 bp across two contigs.

## 3. The references do not generalise past *Cryptococcus*

`high`-confidence calls, by genus of the genome:

```
Cryptococcus  243 genomes   69 with a high call
Kwoniella      24 genomes    0
Dioszegia      14 genomes    0
Tremella       12 genomes    0
Naematelia     10 genomes    0
...everything else            0
```

Both curated records are *C. deneoformans*. This is the same pattern Finding 3
of the Mucoromycota scale note measured — detection quality tracks reference
distance, not genome quality. **The gap value generalises across the order; the
reference set does not.** 128 of 334 genomes have exactly one high call, 3 have
two, 203 have none.

## 4. `high` already separates signal from noise

| conf | n | median best identity | median genes | median span | evidence |
|---|---|---|---|---|---|
| **high** | 134 | **94.7%** | 5 | **44,536 bp** | 100% polished/proteome |
| medium | 1,384 | 43.8% | 3 | **113 bp** | 79% tblastn-only |
| low | 2,283 | 32.8% | 1 | **152 bp** | 100% tblastn-only |

A median medium/low call is **a single ~130 bp tblastn HSP**. These are not
over-merged clusters — the 120 kb gap never chained them, because there is
nothing to chain.

## 5. The over-calling has one specific, fixable cause

**3,210 of 3,801 calls (84%) are a single genomic interval.** The gene-name
sets responsible:

```
 1197  [SXI1]                                  low     1 name  — correct
  960  [SXI2]                                  low     1 name  — correct
  457  [MFa1, MFa2, MFa3]                      MEDIUM  3 names, ONE ORF
  336  [MFalpha1, MFalpha2, MFalpha3]          MEDIUM  3 names, ONE ORF
   68  [MFalpha1, MFalpha2]                    MEDIUM  2 names, ONE ORF
   66  [MFa1..3, MFalpha1..3]                  MEDIUM  6 names, ONE ORF
```

**927 of the 1,384 medium calls are one ORF wearing several gene names.**
`EvidenceFloor.min_hits` counts distinct gene NAMES, and MFa1/MFa2/MFa3 are
three paralogous pheromone-precursor entries in the curated roster that all
match the same tiny ORF — so a single HSP satisfies ">=2 distinct genes".

This is the same bug class as the sexM/sexP cross-match already handled by
`idiomorph.resolve_idiomorph_overlaps`, but that resolver only collapses
MUTUALLY EXCLUSIVE idiomorph pairs. MFa1/2/3 share one idiomorph, so it never
fires. (The 66 mixed MFa/MFalpha cases *are* mutually exclusive and arguably
should already have collapsed — worth a separate look.)

**Every `high` call has >=2 distinct intervals** (median 4). Distinct
non-overlapping intervals looks like a much better admission unit than distinct
gene names. That change is not made here: it is a semantics change to the
evidence floor and belongs to the curator, and a methods review is in progress.

## 6. `mat_locus` is unreachable for 13 of 19 families

`classify_locus` returns `mat_locus` only for a cluster holding a core gene AND
a non-core gene. Families whose roster is entirely `core_MAT` can never satisfy
it:

```
UNREACHABLE  Ascomycota     MTL, PM, mat1, mat2, mat3
UNREACHABLE  Basidiomycota  Aalpha, Abeta, Balpha, Bbeta, HD, MAT, PR, bLocus
reachable    Ascomycota     MAT, MATsc, MATtub, MATyl
reachable    Basidiomycota  aLocus            (rba1 is flanking_variable)
reachable    Mucoromycota   MAT
```

8 of 9 Basidiomycota families. Consequences: the complete JEC21 locus in
Figure A is filed as `idiomorph_gene_only`, the same class as the 95 bp stray
HSP in Figure B; a Basidiomycota `locus_class_tally` reads ~100%
`idiomorph_gene_only` and ~0 loci; and `partial_locus` — which displaces only
`mat_locus` — can never fire in Basidiomycota.

This is pre-existing, not introduced by this branch. It is the same shape as
two bugs already fixed here: `assign_tier` counting `genes_not_searchable`, and
the `pheromone_receptor` alias permanently blocking `core_found`. Both were
"requiring something the roster makes impossible". Not fixed here — the two
candidate fixes (curate real flanking genes into those rosters, or treat a
roster with no non-core gene as satisfying the requirement vacuously) are
curation decisions.

## 7. What this note does NOT establish

* **No ground truth for 332 of the 334 genomes.** Only JEC21 and JEC20 have a
  known answer. Every other statement here is internal consistency or a
  population distribution — not accuracy.
* **Nothing about whether the 134 `high` calls are correct.** They look right
  by identity and structure; that is not the same as being verified.
* **Nothing about `Kwoniella` MAT architecture.** 0 high calls there means the
  references did not reach it. It does NOT mean *Kwoniella* lacks a MAT locus,
  and it does not test whether 120 kb suits *Kwoniella*'s architecture.
* **The 927-call diagnosis is a mechanism, not a fix.** No change has been made
  and the effect of changing the admission unit is unmeasured.

---

# Figures

All rendered from the archived reports.

```
========================================================================================
FIGURE A -- a REAL, COMPLETE locus that is NOT called mat_locus
  C. deneoformans JEC21 (GCF_000091045.1), genome-only, --taxid 214684.
  Curated record 40410_jec21_MAT_alpha spans 103,914 bp; this call is 103,912 bp.
  All 5 alpha-idiomorph genes, every one at 100% identity, idiomorph resolved, high.
  Classified idiomorph_gene_only ONLY because Basidiomycota MAT declares no
  non-core gene, so classify_locus's flanking requirement can never be satisfied.
----------------------------------------------------------------------------------------
  contig NC_006686.1  header says 1,534,710-1,638,622  (span 103,912 bp)
  class=idiomorph_gene_only   confidence=high   idiomorph=alpha   pass=strict
  genes_found=['MFalpha1', 'MFalpha2', 'MFalpha3', 'STE3', 'SXI1']

  #                                  ##                                                ###
  |......................................................................................|
  1,534,710                                                                    1,638,622

  pos (offset)        len  genes at this interval                        ident via
  #            0    113  MFa1,MFa2,MFa3,MFalpha1,MFalpha2             100.0% Xt   (not counted)
  #          727    113  MFalpha3                                     100.0% X
  #       42,762   1317  STE3                                         100.0% X
  #      102,555   1357  SXI1                                         100.0% X

========================================================================================
FIGURE B -- a lone pheromone fragment: idiomorph_gene_only is CORRECT here
  Same genome, different contig. 95 bp, 37.5% identity, tblastn only, nothing near it.
  SAME CLASS AS FIGURE A. One label covers a complete 104 kb locus and a stray HSP.
----------------------------------------------------------------------------------------
  contig NC_006685.1  header says 2,042,184-2,042,279  (span 95 bp)
  class=idiomorph_gene_only   confidence=medium   idiomorph=a   pass=strict
  genes_found=['MFa1', 'MFa2', 'MFa3']
  genes_missing=['SXI2', 'STE3']

  ########################################################################################
  |......................................................................................|
  2,042,184                                                                    2,042,279

  pos (offset)        len  genes at this interval                        ident via
  #            0     95  MFa1,MFa2,MFa3                                37.5% t

========================================================================================
FIGURE C -- the middle ground the classification has to resolve
  Two genes 70 kb apart, 38.2% and 26.7%, both tblastn-only.
  Degenerate remnant, or two unrelated paralogs the 120 kb gap chained together?
----------------------------------------------------------------------------------------
  contig NC_006682.1  header says 390,408-461,327  (span 70,919 bp)
  class=idiomorph_gene_only   confidence=medium   idiomorph=alpha   pass=strict
  genes_found=['STE3', 'SXI1']

  #                                                                                     ##
  |......................................................................................|
  390,408                                                                        461,327

  pos (offset)        len  genes at this interval                        ident via
  #            0    101  SXI1                                          38.2% t
  #       70,626    293  STE3                                          26.7% t

========================================================================================
FIGURE D -- cross-contig merging ON: both idiomorphs, call collapses
  JEC20. Two weak tblastn hits on ANOTHER contig are merged in, bringing SXI1
  (alpha) beside SXI2 (a). Idiomorph -> undetermined, confidence -> low.
  Note the header span (96 kb) against the true evidence extent (1.26 Mb).
----------------------------------------------------------------------------------------
  contig CM152755.1  header says 1,555,944-1,652,483  (span 96,539 bp)
  !! evidence actually spans 390,414-1,652,483 (1,262,069 bp) across 2 segment(s);
  !! fragmented=True -- the header names segment 0 only
  class=idiomorph_gene_only   confidence=low   idiomorph=undetermined   pass=strict
  genes_found=['MFa1', 'MFa2', 'MFa3', 'STE3', 'SXI1', 'SXI2']

  #   #                                                                           ##    ##
  |......................................................................................|
  390,414                                                                      1,652,483

  pos (offset)        len  genes at this interval                        ident via
  #            0    101  SXI1 [CM152764.1]                             38.2% t
  #       70,626    293  STE3 [CM152764.1]                             26.7% t
  #    1,165,530    125  MFa1,MFa3,MFalpha1,MFalpha2,MFalpha3         100.0% Xt   (not counted)
  #    1,174,037    125  MFa2                                         100.0% X
  #    1,177,354   1323  STE3                                         100.0% X
  #    1,259,548   2521  SXI2                                         100.0% X

========================================================================================
FIGURE E -- the SAME genome, merging OFF: correct call
  Identical to D except one flag. SXI1 gone, idiomorph resolves to 'a', high.
  This pair is the evidence for PR #4's default.
----------------------------------------------------------------------------------------
  contig CM152755.1  header says 1,555,944-1,652,483  (span 96,539 bp)
  class=idiomorph_gene_only   confidence=high   idiomorph=a   pass=strict
  genes_found=['MFa1', 'MFa2', 'MFa3', 'STE3', 'SXI2']

  #      #  ##                                                                        ####
  |......................................................................................|
  1,555,944                                                                    1,652,483

  pos (offset)        len  genes at this interval                        ident via
  #            0    125  MFa1,MFa3,MFalpha1,MFalpha2,MFalpha3         100.0% Xt   (not counted)
  #        8,507    125  MFa2                                         100.0% X
  #       11,824   1323  STE3                                         100.0% X
  #       94,018   2521  SXI2                                         100.0% X

========================================================================================
FIGURE F -- what mat_locus looks like when the roster HAS flanking genes
  Absidia cuneospora RSA 623, Mucoromycota. Core sexP in a flanking neighbourhood.
  Note sexM at 36% via tblastn overlapping sexP at 48% via the proteome -- the
  idiomorph cross-match, resolved and marked 'not counted'.
----------------------------------------------------------------------------------------
  contig scaffold_217  header says 12,643-36,677  (span 24,034 bp)
  class=mat_locus   confidence=high   idiomorph=Plus   pass=strict
  genes_found=['algA', 'glrA', 'rnhA', 'sexP', 'tptA']

  ------ ======  ###===================                                           --------
  |......................................................................................|
  12,643                                                                          36,677

  pos (offset)        len  genes at this interval                        ident via
  -            0   1642  algA                                          84.6% P
  =        2,052   1279  tptA                                          75.9% P
  #        4,296    608  sexP                                          48.3% P
  #        4,362    251  sexM                                          36.0% t   (not counted)
  =        5,122   4929  rnhA                                          52.3% P
  -       22,346   1688  glrA                                          77.7% P

========================================================================================
FIGURE G -- partial_locus: the 0.500 tie
  Mucor sp. NRRL A-21230. sexP+sexM+glrA = 3 of 6 = exactly the floor.
  Two ~32% HMG fragments 582 bp apart, plus one real gene 5 kb away.
  8 of the 23 ground-truth genomes carry this same ~6.5 kb element.
----------------------------------------------------------------------------------------
  contig scaffold_131  header says 2,100-8,583  (span 6,483 bp)
  class=partial_locus   confidence=medium   idiomorph=undetermined   pass=strict
  genes_found=['glrA', 'sexM', 'sexP']
  genes_missing=['tptA', 'rnhA', 'algA']

  ###    ###                                                          --------------------
  |......................................................................................|
  2,100                                                                            8,583

  pos (offset)        len  genes at this interval                        ident via
  #            0    182  sexP                                          31.1% t
  #          582    152  sexM                                          32.7% t
  -        5,080   1403  glrA                                          69.0% P

========================================================================================
FIGURE H -- ONE ORF wearing three gene names defeats the evidence floor
  GCA_000442785.1_Cf_30_300r_Split10plusN (Kwoniella/Cryptococcus panel).
  The admission bar is '>=2 DISTINCT genes including >=1 core_MAT'.
  MFa1/MFa2/MFa3 are three roster entries that all match the SAME 95 bp ORF,
  so one tblastn HSP counts as three genes and clears the bar.
  927 of the 1,384 medium calls in the 334-genome panel are exactly this.
----------------------------------------------------------------------------------------
  contig CAUG01000574.1  header says 17,543-17,638  (span 95 bp)
  class=idiomorph_gene_only   confidence=medium   idiomorph=a   pass=strict
  genes_found=['MFa1', 'MFa2', 'MFa3']
  genes_missing=['SXI2', 'STE3']

  ########################################################################################
  |......................................................................................|
  17,543                                                                          17,638

  pos (offset)        len  genes at this interval                        ident via
  #            0     95  MFa1,MFa2,MFa3                                40.6% t

```
