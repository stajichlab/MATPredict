# Pucciniomycotina curation: Sporidiobolales and Pucciniales references

2026-09-26, branch `curation-puccinio` (from `basidio-anchors` 2aa9bf1),
commit `293640d`. Curator ruling: curate Sporidiobolales and Pucciniales now;
check *Wallemia* and add nothing if no MAT locus is described. All records
are pending curator sign-off. Scripts, tables and reports:
`results/2026-09-26_puccinio_curation/` (main checkout).

## Why the 2026-09-20 review found nothing, and what changed

The 2026-09-20 review looked for locus-specific GenBank deposits in
*Puccinia*, *Leucosporidium*, *Rhodotorula* and *Cystobasidium* and found
none it could use. This pass used two sources it did not:

- **Sporidiobolales:** Coelho et al. 2010 (PMID 20700437) and 2011
  (PMID 21880139) deposited complete STE3 receptor CDSs for BOTH MAT alleles
  (A1, A2) in several red-yeast species, some with a complete RibL6 flank.
- **Pucciniales:** no rust HD gene has a locus-specific deposit, but Cuomo et
  al. 2017 (G3 7:361, PMID 27913634) Table S11 names the b mating-type HD
  genes by locus tag in the Pgt, Pt and Pst genome annotations. These are
  admitted as TIER 2 (genome-derived) under the 2026-09-24 ruling.

## Two new families (db/Basidiomycota/order.yml)

| family | scope | vocabulary | genes |
|---|---|---|---|
| `rustHD` | Pucciniales (5258) | pattern `^b[0-9]+$` | bE (HD2), bW (HD1) |
| `redPR` | Sporidiobolales (231213) | enum A1 / A2 | STE3a1 (A1), STE3a2 (A2), RibL6 (optional flank) |

## Records (11; all 17 proteins validate at 100% identity and coverage)

| record | species | source | tier |
|---|---|---|---|
| `5005_cbs-483_redPR_A1` | *Sporobolomyces salmonicolor* | GU474642.1 (Coelho 2010) | 1 |
| `5005_cbs-490_redPR_A2` | *S. salmonicolor* | GU474641.1 (Coelho 2010) | 1 |
| `5003_cbs-491_redPR_A1` | *S. pararoseus* | JN246670.1 (Coelho 2011) | 1 |
| `5003_cbs-499_redPR_A2` | *S. pararoseus* | JN246611.1 (Coelho 2011) | 1 |
| `86836_pycc-4818_redPR_A1` | *Rhodotorula kratochvilovae* | JN246653.1 (Coelho 2011) | 1 |
| `86836_cbs-4785_redPR_A2` | *R. kratochvilovae* | JN246585.1 (Coelho 2011) | 1 |
| `217517_cbs-9654_redPR_A2` | *S. longiusculus* (STE3a2 + RibL6) | JN246612.1 (Coelho 2011) | 1 |
| `418459_crl-75-36-700-3_rustHD_b1` | *Puccinia graminis* f. sp. *tritici* | DS178274.1, PGTG_05143/05144 | 2 |
| `630390_1-1-bbbd-race-1_rustHD_b2` | *P. triticina* | ADAS02000068.1, PTTG_03697/27730 | 2 |
| `1165861_pst-78_rustHD_b1` | *P. striiformis* f. sp. *tritici* | AJIL01000034.1, PSTG_05918/05919 | 2 |
| `747676_98ag31_rustHD_b1` | *Melampsora larici-populina* | GL883124.1, EGG03504.1/EGG03439.1 | 2 |

Annotation problems logged: B5.1 (Pgt bW1 model 579 aa vs 618 aa in the
paper), B5.2 (RefSeq's "MlpbE1" label sits on the bW ortholog; curated by
orthology), B5.3 (Pt b1 models partial), C11 (JN246670.1 title names a RibL6
it has no CDS for).

## Not curated, and why

- **Rust pheromone receptors (STE3.2/STE3.3).** Cuomo 2017 Table S9 names
  them, but a receptor alone cannot pass the 2-modelled-gene bar (the
  pheromone precursors are too short to search) and no flank is defined.
- **Sporidiobolales HD genes.** Only partial HD deposits exist
  (HM143854-HM143857, JN246635.1); HD genes are not near the P/R genes in
  R. toruloides, S. roseus or S. salmonicolor (Coelho 2008).
- **R. araucariae** PYCC 4820 carries both a1 (JN246654.1) and a2
  (JN246588.1) receptors in one strain; left out to keep one idiomorph per
  record and strain.
- **Wallemia.** A putative MAT locus IS described (Padamsee et al. 2012,
  PMID 22326418; Sun et al. 2019 W. mellicola, PMID 31167502; Gostincar et al.
  2019 W. ichthyophaga, PMID 31551960): SXI1, an HMG gene and STE3 in two
  versions. But there is no GenBank deposit of these genes: the sequences are
  in a supplementary file, STE3 was found only by manual search, and the
  papers call the locus "putative". Nothing was added. A tier-2 record from
  the W. mellicola CBS 633.66 annotation is possible and is a curator
  decision.

## Before/after (26 pilot genomes)

Before = the full Basidiomycota run (`run-ad1f865`, old DB, all genomes
phylum-fallback). After = `run-293640d`. Same code; only the DB differs.

| order | genomes | called before | called after |
|---|---:|---:|---:|
| Sporidiobolales | 15 | 4 | 15 |
| Pucciniales | 11 | 0 | 10 |
| Pucciniales, excluding the 3 genomes that supplied a record | 8 | 0 | 7 |

- **Sporidiobolales:** every call is `redPR`, high, `mat_locus`, with RibL6
  as the second gene. The idiomorph is clear: the winning receptor scores
  55-100% identity against 29-41% for the other allele. Species with no
  record are called too: *R. toruloides* (A1 and A2), *R. mucilaginosa* x3,
  *R. graminis*, *R. paludigena*, *Rhodotorula* sp. ZM1.
- **Trade-off:** the 4 "before" calls were Agaricomycete `HD`/`bLocus` hits
  by phylum fallback (HD1/HD2 at 28-40%). With `redPR` scoped to the order,
  Sporidiobolales genomes are routed only to `redPR`, so the HD locus is no
  longer searched there. Those HD hits may be real red-yeast HD genes; no
  red-yeast HD reference exists to test it.
- **Pucciniales:** calls are `rustHD`, `idiomorph_gene_only`, with both bE
  and bW modelled. Non-self genomes: *P. striiformis* GCA_025169535.1 (bE
  95.7%, bW 78.5%), *P. triticina* GCF_026914185.1 and 19NSW04 (78-83%; the
  latter with two HD loci, as a dikaryon should have), and four *Phakopsora
  pachyrhizi* genomes (a different family; bE 49%, bW 31-33%, two loci each).
  The *Phakopsora* calls rest on distant hits and were not checked for
  gene order. *P. striiformis* GCA_002008935.1 stays uncalled (bW fragments
  only, nothing modelled).
- **Runtime:** lineage routing searches one family instead of all.
  Sporidiobolales 254-1,737 s -> 59-94 s; Pucciniales under 750 Mb
  104-444 s -> 21-114 s; the four 1.27 Gb *Phakopsora* genomes
  2,452-2,570 s -> 693-717 s (well inside a 1 h per-genome limit).

## Limits

- n = 26; the *Phakopsora* set is four genomes of one species.
- Several Sporidiobolales pilot genomes are the same species as a record
  (species radius), though not the same strain as far as the metadata show.
- Idiomorph calls in Sporidiobolales are not checked against known mating
  types of the sequenced strains.

## Curator rulings (2026-09-27)

- All 11 records signed off (7 Sporidiobolales tier 1, 4 Pucciniales tier 2).
- Sporidiobolales stay receptor-only (`redPR`) for now; Sporidiobolales HD is
  queued. Next test: extend `rustHD` scope to Sporidiobolales and see whether
  the rust bE/bW genes find their HD locus.
- *Wallemia*: build a tier-2 record from the W. mellicola CBS 633.66
  annotation, labelled putative; report gene adjacency and mating-type mix
  across the 51 genomes.
- Rust receptors (STE3.2/STE3.3) and their pheromone precursors go into the
  queued receptor work; the tiny precursors need a dedicated plan.

## Follow-up results (2026-09-27)

Full note: `results/2026-09-27_puccinio_followup/NOTE.md` (main checkout).

- **rustHD on Sporidiobolales: reverted.** 0/15 rustHD calls; 10/15 genomes
  have a co-located bE+bW cluster at 42-57% identity but nothing models, so all
  were withheld. redPR calls unchanged (22/22); runtime 2.5x. Trial `8a49b00`,
  revert `a6770ca`. Sporidiobolales HD stays queued.
- **Wallemia: putative tier-2 record `671144_cbs-633-66_wallMAT_v1`** (family
  `wallMAT`, `6590929`). BAP31, STE3, CAF1, HMG from JH668224.1 CDS features;
  SXI1 unannotated in the assembly, curated from a miniprot model (278 aa).
  Detect: 51/51 Wallemiales called (0 before); 33 high (v1 type), 18 medium
  (other version, divergent STE3). Check 1: SXI1/HMG/STE3 co-located in 33/33 v1
  genomes. Check 2: W. mellicola 14 v1 : 13 other; W. ichthyophaga 18 : 4, the 4
  being exactly the paper's inverted strains. Pending curator sign-off.
