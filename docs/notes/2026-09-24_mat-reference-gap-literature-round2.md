# MAT reference gap lineages: literature and accession survey, round 2

Date: 2026-09-24. Round 1 is `2026-09-24_mat-reference-gap-literature.md`. This
note does not repeat round 1 findings except where round 2 extends or corrects
them. Nothing in `db/` was changed.

## Method and what the verification codes mean

- Literature: PubMed (E-utilities and the PubMed tool), PMC full text, one web
  search.
- **Y** = nucleotide accession fetched from NCBI nuccore on 2026-09-24 (esummary
  for length/title/organism, efetch GenBank flat file for CDS list). Length,
  definition line, CDS products and protein IDs quoted below come from that fetch.
- **G** = genome-derived coordinate. No locus deposit exists. I measured it with
  local `tblastn` (BLAST+ 2.17.0, `-seg no`) of curated or published proteins
  against the named NCBI sequence. It is a measurement, not a published
  coordinate. E-values and identities are given so the strength is visible.
  Re-check before curation.
- **N** = no verified accession. The paper exists, but the loci are genome-only
  or I did not locate a deposit.
- "No deposit found" means none found by organism + keyword + length searches
  and PubMed-to-nuccore links. It is not proof of absence.
- HMG-box and homeodomain queries hit many non-MAT paralogs. A short (70-100 aa)
  HMG or homeobox hit without synteny is NOT evidence of a MAT gene. I say so
  where it applies.

## Summary table (ranked by expected impact = BFD genomes x distance from any curated record)

| rank | lineage | BFD genomes | pilot call rate | best candidate (species, idiomorph, accession) | verified | architecture flag |
|---|---|---|---|---|---|---|
| 1 | Serinales | 2,368 | 2/75 | *C. albicans* SC5314 MTLa `AF167162.1` + MTLalpha `AF167163.1`; *C. lusitaniae* P1 MTLalpha `CP038489.1`:389,970-403,521 | Y / G | PAP-OBP-PIK inside idiomorph; alpha2 absent in Metschnikowiaceae; D. hansenii and S. stipitis are homothallic-contiguous; *L. elongisporus* and *C. sojae* have NO MAT |
| 2 | Boletales | 511 | not piloted | *Rhizopogon roseolus* Rrsb10 A locus + PR genes `LC662829.1` (121.7 kb) | Y (unpublished) | Bipolar in *R. roseolus* (HD and PR on one 122-kb fragment); HD multiallelic in suilloids |
| 3 | Dipodascales (non-*Yarrowia*) | 398 | not in pilot | none deposited; *Blastobotrys* alpha1 between SLA2 and APN2 (Krassowski 2019) | N | MATalpha2 absent in *Blastobotrys*; *B. proliferans* PHN |
| 4 | Sporidiobolales | 253 | not piloted | none deposited as a locus; Coelho 2008/2010 genome-based | N | HD and PR linked (pseudo-bipolar), HD multiallelic |
| 5 | Polyporales | 198 | not piloted | *Phanerochaete chrysosporium* `HQ188438-43.1` (MIP + HD2 + HD1); *Mycoleptodonoides* `LC730474.1` | Y | Bipolar HD-only (*P. chrysosporium*); MIP and beta-fg flanks |
| 6 | Phaffomycetales | 188 | not piloted | none deposited | N | Flip/flop FF1 in *Cyberlindnera saturnus* clade; FF2 *W. canadensis*; PHC *W. hampshirensis* |
| 7 | Kickxellomycota | 171 | not piloted | none | N | No MAT locus described. Only an ambiguous HMG hit (G) |
| 8 | Pucciniales | 131 | not piloted | Luo et al. 2024 cereal rusts (genome); Ferrarezi 2022 contigs `OM2429xx` | N / Y (STE3/HD fragments only) | Dikaryotic assemblies carry BOTH alleles; HD chr4, PR chr9, STE3.2-1 chr1 |
| 9 | Trichosporonales | 127 | not piloted | Sun et al. 2019, 24 species (genome); *C. oleaginosum* ATCC 20508 | N | Fused HD+PR in all 24 species, one HD gene per mating type |
| 10 | Mortierellomycota | 100 | not piloted | none | N | No sex locus described; no tptA-HMG-rnhA synteny found (G) |
| 11 | Cantharellales | 99 | not piloted | none (Rhizoctonia STE3 mRNA `DQ926703.2` only) | N | - |
| 12 | Russulales | 92 | not piloted | *Heterobasidion parviporum* HR32 `KF280380.1` (MIP + a1/a2 + b1/b2) | Y | Two HD gene pairs at one locus; a2 pseudogene in some strains |
| 13 | Ascoideales | 80 | not piloted | *Ascoidea rubescens* a2 `XP_020044210.1` (round 1); *Saccharomycopsis schoenii* 2026 | N | FF1 in *A. rubescens*; AGE1/STE20 inside idiomorph in *A. asiatica*; PHN *S. capsularis*; expanded MAT in *S. schoenii* |
| 14 | Discinaceae | 72 | not piloted | *Gyromitra venenata* MAT1-1 `JBAMAM010000002.1` (COX13-APN2-MAT1-1-1-SLA2) | G | Mostly heterothallic (Dirks 2025) |
| 15 | Saccharomycodales | 71 | not piloted | none deposited | N | MAT is PRESENT; heterothallic (Krassowski 2019). A no-call is a miss, not loss |
| 16 | Orbiliales | 67 | pilot run | none; APN2/COX13 at `JH119078.1` (A. oligospora) with no MAT gene detected | G (flanks only) | Heterothallic (PCR, Zhou 2021) |
| 17 | Tilletiales | 67 | not piloted | none | N | No MAT locus paper found |
| 18 | Lipomycetales | 59 | not piloted | none deposited | N | Ancestrally PHC, all four genes plus HMGX; SLA2 ... APN2/MLH3 |
| 19 | Hymenochaetales | 55 | not piloted | none usable | N | - |
| 20 | Wallemiales | 51 | not piloted | *W. mellicola* CBS 633.66 putative locus (Padamsee 2012; genome) | N | Sxi1-type HD + HMG gene; two versions differ by inversion |
| 21 | Cystofilobasidiales | 42 | not piloted | *Phaffia rhodozyma* HD1/HD2 single genes `KU3157xx` | Y (single genes) | Homothallism evolved several times |
| 22 | Filobasidiales | 41 | not piloted | none | N | - |
| 23 | Microbotryales | 38 | not piloted | genome assemblies only; `JQ423666.1` is UNVERIFIED with no CDS | N | HD-PR linkage by repeated chromosome fusions; HD loss of function in 3 species |
| 24 | Taphrinales (beyond *T. deformans*) | 22 | 5/22 | *Protomyces lactucae-debilis* `MCFI01000029.1` Mc + Pi 7.5 kb apart | G | Fused M+P locus, like *T. deformans* |
| 25 | Auriculariales | 20 | not piloted | *Auricularia* HD single-gene deposits (`MN2670xx`, `MT8183xx`) | Y (single genes) | Tetrapolar |
| 26 | Basidiobolus | 18 | not piloted | none | N | Two unlinked HMG hits (G), not interpretable |
| 27 | Trigonopsidales | 16 | not piloted | none deposited | N | Heterothallic |
| 28 | Pezizaceae | 8 | not piloted | *Terfezia claveryi* T7 MAT1-1 (JGI Tercla1 scaffold_39:66,965-70,230, paper) | N (APN2 flank G) | Heterothallic; MAT1-1-1 not found by tblastn |
| 29 | *Saitoella* | 2 | 0/2 | weak Pi-like hit `KV454278.1`:167,097-167,699 | G (weak) | Unresolved |
| 30 | *Neolecta* | 1 | 0/1 | Mc `LXFE01003167.1`; Pi `LXFE01003595.1` | G | Different contigs; fragmented assembly |
| - | Glomeromycota | not given | - | *Rhizophagus irregularis* A1 MAT-like HD `KT946687.1` | Y | MAT-like; Lee 2024 argues it is not a mating locus |
| - | Chytridiomycota | not given | - | none | - | No MAT locus is known in any chytrid |

## Correction to an existing record (found during this round)

**`db/Ascomycota/Saccharomycetales/4959_cbs767_MTL_A` says "no MTLalpha1 gene is
present in this strain's assembly." Measurement contradicts this.**

- On `NC_006047.2` (D. hansenii CBS767 chr E), locus_tag DEHA2E19096g codes
  `XP_460134.1`, 212 aa, product "DEHA2E19096p". Its record carries a
  `MATalpha_HMGbox` region (pfam04769).
- Coordinates: complement 1,587,954-1,588,592. This is about 1 kb from MTLa1
  (DEHA2E19118g, 1,589,597-1,590,216) and MTLa2 (DEHA2E19140g, 1,590,508-1,591,344).
- tblastn (G): *C. albicans* alpha1 `AAD51411.1` hits 1,587,966-1,588,574
  (32% id, E=2e-23). *M. guilliermondii* alpha1 `ACS29273.1` hits the same
  place (E=3e-23).
- Krassowski et al. 2019 state that the *D. hansenii* MAT locus "contains
  neighboring MATa1, MATa2 and MATalpha1 genes" and class it as primary
  homothallic contiguous (PHC), with 5 other *Debaryomyces* species the same.
- PIK1/OBP/PAP1 are on the same chromosome but about 700 kb away
  (882,004-889,097). The MTL locus is not flanked by them in this species.
- Consequence: the record's `system: heterothallic` and the absence of alpha1
  are wrong. The alpha1 is real but unnamed in the annotation, so a name search
  misses it. Same failure shape as the annotation-gap lesson.

## Priority A

### 1. Serinales (2,368 genomes; 2/75 called)

**Papers**
- Hull CM, Johnson AD 1999. Science 285:1271-1275. PMID 10455055,
  doi:10.1126/science.285.5431.1271. *C. albicans* MTLa and MTLalpha.
- Butler G et al. 2009. Nature 459:657-662. PMID 19465905,
  doi:10.1038/nature08064. Eight *Candida* genomes. *L. elongisporus* has no MAT.
- Reedy JL et al. 2009. Curr Biol 19:891-899. PMID 19446455 (round 1).
- Lee DK, Hsiang T, Lachance MA 2018. Antonie van Leeuwenhoek 111:1935-1953.
  PMID 29651688, doi:10.1007/s10482-018-1084-y. *Metschnikowia* mating genomics.
  The *C. lusitaniae* MTL organisation is universal in, and exclusive to,
  Metschnikowiaceae.
- Krassowski T et al. 2019. Curr Biol 29:2555-2562. PMID 31353182,
  doi:10.1016/j.cub.2019.06.056 (PMC6692504, read in full for this round).
- Muñoz JF et al. 2018 (*C. auris*) was covered in round 1.

**Verified deposits (Y)**

| accession | bp | taxon / strain | idiomorph | CDS (protein, aa) |
|---|---|---|---|---|
| `AF167162.1` | 11,921 | *C. albicans* SC5314 | MTLa | PAPa c(2019..3686) AAD51407.1 555; OBPa AAD51406.1 433; PIKa AAD51405.1 956; MTLa2 join(8944..9213,9269..9604) AAT08437.1 201; MTLa1 join(9709..9900,9975..10346,10602..10670) AAD51404.1 210 |
| `AF167163.1` | 12,105 | *C. albicans* SC5314 | MTLalpha | alpha2 c(join(1934..2009,2069..2553)) AAD51408.1 186 (product "unknown", no gene name); OBPalpha AAD51409.1 429; PIKalpha AAD51410.1 977; MTLalpha1 7938..8519 AAD51411.1 193; PAPalpha c(8554..10230) AAD51412.1 558 |
| `AY622606.1` | 8,593 | *C. dubliniensis* d1558 | MTLalpha | alpha2 AAU13776.1 187; OBP AAU13914.1 429; PIK AAU13915.1 977; alpha1 AAU13777.1 193; PAP partial AAU13916.1 360. Pujol 2004, PMID 15302834 |
| `HQ696681.1` | 16,160 | *C. orthopsilosis* Cp90-125 | MTLa | GAP1 partial, PAPa, OBPa, PIKa, MTLa2 ADY62690.1 199, MTLa1 ADY62691.1 183, orf19.3202, RCY1. Sai 2011, PMID 21335529 |
| `HQ696679.1` / `HQ696680.1` / `HQ696682.1` | 8,814 / 9,920 / 9,848 | *C. orthopsilosis* J981224 / Cp289 / Cp185 | MTLa | same gene set |
| `HQ696678.1` | 14,971 | *C. metapsilosis* ATCC 96143 | MTLa | round 1 |
| `FJ524850.1` | 14,633 | *C. lusitaniae* CL143 | MTLa | PAP1 ACS29263.1, OBP1 ACS29264.1, PIK1 ACS29265.1, a2 ACS29266.1 176, a1 ACS29267.1 164, RCY1 ACS29268.1 |
| `FJ524851.1` | 14,353 | *M. guilliermondii* NRRL Y-2075 | MTLalpha | HIP1, OBP1 ACS29270.1, PIK1 ACS29271.1, alpha1 ACS29273.1 210, PAP1 ACS29269.1, RCY1 ACS29274.1 |

Note on `AF167163.1`: the alpha2 CDS has product "unknown" and no /gene
qualifier. Any curation step that looks for a gene name will miss it.

**Genome-derived loci (G)**

*C. tropicalis* MYA-3404 (Butler 2009 assembly). The assembly carries both
idiomorphs on separate scaffolds. Krassowski 2019 calls this a diploid assembly
of a heterothallic species.
- MTLa on `GG692408.1` (12,816 bp, whole scaffold). Annotated: PAPa EER30100.1
  c(754..2427), OBP EER30101.1 2826..4106, PIK1a EER30102.1 4383..7229, and a
  148-aa "hypothetical protein" EER30103.1 at 8487..8933 (RefSeq calls it Mtla1p).
- tblastn of *C. albicans* proteins:
  - MTLa2 `AAT08437.1`: two HSPs, 7,360-7,638 (aa 1-92) and 7,692-8,021 (aa 91-201),
    E=4e-41. **MTLa2 is not annotated.**
  - MTLa1 `AAD51404.1`: 8,179-8,367 (aa 1-66) and 8,421-8,912 (aa 64-206),
    E=2e-47. **The annotated a1 (EER30103.1) lacks the first exon.** It is a
    truncated model. Do not curate EER30103.1 as the a1 protein.
- MTLalpha on `GG692402.1` (802,573 bp): RCY1 EER30844.1 c(171190..173538),
  EER30845.1 173707..175893, PAPalpha EER30846.1 177475..179151, alpha1
  EER30847.1 c(179194..179775) (pfam04769), PIK1alpha EER30848.1
  c(180101..183022), OBPalpha EER30849.1 c(183279..184568).
  - alpha2 is not annotated. tblastn of `AAD51408.1` hits 185,474-186,082
    (45% id, E=6e-47). **MTLalpha2 is present but unannotated.**

*Clavispora lusitaniae* MTLalpha, strain P1, `CP038489.1` (chromosome 6; source
publication not checked):
- RCY1 389,970-392,426 (98% id to CL143 RCY1), PAP 396,680-398,344 (70% id),
  alpha1 398,370-398,837 (annotated QFZ29657.1, 155 aa, "putative mating-type
  protein ALPHA1"), PIK 399,025-401,784, OBP 402,445-403,521.
- No alpha2 hit. Consistent with alpha2 absence in Metschnikowiaceae.
- This is the same-species partner for the `FJ524850.1` MTLa record.

*Metschnikowia bicuspidata* NRRL YB-4993 MTLalpha, `NW_017387732.1`
(2,563,417 bp):
- RCY1 1,258,479-1,260,902; PAP 1,267,219-1,268,883; alpha1 1,268,939-1,269,355
  (E=8e-23 vs C. albicans alpha1); PIK 1,269,573-1,272,311; OBP 1,272,726-1,273,832.
- No alpha2 hit.

*Scheffersomyces stipitis* CBS 6054, `NC_009047.1` (chr 7):
- a2 702,296-702,913 (E=8e-22); PAP 703,131-704,759; alpha1 704,856-705,431
  (E=1e-20); PIK 705,688-708,576; OBP 708,916-710,205.
- One locus with a2 AND alpha1 beside PAP/PIK/OBP. No a1 or alpha2 hit at E<1e-5.
- Krassowski 2019: contiguous a1, a2, alpha1, PHC. My search found no a1. The a1
  may be too short or divergent for the *C. albicans* query. **a1 unresolved.**

*C. parapsilosis* FDAARGOS_653 `JABWAC010000006.1`: MTLa2 99% identical to
`AAY33181.1` at 201,579-202,233; PAP/OBP/PIK beside it. No alpha1/alpha2. NCBI
names the PAP protein "PAPalpha" (KAF6061291.1 etc.). That name is a transfer
artifact. Do not read it as evidence of an alpha idiomorph.

**Architecture flags (Serinales)**
- PAP, OBP, PIK are idiomorph-specific paralogs (a and alpha versions differ).
  A detector that uses them as "conserved flanks" will see allele-specific
  divergence.
- The two *C. albicans* idiomorphs have different gene order (queue memory).
- alpha2 is absent in Metschnikowiaceae (*Clavispora*, *Metschnikowia*,
  *C. auris*). alpha2 is present in *C. albicans*, *C. dubliniensis*,
  *C. tropicalis*.
- MTLa1 is a pseudogene in *C. parapsilosis* (Logue 2005).
- **NOMAT species**: *Lodderomyces elongisporus* and *Candida sojae* (Krassowski
  2019; Butler 2009). A no-call there is correct.
- **Homothallic contiguous**: *D. hansenii* and five other *Debaryomyces*,
  *S. stipitis*, *Spathaspora passalidarum* (a1, alpha1, alpha2), *Priceomyces*,
  *Cephaloascus fragrans*.
- **PHN (two non-allelic loci)**: *Wickerhamia fluorescens*, *Yamadazyma
  nakazawae*, *Y. philogaea* (MATa near a telomere, MATalpha at PIK/PAP/OBP).
- Diploid assemblies of heterothallic species show both idiomorphs:
  *C. tropicalis* MYA-3404, *C. albicans* hybrid assemblies. "Both present"
  is not homothallism in these.

### 2. Taphrinomycotina (Taphrinales 5/22, *Saitoella* 0/2, *Neolecta* 0/1)

**Literature**: No MAT analysis was found for *Protomyces*, *Saitoella*,
*Neolecta* or *Taphrina* species other than *T. deformans*. PubMed query
`(Neolecta OR Saitoella OR Protomyces OR Taphrina OR Taphrinomycotina OR
Archaeorhizomyces) AND (mating OR MAT OR homothallic OR heterothallic)` returned
5 hits: 2 on *Pneumocystis* (Hauser 2021 review), 1 *Novakomyces* taxonomy, 2 on
*T. deformans* (already the source of record 5011). Riley et al. 2016 PNAS
(PMID 27535936) is on Saccharomycotina and does not cover Taphrinomycotina.

**Genome-derived (G)**, queries = curated matMc/matPi/matMi/matPc from records
5011, 42068, 4754, 4896, 4899:

| genome | Mc best hit | Pi best hit | reading |
|---|---|---|---|
| *Protomyces lactucae-debilis* `GCA_002105105.1` | `MCFI01000029.1`:30,443-30,988, 41% id over 188 aa to *T. deformans* Mc, E=2e-39 | same contig 38,152-38,538, E=2e-18 | Mc and Pi 7.5 kb apart on one contig: fused locus like *T. deformans*. **Best Taphrinales-proximal candidate** |
| *Taphrina betulina* `GCA_008802775.1` | only a 72-aa HMG-box hit, `VWSI01000026.1`:59,897-60,106, E=4e-13 | `VWSI01000027.1`:16,268-16,750, 35% id over 164 aa, E=2e-26 | Pi convincing. No Mc next to it (search at E<1). Mc hit is probably a non-MAT HMG. Unresolved |
| *Neolecta irregularis* DAH-3 `GCA_001929475.1` | `LXFE01003167.1`:5-574 (contig edge), 42% id over 187 aa, E=6e-41 | `LXFE01003595.1`:11,513-12,163, 34% id over 233 aa, E=5e-35 | Both genes present, different contigs (assembly N50 15.7 kb). Linkage unknown |
| *Saitoella complicata* NRRL Y-17804 `GCA_001661265.1` | only 70-aa HMG-box hit, `KV454279.1`:205,219-205,428, E=6e-14 | `KV454278.1`:167,097-167,699, 31% id over 206 aa, E=1e-22 | Weak. Pi-like hit may be a non-MAT homeobox gene. Unresolved |

The pilot failure for Taphrinomycotina is consistent with these numbers: outside
*Protomyces* and *Neolecta*, the only full-length homology is to Pi.

### 3. Saccharomycotina orders outside Saccharomycetaceae

Source: Krassowski 2019 (PMC6692504), which classifies 332 genomes. I read the
per-family STAR Methods text. No new locus-specific GenBank deposits were found
for these orders. The facts below are from the paper text; coordinates are not
given in the text.

- **Saccharomycodales (Hanseniaspora, 71 genomes): MAT is NOT lost.**
  *H. uvarum*, *H. pseudoguilliermondii*, *H. valbyensis*, *H. osmophila* each
  contain only a MATa or only a MATalpha locus: heterothallic. *Kloeckera
  hatyaiensis*, *H. singularis*, *H. vineae*, *H. clermontiae* assemblies contain
  both, interpreted as diploid assemblies of heterothallic species. Steenwyk
  et al. 2019 (PMID 31112549, PMC6528967, full text checked) discusses loss of
  cell-cycle, DNA-repair and meiosis genes (for example IME1) but does not
  report MAT gene loss. **A no-call in Hanseniaspora is a detection miss.**
- **Phaffomycetales (188)**: *Cyberlindnera saturnus*, *C. mrakii*,
  *C. suaveolens*: FF1 flip/flop, 49-kb invertible region with MATa genes at one
  end and MATalpha at the other; SLA2 and DIC1 lie INSIDE the invertible region.
  NCBI contig named in the paper: `PPNR02000017.1` (DISCOVAR assembly with
  artifactually long IRs). *Starmera quercuum* FF1, contig `PPIB01000006.1`.
  *Wickerhamomyces canadensis* FF2 (MATa mid-scaffold, MATalpha on a 7-kb
  subtelomeric contig). *W. hampshirensis* PHC (all four genes together).
  *Wickerhamomyces* sp. NRRL YB-2243: SLA2-MATalpha1-MATalpha2-ORF-DIC1,
  heterothallic. *C. jadinii*, *C. maclurae*: diploid assemblies. Accessions
  named in the paper were not fetched by me: **unverified**.
- **Ascoideales (80)**: *Ascoidea rubescens* FF1: MATa1-MATa2 and
  MATalpha1-MATalpha2 separated by 44 kb of non-coding DNA, 2-kb IRs, one IR
  beside a telomere. *A. asiatica* heterothallic, MAT flanked by ADR1 and MHP1,
  with **AGE1 and STE20 integrated into the idiomorphs** (a: MATa2-AGE1a-STE20a-
  MATa1; alpha: AGE1alpha-MATalpha1-STE20alpha), no MATalpha2.
  *Saccharomycopsis capsularis* PHN (alpha genes between SLA2 and YJR098C; a
  genes beside an SLA2 pseudogene on a 27-kb contig). *S. malanga* has only
  MATa2 (HET). New: Kriti D et al. 2026, G3 16(5), PMID 41849659,
  doi:10.1093/g3journal/jkag067: chromosome-scale *S. schoenii* genome with an
  "expanded MAT system with multiple active copies that are co-transcribed". Not
  read in full. **Unverified coordinates.**
- **Dipodascales / Trichomonascaceae (398)**: *Blastobotrys proliferans*
  MATalpha1 between full-length SLA2 and APN2, syntenic with *B. adeninivorans*
  (Kunze 2014, Biotechnol Biofuels 7:66). **MATalpha2 is absent throughout
  Blastobotrys.** *B. proliferans* also has MATa2 between SLA2/APN2 pseudogenes
  at a subtelomeric site (PHN). *B. pratensis* lost its MATa genes.
  *Starmerella bombicola*: SLA2-ORF-TFC1 with no identifiable a1/a2; relatives
  *Candida apicola*, *Wickerhamiella domercqiae* have SLA2-MATalpha1-MFalpha-TFC1
  (a pheromone gene at MAT). *Nadsonia fulvescens* var. *elongata* PHN.
  *Geotrichum candidum* single-gene deposits `HF558448.1` / `HF558449.1` (round 1).
  *Sporopachydermia quercuum* PHC. No locus deposit found for *Magnusiomyces*,
  *Galactomyces*, *Sugiyamaella*.
- **Lipomycetales (59)**: ancestral site SLA2 ... APN2 and MLH3. *L. starkeyi*,
  *L. arxii*, *L. mesembrius*, *L. kononenkoae* PHC with all four genes adjacent
  plus an extra HMG gene **HMGX**. *L. japonicus*, *L. oligophaga*,
  *L. suomiensis* PHN. *L. doorenjongii* heterothallic diploid (the only HET
  genome in the genus).
- **Trigonopsidales (16)**: *Trigonopsis variabilis* only MATa, *T. vinaria*
  only MATalpha (HET). *Tortispora caseinolytica* only alpha, *T. starmeri*
  only a, *T. ganteri* diploid assembly. *Botryozyma nematodophila*: partial
  MATalpha2 downstream of SLA2.
- **Pichiales additions**: *Kregervanrija fluxuum* / *K. delftensis* FF1 with a
  12-kb invertible region holding only the four MAT genes. *Saturnispora* FF2.
  *Kuraishia molischiana* PHN. *Pachysolen tannophilus* FF1.
- **PHN pattern**: in 4 of 12 PHN species the second, subtelomeric MAT locus
  sits beside an SLA2 pseudogene. SLA2-based anchoring can land on the silent
  or secondary copy.

### 4. Pezizales and Orbiliales

**Pezizaceae — Andreu-Ardil L et al. 2026.** Mycorrhiza 36(3). PMID 42095936,
doi:10.1007/s00572-026-01266-3 (PMC13152902, read).
- *Terfezia claveryi*: heterothallic. Strains carry TcMAT1-1-1 (Tc1705, genome
  strain T7) or TcMAT1-2-1 (TcLlano). Six spores from one ascus segregated.
- MAT1-1-1 in JGI Tercla1 scaffold 39, coordinates 66,965-70,230 (paper).
  TcMAT1-1-1 protein 254 aa. TcMAT1-2-1 has five introns, three in the HMG box.
- Idiomorph amplicons 2,925 bp (MAT1-1) and 3,275 bp (MAT1-2) with primers in
  conserved flanks. The 3' flank includes an APN-type endonuclease.
- Data availability: "within the paper and its Supplementary Information". No
  nuccore deposit found (searched *Terfezia* and *Tirmania*).
- Also: *Tirmania nivea* G3 genome carries MAT1-1-1 only (Marqués-Gálvez 2021,
  cited in the paper).
- G check on `WHUX01000039.1` (T. claveryi T7 scaffold_39, 393,670 bp): APN2
  hits at 64,942-67,685 (Morchella APN2 AVI60822.1, 58-78% id, E=4e-73). No
  MAT1-1-1 hit at E<1e-3 with *Morchella* (AVI60816.1) or any curated *Tuber*
  MAT1-1-1; best weak hit 72,005-72,157 (E=0.7). **Even with a published
  coordinate, current Pezizales references do not find this MAT1-1-1 by
  tblastn.** The MAT gene itself would have to be extracted from the JGI gene
  model or supplement.

**Discinaceae — Dirks AC et al. 2025.** Mol Phylogenet Evol 205:108286.
PMID 39788220, doi:10.1016/j.ympev.2025.108286. Not in PMC; OSTI and
eScholarship pages gave only the abstract. Predominantly heterothallic; one
colocalized MAT1-1/MAT1-2 case (*G. esculenta* CBS101906) needing confirmation.
No deposits.

G measurement with *Morchella importuna* proteins (`KY782629.1`/`KY782630.1`)
and curated SLA2/APN2/COX13:
- *Gyromitra venenata* `GCA_040804095.1`, contig `JBAMAM010000002.1`:
  COX13 230,823-231,329 (E=2e-24); APN2 232,208-232,768 (74% id to Morchella,
  E=3e-101); **MAT1-1-1 236,696-237,604** (55% id over 119 aa to Morchella
  MAT1-1-1, E=2e-31); SLA2 257,551-260,618 (73% id, E=0). One contig with the
  full COX13-APN2-MAT1-1-1 ... SLA2 arrangement. MAT1-2-1 hits elsewhere are
  short HMG hits (E>=1e-8), not a MAT1-2 idiomorph. Reading: MAT1-1 strain.
- *Gyromitra* sp. 9 `GCA_040803905.1`: MAT1-1-1 `JBAMAC010000027.1`:
  484,441-484,926 (E=2e-62), SLA2 on the same contig 461,733-463,103; APN2 and
  COX13 on `JBAMAC010000057.1`. MAT1-1 strain.
- These are the first concrete Discinaceae coordinates. The MAT1-1-1 exon
  structure was not modelled.

**Morchella beyond Chai 2017**: no new locus deposit found this round. Round 1
records stand.

**Orbiliales (67 genomes).** No MAT locus paper beyond PCR typing (Zhou 2021,
PMID 34576814). No nuccore or protein MAT deposit found for Orbiliaceae. G
measurement with the curated Ascomycota MAT1-1-1/MAT1-2-1 set plus flanks:
- *Arthrobotrys oligospora* ATCC 24927 `GCA_000225545.1`: COX13
  `JH119078.1`:286,280-286,765 (E=2e-26), APN2 287,645-288,439 (E=6e-73), SLA2
  539,158-541,446 (E=0, 250 kb away). **No MAT1-1-1 or MAT1-2-1 hit within
  250-330 kb at E<0.03.** MAT1-2-1 hits elsewhere are 90-125-aa HMG hits
  (E~1e-9) on other scaffolds, not interpretable.
- *Dactylellina haptotyla* `GCA_031310005.1` (chromosome-level): APN2
  `CM062383.1`:1,302,583-1,303,260; COX13 1,304,646-1,304,918; SLA2
  1,088,118-1,089,209. Weak MAT1-2-1 HMG signal at 1,296,550-1,297,172
  (E=1e-4, several curated MAT1-2-1 queries), about 5 kb from APN2.
  **Weak candidate only.**
- *Drechslerella stenobrocha* `GCA_000525045.1`: only short HMG hits.
- Reading: the Orbiliales MAT genes are too divergent for current references.
  Flank anchors (APN2, COX13) are detectable. A reference would need a
  gene model built from RNA-seq or from the *D. haptotyla* candidate.

## Priority B

### 5. Basidiomycota

**Boletales (511).**
- `LC662829.1`, 121,710 bp, *Rhizopogon roseolus* Rrsb10, "A-mating type locus
  and its flanking region 1". Unpublished (Zhang W, Wan JN, Shimomura N,
  Yamaguchi T, Aimi T; title "Evolution of genomic structure of mating type locus
  in a bipolar basidiomycete, Rhizopogon roseolus from its tetrapolar
  ancestors"). CDS: pheromones phb3, phb7, phb2, phb1; pheromone receptors
  Rr-Rcb2 BDD37065.1 (453 aa) and Rr-Rcb1 BDD37067.1 (451 aa); RrA2-Hox2
  BDD37068.1 (576 aa); RrA2-Hox1 BDD37069.1 (651 aa); mip BDD37070.1; beta-fg;
  glycogenin; RPB2; others. **HD and PR genes are within 10 kb on one fragment.**
- `LC662830.1`, 20,246 bp, same strain, "similar with B-mating type region":
  phb6, Rr-Rcb3, Rr-Rcb4 (non-mating-type STE3 paralogs possible; not stated).
- Ke YH et al. 2023. Genetics 224(2):iyad069. PMID 37070772,
  doi:10.1093/genetics/iyad069. Suilloid HD MAT haplotypes. HD locus multiallelic
  with trans-specific polymorphism (*Suillus* + *Rhizopogon*). Deposits
  `ON315855.1`-`ON315867.1` (4.3-4.9 kb, HD1 + HD2 only, no flanks);
  e.g. `ON315855.1` *S. luteus* HD1 WDY60839.1 660 aa, HD2 WDY60840.1 574 aa.
  `ON315861.1` is flagged UNVERIFIED.
- *Serpula lacrymans*, *Coniophora*: no MAT locus deposit found.
- Flag: *R. roseolus* is bipolar by HD-PR linkage (per deposit title and gene
  layout). Suilloids otherwise treated as tetrapolar with multiallelic HD.

**Polyporales (198).**
- James TY, Lee M, van Diepen LTA 2011. Eukaryot Cell 10:249-261. PMID 21131435.
  *Phanerochaete chrysosporium* is bipolar; the single MAT locus is HD genes.
  `HQ188438.1`-`HQ188443.1` (7.7-18.6 kb). `HQ188442.1` (18,639 bp,
  ME-OC-11_c7): MIP partial ADN97186.1, A2 (HD2) ADN97187.1 487 aa, A1 (HD1)
  ADN97188.1 632 aa, HP2.
- `LC730474.1`, 70,496 bp, *Mycoleptodonoides aitchisonii* 50005-7 A locus
  (unpublished, Aimi lab). glydh, sec61, up8, up2, up11, mip BDS00035.1, Mahd2-7
  BDS00036.1 (569 aa), Mahd1-7-like BDS00037.1 ("HD1-like protein without
  homeodomain", 619 aa), beta-fg, glycogenin, RPB2. `LC532152.1` (255,530 bp) is
  a second A-region from the same species.
- Flag: `LC730474.1` HD1 lacks a homeodomain. A domain-required HD1 check will
  reject it.
- *Trametes*, *Ganoderma*: only single-gene deposits (Ganoderma boninense HD1/HD2
  mRNA `ON855036.1`/`ON855037.1`).

**Russulales (92).**
- van Diepen LTA et al. 2013. Mol Biol Evol 30:2286-2301. PMID 23864721.
  *Heterobasidion*. `KF280366.1`-`KF280390.1` (about 8.9-12.4 kb) across
  *H. annosum*, *H. parviporum*, *H. abietinum*, *H. araucariae*.
  `KF280380.1` (11,547 bp, *H. parviporum* HR32): MIP partial AGS09343.1, a1 HD
  AGS09344.1 653 aa, a2 AGS09345.1 609 aa, b1 HD AGS09346.1 643 aa, b2 AGS09347.1
  605 aa, beta-fg fragment.
- Flags: two HD gene pairs (a and b) at the A locus; a2 is a pseudogene in some
  strains (e.g. `KF280389.1`, `KF280385.1` titles); extensive trans-specific
  polymorphism, so alleles are shared across species.
- *Auriscalpium orientale* `MF4405xx`: partial homeodomain genes only.

**Cantharellales (99).** Only *Rhizoctonia solani* STE3-like mRNA `DQ926703.2`
(1,468 bp) and *Thanatephorus* partial ste3 `AY226017.1`. No locus found.

**Hymenochaetales (55).** Only `ON693462.1`/`ON693465.1` (*Inonotus obliquus*,
UNVERIFIED, no annotation) and a partial ste3 (`AY226013.1`). Nothing usable.

**Auriculariales (20).** Single-gene HD and STE3 deposits: *A.
auricula-judae* `MN267026.1`-`MN267033.1` (HD1/HD2 alleles, isolates 14-5 and
18-119), *A. cornea* `MT818361.1` (HD protein 2, 1,270-aa CDS, unpublished),
`MT040090.1`-`MT040093.1`, *A. heimuer* `MN442080.1`. No flanks.

**Sporidiobolales (253).**
- Coelho MA et al. 2008. Eukaryot Cell. PMID 18408057,
  doi:10.1128/EC.00025-08. *R. toruloides* MAT genes. Deposits are short:
  `EU386160.1` (3,420 bp, STE20-A1), `EU386161.1` (3,560 bp, STE20-A2 + RHA2
  pheromone).
- Coelho MA et al. 2010. PLoS Genet 6:e1001052. PMID 20700437,
  doi:10.1371/journal.pgen.1001052. *Sporidiobolus salmonicolor*: a deviation
  from the bipolar-tetrapolar paradigm (PR and HD regions linked in one large
  MAT region; HD multiallelic). The 262 linked nuccore records are all <=2.4 kb
  (STE3, STE20, RibL6 fragments).
- R. toruloides genome scaffolds exist (e.g. CECT1137 `LK052958.1`,
  `LK052966.1`), but they have no CDS annotation. G coordinates not derived.
- **No usable locus record yet.** Building one needs genome coordinates.

**Microbotryales (38).**
- Petit E et al. 2012. Evolution. PMID 23106715. Deposits are short loci
  (`JQ4236xx`, ~600 bp). `JQ423666.1` (14,887 bp, XRN-HD2-HD1) and
  `JQ423663.1`/`JQ423661.1` are flagged UNVERIFIED and carry no CDS. Not usable.
- Lucotte EA et al. 2025. Nat Commun. PMID 40436846,
  doi:10.1038/s41467-025-60222-5. In *M. superbum*, *M. shykoffianum*,
  *M. scorzonerae* the HD genes lost mating function; strains can be homozygous
  for a disrupted HD2. Other *Microbotryum* species link HD and PR by chromosome
  fusions.
- Genome assemblies of *M. lychnidis-dioicae* a1/a2 mating-type chromosomes
  exist (Branco et al. 2017/2018; Hood 2013 PMID 23150606). I did not fetch or
  verify their accessions. **N.**
- Flag: HD pseudogenes are expected. A missing or disrupted HD is not a failure
  in these species.

**Pucciniales (131).**
- Luo Z, McTaggart A, Schwessinger B 2024. PLoS Genet 20:e1011207.
  PMID 38498573, doi:10.1371/journal.pgen.1011207. *P. coronata* f. sp.
  *avenae*, *P. graminis* f. sp. *tritici*, *P. triticina*, *P. striiformis*
  f. sp. *tritici*. HD (bW-HD1, bE-HD2) on chromosome 4; Pra (STE3.2-2,
  STE3.2-3) on chromosome 9; STE3.2-1 on chromosome 1 (non-MAT). Tetrapolar. No
  nuccore links; coordinates are in haplotype-phased assemblies (not fetched).
- Ferrarezi JA et al. 2022. PLoS Pathog 18:e1010439. PMID 35617196. Large
  contigs with partial STE3 or HD1 CDS across many rust genera, e.g.
  `OM242977.1` (72,703 bp, *Puccinia paullula*, STE3.2.3 partial USF89498.1) and
  `OM242968.1` (5,118 bp, *Sphaerophragmium* HD1 partial USF89489.1). Partial
  CDS only.
- Holden S et al. 2023 (PMID 37880702): *P. striiformis* MAT alleles (not read).
- Flag: rust assemblies are dikaryotic. Phased assemblies carry both HD alleles
  and both Pra alleles. Expect two alleles per genome; that is not homothallism.

**Trichosporonales (127).**
- Sun S, Coelho MA, Heitman J, Nowrousian M 2019. PLoS Genet 15:e1008365.
  PMID 31490920, doi:10.1371/journal.pgen.1008365. 24 species. MAT loci fused
  (HD + PR) in all, one HD gene per mating type, highly conserved gene order,
  no extended recombination suppression. Linked nuccore records are the
  *Cutaneotrichosporon oleaginosum* ATCC 20508 assembly (`NW_027072170.1`-
  `NW_027072177.1`) and its mRNAs. Coordinates not extracted. **N.**
- `KM821408.1` (1,928 bp) *C. oleaginosum* ATCC 20509 STE3 only (Kourist 2015).

**Tremellomycetes outside Tremellales.**
- Cystofilobasidiales (42): David-Palma M et al. 2016 PLoS Genet, PMID 27327578
  (*Phaffia rhodozyma* primary homothallic; deposits `KU315762.1`-`KU3157xx`
  are single HD1/HD2 genes, 1.2-1.9 kb). Cabrita A et al. 2021 mBio 12:e03130-20,
  PMID 33593979, doi:10.1128/mBio.03130-20: homothallism arose more than once in
  the order (*Cystofilobasidium*); possible Hd2 homodimer replacing Hd1/Hd2.
  *Cystofilobasidium* STE3 partials `MT5613xx`. No locus deposit.
- Filobasidiales (41): nothing found.

**Tilletiales (67).** No MAT locus paper or deposit found (PubMed and nuccore).

**Wallemiales (51).** Putative MAT locus in *W. mellicola* CBS 633.66
(Padamsee et al. 2012 genome; described in Sun X et al. 2019 Genes 10:427,
PMID 31167502, and Gostinčar C et al. 2019 Front Microbiol, PMID 31551960):
a Sxi1-like HD gene plus an HMG gene; two versions across 25 strains that
differ in some genes and in locus orientation (inversion). Genome-only. **N.**

**Cystobasidiomycetes.** Nothing found.

### 6. Early-diverging lineages

**Mortierellomycota (100).** No sex-locus paper found. NCBI "mating type"
proteins (e.g. GJJ68563.1 *Entomortierella parvispora* "mating-type protein
A1", OAQ36521.1 *Linnemannia elongata* "mating type protein 2, partial") are
pipeline names by similarity, not characterized genes.
G check with the curated Mucorales sexM/sexP/tptA/rnhA/glrA/algA/btbA set:
- *Linnemannia elongata* `GCA_036320915.1`: one HMG gene hit by both sexM and
  sexP queries at `JAXBDG010000013.1`:1,183,306-1,183,545 (80 aa, E=4e-16).
  tptA on contig 8, rnhA on contig 14. No tptA-sex-rnhA synteny.
- *Mortierella alpina* `GCA_977091265.1`: one HMG gene at `CDRYGW010000012.1`:
  3,024,605-3,024,910 (E=8e-16); rnhA on the same contig at 715-716 kb (2.3 Mb
  away); tptA on another contig.
- Reading: no Mucorales-type sex locus is detectable. An 80-aa HMG hit is not
  a sex gene call.

**Kickxellomycota (171).** No sex-locus paper. *Smittium culicis* proteins
named "Silenced mating-type M-specific polypeptide Mc" (OMJ27758.1 etc.) are
pipeline names. *Coemansia reversa* `GCA_002705745.1`: one HMG gene hit by sexM
and sexP at `KZ303509.1`:21,215-21,484 (E=2e-15); glrA and rnhA on other
scaffolds. Not interpretable as a sex locus.

**Basidiobolus (18).** No sex-locus paper. *B. meristosporus* CBS 931.73
`GCA_002104905.1`: sexM best hit `MCFE01000070.1`:23,133-23,378 (E=3e-19);
sexP best hit a different gene `MCFE01000238.1`:97,935-98,162 (E=2e-12); tptA,
rnhA, glrA, algA on other contigs. Not interpretable.

**Glomeromycota.**
- Ropars J et al. 2016. Nat Microbiol 1:16033. PMID 27572831. MAT-like HD locus
  in *Rhizophagus irregularis*. `KT946687.1` (8,001 bp, isolate A1, scaffold
  816): phosphoglycerate mutase family AMM63104.1, hd2 AMM63105.1 509 aa,
  hd1-like AMM63106.1 315 aa, hypothetical AMM63107.1. Also `KT946688.1`-
  `KT946690.1`, `KT954979.1`, `KT962968.1`, `KT962969.1` (other isolates).
- Lee SJ et al. 2024. BMC Genomics 25:888. PMID 39304834,
  doi:10.1186/s12864-024-10770-9: high diversity at this locus, congruent with
  genome-wide loci; they reject a mating role. Riley 2014 (PMID 24033097)
  reports a greatly expanded MATA-HMG family.
- Flag: if used, label "MAT-like HD"; the HMG family expansion makes HMG hits
  uninformative.

**Chytridiomycota.** No MAT locus is known. PubMed search for chytrid or
Blastocladiomycota "mating type", "MAT locus", "sex locus" returned no locus
characterization. A no-call is the expected result.

## Detection-relevant flags added in round 2

1. **Unannotated or truncated core genes in reference assemblies.**
   *C. tropicalis* MTLa2 and MTLalpha2 unannotated; *C. tropicalis* a1 model
   missing exon 1; *D. hansenii* alpha1 annotated only as "DEHA2E19096p";
   *C. albicans* `AF167163.1` alpha2 CDS has no gene name.
2. **True MAT absence**: *L. elongisporus*, *C. sojae* (Serinales). Expected no
   call.
3. **Hanseniaspora retains MAT.** A no-call there is a miss.
4. **Idiomorph-specific "flanks"**: PAP/OBP/PIK (CTG clade), AGE1/STE20
   (*Ascoidea asiatica*). Allelic divergence in these genes is biology.
5. **Two HD pairs per A locus** (*Heterobasidion*); **HD1 without homeodomain**
   (*Mycoleptodonoides*); **HD loss of function** (three *Microbotryum*
   species); **dikaryotic assemblies with both alleles** (rusts).
6. **Linked HD-PR (bipolar)**: *Rhizopogon roseolus*, *Phanerochaete* (HD-only
   bipolar), Trichosporonales (all fused), *Sporidiobolus* (linked, rare
   recombination), some *Microbotryum*.
7. **Tblastn reach**: Orbiliales and *Terfezia* MAT1-1-1 are not found by any
   curated Pezizomycotina protein, although flanks are found. Adding more
   distant Pezizomycotina references will not fix this; lineage-own gene models
   are needed.
