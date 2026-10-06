# Is XP_460134.1 (D. hansenii CBS767) an MTLalpha1 gene?

Question from the curator, 2026-09-24. The curated record `4959_cbs767_MTL_A`
states CBS767 "carries only the A gene configuration (MTLA1, MTLA2); no
MTLalpha1 gene". XP_460134.1 is annotated only as "DEHA2E19096p".

Everything below was measured in this directory on 2026-09-24.

## 1. It carries the MATalpha1 domain, and it is the only protein that does

* NCBI's own protein record annotates residues 103->157 as Region
  `MATalpha_HMGbox`, "Mating-type protein MAT alpha 1 HMG-box; pfam04769".
* `hmmsearch --cut_ga` with Pfam PF04769 (MATalpha_HMGbox) over whole RefSeq
  proteomes (`hmm_*.tbl`). The alpha box occurs only in MATalpha1 proteins:

  | proteome | PF04769 proteins |
  |---|---|
  | **D. hansenii CBS767** | **1: XP_460134.1 (E=2.2e-14)** |
  | C. lusitaniae, C. tropicalis, C. dubliniensis | 1 each |
  | K. lactis, L. thermotolerans | 1 each |
  | S. cerevisiae S288C | 2: NP_009867.1, NP_009869.3 = HMLalpha1 and MATalpha1 |
  | C. albicans SC5314 (reference haplotype), C. parapsilosis, M. guilliermondii, S. stipitis | 0 |

  So a paralog is not an available explanation: D. hansenii has one alpha-box
  protein genome-wide.

## 2. Reciprocal best hits

* XP_460134.1 blastp against 11 yeast proteomes: every hit with E < 1e-3 is a
  PF04769 protein (C. tropicalis 7.5e-28, C. dubliniensis 4.6e-25,
  C. lusitaniae 6.7e-24, K. lactis 3.1e-8, S. cerevisiae 3.2e-5, L.
  thermotolerans 1.1e-4). Nothing else.
* C. albicans MTLalpha1 (AAD51411.1, from the Hull & Johnson locus deposit
  AF167163.1) blastp against the whole D. hansenii proteome: top hit
  **XP_460134.1, 32% identity over 209 aa, E=2.8e-25**. The next hit is
  E=0.55. One clear ortholog.

## 3. Gene tree

MAFFT L-INS-i, ClipKIT kpic-gappy, IQ-TREE 3 (Q.YEAST+I+G4, 1,000 UFBoot),
`alpha1.treefile`. XP_460134.1 falls with the CTG-clade MTLalpha1 proteins
(C. lusitaniae, C. tropicalis, C. albicans, C. dubliniensis), apart from the
Saccharomycetaceae MATalpha1 proteins and the Pezizomycotina MAT1-1-1
outgroup (A. fumigatus, N. crassa). That matches the species tree, as an
ortholog should. Rooted on the Pezizomycotina outgroup, the five CTG-clade
proteins form one group (UFBoot 82) and XP_460134.1 sits with the Candida
MTLalpha1 proteins inside it (93). Support across the tree ranges 48-100 on a
12-sequence alignment of a short, fast-evolving domain; the tree is
corroboration, not the main line.

## 4. It sits in the MTL locus

RefSeq GFF, NC_006047.2:

    1,586,741-1,587,673  +  pirin family protein
    1,587,954-1,588,592  -  XP_460134.1  (alpha box)
    1,590,115-...        -  Mtla1p  (XP_460135.2)
    1,590,991-1,591,344  -  Mtla2p  (XP_460137.2)
    1,592,930-1,598,572  +  myosin 1

alpha1 lies ~1.5 kb from MTLa1, in the same locus. The PAP/OBP/PIK genes that
sit inside the C. albicans idiomorph are not in this window.

## Conclusion

Four independent lines -- the MAT-specific domain (unique in the proteome),
reciprocal best hits with C. albicans MTLalpha1, orthologous placement in the
tree, and position inside the MTL locus next to MTLa1/a2 -- say XP_460134.1 is
the D. hansenii MTLalpha1 gene. CBS767 therefore carries a1, a2 and alpha1 at
one locus. This agrees with Krassowski et al. 2019 (homothallic, contiguous
a+alpha genes). What is NOT shown here is whether the gene is functional or
whether CBS767 is homothallic in practice; that is the literature's claim, not
a measurement from this analysis.
