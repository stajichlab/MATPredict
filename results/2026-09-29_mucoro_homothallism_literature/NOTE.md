# MAT locus structure in homothallic Mucoromycotina: literature note

Date: 2026-09-29. Read-only literature review. No MATPredict data was analysed.
Trigger: 36 genomes have a sexP call and a sexM call on different contigs. The
current `homothallic_candidate` label fires only when both genes sit at one locus.

## Sources read

- Book: Wendland J (ed.) 2016, *The Mycota I: Growth, Differentiation and
  Sexuality*, 3rd ed. (resource/978-3-319-25844-7.epub).
  - Ch. 16, Peraza-Reyes & Malagnac, "Sexual Development in Fungi",
    Sect. II ("Sexual Identity and Mate Recognition"), subsect. 1 (zygomycete
    pheromones). doi:10.1007/978-3-319-25844-7_16
  - Ch. 10, Wöstemeyer et al., "Pheromone Action in ... Zygomycetes ...",
    Sect. II.B.5 "Genetic Control". doi:10.1007/978-3-319-25844-7_10
  - The EPUB has no page numbers, so section numbers are given instead.
- Idnurm A 2011. Sex determination in the first-described sexual fungus.
  Eukaryot Cell 10:1485-91. PMID 21908600, doi:10.1128/EC.05149-11. Full text read.
- Gryganskyi AP et al. 2018. Phylogenetic and phylogenomic definition of
  *Rhizopus* species. G3 8:2007-18. PMID 29674435, doi:10.1534/g3.118.200235.
  Full text read (mating-locus section and Fig. 5 legend).
- Lee SC, Idnurm A 2017. Fungal sex: the Mucoromycota. Microbiol Spectr 5(2).
  PMID 28332467, doi:10.1128/microbiolspec.FUNK-0041-2017. Full text read.
- Lee SC, Heitman J 2014. Sex in the Mucoralean fungi. Mycoses 57 Suppl 3:18-24.
  PMID 25175551, doi:10.1111/myc.12244. Full text scanned.
- Gryganskyi AP et al. 2010. Mating locus in the *R. oryzae* complex.
  PLoS One 5:e15273. PMID 21151560, doi:10.1371/journal.pone.0015273. Abstract.
- Idnurm A et al. 2008. Nature 451:193-6. PMID 18185588, doi:10.1038/nature06453. Abstract.
- Mehta BJ, Cerdá-Olmedo E 2001. Intersexual partial diploids of *Phycomyces*.
  Genetics 158:635-41. PMID 11404328, doi:10.1093/genetics/158.2.635. Abstract.

Sources I could NOT read (these are gaps):
- Schulz E, Wetzel J, Burmester A, Ellenberger S, Siegmund L, Wöstemeyer J 2016.
  "Sex loci of homothallic and heterothallic Mucorales." Endocytobiosis Cell Res
  27:39-57. Not in PubMed. No full text found. It is likely the key paper for
  *Zygorhynchus*/*Mucor* homothallics. A web snippet says homothallic Mucorales
  carry both sexM and sexP. I could not see the locus arrangement. **Get this paper.**
- Schulz E, Wetzel J 2016. Sex-deficient mutants of homothallic *Zygorhynchus
  moelleri*. Mycoscience. doi:10.1016/j.myc.2016.05.002. Paywalled; no abstract.
- Gauger W 1990. The genetics of two homothallic species of the Mucoraceae
  (*M. genevensis*, *Z. exponens*). Sex Plant Reprod. doi:10.1007/BF00189949.
  Paywalled; no abstract.

## What the sources show

### The book (2016)
- Ch. 16 says only that all studied heterothallic zygomycetes have one bipolar MAT
  locus. It says homothallic fungi often carry both idiomorphs in one genome. It
  gives this as a general fungal statement, mostly from ascomycetes. It has no
  data on homothallic Mucorales loci.
- Ch. 10, Sect. II.B.5 lists sexM/sexP homologues in *Rhizopus oryzae*, the
  homothallic *Syzygites megalocarpus* (citing Idnurm 2011), *Mucor mucedo* and
  *M. circinelloides*. It does not describe locus arrangement in homothallics.
- Ch. 10 mentions *Zygorhynchus moelleri* only for polyamine work. It gives no
  sex-locus data for it.

### *Syzygites megalocarpus* (Idnurm 2011): TWO SEPARATE LOCI
- Each strain has two sex loci, one with sexM and one with sexP. The two regions
  were 20.8 kb and 25.6 kb and were assembled separately by inverse PCR.
- Both loci carry copies of the same flanking genes: rnhA (RNA helicase) and glrA
  (glutathione oxidoreductase). So the flank block is duplicated.
- One copy of each duplicated flank gene is a pseudogene. The rnhA next to sexP
  has an inversion and a deletion. The glrA next to sexM has a deletion.
- Transposon remnants sit next to sexP.
- tptA could not be amplified from *Syzygites*. The arbA (BTB) gene sits next to
  sexM. In *Rhizopus* arbA sits inside the sexP allele.
- Idnurm's model: a segmental translocation from a heterothallic, *Rhizopus*-like
  ancestor. This moved one allele to a second position, on another chromosome or
  far away on the same chromosome. The paper does not show which. There is no
  genome assembly or karyotype.
- Caveat in the paper: spores have 20+ nuclei, and no homokaryon was made. So it
  is not proven that both loci are in one nucleus.
- Zygospore progeny were all self-fertile (24 of 24 from CBS 108947).

### *Rhizopus microsporus* var. *azygosporus* CBS 357.93 (Gryganskyi 2018)
- The genome (PJQM00000000) has "two sex loci", one with sexM and one with sexP
  (GenBank MG967659-60).
- The authors say the assembly is poor: low coverage, 16 Mb.
- The authors state that it is not clear whether CBS 357.93 is a true homothallic
  or "an unreduced fusion event" between (+) and (-) strains.
- In *R. microsporus* (+) and (-), the flank is glrA, not tptA. This matches
  *Syzygites*.
- The paper says the tptA-sex-rnhA triplet "is not universally conserved among
  Mucorales". The *R. arrhizus*/*delemar*/*stolonifer* loci have arbA on the
  side opposite rnhA.
- The paper does not state whether the two *azygosporus* loci are on one contig or
  on different contigs. The Fig. 5 legend shows two separate locus diagrams.

### Other taxa named by the curator
- *Zygorhynchus* (now in *Mucor*): I found no locus-structure data in accessible
  sources. The Schulz 2016 papers probably have it.
- *Mucor genevensis*: I found no locus-structure data. Gauger 1990 did genetics
  on it, but I could not read it.
- *Rhizopus homothallicus*: Idnurm 2011 names it as a target for future work. I
  found no published sex-locus description. Recent papers on it are clinical
  case reports.
- Lee & Idnurm 2017 state that the genes behind homothallism in Mucorales other
  than *Syzygites* "have not been characterized". I found no 2018-2026 PubMed
  paper that describes a homothallic Mucorales MAT locus. Gryganskyi 2018 is the
  only one, and its case is ambiguous.

### Mimics: diploids, heterokaryons, fusions
- *Phycomyces*: some cross progeny are diploids or partial diploids that are
  heterozygous for sex. Chromosome loss makes heterokaryons and sectors with
  mixed mating behaviour (Mehta & Cerdá-Olmedo 2001).
- Blakeslee's "homothallic" *Phycomyces* progeny were mitotically unstable. They
  reverted to single-sex strains (discussed in Idnurm 2011).
- Mucorales mycelia are coenocytic and multinucleate. A culture started from a
  multinucleate spore or hyphal tip can carry nuclei of both mating types.
- Mucorales have ancient genome duplications (Corrochano 2016, PMID 27238284).
  Gryganskyi 2018 reports a threefold range in genome size within *R. microsporus*.

### *Umbelopsis*
- None of the sources I read describe sex, zygospores, or a MAT locus in
  *Umbelopsis*. I did not verify the common statement that its zygospores are
  unknown. I found no evidence of a homothallic *Umbelopsis*.

## Answers

1. **Arrangement.** The one well-studied case, *Syzygites*, is (b): two separate
   loci. It is not a tandem or fused locus. *R. azygosporus* also has two loci,
   but its origin (homothallism or fusion) is unresolved. For *Zygorhynchus*,
   *M. genevensis* and *R. homothallicus* there are no accessible data. So (c),
   "variable by species", is neither supported nor excluded. No source reports a
   homothallic Mucorales with sexP and sexM in tandem at one locus.
2. **Synteny.** In *Syzygites* both loci keep an rnhA + glrA flank, with one copy
   of each gene a pseudogene. tptA was not found. In the *R. microsporus* clade,
   glrA replaces tptA. So a strict tptA-sex-rnhA check will fail in exactly these
   lineages. rnhA (and glrA) is the more conserved anchor in the homothallics
   studied.
3. **Umbelopsis and mimics.** There is no evidence for a homothallic *Umbelopsis*.
   There is direct evidence (in *Phycomyces*) that sex-heterozygous diploids and
   heterokaryons exist and are unstable. *R. azygosporus* CBS 357.93 is the
   published case where the authors could not separate homothallism from a (+)/(-)
   fusion.

## Recommendation for detection logic (inference, not from sources)

The facts above suggest that a one-locus rule is the wrong model for Mucorales.
Proposed tiers when sexP and sexM calls sit on different contigs:

Flag as `homothallic_candidate_unlinked` when all of these hold:
- Both calls are high confidence and classifier-typed as sexP and sexM. The two
  proteins are not near-identical, so this is not a split or duplicated single gene.
- Each call has its own flank context: an rnhA hit (or glrA/arbA/tptA) on the same
  contig within the Mucoromycota cluster gap. Do not require tptA.
- Optional support: one of the two flank copies looks degraded (a truncated or
  frameshifted rnhA/glrA). This is the *Syzygites* signature.
- Optional support: the species is described as homothallic. Use this only as a
  prior. Report it as metadata, not as evidence.

Flag as `mixed_or_heterokaryotic_suspect` instead when any of these hold:
- Genome-wide heterozygosity or read-depth evidence of two haplotypes: allele
  ratios near 0.5 at many single-copy genes, or two copies of most BUSCOs.
- Two full, intact, co-linear copies of the whole flank block (tptA/rnhA/glrA)
  with high identity outside the idiomorph, and no pseudogenes. This looks like a
  (+) haplotype and a (-) haplotype side by side.
- sexP and sexM contigs differ in depth (a minority nucleus) or in GC or taxonomy
  (contamination). Check the taxonomic best hit of each contig.
- Species is known to be heterothallic (for example *R. arrhizus*,
  *M. circinelloides*, *Phycomyces*).
- More than one sexP or more than one sexM copy.

Keep the existing same-locus label for true tandem or fused cases.

Open item: none of these thresholds has been tested on MATPredict data.
Testing needs the 36 genomes and ideally reads for depth and heterozygosity.
The flank-degradation test is based on one species (*Syzygites*). It may not
generalise.

---

## ADDENDUM 2026-09-29: Schulz et al. 2016 (read in full)

Source: Schulz E, Wetzel J, Burmester A, Ellenberger S, Siegmund L, Wöstemeyer J
2016. "Sex loci of homothallic and heterothallic Mucorales." Endocytobiosis Cell
Res 27(4):39-57. Local copy: resource/Schulzetal.2016.pdf. Text extracted with
pypdf. This closes the "Get this paper" gap above. Page numbers are journal pages.

### Scope correction
The paper does NOT cover *Mucor genevensis*, *Rhizopus homothallicus* or
*R. azygosporus*. Its homothallic cases are *Syzygites megalocarpus* (from Idnurm
2011), *Zygorhynchus heterogamus* NRRL 1489 (JGI), *Mycotypha africana* NRRL 2978
(JGI) and a partial *Z. moelleri* FSU 531 = CBS 140413 locus (GenBank KX966017).
Heterothallic comparisons: 15+ species (Table 1, p. 42-43; Fig. 2, p. 44).

### Per species
- ***Z. heterogamus*: ONE LOCUS.** sexM and sexP are 5.3 kb apart (p. 46).
  The flank block tptA + algL is duplicated. The two tptA copies are 74% nt
  identical; the two algL copies are 76% nt identical (p. 50; Fig. 6-7). Both
  copies of each look intact; they differ in N/C-terminal length and indels.
  Only one rnhA, intact. Transposable elements flank sexP (Fig. 5). Authors:
  "the sex locus of Z. heterogamus can be interpreted as a chimera from sexM and
  sexP loci" (p. 50). glrA is also near the locus (p. 44).
- ***M. africana*: SAME SCAFFOLD, ~150 kb APART** (p. 46). One rnhA, intact
  (p. 50). glrA near the locus (p. 44).
- ***S. megalocarpus*: "most likely" different chromosomes** (p. 46). This is
  inference from Idnurm 2011, not new data. Two rnhA; the sexP-side copy is a
  pseudogene (inversion + ~1.4 kb deletion, p. 50). The paper contradicts itself
  on which glrA is the pseudogene (p. 43 says sexM side; p. 44 says sexP side).
  Idnurm 2011 is the primary source and says sexM side.
- ***Z. moelleri*:** only sexP + adjacent rnhA found (p. 51). sexM not found.
  The authors predict a *Z. heterogamus*-like locus. This is not tested. (The text
  calls it "heterothallic" once on p. 51; the title and next paragraph say
  homothallic. Treat that as a typo.)

### Heterothallic synteny breaks (relevant to flank rules)
- Lichtheimiaceae: sexM is on a different chromosome from rnhA/algL/tptA (p. 43).
- *A. glauca*, *A. repens*, *B. circina*: rnhA is on a different scaffold from the
  sex gene (p. 43).
- So "no flank next to a sex gene" also occurs in heterothallic genomes.

### Homothallism vs heterokaryosis
- No single-spore, nuclear-state or ploidy data for any homothallic genome here.
- *A. glauca* protoplast fusions of (+) and (-) gave homothallic phenotypes.
  Stable heterokaryons were rare. Most derivatives came from nuclear fusion, then
  chromosome loss (aneuploids with both sex genes in one nucleus) (p. 47). So an
  artificial or natural fusion can look like a homothallic genome.

### Sequence divergence
- HMG domains of SexM and SexP from all species cluster by type (Fig. 3-4).
- SexM looks polyphyletic and species-specific. SexP is more conserved (p. 45).
- Only number for a homothallic: *Z. moelleri* SexP is 62.1% similar to
  *Z. heterogamus* SexP and 56% to *M. mucedo* (p. 51).
- *Umbelopsis rammaniana* HMG protein was ambiguous by BLAST (31% to both). It
  was typed sexP by tree position and a central HMG domain (p. 45).
- Inference: SexP typing should work. SexM typing in homothallics is a higher
  risk. Not tested on MATPredict data.

### Revised criteria (inference; untested)
1. Three arrangements are real, not one. Report them separately:
   - `homothallic_candidate` (same locus, <= cluster gap): *Z. heterogamus*.
   - NEW `homothallic_candidate_same_contig_distant` (same contig, > 50 kb):
     *M. africana* (~150 kb). The current rules call this neither tier.
   - `homothallic_candidate_unlinked` (different contigs): *Syzygites*.
2. Flank rule: accept rnhA, glrA, algL (algA, pfam05426) or tptA; tptA not
   required. Do not require a flank on EACH contig. rnhA can be far from the sex
   gene in *Absidia*/*Backusella*, and sexM is flankless in *Lichtheimia*.
   Require a flank on at least one contig.
3. Degraded flank copy stays optional support. It is not a general signature:
   *Z. heterogamus* and *M. africana* have one intact rnhA.
4. Duplicated flank genes: divergent intact copies (~74-76% nt, as in
   *Z. heterogamus*) support homothallism. Near-identical copies suggest two
   allelic haplotypes (mixed/heterokaryon). The cut-off is not known.
5. Keep GC, depth, duplicated-BUSCO and copy-number checks. The *A. glauca* result
   shows that fusion products can carry both genes in one nucleus, so these
   checks cannot fully exclude fusion.
6. SexM classifier confidence: consider a softer bar for sexM than for sexP.

Expected positives in LCG: *Syzygites* (unlinked); *Zygorhynchus heterogamus*
(same locus); *Mycotypha africana* (same contig, distant). *Z. moelleri*:
unknown (sexM not found). *M. genevensis*, *R. homothallicus*, *R. azygosporus*:
no data in this paper.
