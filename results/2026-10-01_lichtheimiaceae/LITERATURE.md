# Lichtheimiaceae mating type and sex locus: literature review

Date: 2026-10-01. Read-only literature task. Sources were read from full text where noted.
"Source" lines report what a paper shows. "Inference" lines are my interpretation.

## 0. Family delimitation

- Source: Hoffmann et al. 2013, Persoonia 30:57-76, PMID 24027347, doi:10.3767/003158513X666259 (full text read).
  They restrict Lichtheimiaceae to *Lichtheimia* + *Dichotomocladium*.
  They place two neighbouring clades *incertae sedis*: *Rhizomucor*/*Thermomucor* and
  *Fennellomyces*/*Circinella*/*Thamnostylum*/*Zychaea*/*Phascolomyces*.
- Source: Walther et al. 2019, J Fungi 5:106, PMID 31739583, doi:10.3390/jof5040106 (full text read).
  It states: "The Lichtheimiaceae were extended by the genera Circinella, Dichotomocladium,
  Fennellomyces Phascolomyces, Rhizomucor, Thamnostylum, Thermomucor, and Zychaea."
  The main *Circinella* group (type *C. umbellata*) sits in Lichtheimiaceae; the genus is polyphyletic.
- Inference: the family is broad and its internal branches are poorly supported.
  Do not assume one flank set fits all nine genera.

## 1. Mating behaviour

| Taxon | Reported system | Source |
|---|---|---|
| *Lichtheimia* (genus) | Heterothallic; zygospores with equatorial rings, suspensors without appendages | Walther et al. 2019 (PMID 31739583); Alastruey-Izquierdo et al. 2010 (PMID 20357218) |
| *L. corymbifera*, *L. ramosa*, *L. ornata*, *L. hyalospora*, *L. sphaerocystis* | Mating tests used for biological species recognition | Alastruey-Izquierdo et al. 2010, doi:10.1128/JCM.01744-09 |
| *Rhizomucor pusillus* | Mostly heterothallic; rare homothallic strains | zygomycetes.org Rhizomucor page; Vastag et al. 1998 JCM 36:2153 (search snippet only) |
| *R. miehei* | Homothallic | same |
| *R. nainitalensis*, *R. endophyticus* | Homothallic | zygomycetes.org (*R. endophyticus* now in *Mucor*, Walther 2019) |
| *Thermomucor indicae-seudaticae* | Homothallic | zygomycetes.org Thermomucor page (secondary) |
| *Circinella* | Homo- or heterothallic, species-dependent | secondary web sources only |
| *Dichotomocladium*, *Zychaea*, *Fennellomyces* | Not verified. *Fennellomyces heterothallicus* exists (name only, Hoffmann 2013) | — |

- Source: Alastruey-Izquierdo 2010 abstract: "mating tests did not show intrinsic reproductive
  barriers for two pairs of the phylogenetic species." I could only read the abstract.
  Search-engine summaries report 168 crosses (73 intra-, 95 interspecific) with zygospores in 17.
  I did not verify these numbers in the paper.
- Source: Nottebrock et al. 1974 obtained zygospores in *L. corymbifera* x *L. ramosa* crosses.
  Garcia-Hermoso et al. 2009 (PMID 19759217) later separated the two species (as cited by Walther 2019).
  Walther 2019 concludes zygospore presence alone does not prove conspecificity.
- Gap: I found no paper that reports the Plus/Minus type of the sequenced strains
  (*L. corymbifera* FSU 9682, *L. ramosa* FSU 6197).

## 2. MAT (sex) locus structure

- Source: Schulz et al. 2016, Endocytobiosis Cell Res 27(4):39-57 (local PDF, pp. 42-43 read).
  Table 1 lists a SexM for three *Lichtheimia* genomes and no SexP:
  *L. corymbifera* (JGI 12200; strain given as "JMRC FSU 6982", likely a typo for 9682),
  *L. hyalospora* FSU 10163 (JGI 126746), *L. ramosa* (CDS03202.1).
  AlgL, TptA and RnhA are "N/A" for all three in the table.
  Text, p. 43: "Here the sexM gene is located on a different chromosome than the other genes
  from the usual cluster, rnhA, algL and tptA, which have retained synteny."
- Source: CDS03202.1 (NCBI) is a 309 aa "hypothetical protein" LRAMOSA00604 with an
  HMG-box (cd01389, aa 47-118). Strain JMRC FSU:6197. From Linde et al. 2014, PMID 25212617.
- Source: Schwartze et al. 2014, PLoS Genet, PMID 25121733 (not 25502079), doi:10.1371/journal.pgen.1004496
  (full text read). It names no sexM/sexP, tptA or rnhA. It only notes the absence of
  the MAT a1 domain (PF04769) and says Mucorales mating uses sex plus/minus HMG factors.
  Assembly: 209 scaffolds.
- Source: Schwartze et al. 2015 Pathog Dis, PMID 25857734: four *Lichtheimia* drafts
  (1,176-3,968 contigs). No mating content.
- Source: NCBI keyword searches (nuccore, protein) and Europe PMC full-text search found
  no GenBank sexM/sexP deposits and no dedicated sex-locus paper for any Lichtheimiaceae genus.
  Rhizomucor: no sex-locus data found.
- Inference: "different chromosome" in Schulz 2016 rests on draft scaffold assemblies.
  Neither genome was chromosome-level in the cited sources. The data show "different scaffold".
  A real translocation and an assembly break are both possible. No source resolves this.

## 3. Conserved neighbours usable as family flanks

- Source: none found. No paper lists genes adjacent to sexM in Lichtheimiaceae.
  Schulz 2016 gives no protein IDs for the Lichtheimia tptA–algL–rnhA block.

## 4. Implications for MATPredict detection

- Inference: tptA/rnhA/algL anchors will likely not co-locate with sexM in *Lichtheimia*.
  A flank-gated or cluster-gated call may miss the HMG gene. Search sexM/sexP genome-wide.
- Inference: the local neighbours of *Lichtheimia* sexM must be found empirically
  (e.g. compare gene order around sexM across available *Lichtheimia* assemblies).
- Data gap: no source gives a Plus/Minus ratio. Schulz found SexM in 3/3 *Lichtheimia*
  genomes. That is n=3 and says nothing about population ratios.
- Pitfall: Mucorales genomes carry many HMG-box genes. *L. corymbifera* has expanded
  gene families (Schwartze 2014). Use SexM/SexP-specific references, not any HMG hit.
- Pitfall: homothallic *Rhizomucor miehei*, *Thermomucor* and some *R. pusillus* strains
  may carry both sexM and sexP. Do not treat co-occurrence as an error in those taxa.
- Pitfall: highly fragmented drafts (*L. ramosa* N50 about 34 kb) can split any locus.
