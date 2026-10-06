# Agaricomycete A and B subloci: one locus or separate calls?

Date: 2026-09-28. Read-only literature review. No code or data changed.

## Question

In tetrapolar Agaricomycetes, the A (HD) locus has Aα/Aβ subloci and the B (P/R) locus has Bα/Bβ subloci (or groups 1-3 in *Coprinopsis*).
Should MATPredict report each sublocus as its own mating-type call?
Or should it report one locus with sublocus evidence inside it?

## Short answer

The literature treats each sublocus as a separate specificity unit.
Each sublocus is functionally independent: a difference at any one sublocus is enough for compatibility.
Subloci can recombine, so new combinations arise.
But the literature still defines the mating-type **locus** (the "A factor" and the "B factor") as the unit that segregates at meiosis and gives tetrapolar behaviour.
Every paper and review we read says "two mating-type loci, A and B", each made of linked subloci.
Recommendation: report **one A call and one B call** per genome. Record the subloci as structured evidence inside each call. Do not report subloci as separate mating-type calls.

## Sources used and access limits

- Full text read: van Peer et al. 2011 (PMC3140503); Coelho et al. 2017 review (PMC5467461); Chu et al. 2025a/b (PMC12193891, PMC12028200); Simchen 1967 (Genet. Res., rendered from PDF).
- Abstract only (full text not accessible here): Raper et al. 1958; Day 1960; Kües et al. 1994; Pardo et al. 1996; Lukens et al. 1996; O'Shea et al. 1998; Halsall et al. 2000; Riquelme et al. 2005; Specht 1995; Wendland et al. 1995; Vaillancourt et al. 1997; Fowler et al. 2001, 2004; Koltin & Stamberg 1973; Ohm et al. 2010; James et al. 2006, 2013; Lentinula papers.
- Not read: Raper 1966 book; James 2015 and Kües 2015 (Fungal Biol Rev, not in PubMed); Casselton & Kües 2007 (book chapter). Statements attributed to those works below come only from other papers that cite them.

## 1. *Schizophyllum commune*

### A locus (Aα, Aβ)

What the papers show:
- Aα holds one HD1/HD2 pair, Z (HD1) and Y (HD2). They are divergently transcribed (Shen et al. 2001, PMID 11525408, doi:10.1007/s002940100219). The Aα1 allele has only Y (Specht et al. 1994, PMID 8088516, doi:10.1093/genetics/137.3.709).
- Activation needs Y and Z from different Aα specificities (Specht et al. 1992, PMID 1353887, doi:10.1073/pnas.89.15.7174).
- Aβ has several HD genes. Chu et al. 2025 report two HD1/HD2 pairs (S/R and Q/V). The two pairs act independently: "R-S and V-Q interactions independently regulate sexual compatibility" (PMID 40558964, doi:10.3390/jof11060451).
- Aβ proteins alone can drive the A pathway when Aα is deleted (Robertson et al. 1996, PMID 8978032, doi:10.1093/genetics/144.4.1437).
- Specificities: about 9 at Aα and about 32 at Aβ (Raper et al., as cited by Chu et al. 2025, PMID 40278098). Simchen 1967 cites Raper et al. 1960 as "nine and fifty" estimated and "nine and twenty-five" found in the wild (doi:10.1017/S001667230001048X).
- Physical distance: "more than 500 kb apart and frequently recombine during meiosis" (Coelho et al. 2017, PMID 28597825, doi:10.1128/microbiolspec.FUNK-0046-2016). van Peer et al. 2011 give "∼450 Kb" (PMID 21799803, doi:10.1371/journal.pone.0022249). A web search result attributes 550 kb to Ohm et al. 2010 (PMID 20622885, doi:10.1038/nbt.1643). We could not open that paper to confirm the number.
- The *pab-1* gene lies between Aα and Aβ in *S. commune*. Kües et al. 2001 see this as a sign of translocation (PMID 11318102, doi:10.1007/s002940000176). van Peer et al. 2011 link the separation to rearrangement of large gene clusters.
- Recombination frequency (measured): Simchen 1967 found 16/330 = 4.85 ± 1.40% recombinant A factors in one wild isolate. In 36 sibling dikaryons from that isolate it ranged from 0% to 19%, with a mean of 6.58% (Table 1). The paper opens with: "a most striking … puzzling instance of heterogeneity in recombination values was found in different strains … for the linked sub-units of the A incompatibility factor (Raper et al., 1958a, 1960)".
- Tetrad analysis showed non-parental A factors are "the reciprocal products of conventional crossing-over" (Papazian 1951, as cited by Simchen 1967).

### B locus (Bα, Bβ)

What the papers show:
- Bα1 holds one receptor gene (*bar1*) and several pheromone genes (Wendland et al. 1995, PMID 7489716, doi:10.1002/j.1460-2075.1995.tb00211.x).
- Bβ1 holds at least three pheromone genes and one receptor gene. Bβ2 holds one receptor and eight pheromones (Vaillancourt et al. 1997, PMID 9178005, doi:10.1093/genetics/146.2.541; Fowler et al. 2001, PMID 11514441, doi:10.1093/genetics/158.4.1491).
- Vaillancourt et al. 1997: the genes "are clustered in each of two recombinable and independently functioning loci, Bα and Bβ. A difference in specificity at either locus … initiates an identical series of events". Bβ genes function "only within the series of Bβ specificities". Bα and Bβ "arose from a common ancestral sequence".
- Fowler et al. 2001 call them "two redundantly functioning B mating-type loci".
- Specificities: nine at Bα and nine at Bβ (Raper, as cited by Chu et al. 2025, PMID 40278098; Gola & Kothe 2003, PMID 12589467).
- Physical distance: "genes of the … Bα and Bβ mating-type loci are shown to be within a few kilobases of each other". Some pheromones activate both a Bα and a Bβ receptor. Bar8 (Bα8) is "functionally identical" to Bbr1 (Bβ1) (Fowler et al. 2004, PMID 14643262, doi:10.1016/j.fgb.2003.08.009). So the Bα/Bβ boundary is not strictly separate in function.
- Recombination: Bα–Bβ recombination frequency is under genetic control. A modifier gene, *B-rec-1*, maps about 9 map units from Bβ (Koltin & Stamberg 1973, PMID 17248610, doi:10.1093/genetics/74.1.55).
- **Not found:** we could not get a numerical Bα–Bβ recombination frequency from an accessible primary source.

## 2. *Coprinopsis cinerea*

### A locus

What the papers show:
- Each A locus has genes "separated into two functionally independent complexes termed Aα and Aβ" (Pardo et al. 1996, PMID 8878675, doi:10.1093/genetics/144.1.87).
- The archetype holds three paralogous HD1/HD2 pairs: a (in Aα), and b and d (in Aβ). "Different allelic versions of gene pairs are compatible but paralogous genes are incompatible" (Pardo et al. 1996). Extant alleles often lack some genes. For example, A43 carries a non-functional a2-2 and a c1-1 remnant (Kües et al. 1994, PMID 7845358, doi:10.1007/BF00279749).
- Physical distance: Aα and Aβ are "separated … by 7 kb of noncoding sequence" (Kües et al. 1994).
- Recombination: "all recombination events were located in 6 kb of noncoding DNA between the alpha and beta subloci and the rate of recombination in this noncoding region matched that generally observed for this genome. No recombination within gene clusters". The authors "propose that pairs of genes constitute both the sex determining and the hereditary unit of A" (Lukens et al. 1996, PMID 8978036, doi:10.1093/genetics/144.4.1471).
- Coelho et al. 2017 describe the two subloci as able to "recombine at detectable frequencies despite their close proximity".
- **Not found:** the recombination percentage in accessible text. Day 1960 (PMID 17247950, doi:10.1093/genetics/45.5.641) and Lukens 1996 hold the numbers, but only their abstracts were accessible.
- Allele counts: 4, 7 and 3 alleles for the a, b and d pairs give 84 combinations. Raper estimated 160 A specificities (Coelho et al. 2017).

### B locus

What the papers show:
- The B6 locus holds nine genes (three receptors and six pheromones) in 17 kb of mating-type-specific sequence. These are "three functionally independent subfamilies" (O'Shea et al. 1998, PMID 9539426, doi:10.1093/genetics/148.3.1081).
- In B42, "the three genes within each group are kept together as a functional unit". "Different B loci may share alleles of one or two groups". "It is the different combinations of their alleles that generate the multiple B mating specificities" (Halsall et al. 2000, PMID 10757757, doi:10.1093/genetics/154.3.1115).
- Riquelme et al. 2005 found 2, 5 and 7 alleles for groups 1, 2 and 3, giving 70 combinations (PMID 15879506, doi:10.1534/genetics.105.040774). Receptor sequence clades "do not correspond to the groups defined by position".
- Coelho et al. 2017 call the groups "groups or subloci 1 to 3", with "freely interchangeable alleles".
- **Not found:** a direct measurement of meiotic recombination between B groups. That groups reshuffle is inferred from shared group alleles among different B loci, not from counts of meiotic recombinants.

## 3. Other Agaricomycetes

What the papers show:
- *Flammulina velutipes*: the matA subloci (matA3a, matA3b) are 73 kb apart and "linked despite their 73 Kb distance". The matB subloci (matB3a, matB3b) are 177 kb apart. The authors use "loci with redundant function (subloci)" and "Each mating-type locus consists of tightly linked subloci" (van Peer et al. 2011, PMID 21799803). They say typical Agaricales subloci are "closely linked (10–20 kb)".
- *Lentinula edodes*: the B locus is bipartite, and each sublocus has one receptor with one or two pheromones (Wu et al. 2013, PMID 24029079, doi:10.1016/j.gene.2013.08.090). Kim et al. 2020 report five Bα and three Bβ alleles (PMID 32375416, doi:10.3390/genes11050506). Li et al. 2015 saw intralocus recombinants in both A and B, at 28/189 progeny of one strain (PMID 23996277, doi:10.1002/jobm.201300313). This is a count of atypical-mating-type progeny, not a clean recombination frequency. A conserved Bα-like locus 43.3 kb from matB is not mating-type specific (Lee et al. 2021, PMID 35035249). An incomplete HD sublocus lies ~2.8 Mb from matA (Gao et al. 2022, PMID 35205921).
- *Coprinus bilanatus*: HD1/HD2 genes are "distributed over two closely linked subloci, Aα and Aβ" (Kües et al. 2001, PMID 11318102).
- *Coprinellus disseminatus* (bipolar): the A-homologous locus "encodes two tightly linked pairs of homeodomain transcription factor genes" (James et al. 2006, PMID 16461425, doi:10.1534/genetics.105.051128).
- Polyporales: James et al. 2013 describe "a MAT-HD locus" as one unit. The abstract does not use sublocus terms (PMID 23928418, doi:10.3852/13-162).
- *Pleurotus eryngii* and *P. tuoliensis*: the papers found describe one A locus with an HD1/HD2 pair and several B receptors. We did not find sublocus-level functional data (PMID 31448140; PMID 29304732).
- *Heterobasidion*, *Phanerochaete*, *Pholiota*: these are bipolar. The HD locus alone determines mating type (Coelho et al. 2017; James et al. 2011, PMID 21131435).
- *Agaricus bisporus/bitorquis*: **not found** in this search. We have no sublocus data for them.

Inference: the sublocus idea is strongest in *Schizophyllum*, *Coprinopsis*, *Flammulina* and *Lentinula*. Elsewhere, papers usually describe one HD locus and one P/R locus. The number of HD pairs and receptor groups varies by species and by allele.

## 4. Evolutionary interpretation and terminology

What the reviews say:
- Coelho et al. 2017: "multiallelism in both P/R and HD loci is generated by rounds of segmental duplication and diversification of MAT genes resulting in independent, but functionally redundant subloci, upon which recombination can act to give rise to novel allele specificities". The same review still counts two MAT loci (P/R and HD) for tetrapolar species.
- van Peer et al. 2011: "The redundant subloci are a result of doubling during evolution".
- Pardo et al. 1996 call the *C. cinerea* HD pairs "three paralogous pairs". Vaillancourt et al. 1997 derive Bα and Bβ from "a common ancestral sequence".
- Raper's classical terms were "A factor" and "B factor", each with two linked "loci" or "sub-units" (α and β). Simchen 1967 uses exactly these terms.

Inference: the subloci are paralogous specificity units made by tandem or segmental duplication inside one mating-type locus. The curator's "a different type of duplication" fits this reading. They differ from ordinary gene duplicates in three ways. (1) Each sublocus has its own allelic series. (2) Paralogues from different subloci do not cross-activate in most cases (Pardo 1996), though Fowler 2004 reports exceptions in *S. commune* B. (3) Subloci recombine with each other at measurable rates.

## 5. Recommendation for MATPredict

The evidence:
1. Every source defines the mating-type locus as A (HD) and B (P/R). Tetrapolar means two loci, not four or six.
2. Subloci are functionally independent specificity units. Recombination between them is real: 0–19% (mean 6.6%) for *S. commune* Aα–Aβ, which are >450 kb apart. In *C. cinerea* recombination falls in a ~6–7 kb intergenic spacer at about the genome-average rate. But this recombination produces new A or B specificities. It does not produce new loci.
3. How many subloci there are, and where they sit, varies by lineage (7 kb to >450 kb for A; a few kb to 177 kb for B). In some species it varies by allele too, because genes are missing from some alleles.

Recommendation, the same for A and B:
- Emit **one call per locus type** (one HD/A call, one P/R/B call) per haplotype or assembly.
- Inside each call, list the subloci as structured evidence. For each sublocus give: label (α/β or group 1–3), member genes (HD1/HD2 pair, or receptor plus pheromones), coordinates, and gene completeness.
- Group subloci into one locus call when they share the conserved locus neighbourhood (for HD: *mip*/β-fg synteny). This holds even when they are far apart (*S. commune* >450 kb, *Flammulina* 73 kb and 177 kb). Do not use a fixed distance cutoff alone to split them.
- Flag "sublocus separated by >X kb" as information, not as a second locus.
- Keep non-mating-type-specific receptors (e.g. *S. commune brl* genes, *L. edodes* Bα-N) out of the B call, or mark them as non-specific. The literature says they do not set mating type.

What the data do not support: we have no numerical Bα–Bβ recombination rate. We have no sublocus data for *Agaricus*. The *C. cinerea* recombination percentages are only in papers we could not open. None of these gaps changes the recommendation, because the recommendation rests on how the locus is defined, not on the exact rates.
