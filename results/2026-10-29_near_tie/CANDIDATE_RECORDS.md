# Candidate reference records for the SLA2-type genera (draft, 2026-10-29)

Search: PubMed (via the PubMed MCP) and NCBI Protein and Assembly (E-utilities) for MAT loci in Barnettozyma, Petasospora (Cyberlindnera fabianii), Wickerhamomyces, Ambrosiozyma, Saccharomycopsis, Citeromyces, Williopsis (Cyberlindnera jadinii), Hansenula, Kuraishia and Pachysolen. Only three leads are real; nothing here is curated yet.

## 1. *Saccharomycopsis schoenii*: best candidate (published locus, RNA-seq, chromosome-scale assembly)

- **Paper:** Kriti D, Barazandeh M, Ogbede JU, Giaever G, Nislow C. G3 2026;16(5), [10.1093/g3journal/jkag067](https://doi.org/10.1093/g3journal/jkag067) (PMID 41849659, PMC13148404). Chromosome-scale assembly (6 nuclear chromosomes, PacBio HiFi + Hi-C), BRAKER2 annotation, RNA-seq.
- **What it says about MAT:** an expanded MAT system, with SLA2-MAT cassettes (MATa1, MATa2, MATalpha1, MATalpha2) on every chromosome; SLA2 stays tightly linked to the MAT genes on every chromosome while DIC1 is scattered; the cassettes are co-transcribed (189-247 bp splice junctions between MATa1 and SLA2; conserved 58 bp introns between MATa2 and MATa1). The active locus is on chromosome IV: a leftmost alpha-dominant sub-cluster (470,610-475,000: MATalpha2, MATalpha1) and a rightmost palindromic cassette (528,044-537,777; SLA2 anchor flanked by mirrored MATa1 and MATa2 arrays). The species is homothallic (Dost et al., Microbiol Res 2024, [10.1016/j.micres.2024.127691](https://doi.org/10.1016/j.micres.2024.127691), PMID 38492364; *S. lassenensis* is heterothallic).
- **Assembly:** GenBank JBTFFQ000000000, BioProject PRJNA1364152 (data availability). The v0.7.0 genome table has GCA_054953955.1 (NCYC 2975, chromosome level, submitted 2026-02-09), which I expect to be this assembly (not checked against the WGS prefix). Our replay calls six `MTL alpha` loci there (CM143126-CM143131), each with SLA2 (40-45% identity) inside the locus, no DIC1 and no PIK1/OBP1/PAP1: the published SLA2-MAT architecture, not MTL.
- **Record plan:** one record per active locus, `system` homothallic or `unknown`, the gene models of the deposit (annotation plus RNA-seq, tier 2), a family of its own or an SLA2-MAT family (curator).
- **Needs:** the deposit's annotation (GFF) and exact coordinates of the chromosome IV loci; the MATa1/MATa2 intron structure.

## 2. *Barnettozyma californica* and relatives

- **Paper:** Hession C, Byrne KP, Wolfe KH, Butler G. PLoS Genet 2025;21(12):e1011814, [10.1371/journal.pgen.1011814](https://doi.org/10.1371/journal.pgen.1011814) (PMID 41359681, PMC12711049). A centromere formed at the MAT locus of *B. californica* (250 kb from the old one); the same chromosome carries MAT and MRP17-PUT2 about 250 kb apart in all three *Barnettozyma* species studied.
- **MAT structure (S3 Fig, text only read):** a 460 bp centromere-associated sequence lies in a 2.8 kb intergenic region **between MATa1 and MATa2 on one side and MATalpha1 and MATalpha2 on the other**, so one locus carries both gene pairs.
- **Our replay:** five *B. californica* assemblies (UCD09 GCA_003123585.1, UCD89 GCA_003123885.1, NRRL Y-17395 GCA_003706655.x, PB4207 GCA_023627795.1) give one MTL call of 16 kb with SLA2 inside and **DIC1 8.2 kb away**, no PIK1: the SLA2-MAT-DIC1 arrangement.
- **Assembly with the best contiguity:** *B. discipulorum* GCA_049863505.1 (PYCC:10016, complete genome, 2025-04-23; not in our run). The accessions the paper used are in its S1 Table, which was not read.
- **Record plan:** a record from the paper's locus (both idiomorph gene pairs at one locus), pseudohomothallic or `unknown`, tier 2; family decision for the curator (a Phaffomycetales/Barnettozyma family, or MATog-like).
- **Needs:** S1 Table and S3 Fig (supplement) for the assembly and coordinates.

## 3. *Petasospora (Cyberlindnera) fabianii*: annotation-only

- Protein ONH67106 "Mating-type-like protein ALPHA1" (strain 65, MPUK01000005.1, locus tag BON22_3092, MATalpha_HMGbox region): a draft-genome annotation with no paper on the locus (PubMed finds none), tier 3 at best.
- Chromosome-level assemblies: GCA_022641835.1 ("Missouri", CP048715.1): our MATsc call is at 630,345-685,843 with DIC1 1.1 kb from the locus and SLA2 inside; GCA_050575385.1 (TBRC4498, chromosome) is uncalled in v0.7.0.
- **Record plan:** a manual annotation of the CP048715.1 locus (the same procedure as for CBS 573) would make this a real record; needs the curator's go-ahead on the evidence bar (no publication, no RNA-seq).

## Leads not followed
- *Williopsis (Cyberlindnera) jadinii* XP_020071302: hypothetical HMG-box protein, nothing on MAT.
- *Wickerhamomyces anomalus* (isolate KG16, JAHTLX01) is one of the seven Phaffomycetales genomes of Hession et al.; whether the supplement describes its MAT locus was not checked.
- *Ambrosiozyma*, *Citeromyces*, *Hansenula*, *Kuraishia*, *Pachysolen*: no MAT-named protein records or papers found.

## Next steps
1. Fetch the S. schoenii annotation (GCA_054953955.1 or JBTFFQ) and the Hession et al. supplement (S1 Table, S3 Fig) for the coordinates.
2. Draft records for S. schoenii (chromosome IV) and *B. californica* with the existing `curate-db propose` workflow, plus a family decision (separate or MATog-like).
3. Curator: the family structure (one SLA2-MAT family for these genera, or per genus) and the evidence bar for the annotation-only *P. fabianii* lead.

Sources (PubMed): 41849659, 41359681, 38492364.
