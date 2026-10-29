# Ogataea and Pichia reference check (#189)

Scripts: `ogataea_loci_dump.py`, `ogataea_pairing.py`, `pichia_kudriavzevii_calls.py` (read the v2 and v0.7.0 reports on HPCC).

## Ogataea: two MATog loci are two cassettes, not duplicates

The MATog records in the database each carry one locus: CBS 4732 (AJ617305.1) is DIC1, MATalpha2, MATalpha1,
MATa1, SLA2 in one 6.2 kb deposit, and DL-1 (NC_027862.1:422526-429249) has the same five genes in one 6.7 kb
region. In the v2 run (Ogataea routed to MATog), the Ogataea and cross-called genomes (134) fall into these patterns:

| pattern in the genome | genomes |
|---|---|
| two MATog loci, one MATa and one MATalpha, on different contigs | 71 |
| one MATog locus | 54 |
| two MATa, or four MATalpha (a hybrid assembly) | 1 and 1 |
| no MATog call | 7 |

- In the four genomes inspected by hand (`ogataea_loci_dump.py`: O. polymorpha 4329, and three others) the MATa locus is SLA2 with MATa1 and the
  MATalpha locus is DIC1 with MATalpha1 and MATalpha2, each flanked by SLA2 or DIC1 as in the reference records. The gene content of the
  other 67 pairs was not tabulated. Read this way the pair is the two cassettes of the reference arrangement split across contigs, not a
  second family and not a duplicate call.
- Of the 54 one-locus genomes, 9 carry MATa1, MATalpha1 and MATalpha2 together with DIC1 and SLA2 (the
  reference layout in one piece); the others carry one cassette.
- The split looks like an assembly effect at least as much as biology: in 69 of the 71 two-cassette genomes at least
  one of the loci starts within 500 bp of its contig start, and in 37 both do. A possible cause is that the two cassettes
  share flanking genes, which would make them hard to assemble; that was not tested. This data cannot say which cassette is the expressed one, or whether the
  two sit apart in the chromosome or are one locus broken by the assembly (DL-1, which is assembled to a chromosome, has them in one 6.7 kb region).
- Consequence for the routing rule: for MATog the output rule must not merge two loci of different idiomorphs. The
  overlap dedupe (#209) and the cross-call rule only merge overlapping calls of different families, so two
  MATog loci that do not overlap are untouched. Labelling them as MATa and MATalpha cassettes of one pseudohomothallic
  genome is a naming decision for the curator.

## Pichia: heterothallic, no silent loci

Douglass et al. 2018 (PLoS Pathog, PMID 30024981, doi 10.1371/journal.ppat.1007138), read from the PMC full text:
- *P. kudriavzevii* CBS573 is heterozygous at the MAT locus, about 310 kb from the left end of chromosome 2. The
  assembly of the chromosome holds MATa1 and MATa2; the alternative allele with MATalpha1 and MATalpha2 is on a
  short PacBio contig, deposited as **MH260578** (11,807 bp: SUI1, SLA2, MATalpha1, MATalpha2, TGL1, SEC3, all annotated).
- The MAT genes lie between SUI1-SLA2 on one side and TGL1-SEC3 on the other. "The genome contains no silent loci,
  indicating this species is incapable of mating type switching and therefore heterothallic." CBS5147 is also heterozygous.

So for Pichia the second-locus calls are not silent cassettes: in a heterozygous diploid assembly the other
allele can sit on a separate contig. Nothing here tests that for the 27 Pichia genomes in #189; it is the explanation to check first.

### What v0.7.0 does on the four complete *P. kudriavzevii* assemblies

| assembly (strain) | locus found | call |
|---|---|---|
| GCF_003054445.1 (CBS573) | NC_042507.1:307258-310766 | MATsc, medium: sla2 (57%) and a 141 bp piece of MATA1 (51%), MATa idiomorph |
| GCA_003054405.1 (CBS5147) | CP028532.1:302754-306262 | same two genes |
| GCA_014873065.1 (CY902) | CP039615.1:317913-321421 | same two genes |
| GCA_003033855.1 (SJP) | CP021089.1:305747-330409 | MATsc, medium: sla2, MATALPHA1 (36%, unpolished), MATALPHA2 (37%), MATalpha idiomorph; a 24.7 kb span |

The locus is found at the right place, but only through the Saccharomyces references and with one partial gene
(MATA2 missing in the three MATa assemblies), which is the case for a Pichia record.

### Candidate record (not added; needs the curator)

- **MATalpha**: MH260578 as deposited, complete and annotated (SUI1, SLA2, MATalpha1 complement(4106..4645),
  MATalpha2 complement(5119..5918, three exons), TGL1, SEC3).
- **MATa**: CBS573 chromosome 2 (CP028774.1, RefSeq NC_042507.1), same strain, the region opposite MH260578. The RefSeq genes there are all
  "uncharacterized protein", so the check below identifies them against MH260578 and the database (`pichia_cbs573_flank_miniprot.py`,
  `pichia_cbs573_hmm_scan.py`, `pichia_reference_proteins.faa`):
  - **Flanks, same strain and synteny confirmed**: miniprot of the MH260578 proteins on NC_042507.1 gives SUI1 306,404-306,947 (-, 100%), SLA2 307,183-310,401
    (+, 100%), TGL1 312,365-313,684 (-, 98%) and SEC3 313,903-317,976 (-, 100%). So SLA2 and TGL1 bracket a 1.96 kb gap, against 2.2 kb between SLA2 and TGL1 in the alpha allele.
  - **MATa1**: the gene at 310,491-311,300 (-, 810 bp). The O. polymorpha MATa1 protein from the MATog record hits it at 55% over 310,497-310,670 (partial). K. phaffii MATa1 and MATa2 and the C. albicans MTLA1 and MTLA2 give no hit.
  - **MATa2**: the Saccharomycotina MATa2 profile (`db/Ascomycota/classifiers/MATa2_saccharomycotina`) hits a 125 aa open reading frame at about 311,868-312,243 (+; E 2.9e-10, score 31.7)
    and, more weakly, a 110 aa one at about 311,518-311,848 (+; E 1.2e-5), which is the 267 bp ORF I had listed. The two are 19 bp apart, so they are probably parts of one gene
    with an intron, but this scan reads unspliced six-frame ORFs and does not give a gene model. **The MATa2 model needs RNA-seq (SRA SRP139056) or a splice-aware prediction before it goes into a record.**
  - Positive control: the same scan of the MH260578 region finds the MAT1-1-1 (alpha-box) profile on MATalpha1 (190 aa, E 3.2e-7) and nothing else.
  - The earlier statement that 311,575-311,841 is a MATa2 candidate is replaced by the above.
- Both records are the same strain (CBS573), which suits the database rule of same-strain records.
- Scope: Pichia has no family (157 of 158 genomes reach it by phylum fallback). A `Pichia` family scoped to the genus
  (the taxid of *Pichia* itself still has to be looked up), heterothallic, with flank genes SUI1, SLA2, TGL1, SEC3.
  Whether to add it is the curator's decision; adding it would also need the routing rule of `RULE_DRAFT.md` so a Pichia genome is not called by MATsc and MTL as well.

## Not checked

- The 27 Pichia genomes with two systems were not examined individually.
- The MATa2 gene model in CBS573: only an unspliced profile hit; no gene model.
- No read or expression data was used for Ogataea; which cassette is active is unknown.

## CBS 573 MATa2 gene model from RNA-seq (2026-10-29)

Reads: 72.1 M unpaired, hisat2 against NC_042507.1, 96.2% aligned (job 31372344). Scripts: `pichia_mata2_model.py`, `pichia_mata2_hmm.py`; candidates in `pichia_mata2_candidates.faa`. Window coordinates: `win.fa` index 0 is genome position 295,000.

| Item | Result |
|---|---|
| Junction 311,842–311,903 (+ strand) | 27 spliced reads, 62 bp intron, GT…AG |
| Junction 306,703–306,916 (− strand) | 53 reads, 214 bp, GT…AG on the reverse strand; a different gene, not MATa2 |
| Other junctions in 305–320 kb | none |
| Stop codon | spliced ORF ends at 312,245 (two exons) |
| Depth in 311,450–311,640 | 0–2; rises to 23+ at ~311,670 |
| MATa2 profile (supporting evidence only) | 52.4 bits (E 8e-18) on the 202 aa ORF from ATG 311,575; 35.9 bits on the 131 aa ORF from ATG 311,788 |

Settled: a two-exon gene with a 62 bp intron after 311,841 and the stop at 312,245. The HMG box (profile hit) spans exon 1 and the junction.

Not settled: the N-terminus. ATG 311,575 gives the better profile score but lies 95 bp upstream of where RNA coverage begins (~311,670); ATG 311,788 is covered but loses the first part of the HMG box. No junction links the upstream region to exon 1. Needs a stranded sample, or a 5' RACE or protein evidence, to choose.

Depth stays high (100–900) from 312,100 to 313,450, past the stop; it is probably the neighbouring gene or UTR, and was not interpreted.
