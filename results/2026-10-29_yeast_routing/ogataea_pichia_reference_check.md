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
- **MATa**: CBS573 chromosome 2 (CP028774.1, RefSeq NC_042507.1), same strain, the region opposite MH260578. From the RefSeq
  feature table (genes are all "uncharacterized protein", so these are identified by position and size against MH260578, not named in the annotation):
  SLA2 at 307,183-310,401 (+, 3,219 bp, the same length as MH260578 SLA2); the next gene at 310,491-311,300 (-, 810 bp) is the
  MATa1 candidate; a 267 bp ORF at 311,575-311,841 (+) is a short candidate for MATa2 and needs a gene-model check;
  TGL1 candidate at 312,347-313,684 (-, 1,338 bp). The MATa1/MATa2 names need confirming by sequence (the paper lists them but the annotation does not).
- Both records are the same strain (CBS573), which suits the database rule of same-strain records.
- Scope: Pichia has no family (157 of 158 genomes reach it by phylum fallback). A `Pichia` family scoped to the genus
  (the taxid of *Pichia* itself still has to be looked up), heterothallic, with flank genes SUI1, SLA2, TGL1, SEC3.
  Whether to add it is the curator's decision; adding it would also need the routing rule of `RULE_DRAFT.md` so a Pichia genome is not called by MATsc and MTL as well.

## Not checked

- The 27 Pichia genomes with two systems were not examined individually.
- MATa2 in CBS573 (above) is a candidate by position and size only.
- No read or expression data was used for Ogataea; which cassette is active is unknown.
