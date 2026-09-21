# Pucciniomycotina MAT-locus literature review — *Puccinia*, *Leucosporidium*, *Cystobasidium*, *Rhodotorula*

Date: 2026-09-20. Requested: broaden curation "further out" than Ustilaginomycotina,
into Pucciniomycotina, starting with these four genera. This is a **literature
review**, not a curation pass — nothing here has been added to `db/`. Every claim
below is sourced; where evidence is thin or absent, that is stated plainly rather
than filled in.

## Why these four genera are not one story

Unlike the Ustilaginomycotina batch curated earlier this session (which shares one
architecture — tetrapolar HD + PR, `bE`/`bW`/`mfa1`/`pra1` naming — across every
genus), these four span **three different classes** within Pucciniomycotina:

| genus | class | order | family |
|---|---|---|---|
| *Puccinia* | Pucciniomycetes | Pucciniales | Pucciniaceae |
| *Leucosporidium* | Microbotryomycetes | Leucosporidiales | Leucosporidiaceae |
| *Rhodotorula* | Microbotryomycetes | Sporidiobolales | Sporidiobolaceae |
| *Cystobasidium* | Cystobasidiomycetes | Cystobasidiales | Cystobasidiaceae |

Pucciniomycotina is the earliest-branching Basidiomycota lineage, and its MAT-locus
architecture (tetrapolar vs. bipolar, linked vs. unlinked) is exactly the subject of
active research (see Maia et al. 2015 below) precisely *because* it is not uniform.
Treating these four genera as one curation batch the way the smuts were would be
unjustified.

## *Leucosporidium* — real primary paper, but only partial-CDS deposits

**Maia TM, Lopes S, Almeida JMGCF, Rosa LH, Sampaio JP, Gonçalves P, Coelho MA
(2015) Evolution of Mating Systems in Basidiomycetes and the Genetic Architecture
Underlying Mating-Type Determination in the Yeast *Leucosporidium scottii*.
*Genetics* 201(1):75-89. doi:10.1534/genetics.115.177717. PMID 26178967,
PMC4566278 (open access).**

Directly characterizes *L. scottii*'s MAT system: **tetrapolar**, two physically
unlinked loci — a **multiallelic homeodomain (HD) locus** and a **biallelic
pheromone/receptor (P/R) locus** — the same architecture already curated for the
Ustilaginomycotina smuts (HD1/HD2 + pheromone/receptor), just independently
evolved/retained in a different subphylum. The paper argues this supports
tetrapolarity as the *ancestral* state for all Basidiomycota.

**GenBank deposits exist (`KR229960`-`KR229978`, ≥19 strains) but every one
checked is a partial-CDS population-survey amplicon**, e.g. `KR229978.1`
(*L. scottii* CBS 7673): `HD2` at `complement(<1..491)`, `HD1` at `733..>1323` —
both truncated at the sequenced fragment's own edges, Sanger-sequenced PCR
amplicons from an allele-frequency survey, not full-locus deposits. None seen is a
complete, independently-retranslatable CDS the way every record curated this
session has been. A companion whole-genome assembly (`MWVB00000000`) exists and
may carry complete gene models, but that was not checked in this pass — that
would be the productive next step, not the amplicon series.

**Not curated. Real primary source, but the readily available deposits are the
wrong shape for this project's curation pipeline** (which validates a record by
independently re-translating its own claimed CDS span — a `<1`/`>1323` partial
amplicon has no defined span to re-translate against).

## *Rhodotorula* (= *Rhodosporidium*) *toruloides* — richest material, two separate threads, neither immediately curatable

Two genuinely different data threads, both real:

**1. The classic pheromone-precursor genes.** Akada R, Minomi K, Kai J,
Yamashita I, Miyakawa T, Fukui S (1989) Multiple genes coding for precursors of
rhodotorucine A, a farnesyl peptide mating pheromone of the basidiomycetous yeast
*Rhodosporidium toruloides*. *Mol Cell Biol* 9(8):3491-3498.
doi:10.1128/mcb.9.8.3491-3498.1989. PMID 2571924. Three genes, `RHA1`/`RHA2`/`RHA3`,
each encoding a precursor with 3-5 tandem repeats of the rhodotorucine A pheromone
peptide (mating type A) — a different precursor architecture than the single-copy
`mfa1` pattern already curated (tandem-repeat polyprotein release, not one peptide
per gene). Real GenBank accessions exist (11 nuccore hits for "rhodotorucine",
including what are very likely the original 1989 deposits, e.g. `218048`,
`1004344`-`1004348`) but **their gene-level structure was not verified this pass**
— 1989-era deposits predate modern annotation conventions and would need the same
careful CDS/exon check every other record this session received before curating.

**2. Modern chromosome-scale genomes with a repeated MAT-locus signature.**
Searching `Rhodotorula toruloides` + mating/pheromone/STE3/homeodomain returned 91
nuccore hits, dominated by complete chromosome sequences (`AP0417xx` series) across
5 strains (NBRC10513, NBRC10512, JCM10296, JCM10295, JCM10049, JCM10021, JCM10020).
**Chromosome 3 (~1.3-1.6 Mb) and chromosome 14 (~0.77-1.06 Mb) recur across every
strain in this list** — a real, unverified-but-suggestive signal that one or both
carries the MAT locus, consistent with the receptor-precursor-only literature above
mentioning "mating type A" as if biallelic/bipolar-like, which would fit a
Microbotryomycetes-typical linked or partially-linked system. **This is an
inference from accession patterns, not a stated fact from any paper read this
pass** — confirming it would mean reading the annotation of one of those
chromosomes for HD/pheromone/receptor gene models directly.

A 2023 paper (Lopes DD et al., *J Ind Microbiol Biotechnol* 50(1):kuad040,
doi:10.1093/jimb/kuad040, PMID 37989723) develops a PCR mating-type marker across
19 strains but is a biotech/industrial-strain-characterization paper, not a
locus-architecture paper — useful corroborating context, not a curation source.

**Not curated. Two real, promising threads, neither run to ground: the 1989
pheromone genes need modern-annotation-standard verification of their GenBank
CDS structure, and the chromosome 3/14 signal needs someone to actually open one
chromosome's annotation and confirm/deny it carries HD or P/R genes.**

## *Cystobasidium* — genuine absence, not a search gap

Broadened the search from the specific gene-name query (0 hits) to the whole class
(`Cystobasidium OR Cystobasidiomycetes OR Cystobasidiales`, 111 hits) and read
every title in the top 20 by relevance. **None concern mating-type or MAT-locus
biology.** The literature that exists is species descriptions (*C. psychroaquaticum*,
*C. halotolerans*, *C. alpinum*, *C. tubakii*, etc.), general genome sequencing,
biotechnology applications (lipid production, bioremediation), and one clinical
case report (bloodstream infection). This is the same class of finding as `Abeta`
earlier this session: a real absence, not a failure to search hard enough.

**Not curated, and no obvious next step exists** — there is no primary paper to
follow up on. If this genus matters enough to pursue, the next move would be a
targeted search for *any* published MAT-locus work in Cystobasidiomycetes more
broadly (the class, not just this one genus), which was not attempted in this pass.

## *Puccinia* — real paper, real supplementary data, wrong shape for this pipeline

Confirms and extends what was found earlier this session: **Luo Z, McTaggart A,
Schwessinger B (2024) Genome biology and evolution of mating-type loci in four
cereal rust fungi. *PLoS Genet* 20(3):e1011207. doi:10.1371/journal.pgen.1011207.
PMID 38498573, PMC10977897 (open access).** Characterizes tetrapolar HD + PR loci
(on separate chromosomes) across *P. coronata* f. sp. *avenae*, *P. graminis* f. sp.
*tritici*, *P. triticina*, *P. striiformis* f. sp. *tritici* — HD multiallelic in
all four, PR biallelic in three and possibly multiallelic in *P. graminis*.

**Read the paper's own Data Availability statement this pass** (not available
previously): the gene-level data — "presumed CDS of reconstructed HD alleles, Pra
alleles" — is hosted at **Dryad, doi:10.5061/dryad.w0vt4b8zm**, alongside a GitHub
analysis-code repository (`github.com/ZhenyanLuo/codes-used-for-mating-type`), NOT
as individual GenBank gene accessions. The underlying genome assemblies are cited
under several BioProject accessions (`PRJNA39437`, `PRJNA396589`, `PRJNA39801`,
`PRJNA39803`, `PRJNA398546`, `PRJNA415866`, `PRJNA60743` appear in the full text).

**Not curated.** This is real, usable, well-characterized data — genuinely better
characterized than either *Leucosporidium* or *Rhodotorula*'s current state — but
curating from it means: downloading the Dryad dataset, understanding its CDS
coordinate system, and mapping those coordinates onto the correct BioProject genome
assembly accession per species, which is a different and larger task than every
other curation done this session (all of which started from one clean,
individually-deposited GenBank accession).

## Summary and recommendation

| genus | primary literature | gene-level data | curatable now? |
|---|---|---|---|
| *Leucosporidium* | yes, strong (Maia et al. 2015) | partial-CDS amplicons only; a WGS assembly exists, unchecked | no |
| *Rhodotorula* | yes, two threads (1989 classic + 2023 modern) | old accessions unverified; chromosome-scale signal unconfirmed | no |
| *Cystobasidium* | **none found** | none | no — nothing to build on |
| *Puccinia* | yes, strong (Luo et al. 2024) | yes, but Dryad-hosted, needs coordinate mapping | no, not without that extra step |

None of the four is a clean "propose one GenBank accession" case the way every
Ustilaginomycotina record this session was. All four need one more concrete
step before a curation pass makes sense, and the four steps are different:

1. **Leucosporidium**: check whether `MWVB00000000` (the *L. scottii* WGS
   assembly) has an annotated, complete gene model for HD1/HD2 at the same locus
   the amplicon series types.
2. **Rhodotorula**: open one strain's chromosome 3 and chromosome 14 annotation
   directly and check for HD/pheromone/receptor gene models; separately verify
   the 1989 `RHA1`/`RHA2`/`RHA3` deposits' CDS structure.
3. **Cystobasidium**: no path forward identified. Would need a from-scratch
   literature or genome search, not a follow-up on anything found here.
4. **Puccinia**: download the Dryad dataset (doi:10.5061/dryad.w0vt4b8zm) and
   determine how its CDS coordinates relate to the cited BioProject genome
   assemblies.

Recommend deciding, per genus, whether any of these four follow-ups is worth the
additional effort before curating — the same go/no-go this session used for
*Abeta* and the assembly-derived *Cunninghamella*/*Chaetocladium* loci.
