# Orbiliomycetes and Chaetothyriales: published MAT loci (literature pass, 2026-11-08)

Curator ruling 2026-11-07: curate Orbiliomycetes and Chaetothyriales first. This note records what the literature and NCBI actually support. Result: two verified candidate records, both Exophiala dermatitidis (Chaetothyriales); nothing verifiable for Orbiliomycetes.

Method: PubMed and nuccore/gene/protein/assembly lookups via NCBI E-utilities from the HPCC login node (the PubMed MCP connector needed authentication and was not used), plus the full text of PMC6417501 and PMC8465350. Search terms covered Arthrobotrys/Orbilia/Drechslerella/Dactylellina/nematode-trapping with mating/MAT1/MAT1-1-1/heterothallic, and Exophiala/Chaetothyriales/black yeast/Cladophialophora/Fonsecaea/Capronia/Knufia/Rhinocladiella with mating type and MAT1. The PubMed search tool is keyword-based; absence of hits is not proof that no paper exists.

## Papers reviewed

| PMID | DOI | What it establishes | What it does not |
|---|---|---|---|
| 30611834 (Metin et al. 2019, Fungal Genet Biol 124:29-38) | 10.1016/j.fgb.2018.12.011 | E. dermatitidis MAT locus. HMG gene MAT1-2-1 near SLA2 and APN2 in the genome of CBS 525.76 (= NIH/UT8656); MAT1-1 idiomorph 3544 bp with MAT1-1-1 (alpha-box) and MAT1-1-4; MAT1-2 idiomorph 3771 bp with MAT1-2-1 only; 74 isolates typed (35 MAT1-1 : 39 MAT1-2); MAT genes transcribed (RT-PCR); no virulence difference. One GenBank deposit, MH341450 (MAT1-1, strain CBS 132752). | No sexual cycle observed (heterothallism is inferred). No coordinates for the idiomorph boundaries. No accession for MAT1-2 sequence apart from the genome; no per-isolate accessions. Locus tags of the MAT1-2-1 gene are not given. |
| 28348446 (Teixeira et al. 2017, Stud Mycol 86:1-28) | 10.1016/j.simyco.2017.01.001 | Chaetothyriales comparative genomics (source of several genomes). Cited by Metin et al. for genome data. | Its MAT content could not be checked: PMC refused automated access and the abstract does not mention MAT. Not used as evidence for any record. |
| 34576814 (Zhou et al. 2021, Microorganisms 9:1919) | 10.3390/microorganisms9091919 | A. oligospora (Orbiliomycetes) isolates from China were typed for mat1-1 / mat1-2: 197 strains, 84 mat1-1 and 113 mat1-2 (about 1:1.3), taken by the authors as consistent with heterothallism. Primers MAT1F/R, MAT2F/R given. | No GenBank accession for any MAT sequence (accessions only for its, mapk, rpb2, sp, tef-1, tub). No gene structure, flank genes, coordinates or genome position. Uses mat1-1 / mat1-2 names, not MAT1-1-1 / MAT1-2-1. Not enough to build a record. |
| 21909256 (Yang et al. 2011, PLoS Pathog) and 24244185 (Meerupati et al. 2013, PLoS Genet) | 10.1371/journal.ppat.1002179 ; 10.1371/journal.pgen.1003909 | Genome papers for A. oligospora and nematode-trapping fungi (listed by search only). | Abstracts and full text were not read for MAT content; not used. |
| 24507587 (Drechslerella stenobrocha genome) | 10.1186/1471-2164-15-114 | Listed by search only. | Not read for MAT content. |
| 21148974 (Mostert et al. 2003, Mycologia) | 10.1080/15572536.2004.11833069 | Search hit; concerns Togninia/Phaeoacremonium (Calosphaeriales); Phaeomoniella is a Chaetothyriales name only in the abstract. | Not relevant to a Chaetothyriales MAT locus. |
| 29102684, 28098555, 34536100 | 10.1016/j.fgb.2017.10.008 ; 10.7554/eLife.20023 ; 10.1007/s00253-021-11455-z | Search hits (virulence review; A. oligospora odor mimicry; DdaSTE12 in D. dactyloides). | No MAT locus content. |

Also checked by nuccore text search: "Arthrobotrys/Orbiliomycetes MAT1", "Cladophialophora MAT1", "Fonsecaea MAT1", "Exophiala MAT1-1-1 / MAT1-2-1". Hits were MH341450 (Metin) and unrelated "CDK-activating kinase assembly factor MAT1" entries; no Orbiliomycetes or other Chaetothyriales MAT deposit was found.

## Records proposed (db/candidates/Ascomycota/)

1. `5970_cbs-132752_MAT_MAT1-1` : E. dermatitidis CBS 132752, GenBank MH341450.1. MAT1-1-4 (AXZ62320.1, 6 exons, complement 1761-2935) and MAT1-1-1 (AXZ62321.1, 3 exons, 3595-4885), coordinates read from the deposit's own annotation by efetch. `curate-db validate`: accession resolved, both genes 100% identity/100% coverage by independent translation.
2. `858893_nih-ut8656_MAT_MAT1-2` : E. dermatitidis NIH/UT8656 (= CBS 525.76), RefSeq GCF_000230625.1, scaffold NW_008751653.1. SLA2 (HMPREF1120_08859, XP_009161377.1, 289121-292377 +), MAT1-2-1 (HMPREF1120_08862, XP_009161380.1, 303111-304250 +, 379 aa, HMG-box_ROX1-like), APN2 (HMPREF1120_08863, XP_009161381.1, 305794-308458 -). Gene intervals from NCBI Gene (converted from 0-based). MAT1-2-1 validates at 100%. SLA2 (97.8% identity) and APN2 (3.6% coverage) give `warn` only because no exon models were entered and both genes have introns; this is a record limitation, not a mismatch. Identification of HMPREF1120_08862 as the MAT1-2-1 gene is by the HMG-box domain and the SLA2-APN2 position; the paper does not name the locus tag (the paper's text describes the locus as lying between APN2 and HMPREF1120_08861, consistent).

Tiers: locus_existence and idiomorph_assignment tier 1 (sequencing, 74-isolate PCR typing, RT-PCR). boundaries tier 2 (paper gives lengths only; record spans are annotation-derived). `mating_type.system` is `unknown` for both: nobody has observed mating.

All records are `needs_review`; nothing was accepted. Note the MAT1-1 record has no flank genes (SLA2 and APN2 are not annotated on MH341450). Between SLA2 and MAT1-2-1 in the genome are HMPREF1120_08860 and HMPREF1120_08861 (non-MAT, similar proteins occur on the MAT1-1 amplicon).

## Found but not verifiable

- A. oligospora mat1-1 / mat1-2 (Zhou 2021): typing data exist but no accession, structure or flank. NOT proposed.
- Chaetothyriales other than E. dermatitidis (Capronia, Cladophialophora, Fonsecaea, Phialophora, Knufia, Hortaea): no published MAT locus found. Metin et al. refer to other Chaetothyriales (C. epimyces and C. coronata, with the alpha-box gene back to back with MAT1-1-4) from the literature; that source was not identified or checked here.

## What is NOT known

- Whether any Orbiliomycetes species has a described MAT idiomorph structure (genes, orientation, flanks). The only published data found are mat1-1/mat1-2 typing in A. oligospora. A sequence-based proposal would need annotation work on the A. oligospora, A. flagrans and Drechslerella genomes (for example with the E. dermatitidis proteins as seeds) and would be tier 2 at best; that is a detection exercise, not a literature record.
- Whether the Ascomycota `MAT` family can describe these groups. The family `Ascomycota:MAT` (db/Ascomycota/order.yml) has `taxonomic_scope` leotiomyceta (taxid 716546), which covers Eurotiomycetes including Chaetothyriales, and its gene vocabulary already holds MAT1-1-1, MAT1-2-1 and MAT1-1-4 (`optional`, MAT1-1). The two E. dermatitidis records fit that vocabulary and `validate` accepted them. Orbiliomycetes (Orbiliomycetes is outside leotiomyceta, which excludes Orbiliomycetes) would be routed only by the Pezizomycotina fallback_subphylum rule; whether that routes correctly was not tested here.
- Whether the E. dermatitidis idiomorph boundaries match any detection span: the paper lists only lengths, so these records should stay out of any coordinate benchmark (set `excluded_from_coordinate_benchmark: true`).
- Teixeira 2017 MAT content, and the other genomes' MAT loci (E. dermatitidis is the only Chaetothyriales species here with a published structure).

## Detection implication

The Chaetothyriales uncalled rate (Exophiala 19/114 called) can now be tested against a real MAT1-2 and MAT1-1 reference once the curator accepts these two records; the alpha-box protein and the HMG protein are the seeds. Expect other Chaetothyriales genera to carry MAT1-1-4 near MAT1-1-1 if the E. dermatitidis structure is typical, but that is not documented here.
