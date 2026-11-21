# Chaetothyriales before/after with the two *E. dermatitidis* records accepted (#344)

Date 2026-11-20. The before/after the curator asked for in #344 before accepting `5970_cbs-132752_MAT_MAT1-1` (MAT1-1-4, MAT1-1-1; MH341450.1) and `858893_nih-ut8656_MAT_MAT1-2` (SLA2, MAT1-2-1, APN2; GCF_000230625.1). Nothing is accepted in the repository; this is a measurement on a copy of the db. Run record: `results/2026-11-20_exophiala_344/RUN.txt` (frozen v0.9.0 code, `--preset v0.9.1`, db = the v0.9.0 db plus the two records, 223 genomes of the Chaetothyriales order, job 32334032).

## Result
| | Genomes |
|---|---|
| Chaetothyriales genomes of the campaign | 223 (222 with a baseline report) |
| called before and after | 147 |
| **newly called** (uncalled in v0.9.0) | **58** |
| called before, uncalled after | 0 |
| uncalled before and after | 17 |

Newly called by genus: *Phialophora* 20 of 21 genomes, *Exophiala* 13, *Cladophialophora* 10, *Capronia* 7, *Fonsecaea* 5, *Rhinocladiella* 2, one *Herpotrichiellaceae* sp. Still uncalled: *Exophiala* 6, *Cyphellophora* 4, *Fonsecaea* 2, *Herpotrichiellaceae* 2, *Phialophora* 1, *Knufia* 1. For *Exophiala* (114 genomes) that is 108 called after and, from the 13 new ones, 95 before; the issue text says 19 of 114 are called today, which does not match this baseline (I did not find where 19 comes from).

Every other difference, listed in `exo344.diff.tsv` (593 rows):
- **91 loci appear** that were not in the v0.9.0 report and 33 loci disappear; the difference is loci of the 58 newly called genomes plus loci whose outer span changed (the report key includes the span, so a span change shows as one loci row removed and one added); the new loci are 25 `high` and 66 `medium`; idiomorph MAT1-1 50, MAT1-2 39, `undetermined` 2.
- **14 tier changes, all `medium` to `high`**, in already-called genomes (`exo344.diff.tsv`, rows `tier`).
- **2 idiomorph changes**:
  - `GCA_030247915.1_UCR_JES112` (*Neophaeococcomyces*): MAT1-2 to MAT1-1. After the change the MAT1-1-1 (35.06% identity) and MAT1-1-4 (39.78%, now `polished_single`; before 26.6%, `unpolished`) genes are better modelled and carry the call.
  - `GCA_056513255.1_ASM5651325v1` (*Epibryon*): `undetermined` to MAT1-2, also with tier `medium` to `high`.
- **Gene-model changes** in called loci: APN2 (114 rows), MAT1-1-1 (95), SLA2 (95), MAT1-1-4 (83), MAT1-2-1 (65); typically identity rises (for example SLA2 82% to 91%) and exon boundaries move a few bases to a few hundred, (the two records add references for SLA2, APN2, MAT1-1-1, MAT1-1-4 and MAT1-2-1; why each model moved was not traced).

## What this does and does not show
- The records mostly add calls: no genome lost a call, 58 gained one, 14 loci gained a tier, and two idiomorph labels changed. It does not show the new calls are correct: there is no independent evidence for the 58 (no FungiDB annotation of these genomes was compared, no RNA-seq, no read typing), and the single sequenced locus of one species (*E. dermatitidis*) is being used as a reference for a whole order. The metadata says nobody has observed mating (`system: unknown`).
- The medium-confidence calls dominate the new loci (66 of 91) and many have thin identity (30 to 45%); the two idiomorph changes both rest on genes at 35 to 47% identity.
- The 14 tier increases are in already-called genomes; I did not trace why each rose (the gene models changed in the same loci), so they should not be read as new evidence for the call.
- Not run: the 313 FungiDB genomes, the 510 panel genomes, and the other orders (the new records belong to `Ascomycota:MAT`, so a genome of another order that uses them could change). tblastn SEG was left at its default, while v0.9.0 used it for part of the campaign; a few span differences may come from that, not from the records.
- No `locus.gbk` or `locus.gff3` was built for the two records, because that needs an NCBI contact address. The detector needs only `metadata.yaml` and `proteins.faa`, but the accepted tree would carry those files.

## For the curator
Accepting is the curator's decision (`curate-db accept`, then `build-gff`). What I would look at before: the idiomorph change in `GCA_030247915.1_UCR_JES112` and the one in `GCA_056513255.1_ASM5651325v1`, and a sample of the newly called `medium` loci by an independent test (flank synteny and HMG-box or alpha-box content, as in the #341 check). Reject or defer is a defensible choice if the new calls cannot be checked.

## Addendum, 2026-11-20: domain and flank check of the new calls and the two idiomorph changes
Run with the #341 domain and flank code (`results/2026-11-20_exophiala_344/exo344_domain.*.tsv`, RUN.txt) on 93 loci: the 91 newly called loci and the two idiomorph-change loci. For each MAT1-1-1, MAT1-2-1 and MAT1-1-3 model the genome around it (+/-600 bp) was translated in six frames and searched with the alpha-box (PF04769) and HMG-box (PF00505, PF09011, PF14887) models at the gathering cutoffs. A domain hit shows the region encodes an alpha-box or an HMG-box protein, not that it is that MAT gene: an HMG-box paralog has the same domain, and the alpha-box is the more mating-specific of the two.

- **The two idiomorph changes: the call is not confirmed in either; each locus encodes both domains.**
  - `GCA_030247915.1_UCR_JES112` (now MAT1-1, high): MAT1-1-1 (35.1% identity) carries the alpha-box (PF04769, 83.4 bits) and the HMG-box ORF (PF00505, 74.1 bits; the MAT1-2-1 and MAT1-1-3 models sit on the same ORF) lies 1.8 kb away; APN2 and COX13 are 2 to 8 kb to the left. The earlier MAT1-2 call and the new MAT1-1 call both have a domain behind them.
  - `GCA_056513255.1_ASM5651325v1` (now MAT1-2, high; was undetermined): MAT1-2-1 (39.9%, `polished_agree`) carries the HMG-box (65.9 bits) and MAT1-1-1 (47.4%, `polished_disagree`) carries the alpha-box (82.4 bits); SLA2 is 2.5 kb from the MAT1-2-1 gene.
  - Reading: both loci contain an alpha-box gene and an HMG-box gene within a few kb, so a single MAT1-1 or MAT1-2 label cannot be called right or wrong from this evidence. As the alpha-box is the more specific domain and an HMG-box paralog is common, MAT1-1 is the more defensible reading for both, which would make the new `high` MAT1-2 call of the second genome the weak one; this is a reading, not a test (no RNA-seq, no read typing, no other isolate). Neither locus is labelled homothallic by the preset's `both_genes` policy (not examined why).
- **All 93 loci: 4 contain both domains** (the two above and two *Capronia* loci, `GCF_000585565.1` and `GCF_000585585.1`, which are `undetermined`). The other 89 have only the domain of the called idiomorph (alpha-box in 50, HMG-box in 39), which agrees with their label; a weak check, since the HMG-box in particular is not specific.
- **Flank genes:** 48 of the 93 loci have SLA2, APN2 or COX13 within 30 kb (26 of the 27 `high`, 22 of the 66 `medium`); 45 have none (44 `medium`, 1 `high`). The flank genes are the pipeline's own hits, so this is a consistency check, not independent evidence.
- **What it does not show:** that any of the 91 calls is the MAT gene (the HMG-box result is non-specific), or that the 45 loci without a flank are wrong; they are the weakest and the first to look at if a sample is wanted.
