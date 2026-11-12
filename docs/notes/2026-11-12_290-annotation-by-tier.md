# #290: FungiDB annotation check by tier (v0.9.0 run against v0.8.1)

Date 2026-11-12. The annotation is an independent cross-check, not the truth: `novel` means no FungiDB gene overlaps the call (our model may be better, or the annotation may miss the gene), not that the call is wrong.

Inputs: the 313 FungiDB genomes (release 71) run on the frozen v0.8.1 code (`FungiDB_MATpredict`, `--preset v0.8.1`) and on the frozen v0.9.0 code (`7d73be5`, `--preset v0.9.1 --withheld-detail-policy on`, `results/2026-11-09_fdb_v090`); verdicts from `FungiDB_MATpredict/scripts/compare_matpredict_to_gff.py` and `annotation_view.py` run over both; join with `scripts/annotation_by_tier.py`. Results: `results/2026-11-12_290_annotation_check/` (HPCC).

## Agreement by tier (MAT genes in called loci)

| Run | Tier | Genes | Annotated | Partial | Novel |
|---|---|---|---|---|---|
| v0.8.1 | high | 905 | 833 | 15 | 57 |
| v0.8.1 | medium | 329 | 293 | 10 | 26 |
| v0.9.0 | high | 905 | 832 | 15 | 58 |
| v0.9.0 | medium | 344 | 306 | 11 | 27 |

No `low` locus has a gene in these genomes. Every verdict row joins to a report locus.

## v0.9.0 against v0.8.1

- **No verdict changed** on the 1,234 genes present in both runs; none was lost.
- **15 genes only in v0.9.0**, from calls that are new in v0.9.0: 12 in *Bipolaris*, *Curvularia*, *Alternaria* and *Cochliobolus*-type genomes (the `--preset v0.8.1` fix; GAP1, BGL1, MAT1-1-1, MAT1-1-2, MAT1-2-1; all `annotated`) and 3 *Pneumocystis* `matMi` (the corrected record, #302): one `partial`, two `novel`, because the FungiDB models do not cover that ORF.
- **7 genes change tier**: 5 `high` to `medium` (the weak-flank and single-region caps) and 2 `medium` to `high` (*Pneumocystis*, #302). The weak-fragment cap is not in the v0.9.0 preset, so its trigger loci (plan step 3) are not examined here.

## Plan step 3: the weak-fragment cap's trigger loci (not in any preset)

The repo's own rule (`weak_fragment_cap.evaluate`, unchanged) run over the 193 `Ascomycota:MAT` loci of the v0.9.0 FungiDB run with a called idiomorph matches **3** (a one-off script in `results/2026-11-12_290_annotation_check/cap_triggers.py`, not in the repo):

| Genome | Called | Tier | Strongest core gene (the cap's pick) | Annotation at the called idiomorph's core genes |
|---|---|---|---|---|
| *A. terreus* NIH2624 | MAT1-1 | high | MAT1-2-1, polished_single, 44.3 bits: `novel` | MAT1-1-1 `annotated` |
| *C. geophilum* 1.58 | MAT1-1 | medium | MAT1-2-1, polished_single, 45.4 bits: `novel` | MAT1-1-1 `annotated` |
| *S. schenckii* 1099-18 | MAT1-1 | high | MAT1-2-1, polished_single, 37.7 bits: `novel` | MAT1-1-1, MAT1-1-2, MAT1-1-3 all `annotated` |

In all three the cap would pick a thin MAT1-2-1 fragment, a gene of the *other* idiomorph, as the "strongest core gene", while the called idiomorph (MAT1-1) is supported by annotated MAT1-1 genes. So on these genomes the cap would lower the tier of a call the annotation supports. This is three loci in one annotated set, not a measurement of the rule, but it is a reason to read the rule's "strongest core gene" as the strongest gene of the called idiomorph before it is considered for a preset.

## The `novel` and `partial` genes (v0.9.0)

111 genes (85 novel, 26 partial): 107 are `core_MAT` genes (of 689 core genes, 16%) and 4 are flanking genes. The commonest are short ORFs that annotation pipelines often miss: MAT1-2-1 (11), MAT1-1-3 (9), MFalpha (8), MTLA2 (7), RPL39 (7), *Pneumocystis* matMi (6), MAT1-1-1 (5), MAT1-1-2 (4), MAT1-2-4 (4). By family: Ascomycota:MAT 47, Basidiomycota:MAT 23, Ascomycota:MTL 13, Ascomycota:PM 7. They have not been examined genome by genome.

### The 13 long, well-modelled genes with no annotation (candidates for annotation gaps)

Of the 111, 13 are at least 150 aa, polished by both tools (`polished_agree`) and have identity 50-100% at about 100% coverage, yet FungiDB has no overlapping gene (11) or only one on the opposite strand (2). These are the cases most likely to be gaps in the FungiDB annotation, not errors in our model. They are candidates for a curator's look, not a ruling:

| Genome | Gene (tier) | aa, identity | FungiDB |
|---|---|---|---|
| *F. odoratissimum* race 4 | MAT1-1-1 (medium), MAT1-1-2 (medium) | 373 aa 98.9%; 407 aa 81.7% | none |
| *A. lentulus* IFM54703 | MAT1-2-4 (high) | 278 aa 95.5% | none |
| *C. auris* B11245 | MTLA1 (high) | 202 aa 98.8% | none |
| *D. fabryi* CBS789 | MTLA1, MTLA2 (high) | 203 aa 79.1%; 273 aa 88.4% | none |
| *C. tropicalis* MYA3404 | MTLA2 (high) | 220 aa 50.0% | none |
| *P. jirovecii* SE8 | matMc (high) | 251 aa 100% | none |
| *H. vastatrix* race XXXIII | bE, bW (high) | 556 aa 100%; 530 aa 99.2% | none |
| *R. delemar* RA99-880 | sexP (high) | 313 aa 93.3% | none |
| *C. auris* B8441 | MTLA1 (high) | 206 aa 100% | B9J08_001456, 8% overlap, opposite strand |
| *C. auris* B11220 | MTLalpha1 (high) | 172 aa 99.4% | CJI96_0001190, 0% overlap, opposite strand |

*C. tropicalis* MTLA2 at 50% identity is the weakest of these and the one to read first.

### Do the 47 short, weak genes (class D) carry a call?

Class D is the largest group of the 111: genes under about 150 aa that are single-tool, unpolished or tools-disagree (rule-based class, `novel_partial_review_v090_classed.tsv`). For each I looked at the other core genes of the same locus (`classD_support_v090.tsv`, HPCC):

| Support in the locus | Genes | Tier of the locus |
|---|---|---|
| another core gene is annotated by FungiDB | 38 | 30 high, 8 medium |
| other core genes exist, none annotated | 8 | medium |
| the gene is the locus's only core gene | 1 | medium |

So 38 of 47 sit in loci whose call an annotated core gene already supports; the short gene is a secondary piece. The other 9 genes are in 5 loci, **all `medium`**, where nothing independent supports the call: *B. graminis tritici* THUN-12 (MAT1-1; MAT1-1-1 128 aa 40%, MAT1-1-3 96 aa 25%, MAT1-2-1 45 aa), *C. immunda* (sample `CimmundaCBS83496`) (called MAT1-2, with a MAT1-1-3 piece), *F. odoratissimum* race 4 (MAT1-1, with a MAT1-2-7 piece), *S. ludwigii* UTAD17 (MATa; MATA1 86 aa, MATA2 86 aa) and *M. guilliermondii* ATCC 6260 (MTLA2 84 aa, the only core gene, tools disagree). No `high` locus in class D lacks an annotated core gene. This is a triage of where the evidence is thin, not a ruling on any call; the genes of the other idiomorph inside a call (for example MAT1-1-3 in a MAT1-2 call) are worth a curator's look.

## What this does and does not show

It supports the v0.9.0 changes at gene level on these 313 genomes: the tier and label rules moved tiers but no gene-level agreement. It says nothing about genomes outside FungiDB, and the disagreements have not been looked at one by one; the list is `results/2026-11-12_290_annotation_check/v090_rows.tsv` (HPCC). The one-by-one look at the 111 disagreements and step 5 stay open on #290.
