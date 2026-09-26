# Saccharomyces MAT/HML/HMR under-calling: diagnosis

Data: `results/2026-09-23_saccharomyces_bar` (1,553 genomes, 3,101 reported loci). Code: commit 5b2fb12.
Scripts and intermediates are in `results/2026-09-24_sacc_diagnosis/` (`parse.py`, `a*.py`, `flank.py`, `asmstats.py`, `nobar.py`).

## Method

1. I parsed every `detection_report.yaml` and `evidence_diagnostics.jsonl` into `parsed.json`.
2. I computed contig count, N50 and total length from each input FASTA (`asmstats.tsv`).
3. **Cassette-site assay (`flank.py`, `flank_all.tsv`, `sites.json`).** I took S288C chrIII (NC_001135.5, from GCF_000146045.2_R64) and cut three windows per cassette: a 9 kb left flank, the cassette window (HML 11,000-15,500, MAT 198,000-202,500, HMR 292,000-295,500), and a 9 kb right flank. For each genome, I used blastn (>=90% identity) to find where the inner ends of the unique flanks land. When both flanks land on one contig, I measured the interval between them and counted the N bases in it. Status values:
   - `intact`: N = 0 and length within -600/+800 bp of S288C.
   - `N_gap`: the interval contains N.
   - `unresolvable`: a flank is missing, the flanks are on different contigs, or they are more than 50 kb apart.
   The assay only works for *S. cerevisiae*, because other species' flanks are below 90% identity. I ran it on all 1,553 genomes.
4. **Reproduction.** I made a frozen worktree at 5b2fb12 on /bigdata and reran GCA_001671115.1 (taxid 4932). The report has `suppressed_unpolished: 32`, and the result matches the bar run. I then ran it again with `min_polished_genes=0`, patched at runtime with no source edit (`nobar.py`). I removed the worktree afterwards.

## Headline: the pipeline calls every intact cassette

*S. cerevisiae* only: 1,309 genomes, 3,927 cassette sites.

| site | intact, called | intact, not called | N_gap, called | N_gap, not called | unresolvable |
|---|---|---|---|---|---|
| HML | 237 | 0 | 235 | 164 | 672 |
| MAT | 375 | 0 | 465 | 93 | 376 |
| HMR | 248 | 0 | 477 | 39 | 506 |
| total | **860** | **0** | 1,177 | **296** | **1,554** |

Labels at intact sites:

| site | MATalpha | MATa |
|---|---|---|
| HML | 231 | 6 |
| HMR | 2 | 246 |
| MAT | 262 | 113 |

The 8 silent-cassette labels that run against expectation are not verified either way.

Genomes with all three sites intact: 184. Of these, 184/184 have 3 or more loci, and 181/184 have both a MATa and a MATalpha call.
Loci per genome by number of intact sites:

| intact sites | loci per genome |
|---|---|
| 3 | 3:179, 4:4, 5:1 |
| 2 | 2:18, 3:52, 4:1 |
| 1 | 1:44, 2:97, 3:25 |
| 0 | 0:109, 1:153, 2:410, 3:208, 4:7, 5:1 |

When all three cassettes sit on one contig (470 genomes, all species), the label pattern ordered HML-MAT-HMR is:

| pattern | genomes |
|---|---|
| alpha-alpha-a | 297 |
| alpha-a-a | 88 |
| alpha-undetermined-a | 33 |
| undetermined-alpha-a | 20 |
| a-alpha-a | 12 |
| other patterns | 20 |

HMR is labelled `a` in 454 of these 470 genomes.

## H1: the >=2 modelled-gene bar withholds real cassettes. Verdict: NOT SUPPORTED as a major cause.

- `suppressed_unpolished` totals 48,179 (median about 30 per genome). These are almost all background homeodomain clusters. Of the 4,654 admitted clusters that were not reported, 4,629 have a best identity below 50%.
- Clusters at >=90% identity with 2 or more genes that were not reported: 11 across 1,553 genomes (6 admitted, 5 not admitted).
- Uncalled N-gap sites (296): 1 is a 2-gene cluster withheld by the bar, 90 are single-gene >=90% clusters that fail the min_hits=2 evidence floor (not the bar), and 205 have no >=90% hit at all.
- All 169 zero-loci genomes carry a "best cluster carried 0/1 modelled gene(s)" reason (119 with 1 gene, 50 with 0). The reason is literally true, but in 141/169 of these genomes no cluster reaches 90% identity. The withheld "best cluster" is background.
- a1 being the only a-specific gene does not cause the loss. At every a cassette, MATA2 and MATALPHA2 cross-hit the shared X region, so the cassette has 3-4 gene names.
- Example: GCA_001671115.1 CP008442.1:199,580-200,484 (MAT). MATALPHA1 is polished at 99.2%. MATALPHA2 is unpolished because the window contains 902 N. The locus is withheld with 1 modelled gene. With the bar off it is reported with `polished_genes: 1`. The other 31 withheld loci in that genome are all below 90% identity.

## H2: silent-cassette confounding or cross-match. Verdict: NOT SUPPORTED for detection; PARTLY SUPPORTED for labels and confidence.

- **No merging.** The largest reported span is 27.7 kb. Cassettes are 94-187 kb apart and the gap setting is 25 kb. Spans over 5 kb (43 loci) come from low-identity background hits, not from two cassettes.
- **Labels at intact sites are correct.** See the table above.
- **High-identity idiomorph resolutions are only MATA2 vs MATALPHA2:** 416 events in 193 genomes. That pair is non-informative, so the resolutions do not change the label.
- **`undetermined` loci (188, MATA2+MATALPHA2 only)** occur at 0 intact sites, 70 N-gap sites and 118 unresolvable sites. The Y-region gene is missing from the assembly. This is not a cross-match.
- **Confidence artifact (code).** `expected_genes_for_idiomorph` narrows the roster using every found gene, including the `idiomorph_informative: false` genes MATA2 and MATALPHA2. MATA2 hits every alpha cassette, so the roster is never narrowed. Measured: alpha {ALPHA1, ALPHA2} scores high (355), but alpha {MATA2, ALPHA1, ALPHA2} scores medium (1,028). MATa cassettes score high only with all 4 genes (927); with 3 genes they are medium (347). This changes confidence, not detection.
- **Holdout mislabel (curation).** Measured in the holdout record run: with S288C HMRa withheld, the only a1 reference is *K. lactis* (40% identity, vote 47). The 95%-identity MATALPHA1 fragment at HMR outvotes it (vote 82). At species and genus radius, the remaining alpha references match at 52-65% and HMR is not found. The cause is that only one Saccharomyces MATa reference exists.

## H3: no flanking genes. Verdict: NOT SUPPORTED as a cause of missed calls.

- sla2 is `optional` and cha1 is excluded from search. sla2 appears in 9 of 3,101 loci. Flanks are not in the fraction denominator and do not cap the tier.
- `idiomorph_gene_only` (3,092 loci) is a label that follows from this design. It does not signal a failure.
- All 860 intact sites are called with no flank.
- Inference: flanking DNA, not flanking genes, would be needed to name which cassette (HML/MAT/HMR) a call is. The order.yml comment already records this.

## H4: assembly limits. Verdict: SUPPORTED. This is the dominant cause.

Only 184/1,309 *S. cerevisiae* assemblies (14%) have all three cassettes intact.

| N50 bin | genomes | mean loci | frac 0 loci | frac 3+ loci | frac both a and alpha |
|---|---|---|---|---|---|
| <20 kb | 47 | 0.96 | 0.34 | 0.09 | 0.15 |
| 20-50 kb | 43 | 0.81 | 0.49 | 0.07 | 0.19 |
| 50-100 kb | 115 | 1.23 | 0.35 | 0.06 | 0.37 |
| 100-300 kb | 253 | 1.42 | 0.25 | 0.08 | 0.44 |
| 300-700 kb | 61 | 1.85 | 0.20 | 0.34 | 0.54 |
| >700 kb | 1,034 | 2.33 | 0.02 | 0.47 | 0.80 |

Median contig count by loci per genome:

| loci | median contigs |
|---|---|
| 0 | 244 |
| 1 | 50 |
| 2 | 33 |
| 3 | 17 |

Three failure modes, with examples:

1. **N-gap inside a cassette** (1,473 sites; 296 not called). Example: GCA_001592655.1 CP014734.1. At MAT, 198,594-198,895 is exactly 300 N in place of the Y region, and MAT is not called. HML (14,160) and HMR (293,146) each also carry a 300-N gap but keep enough sequence to be called.
2. **Cassette sequence absent, flanks joined directly.** Example: GCA_000769245.1 CM002937.1. The MAT flanks meet at 181,920/181,922 with 1 N, and about 2.3 kb of the cassette is missing. The genome has 0 loci and no hit at or above 60%. The same holds for NCYC1159 (GCA_947343775.1). In 141 of 169 zero-loci genomes, no MAT protein matches at 90% or more.
3. **Contig break at a cassette** (1,554 unresolvable sites; 616 of them still carry a call on a fragment). There are 252 unreported single-gene clusters at >=90% identity in 200 genomes, often at a contig end. Example: GCA_000166955.1 AABY01000077.1:2-394 (93.9%).

Two-locus genomes on one contig, by the spacing between the two loci:

| missing cassette | spacing | genomes |
|---|---|---|
| HML | MAT-HMR, ~94 kb | 300 |
| MAT | HML-HMR, ~280 kb | 64 |
| HMR | HML-MAT, ~187 kb | 32 |

Of the 64 genomes missing MAT, 44 have no evidence cluster at the MAT position.

**Diploid a/alpha collapse: NOT TESTED.** This needs reads or phased assemblies. Inference, consistent with the data but not proven: the MAT Y region is a 2-copy repeat with HML or HMR, and a heterozygous MAT adds a bubble. Either could explain why N-gaps land on the cassettes.

## Fixes ranked by expected recovered calls

| rank | fix | type | measured ceiling |
|---|---|---|---|
| 1 | Accept that short-read assemblies lose cassettes; report per-genome site status (intact / N-gap / broken) so a missing call is not read as absence | inherent / reporting | about 205 N-gap sites plus about 938 broken sites without calls cannot be recovered from the sequence |
| 2 | Report single-gene >=90%-identity polished clusters as a partial-cassette class, instead of dropping them at min_hits=2 | code | up to 252 clusters in 200 genomes (90 at resolvable S. cerevisiae N-gap sites); not deduplicated; which gene each carries is not measured |
| 3 | Exclude `idiomorph_informative: false` genes from roster narrowing | code | 0 new calls; changes confidence for about 1,400 loci (inference: 3-gene alpha cassettes become high) |
| 4 | Let the bar accept 1 modelled gene when it is >=90% identity | code | at most 11 clusters |
| 5 | Add a second Saccharomyces MATa reference and non-cerevisiae Saccharomyces records | curation | fixes the holdout HMR mislabel and miss; no measured effect on the full-db run |

Not measured: whether non-cerevisiae genomes (244; 60 with zero loci) lose cassettes to assembly or to reference distance. The flank assay cannot resolve their sites. The median best-gene identity of their reported loci is 98-100%.
