# Pheromone-receptor arrays in Agaricomycetes: what they are, and whether they change what is called
Status: open (measurement and proposals only; no pipeline change, nothing merged)

## Question
1. How are STE3-like (PF02076) receptor loci organised in Agaricomycete genomes (singletons, tandem arrays), per order?
2. How do the 1,276 v0.6.0 `PR` calls relate to the arrays? Which arrays carry a strict-CAAX ORF but no call?
3. How much of the CAAX flag and of the STE20 co-location is chance, contig length or genome effects?
4. Should arrays change what is reported or called, and how could that be validated?

## Data and method
- Genomes: the v0.6.0 Basidiomycota run (`results/2026-10-03_basidiomycota_v060/`), class Agaricomycetes, 1,854 genomes. **The HPCC scan (`array_scan.py`, 62 array tasks, about 12 min each on 8 CPUs) had finished for 1,853 genomes (one of 500 Mb or more was skipped); the earlier tables used a partial pull of 248 genomes (184 quality-pass). All numbers here use the complete scan.** Archive of the per-genome scan output: `results/2026-10-06_agaricomycetes_pr_arrays/scan_out_full.tar.xz`.
- Quality filter: BUSCO complete >= 70 (BFD `busco_genome`), N50 >= 20 kb, contigs <= 5,000 (removes 5 more): 1,287 genomes, 705 species (binomial). Genomes of one species are not independent: intervals below are percentile bootstraps resampling species (2,000 draws), and key tests are repeated on one genome per species (best BUSCO then N50, 705 genomes).
- STE3-like loci: miniprot of 1,094 STE3 queries (BFD hits and curated receptors, `ste3_all.faa`), alignments covering >= 50% of the query and <= 8 kb, merged per strand (PR #33 `scan_genome.py`). 12,518 loci in the 1,287 genomes; 3.1% overlap their neighbour on the opposite strand.
- Array: single linkage, same contig, gap <= 50 kb between locus ends. Flag: a strict-CAAX ORF (`C[VI][IV][AVMG]` before a stop, Met 20-130 codons upstream, six-frame scan, no annotation) within 10 kb of any member. No label, PR-only call or CAAX status is used as ground truth anywhere in this note.
- Region genes: miniprot of 24 region queries (`region_queries.faa`: 7 STE20 from Cryptococcus and Sporidiobolales MAT records, none from Agaricomycetes, plus HD and MIPBF queries); the top-scoring hit per genome is "the" STE20 (or MIPBF); HD is the pipeline's called HD locus.
- Code (`results/2026-10-06_agaricomycetes_pr_arrays/`): `analyze.py` (previous agent, bug fixed: `x.array` on a Series), `array_verify.py`, `array_verify2.py`, `array_verify3.py` (my independent re-derivation and controls); outputs `v_*.tsv`, `v2_*.tsv`, `v_summary.txt`. Contig lengths (`contig_lengths.sh`) were needed for contig-matched nulls.

## Verification of the earlier numbers
| Earlier figure | Result |
|---|---|
| 184 of 1,854 genomes scanned | True of the partial pull only. Scan complete: 1,853 scanned, 1,287 pass quality. All conclusions below use the full set. |
| Median 9 STE3-like loci per genome, 94.6% with >= 5 | Reproduced on the full set: median 9 (IQR 7-11), 95.5% (93.7-96.9) with >= 5. Max 136 (a Russulales genome). |
| 1,242 arrays, 81% singletons | Subset values. Full: 9,436 arrays, 80.2% singletons. My independent clustering gives the identical partition to `analyze.py` (12,518 loci). |
| 129 PR calls (103 unverified), median call covers 2 loci | Subset. Full: 1,039 calls (944 unverified, 95 verified by other evidence), median 2 loci per call. |
| 213 flagged arrays, 91 uncovered | Subset. Full: 1,633 flagged, 627 uncovered. |
| 1,271 withheld PR loci, 1,007 by modelled-gene bar | Right as a count of withheld loci, misleading for the question. Among the 627 flagged-uncovered arrays, 437 hold a locus withheld at the **fraction floor** and only 15 at the modelled-gene bar; the modelled-gene bar (907 of 910 withheld unflagged arrays) acts on arrays with no CAAX ORF. |
| "13 expected by chance" among the 91 uncovered | **Wrong method.** It multiplied the order's locus-level chance share of all flagged loci (13-14%) by the uncovered arrays. The uncovered flagged arrays are the residual after the calls took the real receptors, so chance is concentrated there. Same shortcut on the full data gives 85; the array-level null below gives **230 of 627 (37%, 32-43)**. |
| STE20 within 250 kb of 83/213 flagged (39%) vs 53/1,029 unflagged (5%) | Same pattern on the full set, 550/1,631 (33.7%) vs 328/7,763 (4.2%); survives the controls (section 4). |
| HD and receptors almost never on the same contig | Confirmed and stronger: 60 of 8,057 arrays share a contig with an HD call where contig length predicts 315 (flagged: 7 against 50). Receptor arrays are depleted on HD contigs, as expected for unlinked A and B loci. |

How the chance figures are computed (new, `array_verify.py` section C): for each array, the probability that a window of the array's span, placed uniformly at random **on the same contig**, has a strict-CAAX ORF within 10 kb (analytic union of intervals; CAAX ORFs within 10 kb of any STE3 locus are removed from the background). Summing over arrays gives the expected number of flagged arrays; share = expected/observed. This fixes the unit (array, not locus) and controls for contig length and contig-level CAAX density. The pipeline's own per-genome figure (`p_random`, 500 locus-sized windows at random over contigs longer than 21 kb) gives 13.0% at locus level; the array-level figure is 18.2% (15.8-20.9). Limit: it does not capture local CAAX density (gene-dense or repeat-poor regions), which would need a model of the genome's gene structure.

## Results
### 1. Sample (`v_sample.tsv`)
| Order | Genomes (all) | Species (all) | Quality + contigs | Species | BUSCO median | N50 median kb |
|---|---|---|---|---|---|---|
| Agaricales | 824 | 385 | 676 | 313 | 99.0 | 273 |
| Boletales | 509 | 347 | 182 | 143 | 98.3 | 60 |
| Polyporales | 198 | 113 | 175 | 107 | 99.2 | 858 |
| Cantharellales | 99 | 18 | 77 | 14 | 97.1 | 346 |
| Russulales | 92 | 72 | 76 | 60 | 98.7 | 479 |
| Hymenochaetales | 55 | 30 | 47 | 29 | 99.0 | 2,338 |
| Auriculariales | 20 | 8 | 14 | 8 | 97.7 | 1,858 |
| other orders (12) | 57 | 38 | 40 | 31 | 98.7 | 640 |
| All | 1,854 | 1,011 | 1,287 | 705 | 98.8 | 327 |

Widening the scan was possible and is done: it is the whole BFD Agaricomycete set. Limits remain: 14 species in Cantharellales (mostly *Rhizoctonia/Ceratobasidium*), 8 in Auriculariales, strain-heavy sampling (1,287 genomes are 705 species), and the 12 other orders hold 40 genomes. The 1,287 genomes are 705 species (1.8 genomes per species on average), which is why the species bootstrap is used.

### 2. Organisation of STE3-like loci (`v2_order_organisation.tsv`, `v2_array_sizes_by_order.tsv`, 95% species bootstrap)
| Order | Species | Loci per genome, median (IQR) | Arrays per genome | Genomes with an array of >= 2 loci | Genomes with a flagged array | Genomes with a PR call |
|---|---|---|---|---|---|---|
| Agaricales | 313 | 10 (8-12) | 7.7 (7.1-8.3) | 91% (85-95) | 89% (82-93) | 63% (49-74) |
| Boletales | 143 | 7 (6-9) | 6.6 (6.1-7.2) | 88% (83-93) | 49% (40-58) | 27% (18-36) |
| Polyporales | 107 | 8 (7-10) | 5.5 (5.2-5.8) | 96% (93-99) | 93% (88-96) | 89% (82-94) |
| Cantharellales | 14 | 7 (6-10) | 7.9 (6.2-12.0) | 45% (24-79) | 86% (67-94) | 87% (68-95) |
| Russulales | 60 | 9 (6-11) | 9.4 (7.0-12.8) | 71% (58-82) | 57% (44-68) | 34% (21-47) |
| Hymenochaetales | 29 | 6 (5-7) | 7.4 (5.9-9.4) | 34% (17-50) | 13% (4-26) | 11% (3-24) |
| Auriculariales | 8 | 6 (5-8) | 4.6 (3.7-5.5) | 86% (58-100) | 86% (58-100) | 86% (58-100) |
| other orders | 31 | 8 (7-11) | 7.6 (6.5-9.1) | 83% (67-95) | 73% (56-86) | 73% (58-86) |
| All | 705 | 9 (7-11) | 7.3 (6.9-7.8) | 85% (80-89) | 78% (73-82) | 60% (52-67) |

Array sizes (all orders): 7,572 singletons (80.2%), 1,164 of size 2, 389 of size 3, 311 of size 4 to 10. 595 of 1,287 genomes have at least one array of >= 3 loci. Median span of arrays with >= 2 loci 20 kb.
- Sensitivity to the gap (`v_gap_sensitivity.tsv`): singletons 84.8% (10 kb), 81.6% (25), 80.2% (50), 79.3% (100), 78.0% (200 kb); arrays per genome 8.0 to 7.0. The headline (about 80% singletons, about 7 arrays per genome) does not depend on the choice.
- Is clustering real? Placing the same loci at random on their own contigs (`v_clustering_vs_null.tsv`), the share of loci in arrays of >= 2 is 22.6% at 50 kb (observed 39.5%) and 12.4% at 10 kb (observed 30.3%). So **excess clustering is real (2.4 times at 10 kb) but 57% of the loci in 50-kb arrays are expected from contig fragmentation alone**. Arrays of size 2 are weak evidence of tandem organisation; arrays of 3 or more are not explained by it. Hymenochaetales (long contigs) and Cantharellales (short arrays) are the extremes.
- 60% of loci (7,572 of 12,518) are singletons, so "arrays" describe a minority of loci (about 40%) that is enriched for CAAX-flagged and called receptors.

### 3. CAAX flag, chance and sensitivity (`v_array_chance_by_order.tsv`, `v_flag_definition_sensitivity.tsv`, `v_flag_by_array_size.tsv`)
| Order | Arrays | Flagged | Expected by chance | Chance share (95% species bootstrap) |
|---|---|---|---|---|
| Agaricales | 5,232 | 1,035 | 144.6 | 0.14 (0.12-0.17) |
| Boletales | 1,201 | 136 | 31.6 | 0.23 (0.18-0.31) |
| Polyporales | 964 | 227 | 47.8 | 0.21 (0.17-0.25) |
| Cantharellales | 610 | 100 | 27.1 | 0.27 (0.20-0.62) |
| Russulales | 716 | 65 | 21.8 | 0.34 (0.25-0.46) |
| Hymenochaetales | 346 | 13 | 10.8 | 0.83 (0.49-2.3): no detectable excess |
| Auriculariales | 64 | 13 | 3.2 | 0.25 (0.18-0.36) |
| other orders | 303 | 44 | 11.0 | 0.25 (0.17-0.35) |
| All | 9,436 | 1,633 | 298.0 | 0.18 (0.16-0.21); one genome per species 0.19 (0.18-0.21) |

- By array size: flagged 628 of 7,572 singletons (8.3%; 184 expected, share 0.29), 454 of 1,164 size 2 (39%; 54 expected), 294 of 389 size 3 (76%), 257 of 311 size >= 4 (83%). Singleton flags are 29% chance-expected (628 flagged, 184 expected, 444 above chance); larger arrays are flagged far above chance.
- Flag definition (all arrays; chance share): CAAX within 5 kb 0.16 (0.13-0.21), 10 kb 0.18, 20 kb 0.22 (0.20-0.25), 50 kb 0.27 (0.25-0.30). **At least two distinct strict-CAAX ORFs within 10 kb: 787 arrays flagged, 33 expected, chance share 0.04 (0.03-0.06)**, and this holds in singletons (179 flagged, 10.7 expected). A wider window adds more chance than signal. A tblastn precursor-homology hit (E <= 1) within 10 kb alone is weak (764 arrays, 420 expected); with a CAAX ORF the joint chance falls to about 7% of uncovered flagged arrays (assumes independence of the two evidences, which probably understates chance because homology hits are often the same ORF).
- Hymenochaetales: the strict-CAAX flag shows no excess in this order (13 flagged, 10.8 expected). The rule cannot be relied on there; the 13 calls rest on it.

### 4. STE20 co-location survives the controls, with a stated limit (`v_STE20_*.tsv`, `v_MIPBF_*.tsv`, `v_region_by_threshold.tsv`)
Raw: within 250 kb of the top STE20 hit, 550/1,631 flagged arrays (33.7%) and 328/7,763 unflagged (4.2%).
- Contig length: flagged arrays lie on shorter contigs (median 0.23 Mb against 0.32 Mb), so contig length is not what produces the excess. A random position in the genome would be within 250 kb of STE20 for 0.6% of arrays; observed 33.7% (flagged) and 4.2% (unflagged). A random position on the same contig: expected 329 against 550 observed (flagged) and 206 against 328 (unflagged). Even conditional on being on the STE20 contig, flagged arrays are nearer (91% within 250 kb against 63%).
- Within genome (Mantel-Haenszel over the 1,003 genomes with both kinds of array and an STE20 hit): flagged 550/1,628 near, unflagged 147/6,117; odds ratio 23.9, species-bootstrap 95% 17.6-33.1; permutation of the flag labels within genome p = 0.0002 (floor of 5,000 permutations; the permutation mean is 163 against 550 observed). One genome per species: 19.5 (14.8-26.9). Stratified by array size as well: 11.0; by array size and contig-length quartile: 13.4 (537 informative strata, flagged 148/684, unflagged 42/1,881). Distance 100 or 500 kb: 22.4 and 21.2.
- By array size (flagged against unflagged, within 250 kb): size 1 19.0% against 3.3%; size 2 30.5% against 10.1%; size 3 46.6% against 20.4%; size >= 4 60.7% against 18.5%. Among flagged arrays at least two CAAX ORFs go with a higher rate within every size class (26.8% against 15.8%, 39.4% against 23.2%, 56.4% against 44.2%).
- By order the flagged rate exceeds the unflagged rate everywhere (Agaricales 27.1 against 3.9, Boletales 19.1 against 1.6, Polyporales 59.9 against 4.2, Cantharellales 56.0 against 3.7, Russulales 33.8 against 5.6, Hymenochaetales 38.5 against 13.7 on 13 arrays, other orders 29.5 against 6.6). Leaving out any one order leaves the odds ratio between 19.3 and 28.9. 532 genomes, 296 species and every order group contribute, one flagged array per genome on average, so no genome or species drives it.
- A second MAT-linked gene is not enriched: MIPBF (the A-locus flank) lies within 250 kb of 0 of 1,630 flagged and 11 of 7,782 unflagged arrays where the genome-position null expects 13.5 and 56. Depletion, not enrichment, so the STE20 effect is not a generic proximity-to-any-query artefact, but the null is also imperfect (arrays avoid MIPBF contigs, as they avoid HD contigs), so absolute excess figures should be read as upper bounds.
- Unflagged arrays still sit near STE20 more often than a random position (4.2% against 0.6%, 7 times), consistent with either unflagged MAT-linked receptors (CAAX misses) or with non-uniform placement of STE3 loci; the data do not separate the two.
- Not independent evidence of function. STE20 was chosen from known MAT linkage (Cryptococcus and Sporidiobolales PR regions); the result is one gene; STE20 has paralogs and the "top hit" per genome is a heuristic. It shows that flagged arrays are non-randomly placed relative to a conserved kinase gene in the same way across orders, which is hard to produce by chance CAAX hits. It does not show that any single array is a mating receptor.

### 5. Calls against arrays (`v_scenarios.tsv`, `v_calls_and_uncovered_by_order.tsv`, `v_genome_rates_by_order.tsv`)
| Order | Genomes | PR calls | Unverified | Unverified, single-locus array | Flagged arrays | Flagged, no call | Chance among all no-call arrays | Genomes with a no-call flagged array |
|---|---|---|---|---|---|---|---|---|
| Agaricales | 676 | 580 | 577 | 141 | 1,035 | 463 | 122.1 | 51% (42-61) |
| Boletales | 182 | 60 | 58 | 12 | 136 | 80 | 27.6 | 34% (26-41) |
| Polyporales | 175 | 198 | 116 | 22 | 227 | 40 | 22.1 | 18% (12-24) |
| Cantharellales | 77 | 99 | 99 | 75 | 100 | 5 | 19.3 | 7% (1-26) |
| Russulales | 76 | 31 | 23 | 6 | 65 | 36 | 18.9 | 32% (21-43) |
| Hymenochaetales | 47 | 13 | 13 | 1 | 13 | 2 | 9.8 | 4% (0-12) |
| Auriculariales | 14 | 13 | 13 | 1 | 13 | 0 | 1.2 | 0% |
| other orders | 40 | 45 | 45 | 26 | 44 | 1 | 9.5 | 3% (0-9) |
| All | 1,287 | 1,039 | 944 | 284 | 1,633 | 627 | 230.4 | 36% (31-43) |

- Calls: 1,039 in 766 genomes (60%); 1,018 arrays carry a call; 14 arrays carry more than one call; 2,475 loci lie in called arrays, 2,311 inside a call and 164 are siblings outside any call; 885 of the 1,018 called arrays are covered entirely. Median call covers 2 loci. One call per array: 1,028 records (gap 25 to 200 kb: 1,057 to 1,004; at 10 kb a call spans several arrays and the count becomes 1,237, so the array gap must be fixed before this is used).
- 95 calls are verified by other evidence; 944 are CAAX-admitted (`unverified`). Composition of the 944: array of >= 2 loci 651; precursor homology (tblastn) 414; >= 2 strict-CAAX ORFs 501; any of the three 804; none 140 (72 of them receptor + CAAX only). By order the unsupported ones are mostly Cantharellales: 73 of 99.
- Flagged arrays no call covers: 627 in 466 genomes; **230 expected by chance among all uncovered arrays (37%, 32-43)**; by order (observed against expected): Agaricales 463 against 122 (26%), Boletales 80 against 28 (34%), Polyporales 40 against 22, Russulales 36 against 19, Cantharellales 5 against 19, Hymenochaetales 2 against 10, other 1 against 9 (the last three at or below chance). By size: singletons 338 against 169 (50% chance), size >= 2 289 against 61 (21%), size >= 3 149 against 25 (17%). With precursor homology as well: 156 against 11 expected (7%); size >= 2 and homology 123 against 9.
- What the pipeline did with them: of the 627, 437 contain a locus withheld at the **fraction floor** (found fraction 0.25 against floor 0.5; every one of the 2,411 floor-withheld PR loci of the run shows exactly 0.25 and the gene content receptor + CAAX scan), 15 at the modelled-gene bar, and 182 have no withheld PR locus (445 have one of either kind). Median receptor identity of the floor-withheld locus is 66.7%, 94% are >= 50%, higher than that of the called unverified loci (62.0%, 84%). Meanwhile 195 called unverified loci have the same gene content (receptor + CAAX scan only) and were admitted. I did not trace why identical gene content passes in one place and sits at fraction 0.25 in another (the admission path through the CAAX-dependent label versus the cluster score); this needs a code reading before any rule is proposed.
- STE20 co-location supports the chance model: flagged arrays covered by a call 42.7% within 250 kb; flagged uncovered of size >= 2 34.6%; flagged uncovered singletons 6.2% (unflagged singletons 3.3%); unflagged uncovered 4.2%. Uncovered singletons look like chance flags; uncovered arrays of >= 2 loci look like the called ones.
- Panel (`panel_agaricomycetes.tsv`, 6 Agaricomycete genomes; only *Coprinopsis cinerea* and *Schizophyllum commune* carry genetically mapped receptors): in both, the mating receptors (3 and 2 loci) lie in the largest array (size 4), which is flagged; the array also holds 1 and 2 other copies of 12 "other" loci, so array membership does not separate a mating copy from a paralog inside the array. No PR call covers the *S. commune* array (the pipeline reports that B locus through its Balpha and Bbeta families, not checked here). Four more genomes (*Trametes*, *Heterobasidion*, *Russula*, *Grifola*) have labels chosen from CAAX position and cannot test anything.

## Solid and uncertain
**Solid (large n, survives controls or checked two ways)**
- Agaricomycetes carry a median of 9 STE3-like loci per genome (95.5% >= 5); 80% of arrays (60% of loci) are isolated singletons; the other 40% of loci sit in arrays of 2 to 10 loci at 50 kb.
- Clustering is above the within-contig null (2.4 times at 10 kb), but about 57% of loci in 50-kb arrays would be expected by contig fragmentation. Arrays of >= 3 loci are the reliable ones.
- The strict-CAAX flag is far above chance in most orders at array level (chance share 0.14 to 0.34, 0.18 overall) and almost entirely chance-free when two distinct CAAX ORFs are required (0.04). Singleton flags are about 30% chance; larger arrays 11 to 12%.
- The 627 flagged arrays with no call are 37% (32-43) chance-expected, not 14%; the earlier "13" was wrong. 86% of the excess (341 of 397) is in Agaricales.
- STE20 lies near flagged arrays 20 to 24 times more often (odds) than near unflagged arrays within the same genome, stable under array size, contig length, species resampling and leave-one-order-out; MIPBF shows no such pattern; HD and receptor arrays are rarely on the same contig.
- Calls and arrays: median call covers 2 loci; 1,039 calls fall in 1,018 arrays, so reporting one record per array would change the count by about 1%.

**Uncertain (say so in any use)**
- Whether any array, or flagged array, is a mating receptor: only two Agaricomycete genomes carry independent labels. The note shows that flagged arrays are non-random and mostly beyond chance, not that they mate.
- The chance model is contig-uniform: local CAAX density or gene-structure effects could move chance by an unknown amount. The pipeline's per-genome figure and the array-level figure differ by 5 points (13% against 18%).
- The STE20 effect: one gene chosen from known MAT linkage, heuristic top hit, no independent labels. Unflagged near-STE20 arrays (4.2%) may be real receptors the CAAX scan misses.
- Hymenochaetales (47 genomes, 29 species; 13 flagged arrays, no excess) and Cantharellales (14 species, 73 of 99 unverified calls unsupported by array, homology or a second CAAX ORF): no conclusion about the CAAX rule is supported there. Auriculariales (8 species) and the 12 other orders (31 species): too few.
- Why identical gene content (receptor + CAAX) is called in 195 places and withheld at fraction 0.25 in 437 flagged arrays: not traced.
- Array size 2 at 50 kb: roughly half or more are chance neighbours on short contigs.
- Species sampling: 705 species out of an order of magnitude more; BFD is dominated by some genera (Agaricales 313 species).

## Does this change what should be called?
For the 944 unverified calls, mostly no: 804 (85%) have array support, which agrees with the chance analysis that the large majority of CAAX calls are above chance, and an array view would restate rather than overturn them. It would flag the 140 with no array support (Cantharellales 73) and the Hymenochaetales calls as the weakest. For the 627 uncovered flagged arrays, partly yes: about 397 of them are above chance, concentrated in arrays of >= 2 loci (228 above chance), in Agaricales and Boletales, and these would add records in 245 genomes that have no PR call now (the genomes with a call would go from 766 to 1,011 of 1,287). The question is the label, not the evidence: these have the same gene content as 195 admitted calls.

## Options (no code change; counts on the 1,287-genome scanned set, 944 unverified calls, 1,039 calls)
| # | Option | Purpose | Before to after | Risk | Validation needed |
|---|---|---|---|---|---|
| 1 | Array-level reporting: one record per array with member count, span, which member carries the precursor or CAAX ORF, number of CAAX ORFs, homology | Report what is found; stop the same array appearing as separate calls; make singletons and 5-locus arrays distinguishable | 1,039 calls to 1,028 array records (1,057 to 1,004 for gaps 25 to 200 kb); 14 merged arrays, 164 sibling loci named; no call added or removed | Low (report layer). The array gap becomes part of the report (10 kb gives 1,237 records) | Fix the gap (50 kb suggested) and state that size 2 is weak; check the 14 merged arrays by hand; unit tests on records |
| 2 | Array-informed confirmation column (`array_support`: array of >= 2 loci, or tblastn precursor homology, or >= 2 distinct strict-CAAX ORFs), never used to lift or drop a call | Tell users which CAAX-admitted calls have more than one CAAX hit behind them | 944 unverified: 804 supported, 140 unsupported (Cantharellales 73, other orders 14, Polyporales 11); calls unchanged | Low if informational. If later used to lift `unverified`, 501 calls with >= 2 CAAX ORFs have 4% chance share (all arrays) | Compare with STE20 co-location by tier (done in this note: more CAAX ORFs go with higher co-location); independent labels before any label change |
| 3 | Admit uncovered flagged arrays as a separate low tier (`array_candidate`, not a PR call), limited to arrays of >= 2 loci, or >= 2 loci with homology, in orders where excess is shown | Recover receptors the call rules withhold without raising the call count | All 627: +60% records, 245 genomes get a first record, 230 (37%) expected chance. Size >= 2: +289 (21% chance), 178 first records. Size >= 2 + homology: +123 (about 7% chance), 89 first records. Agaricales 225 and 111 of those | Medium: chance records at 7 to 37%; Cantharellales, Hymenochaetales and other orders show no excess; the tier must not feed labels or counts of PR calls | Label some genomes independently (see 6); hold-out by order; per-order chance thresholds; STE20 co-location by tier as an independent check (size >= 2 uncovered: 34.6%) |
| 4 | Array-aware polish or fraction floor: let a withheld receptor + CAAX locus pass the floor when its array has >= 2 loci or a second CAAX ORF | Same loci as option 3 but with the pipeline's own gene models | 437 flagged-uncovered arrays hold a floor-withheld locus; 230 of the 289 size >= 2 ones; 225 have >= 50% identity; the modelled-gene bar (907 arrays) is not touched, because those arrays have no CAAX ORF | Needs a pipeline change; unexplained inconsistency with 195 admitted calls of the same content | Trace the fraction path first; full regression panel (below); replay on the 1,854 genomes |
| 5 | Use array structure to separate mating receptors from non-mating paralogs | Choose the receptor copy or locus per genome | A unique flagged array of >= 3 loci exists in 480 of 1,287 genomes; none in 772; two or more in 35 | High: only 2 independent labelled genomes, and the mating array also holds paralogs | Not supportable on this data; needs labels from >= 10 orders |
| 6 | Dedicated Agaricomycete PR family or focused run: receptor array of >= 2 STE3 loci, >= 2 CAAX ORFs, references from Agaricomycete B loci, optional STE20 flank | Test whether an array-based definition reproduces known B loci | Focused re-detect of the 1,287 quality genomes about 64 genome-hours, about an hour elapsed at 64-way (v0.6.0 Agaricomycetes: 92 genome-hours for 1,854 genomes, median 141 s each); one genome per species (705) about 35; the array scan costs 12 min per 30 genomes on 8 CPUs. Development, not compute, is the cost: curated Agaricomycete B-locus records, calls to calibrate against | High: the family needs labels that do not exist (2 independent genomes) | Gate: `testset/regression_panel.tsv` (163 genomes, of which 33 Basidiomycota) plus the 23 Zygo genomes, no changed call outside the new family; the 6 Agaricomycete panel genomes must keep their calls; a hold-out by order of the labelled set |

## Recommendation (ranked)
1. **Option 1 and option 2 together, as reporting only.** They cost nothing in calls, fix a counting artefact (the same array as several calls) and expose the weakest 140 calls and the Cantharellales and Hymenochaetales weakness. Do not use the support column to lift `unverified` yet.
2. **Option 3 restricted to arrays of >= 2 loci with precursor homology (123 records, 118 of them in Agaricales and Boletales) or with a second CAAX ORF (150)**, as a separate tier outside the PR call count, after one labelled check. For the 123 about 9 are expected by chance (assumes independence of CAAX and homology, so a lower bound), against 230 of 627 for the unrestricted rule. Do not admit singletons (50% chance) or any array in Cantharellales, Hymenochaetales or the 12 minor orders.
3. **Option 4 only after tracing the fraction-floor behaviour** (the same gene content is called in 195 places and withheld at 0.25 in 437), because one explanation would make options 3 and 4 the same fix and a second would leave a latent inconsistency in current calls.
4. **Option 6 after labels exist**: at least 10 orders with several genetically mapped receptors each; only 2 genomes have them now. The compute is small; the labelled set is the gate.
5. Option 5: do not pursue on this data.

## Next steps (not done here)
- Trace the fraction-floor path for receptor + CAAX loci and explain the 195 admitted against 437 withheld.
- Label Agaricomycete genomes independently (mapped B-locus receptors: *Coprinopsis*, *Schizophyllum*, *Pleurotus*, *Flammulina*, *Lentinula*, *Ustilago*-like outgroups are not enough), then re-run the STE20 and chance analyses against labels rather than against the flag.
- Repeat the STE20 test with a query set built from Agaricomycete genomes only, and with a second MAT-linked gene chosen blind to this data.
- A local-density chance model (shuffle CAAX ORFs within gene-matched windows) to bound the chance share.

## Files
`results/2026-10-06_agaricomycetes_pr_arrays/`: `array_verify.py`, `array_verify2.py`, `array_verify3.py`, `contig_lengths.sh`, `scan_out_full.tar.xz`, `v_*.tsv`, `v2_*.tsv`, `v3_summary.txt`, `v_summary.txt`, `v_arrays.tsv.gz`, `v_calls.tsv.gz`, plus the regenerated `order_summary.tsv`, `call_coverage_*.tsv`, `scenarios.tsv`, `array_organisation.tsv`, `flagged_vs_unflagged_arrays.tsv` (full scan, `analyze.py`).
