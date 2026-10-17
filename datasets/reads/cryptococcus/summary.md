# Cryptococcus Illumina paired-end WGS read inventory (ENA metadata only)

Source: ENA portal API, `result=read_run`, query
`tax_tree(5206) AND library_strategy="WGS" AND library_source="GENOMIC" AND instrument_platform="ILLUMINA" AND library_layout="PAIRED"`
(queried 2026-10-17; 7,135 runs returned for the whole genus). No reads were downloaded.
All numbers below were computed from the returned table. Genome size for coverage = 19 Mb (est_coverage = base_count / 19e6; ENA base_count is raw bases, so this over-estimates usable coverage).

## Files
- `manifest.tsv` - 6,991 runs in the species complex (one row per run). Extra columns beyond the request: `var_from_text`, `lineage_source`, `species_from_lineage`, `reference_strain`, `ref_text_mention`, `lab_flags`, `species_norm`, `experiment_accession`, `study_title`.
- `selected.tsv` - 312 runs, same columns plus `selection_reason`.
- `download.sh` - Slurm-array downloader for `selected.tsv` (syntax-checked with `bash -n`, never run).

## Scope decisions
- Kept: names C. neoformans (incl. strain-level taxids such as "C. neoformans H99"), C. deneoformans, C. gattii (umbrella name), C. bacillisporus, C. deuterogattii, C. tetragattii, C. decagattii (taxid 1859122), and "C. neoformans x C. deneoformans" (taxid 192011, AD hybrids).
- Excluded as outside the complex (144 runs, not in manifest): C. amylolentus 100, "Cryptococcus sp." 20, C. depauperatus 12, C. hyracis 8, C. wingfieldii 4.
- No C. gattii x neoformans hybrid taxon was returned. No C. neoformans var. name appears as a taxon (all species-level).
- Not filtered out by me: `library_selection` oddities (RT-PCR 83, DNase 119, size fractionation 30, PCR 47, "5-methylcytidine antibody" 4 runs). These are all submitted as WGS/GENOMIC and look like mis-set selection fields; they are noted in `notes`. The 4 MeDIP-labelled runs (PRJNA1281806, study title is garbled ASCII codes in ENA) deserve a manual look.

## Counts (manifest, 6,991 runs; 4,879 distinct BioSamples; 471 BioProjects)
| species (norm.) | runs | samples | low_cov (<20x) | lab-derived flag | est. size (GB, both mates + any unpaired) | selected |
|---|---|---|---|---|---|---|
| C. neoformans | 5,228 | 3,400 | 401 | 961 | 4,483 | 50 |
| C. deneoformans | 731 | 635 | 20 | 127 | 620 | 50 |
| C. deuterogattii | 294 | 234 | 6 | 51 | 516 | 50 |
| C. gattii (umbrella name) | 271 | 198 | 15 | 0 | 388 | 50 |
| C. bacillisporus | 231 | 183 | 1 | 0 | 315 | 50 |
| C. tetragattii | 155 | 149 | 3 | 0 | 162 | 50 |
| C. neoformans x deneoformans (AD hybrid) | 79 | 78 | 2 | 67 | 151 | 10 |
| C. decagattii | 2 | 2 | 0 | 0 | 1 | 2 |

Size column is the total fastq bytes of the _1/_2 files of all runs in that species (rounded). The whole manifest is 6.64 TB.

### Lineage (explicit in metadata only)
VGIII 188, VGII 141, VGIV 120, VNB 36, VGI 28, VNII 6, VNI 3, **none found 6,469 (92.5%)**.
264 of the 522 assignments come from sample-level fields, 258 from a single-lineage study title. Multi-lineage study titles (e.g. "VGI, VGII, VGIII and VGIV") are deliberately NOT assigned to samples. Serotype field: A 1,051; D 437; B 20; AD hybrid 10; C 3; blank 5,468. `var_from_text` (var. grubii / var. neoformans in titles) gives C. neoformans 2,592 and C. deneoformans 572 runs. **VNI/VNII/VNIII/VNIV are essentially not recoverable from metadata**; most C. neoformans lineage must come from genotyping the reads.

### Top BioProjects (runs)
PRJNA269039 485 ("Saccharomyces neoformans Genome sequencing", old taxon name in title); PRJEB27222 414; PRJNA967729 390; PRJNA932005 384; PRJNA420966 320 (XL280 progeny); PRJNA764746 284; PRJNA999197 211 (hybrid speciation); PRJNA1328773 152; PRJEB76062 150 (clinical, China); PRJNA804139 120 (VGIV); PRJNA1054030 117; PRJNA244927 117 (VGIII).

### Coverage distribution (est. x, runs)
<5: 90; 5-10: 26; 10-20: 332; 20-50: 2,255; 50-100: 1,954; 100-200: 1,489; 200-500: 786; >=500: 59. Low-coverage (<20x): 448 runs (6.4%).

### Duplicates / multi-run samples
526 BioSamples have more than one run (max 36 runs for one sample); 2,638 runs belong to such samples. Pick one run per sample unless you want to merge them.

## Flags and data-quality concerns
- **Lab-derived runs** (regex on titles/strain, heuristic, expect false positives/negatives): cross_progeny 923, exp_evolution_or_passage 119, both 27, mutant_or_engineered_possible 137. Kept in manifest, excluded from selection. Many are valuable for MAT validation (progeny with known MAT alleles, e.g. PRJNA420966 XL280 crosses, PRJNA747819 B3502 crosses) - consider a separate "known-segregation" set.
- **AD hybrids**: 79 runs, but 67 are flagged lab-derived (hybrid evolution experiments, PRJNA999197 etc.); only 10 natural/unflagged and eligible. Hybrids carry two MAT haplotypes; read typing must expect both.
- **Species label suspects: only 3 runs flagged** by my rules: SRR10233430 (labelled C. neoformans, strain "Fungi", lineage VGI in metadata), SRR13663995 and SRR13663996 (labelled AD hybrid, serotype D / A). Checks used lineage-vs-species, var. text-vs-species, serotype-vs-species. Limits: the check cannot detect mislabels without lineage/serotype/var text (most runs).
- **C. gattii umbrella name**: 271 runs use the pre-2015 name. 63 of them carry lineage VGIII/VGI text; per the 2015 split VGIII = C. bacillisporus, VGI = C. gattii s.s. `species_from_lineage` records this. These were not treated as errors. `species_norm` stays the ENA name.
- **Missing metadata**: country blank/missing 2,210 runs; host blank/missing 4,803 runs; strain/isolate often only in sample title.
- **Unpaired files**: 177 runs list 3 fastq files (an extra `<run>.fastq.gz` of orphan reads); download uses only _1/_2. 5 runs (ERR1395164, ERR1395165, SRR1045115-17) are layout=PAIRED but have no _1/_2 pair on ENA (not downloadable by this script; excluded from selection).
- **Garbled study title** on 5 runs (ASCII character codes) - ENA metadata oddity, kept.
- Coverage is a raw-base estimate; actual mapped coverage will be lower.

## Reference strains (strict match on strain/isolate, or ENA taxon name for runs with blank strain)
| strain | runs found | notes |
|---|---|---|
| H99 | 50 | many are in a mutant-flagged study (PRJNA791949, SRR17322558 onward); selected SRR33889912 (378x) |
| KN99 | 75 | KN99alpha/a only; selected SRR22806371 |
| JEC21 | 27 | selected SRR33159723 and SRR25436383 |
| JEC20 | 38 | selected SRR11576540 |
| R265 | 29 | selected SRR1520841, SRR8653984 |
| WM276 | 1 | SRR643918, 115x (species name "C. gattii WM276") - selected |
| XL280 | 5 strict (+322 progeny/derivatives by text) | selected SRR1505176 |
| B-3501 | **0** | absent; B-3502 appears only in text of progeny samples (PRJNA747819 "B_xx", flagged cross_progeny); no B-3502 parent run found by strict match |
| E566, IND107, CDC15 | 6 / 7 / 0 strict | CDC15 only in text of one China run |

## Selection (selected.tsv, 312 runs, ~539 GB, 103 BioProjects, median est. 116x, min 29.7x)
Rule (implemented in the build script; reproduce by re-filtering manifest.tsv):
1. Eligible: _1 and _2 fastq present; est_coverage >= 20x; species_label_suspect = no; lab_flags empty.
2. One run per BioSample: highest min(coverage,100x), ties broken by smaller download.
3. Per species cap of **50** (recommended; matches your suggestion). Species with fewer eligible runs keep all (AD hybrid 10, decagattii 2).
4. Within a species: first the reference strains present (one run per reference name), then a two-level round-robin: across lineage groups (explicit lineage or "unknown"), and within each group across BioProjects (least-used project first, then highest capped coverage). This spreads over lineages and studies rather than taking the 50 deepest runs of one project.
Coverage is capped at 100x only for ranking; very deep runs (up to ~1,130x) can still be picked and are expensive. Consider subsampling after download.
Eligible samples available per species (before cap): C. neoformans 2,500; C. deneoformans 508; C. deuterogattii 199; C. gattii 197; C. bacillisporus 182; C. tetragattii 149.

## Could not determine / open questions
- VNI vs VNII vs VNB for most C. neoformans (lineage absent from metadata); the selected 50 neoformans runs are mostly "unknown" lineage and may be skewed toward whichever BioProjects round-robin favoured (VNB 23 explicit).
- Whether lab-derived/progeny runs should be a separate validation set.
- Whether to accept the 83 RT-PCR / 119 DNase / 4 MeDIP selection-labelled runs (kept, flagged in notes, unvetted).
- Whether AD hybrids deserve a cap exception and inclusion of lab hybrids.
- SRA-only runs not yet mirrored to ENA are not captured; counts are ENA-as-of-today.
