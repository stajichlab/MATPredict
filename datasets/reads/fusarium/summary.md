# Fusarium public Illumina WGS read inventory (ENA metadata, queried 2026-10-17)

Everything here was computed from ENA portal API metadata. No read data was downloaded.
Files: `manifest.tsv` (all runs), `selected.tsv` (proposed subset), `download.sh`, `build_manifest.py` (regenerates the tables from ENA).

## Query
`result=read_run`, `tax_tree(5506) AND library_strategy="WGS" AND library_source="GENOMIC" AND instrument_platform="ILLUMINA" AND library_layout="PAIRED"`
returned **3310 runs** (3,310 unique run accessions, 2600 BioSamples, 253 BioProjects, 128 species labels).
- Synonyms: the ENA taxonomy has **no** taxon named Neocosmospora, Gibberella, Haematonectria, Cyanonectria or Albonectria (name lookups return nothing, read_run name searches for Neocosmospora/Haematonectria/Cyanonectria/Albonectria return 0 runs; the Gibberella wildcard search returned an HTML error page, so only its taxonomy lookup, which is empty, was confirmed). The *F. solani* species complex sits inside taxid 5506 (lineage "Fusarium solani species complex"), so tax_tree(5506) already covers it. *Gibberella zeae* is a synonym of taxid 5518 (*F. graminearum*).
- Same filter on other strategies, for context (Illumina PE, tax_tree 5506, not in manifest): RNA-Seq 3,398 transcriptomic + 101 genomic, AMPLICON 362, ChIP-Seq 138, Hi-C 26, ATAC 14, RAD-Seq 9; plus **OTHER/GENOMIC 307, WGA 27, WCS 12** that could include some real genome sequencing but were excluded by the WGS rule (not inspected). Single-end Illumina WGS: 169 runs, excluded.
- species_group comes from the NCBI lineage of the run's taxid (oxysporum / fujikuroi / solani species-complex nodes). There is no FGSC node in NCBI, so graminearum_complex is assigned by epithet (graminearum, asiaticum, boothii, meridionale, gerlachii, aethiopicum, louisianense, cortaderiae, vorax, ussurianum, mesoamericanum, austroamericanum, brasilicum, nepalense, acaciae-mearnsii; only 7 of these occur). *F. pseudograminearum*, *F. culmorum*, *F. sambucinum*, *F. poae* are in group `other`. `oxysporum` group includes *F. odoratissimum* (TR4) and other species NCBI placed in the oxysporum complex; species is kept in `species`.

## Counts by species_group
| group | runs | species labels | selected | all-runs GB | selected GB |
|---|--:|--:|--:|--:|--:|
| fujikuroi_complex | 380 | 34 | 257 | 741 | 575 |
| graminearum_complex | 431 | 7 | 178 | 663 | 257 |
| other | 799 | 37 | 348 | 1025 | 756 |
| oxysporum | 1505 | 13 | 200 | 2727 | 572 |
| solani | 195 | 40 | 168 | 540 | 503 |

Total download size: all 3310 runs = **5.70 TB** (ENA fastq_bytes, both mates, where listed); selected 1151 runs = **2.66 TB** (median 1.7 GB/run).

## F. graminearum
- **403 runs / 401 BioSamples** labelled *Fusarium graminearum* (taxid 5518); 431 runs in the whole graminearum_complex (adds boothii 9, meridionale 7, asiaticum 5, gerlachii 3, aethiopicum 2, louisianense 2).
- Of the 403: 109 are flagged as engineered or lab-mutant lines (CRISPR, deletion, suppressor, reporter strains; keyword match on study/sample titles, so this is a lower bound), 8 are <20x, 3 have no usable mate pair. **285 BioSamples (286 runs) remain eligible.**
- 248 of those 286 come from three BioProjects (PRJNA1036432 French isolates 102, PRJNA677929 87, PRJNA397890 population genomics 59); the rest are spread over 18 other BioProjects (21 BioProjects in total).

## Species table (runs / BioSamples / eligible BioSamples / selected)
Eligible = paired files present, not pooled, not engineered-flagged, >=20x, not already used, not the Fola cohort project; counted as one per BioSample.

| species | group | runs | biosamples | eligible biosamples | selected |
|---|---|--:|--:|--:|--:|
| Fusarium oxysporum | oxysporum | 1439 | 819 | 691 | 141 |
| Fusarium graminearum | graminearum_complex | 403 | 401 | 285 | 150 |
| Fusarium pseudograminearum | other | 276 | 238 | 51 | 50 |
| Fusarium sambucinum | other | 159 | 159 | 157 | 50 |
| Fusarium proliferatum | fujikuroi_complex | 111 | 109 | 109 | 50 |
| Fusarium culmorum | other | 106 | 106 | 105 | 50 |
| Fusarium fujikuroi | fujikuroi_complex | 100 | 100 | 100 | 50 |
| Fusarium poae | other | 90 | 90 | 89 | 50 |
| Fusarium andiyazi | fujikuroi_complex | 36 | 36 | 36 | 36 |
| Fusarium sacchari | fujikuroi_complex | 34 | 34 | 34 | 34 |
| Fusarium sp. | solani | 33 | 32 | 31 | 31 |
| Fusarium cf. solani | solani | 30 | 24 | 22 | 22 |
| Fusarium venenatum | other | 28 | 28 | 27 | 27 |
| Fusarium odoratissimum | oxysporum | 26 | 19 | 19 | 19 |
| Fusarium falciforme | solani | 26 | 25 | 25 | 25 |
| Fusarium verticillioides | fujikuroi_complex | 25 | 23 | 21 | 21 |
| Fusarium avenaceum | other | 21 | 21 | 21 | 21 |
| Fusarium solani | solani | 20 | 20 | 18 | 18 |
| Fusarium cf. falciforme | solani | 18 | 10 | 10 | 10 |
| Fusarium keratoplasticum | solani | 17 | 17 | 17 | 17 |
| Fusarium lateritium | other | 16 | 16 | 16 | 16 |
| Fusarium circinatum | fujikuroi_complex | 13 | 13 | 13 | 13 |
| Fusarium petroliphilum | solani | 13 | 13 | 13 | 13 |
| Fusarium xylarioides | fujikuroi_complex | 11 | 11 | 11 | 11 |
| Fusarium euwallaceae | solani | 10 | 9 | 8 | 8 |
| Fusarium boothii | graminearum_complex | 9 | 9 | 9 | 9 |
| Fusarium purpurascens | oxysporum | 9 | 9 | 9 | 9 |
| Fusarium annulatum | fujikuroi_complex | 7 | 7 | 7 | 7 |
| Fusarium meridionale | graminearum_complex | 7 | 7 | 7 | 7 |
| Fusarium commune | other | 7 | 6 | 3 | 3 |
| Fusarium tardichlamydosporum | oxysporum | 7 | 7 | 7 | 7 |
| Fusarium cf. verticillioides | fujikuroi_complex | 6 | 4 | 4 | 4 |
| Fusarium phialophorum | oxysporum | 6 | 6 | 6 | 6 |
| Fusarium kuroshium | solani | 6 | 6 | 6 | 6 |
| Fusarium nygamai | fujikuroi_complex | 5 | 5 | 5 | 5 |
| Fusarium asiaticum | graminearum_complex | 5 | 5 | 5 | 5 |
| Fusarium duoseptatum | oxysporum | 5 | 5 | 5 | 5 |
| Fusarium ambrosium | solani | 5 | 4 | 4 | 4 |
| Fusarium cf. nygamai | fujikuroi_complex | 4 | 2 | 1 | 1 |
| Fusarium brachygibbosum | other | 4 | 4 | 4 | 4 |
| Fusarium cf. hostae | other | 4 | 3 | 3 | 3 |
| Fusarium equiseti | other | 4 | 4 | 4 | 4 |
| Fusarium hostae | other | 4 | 2 | 2 | 2 |
| Fusarium langsethiae | other | 4 | 2 | 2 | 2 |
| Fusarium musae | other | 4 | 3 | 3 | 3 |
| Fusarium tricinctum | other | 4 | 2 | 2 | 2 |
| Fusarium vanettenii | solani | 4 | 3 | 3 | 3 |
| Fusarium virguliforme | solani | 4 | 4 | 4 | 4 |
| Fusarium temperatum | fujikuroi_complex | 3 | 2 | 2 | 2 |
| Fusarium gerlachii | graminearum_complex | 3 | 3 | 3 | 3 |
| Fusarium dimerum | other | 3 | 3 | 2 | 2 |
| Fusarium nematophilum | other | 3 | 3 | 3 | 3 |
| Fusarium redolens | other | 3 | 2 | 2 | 2 |
| Fusarium grosmichelii | oxysporum | 3 | 3 | 3 | 3 |
| Fusarium veterinarium | oxysporum | 3 | 3 | 3 | 3 |
| Fusarium tucumaniae | solani | 3 | 3 | 3 | 3 |
| Fusarium chuoi | fujikuroi_complex | 2 | 2 | 2 | 2 |
| Fusarium pseudonygamai | fujikuroi_complex | 2 | 1 | 1 | 1 |
| Fusarium subglutinans | fujikuroi_complex | 2 | 2 | 2 | 2 |
| Fusarium aethiopicum | graminearum_complex | 2 | 2 | 2 | 2 |
| Fusarium louisianense | graminearum_complex | 2 | 2 | 2 | 2 |
| Fusarium acuminatum | other | 2 | 2 | 2 | 2 |
| Fusarium cerealis | other | 2 | 2 | 1 | 1 |
| Fusarium cf. aywerte | other | 2 | 1 | 1 | 1 |
| Fusarium citri | other | 2 | 2 | 2 | 2 |
| Fusarium clavum | other | 2 | 2 | 2 | 2 |
| Fusarium delphinoides | other | 2 | 2 | 2 | 2 |
| Fusarium irregulare | other | 2 | 2 | 2 | 2 |
| Fusarium sinense | other | 2 | 1 | 1 | 1 |
| Fusarium sporotrichioides | other | 2 | 2 | 2 | 2 |
| Fusarium nirenbergiae | oxysporum | 2 | 2 | 2 | 2 |
| Fusarium tardicrescens | oxysporum | 2 | 2 | 2 | 2 |
| Fusarium akasia | solani | 2 | 2 | 1 | 1 |
| Fusarium awan | solani | 2 | 2 | 2 | 2 |
| Fusarium bataticola | solani | 2 | 1 | 1 | 1 |
| Fusarium cucurbiticola | solani | 2 | 1 | 1 | 1 |
| Fusarium duplospermum | solani | 2 | 2 | 2 | 2 |
| Fusarium mekan | solani | 2 | 2 | 2 | 2 |
| Fusarium metavorans | solani | 2 | 2 | 1 | 1 |
| Fusarium obliquiseptatum | solani | 2 | 2 | 2 | 2 |
| Fusarium papillatum | solani | 2 | 2 | 2 | 2 |
| Fusarium pseudensiforme | solani | 2 | 2 | 2 | 2 |
| Fusarium rekanum | solani | 2 | 2 | 2 | 2 |
| Fusarium tuaranense | solani | 2 | 2 | 2 | 2 |
| Fusarium variasi | solani | 2 | 2 | 2 | 2 |
| Fusarium acutatum | fujikuroi_complex | 1 | 1 | 1 | 1 |
| Fusarium ananatum | fujikuroi_complex | 1 | 1 | 1 | 1 |
| Fusarium anthophilum | fujikuroi_complex | 1 | 1 | 1 | 1 |
| Fusarium bactridioides | fujikuroi_complex | 1 | 1 | 1 | 1 |
| Fusarium begoniae | fujikuroi_complex | 1 | 1 | 1 | 1 |
| Fusarium brevicatenulatum | fujikuroi_complex | 1 | 1 | 1 | 1 |
| Fusarium concentricum | fujikuroi_complex | 1 | 1 | 1 | 1 |
| Fusarium denticulatum | fujikuroi_complex | 1 | 1 | 1 | 1 |
| Fusarium ficicrescens | fujikuroi_complex | 1 | 1 | 1 | 1 |
| Fusarium fracticaudum | fujikuroi_complex | 1 | 1 | 1 | 1 |
| Fusarium globosum | fujikuroi_complex | 1 | 1 | 1 | 1 |
| Fusarium guttiforme | fujikuroi_complex | 1 | 1 | 1 | 1 |
| Fusarium lactis | fujikuroi_complex | 1 | 1 | 0 | 0 |
| Fusarium napiforme | fujikuroi_complex | 1 | 1 | 1 | 1 |
| Fusarium pseudoanthophilum | fujikuroi_complex | 1 | 1 | 1 | 1 |
| Fusarium pseudocircinatum | fujikuroi_complex | 1 | 1 | 1 | 1 |
| Fusarium ramigenum | fujikuroi_complex | 1 | 1 | 1 | 1 |
| Fusarium udum | fujikuroi_complex | 1 | 1 | 1 | 1 |
| Fusarium azukiicola | other | 1 | 1 | 1 | 1 |
| Fusarium chlamydosporum | other | 1 | 1 | 1 | 1 |
| Fusarium coffeatum | other | 1 | 1 | 1 | 1 |
| Fusarium incarnatum | other | 1 | 1 | 1 | 1 |
| Fusarium kistleri | other | 1 | 1 | 1 | 1 |
| Fusarium rubellum | other | 1 | 1 | 1 | 1 |
| Fusarium terricola | other | 1 | 1 | 1 | 1 |
| Fusarium torreyae | other | 1 | 1 | 1 | 1 |
| Fusarium cugenangense | oxysporum | 1 | 1 | 1 | 1 |
| Fusarium foetens | oxysporum | 1 | 1 | 1 | 1 |
| Fusarium kalimantanense | oxysporum | 1 | 1 | 1 | 1 |
| Fusarium brasiliense | solani | 1 | 1 | 1 | 1 |
| Fusarium cuneirostrum | solani | 1 | 1 | 1 | 1 |
| Fusarium drepaniforme | solani | 1 | 1 | 1 | 1 |
| Fusarium floridanum | solani | 1 | 1 | 1 | 1 |
| Fusarium lichenicola | solani | 1 | 1 | 1 | 1 |
| Fusarium neocosmosporiellum | solani | 1 | 1 | 1 | 1 |
| Fusarium oligoseptatum | solani | 1 | 1 | 1 | 1 |
| Fusarium phaseoli | solani | 1 | 1 | 1 | 1 |
| Fusarium piperis | solani | 1 | 1 | 1 | 1 |
| Fusarium suttonianum | solani | 1 | 1 | 1 | 1 |
| Fusarium tenuicristatum | solani | 1 | 1 | 1 | 1 |
| Fusarium tonkinense | solani | 1 | 1 | 1 | 1 |
| Fusarium waltergamsii | solani | 1 | 1 | 1 | 1 |
| Fusarium warna | solani | 1 | 1 | 1 | 1 |

## Top BioProjects (by runs in manifest)
| BioProject | runs | selected | top species | center |
|---|--:|--:|---|---|
| PRJNA38293 | 414 | 2 | F. oxysporum (414) | BI |
| PRJNA412392 | 393 | 28 | F. oxysporum (343), F. cf. falciforme (18) |  |
| PRJNA306180 | 246 | 42 | F. pseudograminearum (246) |  |
| PRJNA1398885 | 160 | 50 | F. sambucinum (159), F. venenatum (1) | SUB15834525 |
| PRJNA1036432 | 103 | 37 | F. graminearum (103) | SUB13730439 |
| PRJNA1367988 | 102 | 1 | F. oxysporum (102) | SUB16240615 |
| PRJNA1027381 | 96 | 46 | F. fujikuroi (90), F. proliferatum (6) | SUB16215754 |
| PRJNA749912 | 95 | 39 | F. culmorum (95) | SUB10063328 |
| PRJNA677929 | 93 | 38 | F. graminearum (93) | SUB12014841 |
| PRJNA895867 | 80 | 77 | F. andiyazi (35), F. sacchari (29) | SUB12231371 |
| PRJNA1227444 | 78 | 38 | F. poae (78) | SUB15401034 |
| PRJNA397890 | 71 | 49 | F. graminearum (60), F. boothii (5) | SUB2959261 |
| PRJNA517364 | 67 | 11 | F. proliferatum (67) | SUB5439097 |
| PRJNA1391476 | 64 | 62 | F. falciforme (16), F. keratoplasticum (13) | SUB15834300 |
| PRJNA1098703 | 62 | 0 | F. oxysporum (61), F. commune (1) | SUB14368275 |

## Coverage distribution (est_coverage = base_count / 40 Mb, all runs)
| x | runs |
|---|--:|
| 0-5 | 120 |
| 5-10 | 180 |
| 10-20 | 148 |
| 20-50 | 1126 |
| 50-100 | 886 |
| 100-200 | 640 |
| 200-+ | 210 |

Median 53.7x. 448 runs (<20x) carry `low_coverage=Y`. The 40 Mb divisor is nominal: *F. graminearum* is about 36-38 Mb, *F. oxysporum*/*solani* 50-60+ Mb (true depth of those is roughly 25-35% lower than shown), and accessory chromosomes in ff.sp. add more. The MAT locus is single-copy so this is a rough guide; the earlier Fola reads-typing timing was quoted at 8M reads per strain.

## Already used in MATPredict (read-only check of docs/, results/, panels/, analysis/, db/)
- **Marked `yes`**: SRR28734943 (strain 20a) and SRR28734937 (JCP043), both PRJNA1098703, *F. oxysporum* f. sp. lactucae.
- **Marked `likely`**: the other 60 runs of PRJNA1098703. The Fola read-typing set (148 strains, results/2026-10-06_fola_reads_type/samples.tsv) was built from local trimmed fastq under the curator's `nicolel/Fola/` directory; the repo only names two SRR accessions, so which SRA runs went into the 148 is **not determinable from the repo**. PRJNA1098703 (49 lactucae runs plus other oxysporum/commune) matches the SRR287349xx accession range and is the best guess. The Fola set also includes an in-house AVITI batch (not public).
- No other Fusarium SRA accession appears in the repo. Fusarium MAT reference records in the database are *F. fujikuroi* (5127_mo44 / mo45) and *F. graminearum* (5518_3639); none for *F. oxysporum* at the time of the 2026-10-05 spec.
- Both `yes` and `likely` runs are excluded from selected.tsv.

## Selection (selected.tsv, 1151 runs, 127 species)
Rule, implemented in `build_manifest.py`:
1. Candidate runs: mate pair listed in ENA, not a pooled run, not engineered-flagged, est_coverage >= 20x, not already used, not in PRJNA1098703.
2. One run per BioSample: the highest-coverage run.
3. Per species cap: **50**, except *F. graminearum* **150** and *F. oxysporum* **150**.
4. Greedy pick within a species: fewest already-picked from the same forma specialis (oxysporum only; max 6 per named f. sp., 30 with no f. sp. given), then fewest from the same BioProject, then fewest from the same country (missing values grouped as unknown), then highest coverage (coverage scored as min(x, 60) so 60x and 300x tie), then smaller download.
5. `selection_rank` records the pick order within each species, so you can truncate (e.g. keep rank <= 25) to shrink the set.

Why the exceptions: *F. graminearum* is the best sampled species here (285 eligible BioSamples) and the user's main population target; 150 reaches about half of it while capped by project, and beyond ~150 additions are mostly more isolates of the same three populations. *F. oxysporum* is a species complex with 31 formae speciales named in the metadata, so 50 would cover under 10 ff.sp.; 150 allows 6 each across the named ff.sp. In practice the rule selects 141 *F. oxysporum* runs (all 31 named ff.sp. plus 30 with no f. sp. recorded), below the cap because the per-f.sp. limits bind. Species with fewer eligible BioSamples than the cap are taken whole (93 of the 128 species labels have fewer than 5 eligible BioSamples). *F. lactis* has 0 selectable runs.

Selected coverage: median 76x; 218 at 20-50x, 520 at 50-100x, 413 at >=100x. Spread over 220 BioProjects.

## Could not determine / data-quality concerns
- **Species labels are the submitter's, unverified.** `species_label_suspect` flags only 104 runs: `cf.`/`aff.`/`sp.`/`unidentified` labels (97) and sample-level text naming a different species (7, e.g. *F. odoratissimum* samples titled F. oxysporum f. sp. cubense TR4, which is a reclassification). Genuine mislabels (known problem in oxysporum and the FSSC) will not be caught without sequence-based checks; the MAT reads-typing result could itself help.
- **Forma specialis** is given in the ENA scientific_name for 633 of 1,439 *F. oxysporum* runs; where it was recovered from strain/sample text the note says `f.sp._from_<field>`. 806 oxysporum runs have none recorded.
- **Pooled runs**: 18 runs of PRJNA38293 (Broad *F. oxysporum* comparative project) list a run-level taxid list for several isolates at once; flagged, excluded from selection.
- **Files**: 10 runs have no _1/_2 pair in ENA (a single fastq or none; notes column); 82 runs also list an extra unpaired file, which manifest and download.sh ignore. All 3,300 usable runs have md5 and byte counts.
- **Multiple runs per BioSample**: 293 BioSamples have >1 run (1,003 rows flagged `duplicate_sample`), mostly PRJNA38293 and PRJNA412392. The same strain name also appears under several BioSamples (note `same_strain_name_in_N_biosamples`), so one-per-BioSample does not guarantee one-per-isolate.
- **Lab mutants** (139 runs, 109 *F. graminearum*, 23 *F. pseudograminearum*, 7 *F. odoratissimum*) are mostly PH-1/CS3005/RBG5266 derivatives; keyword-based detection is incomplete, and untitled mutant lines in other projects will remain in selected.tsv.
- **Metadata gaps**: country missing for 132 runs (plus placeholder strings like "missing: data agreement established pre-2023"), host missing for 1,336, collection date inconsistent. Strain/isolate field is filled (falls back to sample alias).
- **Coverage** is raw base_count over a nominal genome; no adapter/contamination/quality checks were made. ENA mirrors NCBI SRA with a short lag and omits controlled-access data.
- **Strategy filter**: OTHER/GENOMIC (307), WGA (27) and single-end WGS (169) runs were not examined and may hold usable data.
