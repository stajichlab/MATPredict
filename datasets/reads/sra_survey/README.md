# Survey: fungal population read sets for `matpredict reads-type`

Date of the ENA query: 2026-10-26. Nothing is downloaded: this counts what exists. Plan and tracking: the issue linked from `datasets/reads/README.md`.

## Question
Which fungal species have read sets from at least 4 strains (pooled across studies) where a raw-read MAT call could be run, and which of them MATPredict already has references for?

## How the table was made
ENA mirrors SRA. Two read_run queries with the ENA portal API (`https://www.ebi.ac.uk/ena/portal/api/search`), both with
`query=tax_tree(4751) AND library_strategy="WGS" AND library_source="GENOMIC" AND instrument_platform="ILLUMINA"` and `format=tsv&limit=0`
(245,243 runs; 599,027 fungal runs of any kind):
- `fields=run_accession,study_accession,sample_accession,secondary_sample_accession,scientific_name,tax_id,library_layout,base_count,read_count,center_name,collection_date`
- `fields=run_accession,study_title,sample_title,sample_alias,strain,isolate,experiment_title,first_public,country,host`

Then `python3 survey_ena.py --runs RUNS.tsv --meta META.tsv --out . --db ../../../db` (about a minute). The two ENA tables (31 MB and 58 MB) are not committed; they are reproduced by the queries above.

## What is counted
- **Species**: the first two words of the scientific name, with a genus and a lower-case species epithet. Metagenomes, unidentified and `sp.` entries are left out. ***Saccharomyces cerevisiae* is skipped** (`--skip-species`, curator ruling 2026-10-26: the focus is outside the lab-model yeast; its public sets are mostly mutant collections, evolution lines and crosses, and it would be 42,000 of the strains). Other *Saccharomyces* species are kept.
- **Runs**: paired-end, at least 300 Mb (about 10x for a 30 Mb genome; a per-species genome size is an open improvement).
- **Strains**: distinct strain, isolate or sample alias (the first that is given) per study; a strain with several runs counts once and keeps its largest run.
- **Pooled across studies**: a strain name seen in several studies of one species counts once. Different studies often use different names for one strain, so this still over-counts somewhat.
- **Left out of the pooled count** (kept in `studies.tsv` with `excluded_as`): studies whose title, or at least half of whose strain names and sample titles, match
  - `mutant`: knockouts, deletion or insertion collections, CRISPR or transposon libraries, mutant strains, and engineered or tagged strains (`engineered`, `constructed`, `reporter`, `GFP`/`BFP`, `tagged`, `edited`, `suppressor`, `allele`, `lacking`, `absence of`, `null`, `point mutation`, `chimeric`);
  - `evolution`: experimental or laboratory evolution, serial passage, mutation accumulation lines, time courses;
  - `cross`: progeny, segregants, backcrosses, tetrads, recombinant inbred lines (these are not natural isolates, but a cross has an expected 1:1 idiomorph ratio and could be a useful test later).
  A title that names a natural or clinical sampling (`surveillance`, `point prevalence`, `survey`, `clinical isolates`, `outbreak`, `epidemiolog`, `GWAS`, `natural isolates`, `wild isolates`, `bloodstream`, `candidaemia`, `genomes project`, `population genom`) is not left out as `evolution` or `cross` on a stray keyword; it never rescues a `mutant` study.
  This is a keyword filter on titles and names, not a review. It has false positives and misses (for example a study named only by a project code). Treat `excluded_as` as a first pass.
- **Reference tier** (`species_summary.tsv`): A = MATPredict has accepted records for two or more idiomorphs of this species; A1 = one idiomorph; B = records for another species of the genus only; C = none. Matching is by deposited name, so a synonym (*Pyricularia oryzae* / *Magnaporthe oryzae*, *Nakaseomyces glabratus* / *Candida glabrata*) reads as B or C. Candidate records are not counted.

## Result (2026-10-26)
245,243 fungal runs, 105,479 paired runs of at least 300 Mb with a species name (without *S. cerevisiae*). 681 species have at least 4 strains before the exclusions and **653 species (73,949 strains) after them**: 2,040 strains were left out as mutants or engineered strains, 823 as experimental evolution and 2,335 as crosses.

| Reference tier | Species | Strains |
|---|---|---|
| A (two or more idiomorphs) | 44 | 41,448 |
| A1 (one idiomorph) | 19 | 1,683 |
| B (genus only) | 168 | 11,025 |
| C (none) | 422 | 19,793 |

The strain counts are dominated by a few species: *Candidozyma auris* (24,888, outbreak surveillance, largely clonal), *Zymoseptoria tritici*, *Aspergillus fumigatus*, *Cryptococcus neoformans*, *Candida albicans*.

Spot checks (2026-10-26): *Candida albicans* (77 studies): the evolved isolates, crosses, mutants and the GFP/BFP-engineered set are left out and the clinical and surveillance sets (the 182-isolate set, bloodstream isolates, "One Day in Denmark") are kept; *Candida auris* reads sit under a *C. albicans* study label and are counted as the label says. *S. paradoxus* and *S. uvarum*: natural and wild-isolate sets are kept, "evolved" and "resistant mutants" sets are left out; a few experimental sets without such wording ("experimental hybrid populations", a recombination-landscape study) are still kept.

## Not checked
- The strains are counted from metadata; no read set was opened, so a mislabelled species or a mixed sample is not caught (the Coccidioides set had two relabelled strains).
- Clonal sets (outbreaks, one lineage sampled many times) count as strains here; their idiomorph counts are not a population frequency.
- Restricted-access or withdrawn runs are not distinguished.
- Non-Illumina data (long reads, Ion Torrent) are left out of this survey.
