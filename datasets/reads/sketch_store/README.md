# Sketch store: one reusable read sketch per isolate

A sketch store keeps, for every public and local isolate of one species, a small FracMinHash sketch of a 30x subsample of its reads.
The same sketches answer several questions without downloading the reads again:

- **MAT**: idiomorph-specific k-mers (the MATPredict reads panel), beside the MATPredict `reads-type` call made on the same subsample
- **target screens**: presence of any sequence (a gene, a Starship, a gene cluster) and its relative copy number
- **population structure**: all-against-all distances, PCoA, clone groups and a FastME tree
- **QC**: mixtures, contamination and mislabelled samples (outlying distances, too many k-mers)

## What is kept per isolate

`sigs/<sample_id>.sig.zip`: sourmash FracMinHash, k=21 and k=31, scaled=100, with abundances; k-mers seen once are dropped.
About 5 MB per *A. fumigatus* isolate (both k), 25 GB for ~4,850 isolates. `mat/` holds the reads-type result and the GmG screen of the
same subsample; `stats/` the bases used. Reads are not kept.

Why these settings (tested 2026-11-07 on 136 *A. fumigatus* strains with SNP calls and a Starship truth set,
`results/2026-11-07_sketch_pcoa_mat` on the HPCC):

| choice | reason |
|---|---|
| 30x subsample, by genome size | a fixed 2M spots gives 50x for 12 Mb yeasts but under 6x for 100 Mb genomes; the >=2 filter needs depth |
| drop k-mers seen once | three quarters of an unfiltered sketch are sequencing errors (4x the size) |
| scaled=100 | an idiomorph gives ~22 specific hashes and a 638 bp gene 3; scaled=1000 cannot see either. Distances are the same at 100, 500 and 1000, so the distance step downsamples to 1000 |
| abundances | the distance needs them (solid containment, below), and screens report copy number |

Results of that test: MAT from the sketch agreed with reads-type in 129 of 129 single calls; PCoA agreed with the SNP PCoA
(Procrustes 0.63 on 2 axes, 0.70 on 3; nearest neighbour among the 3 nearest by SNPs for 89 of 113) but fine-scale distances did not
(rank correlation 0.34); element-specific Starships were called almost without error, the overlapping Janus family and Osiris-h3
were not; *bafA* 17 of 19, and 4 of the 5 disagreements were errors of the truth proxy or the sample, checked by read depth.

**Solid containment.** Plain Jaccard or containment cannot separate *A. fumigatus* strains: missing k-mers at 30x and surviving error
k-mers make two halves of one run as different as two strains. A k-mer counts in sample a's denominator only if its abundance in a is
at least max(3, 0.25 x median), and as shared if b has it at all. Same-strain median 0.99996, between-strain p99 0.999.

## Layout (HPCC)

```
/bigdata/stajichlab/jstajich/projects/MATPredict/sketch_store/<species>/
  manifest.tsv   one row per sample (build_manifest.py)        env.sh    panel, genome size, coverage, local input
  bin/           frozen copies of the scripts used for the build
  sigs/ mat/ stats/ logs/ lists/   jobs.tsv (every submission)
  sketched.tsv   cache index: run_accession, sample_id, params, path (rewritten every planner round; check before downloading)
  targets/       target FASTAs, screens.tsv (the standing screens), truth tables
  query/         qc/, mat/, screen_<name>/, distances/ (ani_k21.npy, pcoa.tsv, samples_qc.tsv with clone groups, tree_fastme.named.nwk)
```

## Build and query

```
python build_manifest.py --runs <PopSeq>/cache/fungi_runs_ena_<date>.tsv.gz --taxmap <PopSeq>/.../taxonomy_species_map.csv.gz \
    --species-taxid 746128 --genome-mb 29.4 --local-samples .../Afumigatus_Global/samples.csv --out STORE/manifest.tsv
cp sketch_one.sh sketch_query.py run_queries.sh submit.sh build_manifest.py STORE/bin/
cp ../sra_survey/pilot/genotype_mat_one_accession.sh ../sra_survey/pilot/screen_reads.py STORE/bin/
bash STORE/bin/submit.sh STORE [JOBS_TO_WAIT_FOR]      # self-feeding planner: rounds sized to the queue limit, retries, then run_queries.sh
python sketch_query.py --store STORE screen --name mygenes --targets genes.fa     # any new question, any time
python mat_frequency.py --store STORE --out OUTDIR      # MAT1-1 fraction per sample and per clone group, by study/country/lineage
python starship_dynamics.py --store STORE --out OUTDIR --truth STORE/targets/truth_store.tsv   # Starship/bafA prevalence, co-occurrence, gain/loss
```

`STORE/exclude.tsv` (sample_id, reason, decided) lists curator exclusions, for example a low-quality sequencing; `mat_frequency.py`
leaves them out, together with samples of another species (mean ANI < 0.95) and QC-flagged mixtures, and says why in its
`samples.tsv`. Its clone-corrected frequency gives each clone group (ANI 0.9995) one vote.

### Starship presence and gain/loss (`starship_dynamics.py`)

Reads the standing screen (`query/screen_starships_bafA/screen.tsv`: fraction of each target's specific hashes found), calls each
target present (fraction >= 0.5), absent (< 0.2) or missing (between; partial hits are never read as loss), takes one vote per clone
group like `mat_frequency.py`, and writes the clone-by-element matrix, prevalence (overall, by lineage and country, Wilson intervals),
pairwise co-occurrence, and Fitch parsimony on the store's distance tree (one tip per clone group; score, and gains/losses for the
root-absent and root-present reconstructions). `--set TARGET=FRACTION` changes the present cut for one target; `--truth` writes
`threshold_sweep.tsv` (TP/FN/FP/TN per cut) to choose it. The tree is a sketch distance tree and Fitch resolves multifurcations left
to right: read the counts as a pattern. Tests: `tests/test_sketch_store_starship_dynamics.py` (Newick parsing, Fitch against a
Sankoff calculation on random trees, clone votes, thresholds, outputs).

## Site report

`export_report.py --store STORE --matfreq MATFREQ_DIR --species "..." --slug <slug> --titles ../sra_survey/pilot/study_titles.tsv
--out reports/<slug>` (on the HPCC) writes `reports/<slug>/population.json` and `runs_public.tsv`, which are committed;
`scripts/build_reads_site.py` copies them for `web/intro/reads/population.html?s=<slug>`, lists the report on the raw-reads overview
and links it from the species report. Local isolates count in every number but are never named (curator decision 2026-11-09):
the run table is public runs only and the PCoA points carry no identifiers.

Public samples: one run per BioSample (paired, most bases), natural isolates only (`survey_ena.classify_study` drops crosses,
mutants, evolution lines). Local samples: every FASTQ sample of the Population_Genomics subproject; a local strain that matches a
public strain name is kept and named in `also_public` (two sequencings of what may be one strain, a useful same-strain check).
Adding samples later: append rows to the manifest and run `submit.sh` again; samples with a sketch are skipped.

Downloads: public runs come through the MATPredict reads worker (AWS open data, then NCBI prefetch). Some ENA runs have an NCBI record
but no loaded data (2026-11-08: 48 of 1,887 ENA *A. fumigatus* runs, deposited with sibling runs that did load, so they will not sync);
for those `sketch_one.sh` streams the first N read pairs from ENA (`ena_fastq`, HTTPS) and types them without the GmG screen.
When ENA refuses HTTPS (403 for 17 runs from the HPCC on 2026-11-09, after earlier retries) it tries the same files over FTP.

Several sessions build sketches on the HPCC; `/bigdata/stajichlab/jstajich/projects/MATPredict/sketch_store/REGISTRY.tsv` records
who builds which species or run list with which parameters. Add a row before starting a build.

`validation/` holds the scripts of the 136-strain test that chose these settings (run from the HPCC results folder named in its README).
