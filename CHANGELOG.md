# Changelog

All notable changes to MATPredict. Versions follow [semantic versioning](https://semver.org/);
release tags may carry a clade suffix naming the lineage whose infrastructure that
release completed.

## [Unreleased]

### Added
- Platforms: `pixi.toml` and `pixi.lock` now cover `osx-arm64` and `linux-aarch64` as well as `linux-64`; a `tests` workflow runs the suite (`pixi run -e test pytest tests`) and a smoke test on Linux and macOS; the Docker image is built for `linux/amd64` and `linux/arm64`. The re-lock also moved the `linux-64` pins, so run the regression check against it.
- Five accepted tier-3 records (curator sign-off 2026-10-30, J. Stajich): *Saccharomycopsis schoenii* GDD12 MATss (MATa-MATalpha left cluster and MATa right cluster, `db/Ascomycota/Ascoideales/`), *Barnettozyma californica* PB4207 MATbz and *Petasospora fabianii* MATpf (JOY008 MATa, TBRC 4498 MATalpha, `db/Ascomycota/Phaffomycetales/`). The gene models are our own (RNA-seq for *S. schoenii*, sequence-based for the others). With them, 20 of the 21 genomes that lose every MATsc/MTL call under the v0.8 cross-family rules are called again by their own family (trial run, `--preset v0.8`); *Ogataea degrootiae* GCA_019232045.1 stays uncalled.
- Raw-reads wave 2 (`datasets/reads/sra_survey/pilot/wave2/`): 55 more tier A, A1 and B species typed (1,468 runs; 13 confirmed, 16 provisional, 7 provisional-homothallic, 49 withheld), yeast idiomorph names (a / alpha, P / M, A / B) carried through, *Botrytis cinerea* typed on all 333 eligible samples. Species reports gain an allele-count table for all studies pooled or study by study, per run or per clone, and clone groups for every species with 10 or more sketched runs (average linkage at 0.998; counts per clone only where the clone number is stable). Output only; the genome tables and the database are unchanged.
- Idiomorph labels beside the raw names (curator rulings of 2026-10-28): `scripts/idiomorph_label_registry.json` (lineage and reference-organism codes, cell tokens, gene maps), `scripts/idiomorph_labels.py`, `scripts/label_loci_table.py`; the site data carry `label_alleles` in `taxa.json`, `call_label` / `alleles_label` in the genome tables and `data/idiomorph_labels.json`, and the explorer shows `Lineage__System_Cell` labels (for example `Taphr__SchPom_M`, `Sacch__CanAlb_alpha`, `Peziz__AspFum_MAT1-1`) with the raw name and the reference organism in the tooltip, and a switch back to the raw names. Output only: the database, the reports and the raw columns are unchanged. All 27,885 v0.7.0 loci are labelled. Duplicate family calls (`scripts/dedupe_family_calls.py`) are not removed from the counts.
- Raw-reads pages (`web/intro/reads/`): an overview of the wave 1 read typing (112 species, 1,629 runs, 599 of 605 single-allele calls agree with an assembly call) with the confidence of each species panel (confirmed / provisional / withheld), and one report per species (calls, agreement with assemblies, calls by study, every run, TSV download). Data written at site-build time by `scripts/build_reads_site.py` from `datasets/reads/sra_survey/pilot/wave1/`; the landing page links to it. The counts are of public runs, not population frequencies. The three curated collections (Fola 148 strains, *A. fumigatus* 331 read sets, *Coccidioides* 558 strains) have their own reports on the same pages, built from their committed tables; the old raw-reads table and Coccidioides chart on the landing page are replaced by a pointer.
- Logo (a lens finding the locus in the explorer's gene-map colours): `web/intro/logo/` with horizontal, stacked, mono and dark variants and a favicon, wordmark outlines traced from Nunito (SIL OFL). In the site header (theme-aware: the dark file for the dark and neuromancer themes), as the tab icon on the landing page, explorer, About page and browser tool, and at the top of the README.
### Changed
- Site wording: "Gating" becomes "Checks and confidence" (nav item "Confidence", table headers Check / What it tests, pipeline step "Verify"); the old `#gates` anchor still resolves. Wording only: no data or code change.
- Locus maps (explorer) mark the outcome of a resolved idiomorph pair: the gene that won the call keeps its solid fill and the gene set aside (`lost_to` in the report) is drawn pale with a solid outline and labelled "(set aside)"; the tooltip gives the basis and scores (HMM bitscore or alignment identity, for example "won over sexM (HMM score 364 vs 71)"), the locus line names the winner, and the legend has a new entry. The locus shards carry the verdict as a seventh item of a gene entry (`[W|L, other gene, basis, scores]`); a gene outside a resolved pair keeps the six-item form, so readers of the old shards are unaffected. Shards rebuilt for v0.7.0 data.
- Mucoromycotina tables: a blank order or family is filled from the genus (`scripts/mucoromycotina_genus_lineage.tsv`, NCBI taxonomy of the 51 genera, option `--genus-lineage`); a given order or family is never overwritten. The ZyGoLife LCG and Jena genomes had no order or family, and none of the Mucoromycotina inputs carried a family, so they sat under "(none)" in the explorer. Only the `order` and `family` columns of `genomes.tsv` and `loci.tsv` change (685 and 974 rows); calls, loci and `by_class.tsv` are identical. `build_site_data.py` now carries the family into the taxa rows and the genome tables (it was hard-coded blank), so the explorer groups Mucoromycotina by order, family and genus.


## [0.7.0] — 2026-10-26 — `v0.7.0`

The clean-run release. It adds the detection policies the curator rulings of 2026-09 to 2026-10 asked for (all **off by default**; the definitive run passes them with `--preset v0.7`, see `docs/notes/2026-10-22_clean-run-plan.md`), a scoped per-taxon roster, many reference-database records and exon structures, the browser page and landing site, and a series of fixes found by the pre-run regression gate. The package version moves from 0.6.3 (never tagged) to 0.7.0.

### Overview by area

#### Detection options (default off)
- `--partial-core-policy` with a companion-gene veto (#95, #102): a partial core hit against a full-length other-idiomorph gene is decided on the profile bitscore and coverage (Botrytis `MAT1-2`).
- `--strong-single-gene-policy` (#97), `--confidence-cap-policy` (#100), and `--homothallic-policy both_genes` (#117, #124: a locus with a complete alpha-box and a complete HMG-box gene is a `homothallic_candidate` with idiomorph `undetermined`, vote winner kept in the note).
- `--scoped-roster-policy` and a per-gene `taxonomic_scope` and `min_identity` in `order.yml` (#132, #136): GAP1 and BGL1 as optional one-member flanks scoped to Pleosporaceae (BGL1 floor 60%); recovers 90 of the 101 previously uncalled *Bipolaris* / *Curvularia* / *Exserohilum* genomes; `detect audit-scope` also lists records carrying a scoped gene from outside its scope.
- `--preset v0.7` (#158): one option sets the five policies above to the values that produced the v0.7.0 site data; an explicit policy flag overrides it; the defaults are unchanged.

#### Detection and reports
- Annotation tracks, gene-model tracks in the report and GenBank as a `detect` input (#69, #88); helixer-lite on locus windows and a flank-only helixer rule (#69, #70); `--polish-strong-identity` and a polished core pair on top of the cap (#40, #104).
- Parallel polishing (#68), an inline executor for one worker (#74), and exonerate aligned against only the requested gene's proteins (#83, about 6 times faster, same calls).
- Profile-HMM alpha-box and HMG-box gene-identity profiles for Ascomycota (#86); a subphylum guard that stops Pezizomycotina labels landing on Saccharomycotina yeasts, with MATyl / MATdc / MATkp / MATog families for Dipodascales and Komagataella / Ogataea (#106, #109).
- Optional gene-level `status: remnant|pseudogene` (#101, #103) and the schema value `mating_type.system: unknown` (#127).
- Remnant proteins leave the search references (#103); the `unisexual_switching` system value (#98).
- Defective gene models flagged (#152): a polished model with a frameshift or an internal stop codon, report-only. A report records the code it ran from (#157, #160): git describe, commit, dirty, and the source checkout's version.

#### Tools and data
- `scripts/deposit_to_record.py` (v2; #128): deterministic generator of core-only candidate records from GenBank deposits, with the curator rulings built in; it writes candidates only.
- Validator rule: every exon lies inside its gene's start/end (#126). Truth sets and held-out tests: BioSample truth (#82, #91), Ascomycota held-out setup (#79), coordinate truth.
- Reads typing datasets and k-mer panels for *Coccidioides*, *Cryptococcus*, *Fusarium* (#110, #112); per-family performance table (#108).
- `defective_genes` in the clean-run tables and the site data (#159).

#### Database records
- Exon structures for 119 reference genes (#119, #126); *Hemileia vastatrix* rustHD (#116); *C. fimbriata* (#98); *V. dahliae* (#62); *Fusarium oxysporum* (#29); *Mycotypha africana* with a classifier rebuild (#25); ten Dothideomycetes references (#36); Dipodascales, *Komagataella* and *Ceratocystis* curation (#98, #106, #109); `MAT1-2-3` in the Ascomycota vocabulary (#63).
- Corrections: *Tuber* `MAT1-1-1` A to S (#119, `record_version` 2 in #123); *Xanthoria* liq80xsp / liq146xsp `MAT1-2-1` (350 aa, 3 exons), *C. auris* B11243 full-length `MTLA1`, *Okayama* B43 notes (#133); *C. cinerea* B43 v2 (#56); the C4 / C5 *C. heterostrophus* records at `record_version` 2 (#132, #136).
- Four candidate records for the *Candida parapsilosis* group and *Aciculoconidium* (#147, #150); they do not affect detection until promoted.

#### Web
- Browser page and landing site: locus downloads, taxon results explorer, per-genome tables, locus maps with NCBI links, About page, theme switcher, *Coccidioides* read-typing results (#64, #67, #72, #73, #77, #78, #80, #84, #85, #87, #90, #111).

#### Changed
- Polish alignments are capped at a 5 kb intron (#92), settable per locus (`max_intron_bp` in `order.yml`) and 10 kb for the Basidiomycota bLocus (#153); a locus verdict sets a pair loser aside only where its winner has a scored protein and, for first-pass-only decisions, is not shorter than 0.8 times the loser's protein (#140; narrows the "winners only" rule of #64 that the regression gate found withheld Mucoromycota `Plus` calls).
- `taxonkit` failures are errors and not "taxonomy not current" (#93, #99).

#### Fixed
- A gene of the called idiomorph is no longer shown as lost (#94); `curate-db build-gff` keeps accession-less genes on assembly segments (#134, #121); the fast-path polish test no longer depends on thread order (#131).


### Added (details)
- `detect --preset v0.7`: one option sets the five policy options to the values that produced the v0.7.0 site data (`--partial-core-policy coverage_profile`, `--homothallic-policy both_genes`, `--scoped-roster-policy on`, `--confidence-cap-policy off`, `--strong-single-gene-policy off`). A policy flag given explicitly overrides the preset, whatever its value; the code defaults are unchanged (every policy off); the report's `run.parameters` records the preset and the resolved policies. A later release adds its own preset instead of editing this one. It applies to `matpredict detect` and to the scripts that call it (`run_clade_panel.slurm`, through `DETECT_ARGS`); the rollout path (`batch_runner.py`, `run_detection_batch.py`) passes none of the five policies, with or without a preset, as before.
- Defective gene models (#152, issue #151): a polished model whose CDS has a frameshift (exonerate) or an internal stop codon (any tool, read in the frame with the fewest) is reported as defective: per gene `defective` with `frameshifts` and/or `internal_stops`, per locus `defective_genes`, and in the HTML report. Report-only (curator ruling 2026-10-25): no tier and no `polished_genes` count reads it. It says the model does not translate cleanly in this assembly; a real pseudogene shows in many isolates of a species (*L. parapsilosis* `MTLA1` in 97 of 97 loci with the intact *L. orthopsilosis* record as reference), an assembly indel in few (conserved *C. auris* genes about 1-1.5%). The stop count is a lower bound for short models. In the pilot of the definitive run it fires on 7.5% of polished Ascomycota entries, mostly low-identity `exonerate` models of the flank genes `APN2`, `COX13` and `SLA2` in distant orders, where an alignment artefact is likelier than a pseudogene.
- The code a report ran from (#157, issue #139): `run.code` carries `git_describe`, `commit`, `dirty` (any change or untracked file under `src`, `db`, `scripts` or `pyproject.toml`) and `source_version`; absent when the package is not in a git checkout of this project; an unreadable checkout is logged.
- `defective_genes` in the clean-run tables and the site data (#159): a `defective_genes` column in `loci.tsv`, `n_defective_loci` in `genomes.tsv`, `defective_loci` in `by_class.tsv` and `summary.json`, a trailing note field in the locus shards and a trailing `n_defective_loci` column in the site genome tables. All trailing, blank or 0 for older tables.
- Four candidate reference records, not used by detection until promoted (#147, #150): *Aciculoconidium aculeatum* `MTLA1` from a BFD annotated genome; *L. orthopsilosis* MTLa (`HQ696681.1`), *L. metapsilosis* MTL alpha (`HQ696678.1`; the deposit's gene names say MTLa, the paper says alpha) and *L. parapsilosis* MTLa with `MTLA1` recorded as a pseudogene (`AY961981.1`), tier 2.
- `scripts/deposit_to_record.py` (v2): deterministic generator of core-only reference records (MAT1-1-1 / MAT1-2-1) from GenBank deposits, with the curator rulings of 2026-10-22 built in (tier by PubMed citation, `system: unknown` by default, two-gene, partial, doubtful-structure, Helotiales `MAT1-1-1` and non-binomial deposits held for hand review, curated records win on conflict). Writes candidates only; no record is added to `db/` by this change. Dry run over the 106 deposits in `results/2026-10-22_curate_deposits_pilot/`.
- `mating_type.system: unknown` in the metadata schema (curator ruling 2026-10-22): the breeding system is not established. Records generated from a single-gene deposit use it instead of a guess. No detection output reads the field.
- `detect --homothallic-policy both_genes` (default off; PR #117): a locus holding a complete alpha-box gene and a complete HMG-box gene is classed `homothallic_candidate` with a `homothallic_note`, and (curator ruling 2026-10-22) its idiomorph is reported `undetermined` with the vote's winner kept in `homothallic_note.vote_idiomorph`. A candidate, not a statement of homothallism. Measured on 3,181 genomes: 99 loci re-classed, 0 changes on the controls (`docs/notes/2026-10-21_homothallic-patch-measurement.md`).
- Optional gene-level `status: remnant|pseudogene` (and `remnant_note`) in the metadata schema; absent means intact. Metadata only: no detection output changes. `validate_gene_status`
  checks it against `present` (no model recorded vs a model recorded); the DuckDB `locus_gene` table gets `status` and `remnant_note` columns; `locus.gff3`/`locus.gbk` writers mark a status
  gene `/pseudo`; the self-consistency benchmark leaves such genes out of the expected set. Seven records are labelled (decision paper `docs/notes/2026-10-17_remnant-status.md`, option B step 1).

### Changed (details)
- A report's `run.matpredict_version` is now the version in the `pyproject.toml` of the source checkout the code runs from, when there is one (a `src/` layout whose project is MATPredict), and the installed package metadata otherwise. The metadata alone is stale on a shared environment whose editable install was never refreshed (issue #139: reports said 0.6.1 for code at 0.6.3). A wheel install is unchanged; `run.code` (git describe, commit, dirty, `source_version`) is unchanged. The HTML report's "rendered by" footer and `<meta generator>`, and `matpredict -V`, use the same version, so they agree with the report.

### Fixed (details)
- `curate-db build-gff` no longer drops the CDS feature and `proteins.faa` entry of a gene with no `protein_accession` on an `assembly`-type segment (#121).
  The protein of such a gene is translated from its recorded coordinates (exons, `codon_start`, `transl_table`); on an `assembly` segment the GCA_/GCF_ accession
  cannot be fetched, so the translation is now fetched from the segment's `seq_region` contig, as the locus sequence already was. A present gene with coordinates
  whose protein still cannot be derived now raises a `UserWarning` naming the record and gene. `validate`'s sequence cross-check is unchanged. One record was exposed
  (`5346_a43-b43-okayama-7_PR_B43`, genes phb3.1, phb1.1, phb2.3, phb2.1, phb3.3).
- The genes shown for a locus (GFF3, the figure, the track report) no longer include a gene of the called idiomorph marked as having lost a pair. A loss was recorded
  by gene name from every resolution event of the locus, although resolution happens per group of overlapping hits, so one weak hit elsewhere (in *Xylaria bambusicola*
  GCA_033171785.1, MAT1-1-5 beating a 15 kb-distant MAT1-2-1 hit) hid the polished MAT1-2-1 model. `lost_to_for_call` keeps only events whose winner belongs to
  the called idiomorph and whose loser does not. Calls, confidence, loci and gene models are unchanged. On 166 curated loci all reports are identical apart from
  `lost_to`; it changes in 37 records, and every one of the 71 changed labels is a called-idiomorph gene that was hidden and is now shown (one winner label changed).
  `idiomorph_resolutions` entries now carry `winner_span` and `loser_span`, where each verdict applies.
- Polish alignments are capped at a 5 kb intron (`miniprot -G`, `exonerate --maxintron`; `MATPREDICT_MAX_INTRON_BP` overrides). With no limit miniprot
  could join a weak gene to distant exons for a few extra points of score: in *Xylaria bambusicola* (GCA_033171785.1) the MAT1-2-1 model was 8.3 kb, with
  5.5 kb and 2.3 kb "introns" that held APN2 and COX13; it is now the 331 bp HMG fragment. The 609 introns in the curated records have a median of 58 bp and a
  maximum of 517 bp. On the 166 curated loci (leave-one-out) 163 reports are identical; the others change a model or a withheld span (one Basidiomycota
  pheromone receptor goes from a 19.6 kb exonerate model to a 1.2 kb one). The campaign counts above, 38 Ascomycota and 20 Basidiomycota loci whose core
  model swallows a flank gene, will change on a re-run.
- The intron cap is now settable per locus (`max_intron_bp` in `order.yml`; precedence: explicit argument, `MATPREDICT_MAX_INTRON_BP`, the locus value, 5 kb). The
  Basidiomycota bLocus (bE/bW) is set to 10 kb: the global 5 kb lost the whole b-locus call in genomes whose genes the aligners still model at 8-10 kb. Of the 7
  Basidiomycota genomes that lost a call to the 5 kb cap, 3 b-locus genomes return with ordinary models (largest intron 2.2-2.8 kb), 3 only return at 50-200 kb with
  26-78 kb chimeric models and stay lost, and the pheromone-receptor one needs 20 kb and stays lost (a known limitation: long-intron pheromone-receptor models are
  chimeric, as found for the 5 kb cap). On 216 Ustilaginales/Malasseziales genomes plus the 7 the only change is the 3 restored calls; on 500 random Basidiomycota
  genomes 490 are identical and the rest are one 10.5 kb bW model (7.7 kb intron, call unchanged), two small bW boundary shifts and one new low-confidence HD locus.
  Evidence: `results/2026-10-25_intron_cap/`.

### The definitive run
The v0.7.0 data come from the freeze tag `cleanrun-2026-10c` (commit `b317ece`), run with `--preset v0.7`, on database digest `285ee309...` (177 records; candidate records excluded). To reproduce the site data run `matpredict detect --preset v0.7 ...`. Run on 2026-10-26 in three throttled campaigns (48 CPUs); every report carries the tag, `dirty: false` and version 0.7.0. Tables, manifests and report archives: `results/2026-10-26_v070_*`.

| Campaign | Genomes listed | Called | Uncalled | No report (suppressed) | Loci | Loci high / medium / low |
|---|---|---|---|---|---|---|
| Ascomycota | 19,415 | 17,592 | 1,751 | 72 | 20,985 | 13,904 / 6,974 / 107 |
| Basidiomycota | 3,275 | 2,977 | 277 | 21 | 5,986 | 3,802 / 2,173 / 11 |
| Mucoromycotina | 974 | 868 | 105 | 1 | 914 | 689 / 212 / 13 |

Against v0.6.0 (the `change_vs_baseline` field of `genomes.tsv`): Ascomycota 150 calls lost and 472 new, 387 confidence up and 588 down, 219 idiomorph flips; Basidiomycota 6 lost and 8 new, 2 confidence down, and 7 genomes that v0.6.0 called now have no report (4 *Austropuccinia psidii*, 3 *Phakopsora pachyrhizi*: the BFD suppress list now excludes them for annotation-pipeline reasons, so MATPredict never ran on them); Mucoromycotina 3 new and 3 idiomorph flips. The reasons for the Ascomycota confidence drops are in the `tier_note` column (see Known issues).

Gates: every listed genome has a report except the same 72, 21 and 1 suppressed genomes as the first clean run (`cleanrun-2026-10b`); no report is empty. 36 of 23,570 reports differ from the first clean run (Mucoromycotina 0). 19 Basidiomycota are reproduced by forcing the per-family intron cap (#153) back to 5 kb, including 4 of the 5 changed calls; for one of them (GCA_011057475.1) a rerun of the first run's own code also gives the new result, so the split may be 18 and 12. The other 17 (5 Ascomycota, 12 Basidiomycota) are run-to-run variance in which of two equal-scoring copies of a gene is polished, reproducible with the first run's own code (#161). A 267-genome pilot run twice gave no differences.

The defective-gene flag (#152) marks at least one gene in 5,160 of 20,985 Ascomycota loci, 1,038 of 5,986 Basidiomycota loci and 156 of 914 Mucoromycotina loci. It is report-only. In Ascomycota most are low-identity flank-gene models; in Basidiomycota most are core genes (`HD1`, `beta_fg`, `HD2`, pheromone receptor), concentrated in Boletales and Wallemiales (see Known issues).

### Known issues
- Detection rules the curator has chosen to report and not change in this release: the idiomorph margin falls back to an identity gap for single-candidate loci (#137, #145), the vote on short fragments and tiny margins (#138), a noise third gene breaking the 0.5 tie (#129). The tables carry `tier_note` and `label_unsupported` for them.
- A better `MTLA1` hit can lower confidence because the unpolished-gene noise floor is a 39-bit cliff (#149); `MTLA1` is a polished model in 1,182 Serinales A loci but only a 20-28% fragment in 43 medium or low ones, and the rule is left as it is.
- The Basidiomycota b-locus intron cap of 10 kb (#153) trades three restored calls for one suspect 10.5 kb bW model in a 500-genome sample (#154); the pheromone-receptor genome GCA_900009495.1 stays lost.
- `defective` marks a model that does not translate cleanly in this assembly, which is not the same as a pseudogene (see above). In Basidiomycota it marks 1,038 of 5,986 loci and 794 of 1,191 flagged gene entries are core genes; 18 of 51 *Wallemia* genomes carry the same `STE3` model with one internal stop. Neither has been checked against the genome sequence.
- The rollout path (`batch_runner.py`, `run_detection_batch.py`) passes none of the policy options and writes reports without a `run` block.
- Rust genomes not searched (#165): 7 genomes that v0.6.0 called (4 *Austropuccinia psidii*, 3 *Phakopsora pachyrhizi*) have no report, and two more large rust genomes (a fifth *A. psidii*, *Hemileia* GCA_004125335.1) were never searched, because the BFD suppress list that the campaign reads excludes them for annotation-pipeline reasons, not MAT reasons. They are counted as "no new report", not as lost calls. Planned for the next release.
- 102 pheromone-receptor loci (94 genomes, 87 in Agaricales) that v0.6.0 withheld are now admitted as medium-confidence calls. Each carries a pheromone tblastn hit from the *Coprinopsis cinerea* Okayama-7 B43 record, which gained its pheromone proteins in the `build-gff` fix (#121). Accepted as known: some may count a pheromone gene and the receptor's own CAAX ORF as two genes. Not verified against the genomes.
- The 10 kb bLocus cap (#153) admits 6 bE/bW models in 5 genomes with introns of 5.1 to 9.6 kb at 12 to 33% coverage (the first clean run has none); they may be chimeric (#154).
- Run-to-run variance (#161): the same code on the same input can polish a different one of two equal-scoring copies of a gene (about 5 of 19,343 Ascomycota and 31 of 3,254 Basidiomycota reports differ in a model; one genome's calls changed).
- Held out of this release and tracked in #155: promoting the candidate records, the noise-floor rule, the margin rules, and the rest of the curator's list.

## [0.6.2] — 2026-10-15 — `v0.6.2`

Browser version (GitHub Pages), a profile-bitscore rule for sexM/sexP-style pairs, and a winners-only locus verdict. See the entries below and the pull request.

### Added
- `supported_span` on every called locus (report only; calls, clustering and the cluster span are unchanged): the extent of the
  locus's modelled genes plus its own hits at or above a bitscore floor (`--supported-min-bitscore`, default 33), with `beyond_supported_bp`. A strong hit of another family, or a weak hit of any family, that chains into the cluster
  stretches the cluster span but not this one; real own-family hits are kept (unlike `core_span`). In `detection_report.yaml`, on
  the GFF3 `MAT_locus` line (`supported_start`, `supported_end`, `beyond_supported_bp`) and in the HTML report ("Supported span",
  shown when the reported span rests partly on weaker hits). For a merged A/B call it is the union of the members'. Validation:
  `results/2026-10-08_supported_span/` (0 call differences on 34 genomes at five floors; every reported gene inside the span).
- `core_span` and `supported_span` also on every withheld locus in `suppressed_loci`, so span changes can be judged on withheld loci
  (`results/2026-10-08_supported_span_withheld/`: of 89 genome/family/contig keys whose cluster span changed with a curation edit,
  86 kept the same supported span).
- `core_span` on every called locus (report only; calls and clustering are unchanged): the extent of the locus's own
  gene models on its contig and `beyond_core_bp`, how much of `start`-`end` lies outside them. The cluster span is built
  from the hits of every family, so one weak hit from another family's short query can stretch a locus (T48-F: a called HD
  locus grows from 12.3 kb to 35.2 kb; see `analysis/2026-10-08_cinerea-b43-trace.md`). Written to `detection_report.yaml`,
  to the `MAT_locus` line of the GFF3 (`core_start`, `core_end`, `beyond_core_bp`) and to the HTML report ("Own genes").

### Fixed
- When two curated reference records tie exactly in alignment score for one gene, the polisher now takes the one with the higher
  identity, then the lower record id, instead of whichever alignment the tool listed first (the tool's order among equal scores is not guaranteed). Found on
  Leppa1, where one gene model was attributed to the A1163 record in most runs and to the Af293 record in some. Where the tied
  alignments differ in coordinates or identity, the reported model, and scores computed from it (idiomorph classifier scores),
  change too. On the 110 Ascomycota and Mucoromycota panel genomes, 8 reports differ from the old order-dependent output: 2 called
  loci (Leppa1: the record named; ASM2413939v1: a gene model 24,889-25,143 at 27.06% identity became 24,895-25,140 at 31.71%) and
  6 withheld loci (gene-model boundaries shifted by 6-36 bp; classifier scores changed, for example 36.8 to 39.4 and 67.3 to 64.4).
  No call, cluster span, gene set or merged_from changed on the 144 panel genomes, but a score near a threshold could change a call.
- PDF report: gene-evidence tables no longer run past their card. In print, header, number and status cells may wrap, and
  table cells break long names (for example `fungal_mating_type_pheromone`). Found on real campaign reports (Serpula
  lacrymans, 78 px past the card). New layout test `tests/report/test_report_layout.py` checks every fixture, including one
  real campaign report.
- `detect --genome` on a gzip or zstd file now stops at once with a message that says to decompress it (it used to fail
  minutes later with a BLAST or index error). Recognised by the first bytes, not the file name.
- `matpredict report table` writes `genomes.tsv` and `loci.tsv` from a directory of detect runs (optional sample-sheet join,
  rows for failed, skipped and never-run genomes). `matpredict report summary` groups the table by any column: recovery,
  idiomorph counts, both-idiomorph and homothallic-candidate signals.

### Changed
- The polish cap now has an identity tier. With the default cap of 6 per family, every admitted cluster whose best
  identity is 50% or more is polished even past the cap, and the remaining slots are filled in the usual rank
  (distinct genes, identity, hits); `--polish-strong-identity PCT` sets the threshold, 0 restores the plain cap. The plain
  rank put a true locus hitting 2 genes at ~100% behind six noise regions hitting 3 genes at 33-43%, so after the
  Dothideomycete records were added 15 genomes (10 *Zymoseptoria*, including IPO323) lost their call to the cap. Replay over
  2,722 Dothideomycete genomes and 2,582 true loci: 0 lost (the plain cap loses 17) at the same polishing work (16,284 against
  16,281 clusters); a cap of 15 would need 2.5 times the work. A genome with more than 6 strong clusters now polishes all of
  them (2 of 2,722). Calibrated on Dothideomycetes; see `analysis/2026-10-05_dothideomycetes-full-run.md`.

### Added
- `matpredict report genome --run DIR [--pdf FILE]`: a self-contained HTML report of one `detect` run (no network,
  no external files), printable to PDF from a browser or written with `--pdf` (WeasyPrint if installed, else headless
  Chrome/Chromium). It leads with the result in plain words (mating type and confidence; or "no call", "not searched",
  "two idiomorphs: needs review"), then one card per called locus: a gene-order figure (inline SVG; role, model status,
  strand, exons, alternate model, called-locus bracket, contig end), the idiomorph evidence and the gene table; then
  what was searched, withheld candidate loci, provenance and a glossary. Light and dark on screen, light in print.
  Revised after an independent web/data-design review (`analysis/2026-10-07_report-design-review.md`).
- `detection_report.yaml` gains a `run` block (first key): sample, organism, taxid, phylum, MATPredict version,
  database content SHA-256 and record count, taxonomy source, genome file name, SHA-256, contigs, length and N50,
  parameters, start time and wall time (`detect.provenance`). Additive; reports without it still render.
  `detect` gains `--sample` and `--organism` (shown in the report; routing still uses `--taxid`).
- `detect` writes `report.html` in `--out-dir` by default (curator's choice 2026-10-07: on by default, opt out).
  `--no-html` or `MATPREDICT_HTML=0` turns it off; `--pdf` adds `report.pdf`. A report or PDF error is logged and
  never fails the run. The batch scripts (`run_clade_panel.slurm`, `run_polish_ab.slurm`, `zygo_regression.py`,
  `run_holdout_benchmark.py`) set `MATPREDICT_HTML=0`: an environment variable rather than the flag, because they may
  run an older frozen worktree that would reject `--no-html`. `batch_runner` (in-process) writes no report.
- WeasyPrint (69.x, conda-forge) joins the pixi environment, `environment.yml` and so the Docker image, so `--pdf`
  works everywhere MATPredict is installed: 25 packages, 10.3 MB download (Pango, Cairo, HarfBuzz, fonts). Chosen
  over a headless Chromium in the image (about 300 MB). No locked version of any existing package changed. The
  Docker CI smoke test renders a PDF with no network.
- Report-only cassette fields on the receptor arrays (stacked on the arrays entry below; assessment
  `analysis/2026-10-06_b-locus-clustering.md`, option 1). Per array, per PR call and in `loci.tsv`:
  `receptor_cassette_loci` (loci with 2 or more strict-CAAX ORFs within 5 kb), `receptor_cassette_class` (`none`, `B`, or `C` when 2 or
  more of the ORFs also carry tblastn precursor homology; best over the array), `receptor_cassette_members` (locus and ORF
  coordinates, `|`-joined in `loci.tsv`) and `receptor_cassette_max_caax_orfs` (the maximum number of strict-CAAX ORFs within the window of any single locus of the array; not a count of cassettes). Reuses the existing
  strict-CAAX and precursor hits; descriptive only (circular for CAAX-admitted calls, class C a self-hit for species with
  curated precursors, tandem receptors may merge into one locus) and no call, tier, label, confidence or count changes. See
  `docs/receptor-arrays.md`.
- Report-only pheromone-receptor arrays (`detect.receptor_arrays`; study
  `analysis/2026-10-06_agaricomycetes-pr-arrays.md`, options 1 and 2). The STE3-like receptor hits of every family
  with a `pheromone_precursor_scan` (Basidiomycota PR) are merged per strand into loci and grouped into arrays (same
  contig, gap of 50 kb or less); a locus needs hits covering 50% or more of one reference receptor and a span of 8 kb or less. Each PR call gains `receptor_array_id`, `receptor_array_size`, `receptor_array_members`, `receptor_array_support`
  (`supported` when the array has 2 or more loci, a pheromone-precursor homology hit, or 2 or more distinct strict-CAAX
  ORFs; else `unsupported`) and `receptor_array_support_reasons`. `detection_report.yaml` gains `receptor_arrays` (each array
  once, with its call count) and `receptor_arrays_note`. A flag only: calls, tiers, confidence, verification labels and
  counts are unchanged. Arrays hold mating and non-mating receptors (paralogs sit beside the mating copies in
  *Coprinopsis cinerea* and *Schizophyllum commune*), so membership does not show that a locus is a mating receptor.
  The per-locus table (`loci.tsv`, written by campaign scripts) takes the same five columns; `receptor_arrays.loci_columns`
  builds them from a report entry. Per-array keys of `receptor_arrays` use the same `receptor_array_*` names. Arrays are for the PR family only (B-alpha/B-beta receptors form arrays only when merged into a PR call; B-locus sublocus structure is a planned assessment); the pheromone receptors are part of the B mating-type locus in tetrapolar Agaricomycetes. See `docs/receptor-arrays.md`.
- Offline NCBI taxonomy. `matpredict curate-db build-taxonomy` writes a slim
  all-taxa table (taxid, parent, rank, genetic code, scientific name; merged
  ids) from an NCBI taxdump directory, `taxdmp_*.zip` or `taxdump.tar.gz`
  (plain, `.gz` or `.zst` files). With `$MATPREDICT_TAXONOMY` set, lineage,
  phylum and genetic code come from it; `$MATPREDICT_OFFLINE=1` never calls
  NCBI. Checked against 8,185 cached efetch answers: phylum and genetic code
  identical for all; lineage identical for 8,179 (6 re-parented by NCBI between
  fetch and snapshot). Reports gain `taxonomy_source`.
- The Docker image builds the table from the dated archive
  `taxdmp_2026-10-01.zip` (SHA-256 pinned) and runs offline by default.

## [0.6.1] — 2026-10-04 — `v0.6.1`

### Fixed
- Genomes that use an NCBI genetic code exonerate has no built-in table for no
  longer fail. exonerate 2.4.0 has tables 1-6, 9-16 and 21-23; by id, any other
  code exited 1 and the whole genome got no report. Code 26 (Alaninales:
  *Pachysolen*, *Nakazawaea*; CUG = Ala) failed 20 of the 19,415 BFD Ascomycota
  genomes in the v0.6.0 run. exonerate now gets such a code as the NCBI table
  string (64 amino acids, TCAG order, from Biopython), which it accepts without
  any change to exonerate. Only a code with no NCBI table skips exonerate.

### Changed
- No default NCBI e-mail. The address comes from `$MATPREDICT_NCBI_EMAIL` or the
  `[ncbi]` table of `~/.config/matpredict/config.toml` (`$MATPREDICT_CONFIG`,
  `$XDG_CONFIG_HOME`). Without one, requests carry no e-mail and one warning is
  logged. Requests now send `tool=MATPredict`.
- The NCBI response cache key leaves out `email`, `api_key` and `tool`. Entries
  cached under the old full-URL key are still found and copied forward.
- README rewritten: goals, workflow, curation process, quickstart, usage,
  validation, citation.

### Added
- `pixi.lock` is now committed, so environments and images are reproducible.
- `environment.yml` for conda/mamba users (exported from `pixi.toml`).
- `Dockerfile` (pixi build stage, Ubuntu 24.04 runtime with the environment,
  the package and `db/`).
- `.github/workflows/docker.yml`: builds and smoke-tests the image on pull
  requests that touch the image inputs; on a `v*` tag, a published release or a
  manual run it also pushes `ghcr.io/stajichlab/matpredict:<version>` and
  `:latest`.

## [0.6.0] — 2026-10-03 — `v0.6.0`

Post-#9 work merged as PR #10: confidence and tier rules, the V3 polish cap,
the flank-carried rule, `not_searched` routing, the MAT-gene gate, the P1
paralog class, deterministic classifier builds, scope-only families,
`idiomorph_class`, the regression check, new curated records, and held-out
validation (LCG 536/621, Jena 61/64, Zygo 23/23). See the `v0.6.0` tag message
and `analysis/INDEX.md`.

## [0.5.0] — 2026-09-20 — `v0.5.0-mucoromycota`

Mucoromycota detection becomes usable end to end: idiomorph calling works, the
scoring denominator is honest, and the reference set covers both idiomorphs.

**Measured on 23 ground-truth Mucoromycota genomes** (`testset/Zygo`, curated
absence/presence tables):

| | before | after |
|---|---|---|
| loci on the correct scaffold | 23/23 | 23/23 |
| **idiomorph called correctly** | **0/23** (all `undetermined`) | **23/23** |
| Plus `fraction_found` | 0.571 | **0.833** |
| Minus `fraction_found` | 0.571 | 0.800 |
| roster genes with no reference | `algA`, `glrA` | **none** |
| confidence | — | 18 high / 5 medium |

**Measured on a 44-genus discovery sweep**: strict locus calls rose from
**6 genera to 14**, with `algA`/`glrA` — previously unfindable — present in 12.

### Added

- **Idiomorph resolution.** `sexM` and `sexP` share an HMG box, so one locus
  gene draws both curated references; this happened in 23/23 genomes and made
  the idiomorph uncallable. Overlapping pairs are now collapsed to one gene,
  with the loser annotated `superseded_by` rather than deleted so the ambiguity
  stays in the report.
- **Search-path tiebreak.** The winner is decided by *which search found the
  gene* — the annotated proteome outranks a genome-wide tblastn rescue —
  falling back to identity when both come from the same search. Identity alone
  scored 19/23 once new references were added; this scores 22/23, and 16/16 for
  both on the smaller reference set.
- **`locus_class`**: `mat_locus`, `homothallic_candidate`, `idiomorph_gene_only`.
  Says *what* was found, orthogonally to `detection_pass` (*how* it was
  admitted). A lone `sexM` or `sexP` with no flank is **kept deliberately** as
  training material for a future per-idiomorph HMM.
- **Homothallic detection.** Both idiomorphs in one locus is real biology —
  *Syzygites megalocarpus* encodes both HMG transcription factors, each flanked
  by its own intact gene with the other pseudogenised (Idnurm 2011). Requires
  both genes proteome-supported and within `max_homothallic_separation_bp`.
- **Relaxed second pass**, run only when the strict pass finds nothing
  genome-wide, admitting on the existing `EvidenceFloor` bar and capped at
  medium confidence.
- **Per-locus curation thresholds** in `order.yml`: `min_idiomorph_margin`
  (5.0) and `max_homothallic_separation_bp` (20 kb), alongside the existing
  `max_cluster_gap_bp`.
- **Report fields**: `idiomorph_margin`, `idiomorph_resolutions` (both members'
  identity, coverage and overlap), `locus_class`, `detection_pass`.
- **Diagnostics corpus**: every row now carries `run_id` and `genome_id`, and
  each idiomorph resolution is logged with the measurements a recalibration
  needs. Without `genome_id` no per-genome statistic could be computed from a
  batch at all.

### Fixed

- **`fraction_found` counted genes that could never be found.** `algA` and
  `glrA` were in the Mucoromycota roster with no reference protein anywhere in
  `db/`, so every genome was scored against an unreachable ceiling. The
  denominator now counts only genes present in the reference set the run
  actually searched with. Without this, resolving the sexM/sexP overlap alone
  would have moved Plus to exactly 0.500 — the rejection boundary.
- **A polish crash cost an entire genome.** `exonerate --refine` segfaults
  deterministically on some windows (verified: `--refine full` also crashes,
  no-refine succeeds, miniprot succeeds). Signal deaths now retry without
  `--refine`, then fall back to miniprot; non-zero exits stay fatal so a
  missing binary is still loud. Observed rate 1 in 44 genomes.
- **Relaxed calls carried no gene evidence** — no coordinates, identities or
  references — making them impossible to verify.
- **Resolution ran only before polishing**, but polishing moves coordinates, so
  pairs that overlap only afterwards escaped collapse. It now runs again after.
- **Alignments under 90 bp are dropped.** A 27 bp "gene" (nine codons) was
  being reported as a locus. The bar is on the alignment, well below the 60 aa
  short-ORF floor, so small MAT genes still pass. It is not an identity bar:
  Minus identities run 25.9–43.5%.
- **Reported evidence could be a superseded hit** while the gene counted as
  found via a different live one.
- **`score_match` passed on identity alone**, ignoring the coverage it
  computed: a record validated as `identity=100.0, coverage=0.86, pass`. Now
  warns below 80% coverage — a warning, not a failure, because 7 accepted
  records legitimately sit at 3.0–86.7%.
- **Evidence floor counted one gene as two** when a cluster's only hits were an
  unresolved `sexM`/`sexP` pair on the same protein.

### Database

Mucoromycota grew from 6 to 15 accepted records (9 Minus / 6 Plus, 2
homothallic); every gene validates at 100% identity **and** 100% coverage.
All are tier-1 published deposits.

| record | accession | contribution |
|---|---|---|
| *Mooraboolomyces wintlei* | `OR965930.1` | the only source of `algA` and `glrA` proteins |
| *Absidia urquhartii* ×2 | `PP971768/9` | first Cunninghamellaceae; first Minus outside *Mucor*/*Phycomyces*/*Rhizopus* |
| *Mucor mucedo* | `JN587498.1` | the alginate lyase the `algL`/`algA` question concerns |
| *Rhizopus azygosporus* | `MG967659.1` | second `glrA` |
| *Blakeslea trispora* | `HG939558.1` | `sexM` only — its `tptA` (38 aa) and `rnhA` (82 aa) are fragments |
| *Syzygites megalocarpus* ×2 | `JN112239/40` | **first homothallic reference**; third `glrA` |
| *Parasitella parasitica* | `KY081664.1` | Parasitellaceae Minus |

Assembly-derived loci for *Cunninghamella* and *Chaetocladium* were considered
and **rejected**: they are tier-2 homology-inferred, which the database design
spec excludes in this phase, and a literature search confirmed no MAT-locus
deposit exists for those genera (nor for *Actinomucor* or 11 others checked).

**Reference balance is not neutral.** Adding two Plus references flipped four
*Cunninghamella* Minus genomes to Plus; check per-idiomorph counts before any
ingest.

### Documentation

- `docs/notes/2026-09-20_algA-glrA-mucorales-literature.md` — establishes that
  `algA` and `algL` are one gene (62.8% identity between the two deposited
  proteins; the same author uses both spellings in one publication), that the
  2017 review carries **no citation** for either gene, and that `glrA` traces
  to Idnurm 2011.
- `docs/HANDOFF-detect-scoring-idiomorph.md` — every threshold with the
  measurement behind it, and the open questions.

### Notes

Every threshold is provisional and documented at its definition with its
evidence. `min_identity` remains `None`: Minus identities run 25.9–43.5%, so
any global cutoff would destroy Minus detection.

---

## [0.4.0] — 2026-09-20 — detection targeting (branch `detect-targeting`)

- Phylum routing: the query set narrows from 181 to 19 proteins for a
  Mucoromycota run; `--phylum` override; `routing_mode` in every report.
- `max_cluster_gap_bp` became per-locus curation data in `order.yml`
  (Mucoromycota `MAT` = 50 kb); a run takes the maximum over routed families.
- `EvidenceFloor` defaults to the curator's ruling: ≥2 distinct genes including
  ≥1 `core_MAT`.
- `run_batch` and a SLURM wrapper, so a batch pays the same narrowed cost.
- `--emit-cds-fasta`, an N+1 genome-parse fix and 60-column FASTA wrapping.

Validated on 23 ground-truth genomes at 100% per-gene recall, mean 41.8 s per
genome; the same pipeline pre-fix had been killed after >24 minutes on one
small yeast genome without completing.

## [0.3.0] — 2026-09 — localize-then-polish search

- Genome-wide `tblastn` localization followed by per-gene refinement with both
  `exonerate --refine region` and `miniprot`, classified into four explicit
  polish statuses; only genuinely unconfirmed genes cap confidence at medium.
- Fast-path rescue for core genes a supplied annotation does not contain — the
  blind spot the pipeline exists to close, since small pheromone-precursor
  genes are routinely missing from whole-genome annotation.
- Fragmented-locus handling: multi-segment calls with per-segment
  `contig_edge_distance`, scoped per contig for valid GFF3 parentage.
- Genome acquisition, SLURM-sized batch orchestration and rollout aggregation.

## [0.2.0] — 2026-09 — detection pipeline

- `matpredict detect`: family routing by taxid, diamond fast path against a
  supplied proteome, gap-based clustering, fractional per-family scoring with
  cross-family ambiguity detection, confidence tiering, idiomorph assignment.
- GFF3 and YAML report writers; families that fall short are reported as "not
  detected" with a reason rather than silently omitted.
- Short-ORF genes distinguished from real absences via `genes_not_searchable`.
- Leave-one-out benchmark suite.

## [0.1.0] — 2026-09 — curated reference database

- `matpredict curate-db`: propose → validate → accept/reject lifecycle, with
  GFF3, GenBank and protein FASTA export and a DuckDB query cache.
- Record schema with per-claim evidence tiers and citations, multi-exon
  `exons`/`codon_start`/`transl_table`, and independent translation
  verification — validation re-derives each protein from the record's own
  claimed coordinates rather than comparing a fetched protein to itself.
- Cached NCBI E-utilities and UniProt clients with retry on 429/5xx.
- `order.yml` controlled vocabulary per phylum, with `taxonomic_scope` routing
  and an optional `gene_class` for cross-family gene grouping.
- Locus and synteny diagram commands (`draw-locus`, `draw-synteny`).

---

## Project state at 0.5.0

- **68 accepted records** — 43 Ascomycota, 10 Basidiomycota, 15 Mucoromycota
- **19 loci** defined across three phyla
- **463 tests** passing, 46 test modules, 35 source modules
