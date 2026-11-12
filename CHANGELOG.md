# Changelog

All notable changes to MATPredict. Versions follow [semantic versioning](https://semver.org/);
release tags may carry a clade suffix naming the lineage whose infrastructure that
release completed.

## [Unreleased]

## [0.9.0] — 2026-11-12 — `v0.9.0`

> The entries below are those of the merged pull requests of the v0.9.0 cycle. The definitive run is described in the next section: it is the frozen commit `7d73be5` with `--preset v0.9.1 --withheld-detail-policy on`, and the 372 genomes that moved under the `tblastn` outlier policy come from reruns with it off. Rules that exist but are in no preset: `--weak-fragment-cap-policy`, `--weak-coverage-noise-policy`, `--strong-single-gene-policy`, `--confidence-cap-policy`.

### The definitive run (the 372 genomes that moved under the `tblastn` outlier policy are taken from reruns with it OFF; curator's choice, 2026-11-12)
Frozen commit `7d73be5`, `--preset v0.9.1 --withheld-detail-policy on`, no weak-fragment cap. 23,579 reports (19,343 Ascomycota, 3,263 Basidiomycota, 973 Mucoromycotina), every one stamped with the same commit, database hash and preset. Against v0.8.1:
- Mucoromycotina: 973 of 973 identical.
- Basidiomycota: 202 loci `high` to `medium` (#167), 10 scope labels (#170), one extra `kept_both` locus (GCA_017654745.1); no coordinate changes.
- Ascomycota: 108 loci become `undetermined` by the short-fragment vote (#138), 23 `high` to `medium` by the weak-flank cap (#129; the approved list), one extra `kept_both` locus (GCA_003706235.2), three *Pneumocystis* genomes gain `matMi` from the corrected record (#302), and GCA_052323795.1 moves its outer span by 16.7 kb from that same record (core and supported spans unchanged; traced by code/database split on 2026-11-12); approved by the curator, 2026-11-12). With the `tblastn` outlier policy off for the 372 genomes it moved, every other span and gene set equals v0.8.1.
- Calls: 17,618 Ascomycota, 2,995 Basidiomycota and 868 Mucoromycotina genomes called; 20,844, 5,994 and 914 loci.
- Determinism: 60 genomes run twice agree as sets; the order of `genes_found` can differ with the Python hash seed (one genome in 60).

### Added
- `--preset v0.9.1` (curator's ruling 2026-11-07): v0.9 plus three rules that lower a tier or withhold an idiomorph label: `--weak-flank-cap-policy on` (#129: 23 loci in the v0.8.1 run go from `high` to `medium`), `--short-fragment-vote-policy undetermined` (#138: 107 loci become `undetermined`) and `--single-region-cap-policy on` (#167: 202 Basidiomycota loci go from `high` to `medium`). The weak-coverage rule (#149: 109 loci `medium` to `high`) is held out for the curator's review. It is the preset of the v0.9.0 full run.
- `--preset v0.9` (curator's ruling 2026-11-07): v0.8.1 plus the three report-only rules of the v0.9.0 cycle: `--cross-family-policy keep_labelled` (#255: both calls kept where the kept call has no idiomorph label and the other has one; one genome in 23,579 v0.8.1 reports is affected), `--idiomorph-margin-policy separate` (#145: a single-candidate identity gap is reported as `idiomorph_identity_gap`; 74 of the loci in a 300-genome panel move the value), and `--scope-label-policy on` (#170: 10 `Basidiomycota:MAT` loci in the v0.8.1 run, on STE3, PAN6 and RPL39 alone outside the family's scope, are marked `family_unconfirmed`). No call, coordinate or tier changes. The three rules that lower a tier or withhold an idiomorph label (`--weak-flank-cap-policy`, `--short-fragment-vote-policy`, `--single-region-cap-policy`) are in `--preset v0.9.1`, not in `v0.9`. `--weak-fragment-cap-policy` and `--weak-coverage-noise-policy` are in no preset.
- `docs/detect-options.md` refreshed: the policy table now includes `--weak-fragment-cap-policy`, the `MATPREDICT_TBLASTN_SEG*` variables are listed, and what is merged but unreleased is separated from what was tried and not built (`--idiomorph-rank-basis`, `--polished-vote-policy`, with `docs/notes/2026-11-07_polished-vote-not-built.md`).
- Opt-in `--weak-fragment-cap-policy off|on` (#275; default `off`, in no preset): an `Ascomycota:MAT` call at `medium` or `high` whose strongest core gene (`MAT1-*`) is not a both-tools-agree model and scores under 50 bits is capped at tier `low`; the idiomorph, the locus and the class are never changed, and the report's `weak_fragment_cap` records the gene, its status and bitscore. On the v0.8.1 run it matches 108 of 12,051 loci and would change the tier of 68 (23 medium, 45 high). Details and the other triggers measured: `docs/notes/2026-11-06_weak-fragment-cap.md`.
- `docs/detect-options.md`: one page listing every `detect` option, policy, preset and environment variable with its default and what it does (tables generated from `--help`), and what is not in a release yet; `tests/test_docs_detect_options.py` fails when a policy, preset or environment variable is added without being listed.
- Opt-in `tblastn` filtering policy for the slowest genomes (#276): `MATPREDICT_TBLASTN_SEG=yes` with `MATPREDICT_TBLASTN_SEG_GENES=SLA2,APN2,MAT1-1-1` (or `MATPREDICT_TBLASTN_SEG_MIN_AA`) runs the genome-wide localization as two `tblastn` calls: `-seg no` for every other query and `-seg yes -soft_masking true` (masked regions used for seeds only) for the listed gene names. On the two genomes that hit the 3,600 s limit in the v0.8.1 run it takes 3,063 s and 1,685 s to 468 s and 288 s with identical loci; on 300 Ascomycota genomes 296 are identical (no change of family, idiomorph, confidence or gene set) and on 120 Basidiomycota genomes all 120. It gives no speed-up on typical genomes. The default is unchanged (`-seg no`). Measurements and what does not work: `docs/notes/2026-11-06_tblastn-seg-policy.md`.
- Opt-in `--weak-coverage-noise-policy off|on` (#149; default `off`, in no preset): an unpolished gene whose every hit covers under 30% of its reference is noise, whatever its bitscore, so it no longer caps the tier at `medium`. It acts only on a locus whose idiomorph is not `undetermined` and whose best core-gene coverage is at least 50% (curator, 2026-11-12). Tier only (`medium` to `high`). A first version without the two limits moved 109 loci from `medium` to `high` in 970 genomes of the v0.8.1 run; the limited version has not been run on a panel, so its count is not measured. It is not in the frozen v0.9.0 run.
- Opt-in `--withheld-detail-policy off|on` (#313; default `off`, in no preset; report only): with `on`, each `suppressed_loci` entry also lists the genes of the withheld locus (name, role, polish status, identity, coverage, bitscore). No call, tier or coordinate changes. The definitive run uses it.
- Label tables: a `nearest_reference` column (#200, #311; curator's ruling 2026-11-07). For each core gene it gives `GENE=SYSTEM`, the system code of the class of the matched record's organism, genes joined by `|`. A matched record whose class has no code is written as its `record_id`. `loci.tsv` gets a `core_reference` column with `gene:record_id:identity:bitscore` per gene. Tables only; no report, call or tier changes.
- Raw-reads gate causes (#195, #196; PR #325). `datasets/reads/sra_survey/pilot/causes.tsv` holds one cause per species: `homothallic` and `silent_cassettes` give `NOT_TYPABLE`, `zygosity` joins the zygosity gate, and a species with no cause is `UNVALIDATED:cause_unchecked`. The 50%-`both` shortcut no longer gives `PROVISIONAL_HOMOTHALLIC`. `cause_evidence.tsv` records the evidence for 42 species. The reads pages show the cause wording of the ruling of 2026-11-09.
- Raw-reads gate, thin evidence (#327; curator's ruling 2026-11-09): a species with one strain or one to three runs is typed from reads and capped at `UNVALIDATED`, with the reason shown, instead of being held. A species is held only when an idiomorph has no locus or no run.
- Raw-reads clone groups (#320) use the error-tolerant solid-containment distance of `sketch_query.py`: 202 species, clone ids at 0.9999, 0.9995, 0.999 and 0.998. The site uses these tables in place of the older plain-containment tables.
- *A. fumigatus* population report from the sketch store (#326): 4,849 isolates (4,736 kept); MAT1-1 in 50.4% of clones (95% CI 48.2-52.5%, 2,077 clones). Also an ENA FTP fallback for runs where ENA answered HTTPS with 403 (#326).
- Sketch store (`datasets/reads/sketch_store/`; #315, #319, #329, #330, #333): one `sourmash` sketch per isolate (k=21 and k=31, scaled 100) answers MAT, target-screen, population-structure and QC questions without downloading reads again; Starship and bafA prevalence, co-occurrence and gain/loss on the tree; a baf family screen. Research tooling; `detect` does not use it.
- Landing page: Fola and *A. fumigatus* read-typing charts, tables and notes (#301). They describe the v0.7.0 reads code. #172 stays open for the Verticillium, Cryptococcus and other Fusarium rows.
- FungiDB annotation check (#290, #336): `scripts/annotation_by_tier.py` and `docs/notes/2026-11-12_290-annotation-by-tier.md`. No gene-level verdict changed on the 1,234 genes present in both the v0.8.1 and v0.9.0 runs; 15 genes appear only in v0.9.0 and none is lost. Annotation is a cross-check, not the truth. #290 stays open for the remaining steps.

### Changed
- Clean-run tables: the Basidiomycota defective flag is a model check (#168, #304; curator's decision 2026-11-07; report only). `defective_genes` lists only a core gene at identity of at least 60% and coverage of at least 80%; every other flagged gene moves to the new column `low_confidence_model_genes`. In the v0.8.1 run, 1,040 loci were flagged; 153 keep a listed gene and 887 moved. Outside Basidiomycota nothing changes. Tiers, calls and reports are untouched.
- Speed: the curated YAML (`order.yml` and the `metadata.yaml` records) is parsed with libyaml's C loader, with the pure-Python loader as fallback (#287, #292). Detection is identical to the frozen v0.8.1 reports on 823 genomes.

### Fixed
- An EMPTY taxonomy answer is now recorded (#324): when the lineage and phylum lookups return nothing (no exception, for example a cache miss on a compute node without network) the report's `routing_error` says so and the not-searched reason reads "taxonomy lookup failed (...)", where it said only "no taxid, or no phylum could be resolved". Seven *C. fimbriata* assemblies (taxid 5158, a valid taxon) were left unsearched this way in the v0.9.0 run. Routing is unchanged. The panel script's taxonomy warm-up now counts an empty answer as a failure, warns, and stops the job when `TAXONOMY_WARMUP_STRICT=1` (default `0`: warn only). Tests: `tests/detect/test_empty_taxonomy_lookup.py`, `tests/test_clade_panel_warmup.py`.
- `genes_found` listed optional genes of different slots (GAP1 and BGL1 in *Cochliobolus*) in hash-seed order, so the same genome could report them in either order from run to run (found 2026-11-12 in the 60-genome determinism check of the v0.9.0 run: 1 genome). They now follow the roster order. Output order only; no call, tier, label or coordinate changes. Test: `tests/detect/test_genes_found_order.py`.

## [0.8.1] — 2026-11-06 — `v0.8.1`

> The entries below are those of the merged pull requests; the definitive run and its counts are in "The definitive run" below. v0.8.1 is a re-run of the definitive v0.8.0 run with: the #228 fallback after withheld hits (preset `v0.8.1`), the *S. pombe* silent-pair rule (#254; the #252 tie-break is reverted), report-only `resolution_flags`, the cross-family `cassette` indication with its shared-flank measurements, `two_idiomorphs` for Basidiomycota MAT and redPR, the rust allow list (#165) and the *Hyphopichia* MTLA1 record (#148; 11 of 389 sampled Serinales genomes drop a tier, tracked under #149). Not in v0.8.1: #255, the MATam family, the Ogataea candidate records.

### Changed
- *S. pombe* silent cassettes (#254; curator ruling 2026-11-05): a merged `mat1` call whose genes hold a P cluster and an M cluster 5-20 kb apart on one contig (the silent `mat2`-P / `mat3`-M pair, about 11 kb apart) is split into a `mat2` P call and a `mat3` M call, each with a `cassette` block (role `silent_pair`, the partner family and the distance). The #252 tie-break by core identity is reverted (it relabelled a lone P cassette `mat2`, probably wrong: the records are cassette-only so identity cannot tell mat1 from mat2); a lone cassette elsewhere stays `mat1`; a pair farther apart (35 kb in GCA_006386235.1) stays one call marked `arrangement_ambiguous`. It applies only where the species-scoped `mat2` and `mat3` families are in the genome's lineage (*S. pombe*), as part of `--cross-family-policy on` (`detect/silent_pair.py`). It changes *S. pombe* locus counts and labels.

- Basidiomycota:MAT and Basidiomycota:redPR now report `two_idiomorphs` (#169; curator, 2026-11-05; `two_idiomorphs_report: true`): a genome that calls both idiomorphs of one of these families gets the report-only statement already used for Mucoromycota:MAT (arrangement, evidence, possible causes; no cause asserted, no call, idiomorph or confidence changes). On the nine v0.8.0 genomes with both idiomorphs (six MAT: *C. neoformans* x2, a *Mrakia* and others; four redPR) every one now carries the statement; one is `same_locus`, the others `unlinked`.

### Added
- Cross-family `cassette` indication (curator ruling 2026-11-04, #189): a genome with calls of two families that have complementary idiomorphs (MATa-type and MATalpha-type) at different positions marks each call `cassette` with `cross_family: true`, `partner_families` and role `cassette_undetermined`, as the within-family MATog pair already is. Report only; no call, family or tier changes. Each block carries `shared_flank_check` measurements (`detect/shared_flank.py`, blastn on each locus +/- 20 kb, every partner checked): status, aligned bp, best identity, `both_sides` (shared sequence left and right of the idiomorph's genes in both loci) and `repeat` (a flank segment found at a third place in the genome). The role stays `cassette_undetermined` whatever is shared, except that it becomes `allelic_pair` (one heterozygous locus, e.g. *Wickerhamomyces edaphicus*) when the sequence is shared on both sides and is not a repeat; real cassette systems share sequence too (*Saccharomycopsis schoenii*), so sharing alone is only a measurement. `not_run` without blastn.

### Added
- A MATPredict-only allow list (#165): `db/allow.txt` removes genomes from the merged suppress set, read by `detect suppress-filter` and the batch scripts. It lists the nine rust genomes (five *Austropuccinia psidii*, one *Hemileia vastatrix*, three *Phakopsora pachyrhizi*; 0.56-1.27 Gb) that the shared BFD list suppresses for annotation-pipeline reasons ("EVM never completes", "too large", "too few models") and that MAT detection can search; run them with a long `GENOME_TIMEOUT`. GCA_943181505.1 ("Too small, 800kb") stays suppressed.

### Added
- Database: the candidate CTG-clade MTLA1 record for *Hyphopichia burtonii* (`717740_makgeolli_MTL_A`, GCA_003856795.1; curator decision 2026-11-05, #148). On a 389-genome Serinales sample it raises 9 genomes' confidence, lowers 11 (the 39-bit noise floor of #149 and the 10-point cross-hit guard acting on weak fragments; see `results/2026-11-04_148_hyphopichia`) and leaves 354 unchanged; the drops are known and tracked in #149.

### Added
- Report-only `resolution_flags` on every detected locus (curator ruling 2026-11-04, #256): `contig_edge` (a core gene within 500 bp of a contig end) and `partial_core` (a core model covering under 50% of its reference at under 90% identity; a high-identity exon piece of an intact gene is not flagged), with the genes in `resolution_flag_genes`. They say a called locus may be a fragment; no call, score or confidence tier changes. New module `detect/resolution_flags.py`; `DetectionOutcome.contig_lengths`.

### Fixed
- Browser tool: a single failed download no longer ends a run. A 429, 5xx or network error on any asset (the database files, the wheels, the WebAssembly tools) is retried three times with a growing wait, or the wait the server asks for in `Retry-After`; each retry appears in the log. A permanent failure such as a 404 stops at once with a message that names the file and the status ("could not fetch matpredict.whl: HTTP 503, tried 4 times"). A first-time failure showed up as `HTTP 503` on `matpredict.whl` during a second run in the same tab, when the browser's 10-minute cache had expired and about 800 files were requested again at once.

### Changed
- Raw-reads typing, species gate (curator, 2026-11-02): a species with a usable panel is no longer hidden for missing a strict limit. A new decision, `UNVALIDATED`, covers a species that types at least half of its finished test runs and whose only failures are soft (typed below 75% or 85%, `both` above 20%, shared depth below 3, or agreement short of the limit with fewer than ten comparison runs). It gets the full run, and the site shows it as "unvalidated" with its reasons. `HOLD` ("withheld") is left for what makes calls unusable: no panel or too few runs or loci, under half the runs typed, too many failed downloads, or agreement short of the limit on ten or more comparison runs. `wave_driver.py` has the rule as one function (`decide`), a `regate` command that decides again from the numbers in a gates table, and `submit-full` includes the unvalidated species. Applied to both waves: 28 of the 79 held species become unvalidated (20 in wave 1, 8 in wave 2); no other decision changes (`wave*/regate_changes.tsv`). 51 stay withheld: 16 panels that did not build, 14 under half typed, 15 with too few runs or loci, 3 panel errors, and 3 with other hard failures. The overview and species pages gain the category, with its definition and, on a species page, why.
- Raw-reads overview (`reads/`): the species from public runs are grouped by taxonomy, phylum then class then order, as collapsible sections that each show their number of species and runs and how many are confirmed, provisional or withheld. Phyla and classes start open and orders closed; a search or filter opens the matches, and Expand all and Collapse all are above the list. A "View" box switches to one flat table, which, like the grouped tables, has a Family column and still sorts on any column. Each species page shows its lineage under the title. The lineage comes from the genome tables of the site's default release (`scripts/build_reads_site.py` looks a species up by its name, else its genus) and is written into the data as `phylum`, `class`, `order` and `family`. All 197 species and the three curated collections are placed, all in Ascomycota; the Basidiomycota and Mucoromycota panels are not done yet. The wording about "wave 1" and "wave 2" is gone from the public pages.
- Raw-reads overview (`reads/`): the species table no longer has columns headed "1st" and "2nd". One column, "Alleles called (runs)", shows each species' own allele names with the runs typed as each (for example `MAT1-1 4 · MAT1-2 5`, or `a 3 · alpha 4` for a yeast) and sorts by their total. The confidence chart's "1st allele" and "2nd allele" are one category, "one allele called", since the species differ in their allele names. The species page headers already named the alleles; their placeholder before the data loads no longer says "1st" and "2nd".
- Site and README: the genomes called "BFD" (the project's internal name) are now named "NCBI public genomes" on the landing, results and About pages, in the dataset labels of `summary.json` and in the README results table. The landing and explorer charts widen their label column to fit the longer name. Internal names are unchanged: the `BFD` source value in the data, the `bfd_*` dataset ids and the `genomes_bfd_*.tsv` download names, so existing links and scripts keep working.
- Browser tool: the database reaches the page as one archive (`db.tar.gz`, 1.5 MB, with its SHA-256 in `db-manifest.json`) instead of 767 separate files. A run makes one request for it, checks the size, checksum and file count, and unpacks it in about 0.1 s; before, the same step made 767 requests in about 14 s, the burst that drew the 503. The site build writes it with `web/spike/page/build_db_manifest.py` and no longer ships the loose files under `_assets/db`. A manifest without an archive (the spike's `assemble.sh`) still mounts file by file.

### Added
- Detection (#228): `--fallback-after-withheld-policy on` (off by default; on in the new preset `v0.8.1`). When a genome routed by taxid or lineage reports no locus because the modelled-gene bar withheld every hit (an *Ogataea* genome with only withheld `MATog` hits), the phylum-fallback families (`MATsc`, `MTL`) are searched as well; their calls follow the usual out-of-scope rules and the report's `fallback_search` says what was added. It does not apply when the routed family reports any locus, and needs `--taxid`.
- Browser tool: the upload dialog also offers `.fsa`, `.seq` and `.dna` (FASTA files; GenBank files are recognised by their content, so a `.seq` or `.dna` GenBank file works too). `detect` on the command line never filtered by extension.
- Raw-reads typing: an optional contamination screen (`SCREEN=1` in `genotype_mat_one_accession.sh`). GmG classifies the same 2M-spot subsample, and `screen_reads.py` flags a read set that is `mixed` (GmG's own mixture flag), a `genus_mismatch` (the whole set is another organism) or `borderline` (a second genus that looks real but is not flagged, which is read deeper, 6M spots, and screened again). The set is still typed. The result is `<run>.screen.json`, the screen columns of `results_all.tsv`, and, on a species page, the Note column and a "runs flagged" card. Off by default, so existing runs are unchanged. The rule and its thresholds come from the spike-in (`results/2026-11-02_contam_spikein/`); on the 558 *Coccidioides* read sets it flags 3 wrong-genus sets, 20 mixed ones (the audit's TX7 and 19 more, four with 31-37% *Histoplasma*) and 4 borderline.
- Raw-reads screen: `no_call` now means GmG made neither a genus nor a species call. The first real wave flagged 16 Clonostachys rosea runs as `no_call` although GmG called the species.
- Raw-reads typing, wave 1 and 2: the 28 `UNVALIDATED` species were run in full with the contamination screen (500 runs; 144 new runs and the screen columns merged into `results_all_final.tsv` by `merge_screen.py`). 28 runs are flagged (mixed 16, genus_mismatch 9, borderline 3). With the added runs four species move up: *Trichoderma viride* to PASS; *Pyrenophora teres*, *Blastomyces dermatitidis* and *Candidozyma vulturna* to PROVISIONAL.
- Raw-reads species pages: every table sorts when you click a column header (click again to reverse; a ▲ or ▼ marks the order). Numbers sort as numbers, run accessions in natural order (SRR2 before SRR10), and empty cells last. The order is kept when the runs table is filtered or the counts are regrouped, and the pooled row stays last. The species overview table on `reads/` shows the same mark. The sorter is `web/intro/sortable.js`, tested in `web/intro/sortable.test.mjs`.
- The pre-publish site check (`web/spike/page/check_site.py`) also checks the site's HTML pages: every stylesheet, script and module import they refer to must be in the built site.
- Pages workflow: a manual run builds and checks the site and keeps it as the `site-preview` artifact without publishing it (the `build_only` input, ticked by default; a `v*` tag still publishes), and every run now checks the assembled site first (`web/spike/page/check_site.py`: required files, every module and file the page refers to, database archive against its manifest).
- Browser tool: a "Guess phylum with GmG" button next to the phylum pulldown. It runs [GmG](https://github.com/stajichlab/GmG) (WebAssembly, in a worker) on the loaded genome and selects the phylum; a weak or non-fungal guess leaves the pulldown unchanged and says why ("Accept weak guesses" is off by default). The 67 MB reference is downloaded once and cached in the browser; the genome is not uploaded. Same rules as `detect --phylum auto`.
- Browser tool: a Reset button (reloads the page, which frees the interpreter and the mounted database of earlier runs; choose the genome again) and a Copy log button (the log with the page address, browser and time, for a bug report). The log is cleared and its timer restarts at each run.
- `matpredict detect --phylum auto` (and `auto-weak`): ask [GmG](https://github.com/stajichlab/GmG) which phylum the genome is in, then search that phylum's curated families (the order or family taxid GmG gives is tried first, so the lineage routing applies, and the genetic code is read from it). A weak call, a genome GmG says is not a fungus, no call, or a phylum with no curated family is NOT searched and the report says why (`parameters.phylum_guess`); `auto-weak` also accepts a weak call; `--exhaustive` still forces a search. An explicit `--phylum NAME` or `--taxid` is unchanged and wins. `--gmg-index FILE|NAME|latest|URL` picks the GmG reference index (downloaded once to `~/.cache/matpredict/gmg`, checked against GmG's published sha256, pinned by default, with a warning when a newer one exists). Needs the `gmg` binary or GmG's Python bindings; needs GmG with `routing_hint` for the full rule (an older GmG is read with the same rule). Design: `docs/notes/2026-11-02_gmg-phylum-routing.md`.

### The definitive run
The v0.8.1 data come from commit `50bc950` (main after the merge of #272; `src`, `db` and `scripts` as released), run with `--preset v0.8.1` (the v0.8 policies plus `--fallback-after-withheld-policy on`) and, for the nine rust genomes of the allow list, `GENOME_TIMEOUT=14400`. The lists are those of the v0.8.0 run: 19,415 Ascomycota, 3,275 Basidiomycota and 974 Mucoromycotina genomes. Two Ascomycota genomes (GCA_016165855.1, GCA_023522825.1) hit the 3,600 s limit on an older node and were re-run alone on an EPYC node with the same code and preset (3,063 s and 1,685 s); the tables take those reports. Folders: `results/2026-11-05_v081_{ascomycota,basidiomycota,mucoromycotina}`.

| Campaign | Genomes listed | Called | Uncalled | No report (suppressed) | Loci | Loci high / medium / low |
|---|---|---|---|---|---|---|
| Ascomycota | 19,415 | 17,618 | 1,725 | 72 | 20,843 | 14,623 / 6,120 / 100 |
| Basidiomycota | 3,275 | 2,995 | 268 | 12 | 5,993 | 3,822 / 2,161 / 10 |
| Mucoromycotina | 974 | 868 | 105 | 1 | 914 | 689 / 212 / 13 |

Against v0.8.0 (23,570 genomes with a report in both runs): 23,214 are identical in calls, confidence, coordinates and gene models and 356 differ; the nine allowed rust genomes have no v0.8.0 report. No call is lost, no idiomorph changes and the Mucoromycotina campaign is unchanged. Every difference was traced with control runs on the frozen code (`results/2026-11-06_attrib_v081` on the HPCC), and nothing is unexplained:
- *S. pombe* silent pair (#254, #274): 37 Ascomycota genomes (loci 57 to 97). The merged `mat1` M call becomes a `mat2` P and a `mat3` M call (40 pairs); 7 `mat1` M and 10 lone `mat1` P calls stay. No other family, idiomorph (outside the split) or confidence changes.
- *Hyphopichia* `MTLA1` record (#148, #273): with the record removed from the database (v0.8.1 code and preset), all 302 affected genomes return exactly to their v0.8.0 report. 190 genomes change a call field: 4 are new calls (*Metschnikowia*, *Cephaloascus*, *Saccharomycopsis*, *Dujonia*), 160 have a changed locus span (mostly Serinales `MTL`, also a few Pichiales and Phaffomycetales `MATsc` genomes searched with the same query set), 36 loci rise and 14 drop a tier (the #149 floor and cross-hit guard: 7 `distant_fragment`, 1 `polished_model_lost`, 6 without a listed mechanism). 112 more genomes keep the same calls and differ only in which reference record supplies a gene model. *Cephaloascus* GCA_030575015.1 is called `MTL` A (medium) through the record; without it the #228 fallback would call `MATsc` MATa instead.
- Fallback after withheld hits (#228): 17 genomes are called that were not before (7 Ascomycota: *Naumovozyma*, *Ogataea*, *Babjeviella*, and one each of *Terfezia*, *Monosporascus*, *Biscogniauxia*, *Ascodesmis*; 10 Basidiomycota, 15 loci: *Tremellales* `PR` and `HD`, one *Puccinia striiformis* `PR` and one *Meruliopsis* `bLocus`). In each, the report's `fallback_search` names the withheld routed family; running the Ascomycota ones with `--preset v0.8` brings back the v0.8.0 report. The fallback searched 1,147 Ascomycota and 69 Basidiomycota genomes and found a call in these 17.
- Rust allow list (#165, #269): the nine genomes are searched; eight are called (`rustHD`, 14 loci) and one (GCA_001443065.1) is uncalled. The 12 remaining suppressed Basidiomycota genomes have no report.
- Report-only changes with no effect on a call: `resolution_flags` is set on 11,857 Ascomycota, 3,489 Basidiomycota and 236 Mucoromycotina loci (#257); `two_idiomorphs` is reported for 9 Basidiomycota genomes (#271); the cross-family `cassette` indication with shared-flank measurements (#264) changes no call or family label in any campaign.

The determinism check (60 Ascomycota genomes drawn with seed 20261106, 58 of them with a report; run twice on this commit with `--preset v0.8.1`) gave no differences in family, idiomorph, confidence, contig, start or end between the two runs or against the campaign; the order of the gene names in `genes_found` differed in one genome between the campaign and the re-runs, which agree with each other (`results/2026-11-06_determinism_v081` on the HPCC).

## [0.8.0] — 2026-10-31 — `v0.8.0`

(Absorbs the planned v0.7.1 site release: the labels and reads pages, the yeast routing and family work (#188, #189) and the browser examples. Unlike v0.7.0, `src` and `db` change: the new `--cross-family-policy` and `--preset v0.8` (off unless asked), a wider MATog scope and a new MATpk family change which families search Ogataea and *P. kudriavzevii* genomes. The site's data and counts are those of the v0.8.0 run (below), made on commit `a2933f7` with `--preset v0.8`.)

### Added
- Idiomorph labels beside the raw names (output only; curator rulings of 2026-10-28, `docs/notes/2026-10-28_idiomorph-naming-proposal.md`): grammar `<Lineage>__<System>_<Cell>` with frozen five-letter lineage codes and reference-organism system codes, in `scripts/idiomorph_label_registry.json` (single source of truth, with cell tokens and gene maps), `scripts/idiomorph_labels.py` and `scripts/label_loci_table.py`; all 27,885 v0.7.0 loci labelled (`results/2026-10-28_idiomorph_labels_v070/`). The site data carry `label_alleles` in `taxa.json`, `call_label` / `alleles_label` in the genome tables and `data/idiomorph_labels.json`. The explorer and the reads pages show labels (for example `Taphr__SchPom_M`, `Sacch__CanAlb_alpha`, `Peziz__AspFum_MAT1-1`) beside the raw names, with the raw name and the reference organism in the tooltip and an "Allele names" switch back to the raw names. Labels are derived from the family called and the raw name; the database, the reports, the calls and the raw columns are unchanged. Duplicate family calls (`scripts/dedupe_family_calls.py`) are not removed from the counts.
- `scripts/dedupe_family_calls.py` marks loci called by two families at identical coordinates (185 loci in 164 genomes of v0.7.0) and `scripts/idiomorph_crosscall_flags.py` lists the 181 yeast genomes with loci from two MAT systems (`results/2026-10-28_idiomorph_crosscall_flags/`); neither changes a call (routing re-run: #188, curation: #189).
- Taphrinomycotina MAT candidates (466 proteins, polishing probe, Pneumocystis synteny; `results/2026-10-28_taphrinomycotina_mat_candidates/`): analysis only, no records added (#202).
- Clone groups from read sketches (`datasets/reads/sra_survey/pilot/clone_pilot/`, `clone_clusters.py`), the raw-reads orchestration scripts archived under `datasets/reads/sra_survey/pilot/orchestration/` (#204) and tests for the reads page builder.
- GenBank input in the browser page: upload a `.gb`, `.gbk`, `.gbff` or `.genbank` file (optionally gzipped) and its sequence is searched like a FASTA upload; contig names are the VERSION accessions, as on the command line. The annotation is not used in the browser (the protein route needs DIAMOND, which has no WebAssembly build). On the command line, a GenBank genome with sequence but no CDS translations (an unannotated deposit) is now searched as a genome, with a message, instead of stopping with an error; an annotated file still takes the proteome route.
- Platforms: `pixi.toml` and `pixi.lock` now cover `osx-arm64` and `linux-aarch64` as well as `linux-64`; a `tests` workflow runs the suite (`pixi run -e test pytest tests`) and a smoke test on Linux and macOS; the Docker image is built for `linux/amd64` and `linux/arm64`. The re-lock also moved the `linux-64` pins, so run the regression check against it.
- Five accepted tier-3 records (curator sign-off 2026-10-30, J. Stajich): *Saccharomycopsis schoenii* GDD12 MATss (MATa-MATalpha left cluster and MATa right cluster, `db/Ascomycota/Ascoideales/`), *Barnettozyma californica* PB4207 MATbz and *Petasospora fabianii* MATpf (JOY008 MATa, TBRC 4498 MATalpha, `db/Ascomycota/Phaffomycetales/`). The gene models are our own (RNA-seq for *S. schoenii*, sequence-based for the others). With them, 20 of the 21 genomes that lose every MATsc/MTL call under the v0.8 cross-family rules are called again by their own family (trial run, `--preset v0.8`); *Ogataea degrootiae* GCA_019232045.1 stays uncalled.
- Raw-reads wave 2 (`datasets/reads/sra_survey/pilot/wave2/`): 55 more tier A, A1 and B species typed (1,468 runs; 13 confirmed, 16 provisional, 7 provisional-homothallic, 49 withheld), yeast idiomorph names (a / alpha, P / M, A / B) carried through, *Botrytis cinerea* typed on all 333 eligible samples. Species reports gain an allele-count table for all studies pooled or study by study, per run or per clone, and clone groups for every species with 10 or more sketched runs (average linkage at 0.998; counts per clone only where the clone number is stable). Output only; the genome tables and the database are unchanged.
- Raw-reads pages (`web/intro/reads/`): an overview of the wave 1 read typing (112 species, 1,629 runs, 599 of 605 single-allele calls agree with an assembly call) with the confidence of each species panel (confirmed / provisional / withheld), and one report per species (calls, agreement with assemblies, calls by study, every run, TSV download). Data written at site-build time by `scripts/build_reads_site.py` from `datasets/reads/sra_survey/pilot/wave1/`; the landing page links to it. The counts are of public runs, not population frequencies. The three curated collections (Fola 148 strains, *A. fumigatus* 331 read sets, *Coccidioides* 558 strains) have their own reports on the same pages, built from their committed tables; the old raw-reads table and Coccidioides chart on the landing page are replaced by a pointer.
- Logo (a lens finding the locus in the explorer's gene-map colours): `web/intro/logo/` with horizontal, stacked, mono and dark variants and a favicon, wordmark outlines traced from Nunito (SIL OFL). In the site header (theme-aware: the dark file for the dark and neuromancer themes), as the tab icon on the landing page, explorer, About page and browser tool, and at the top of the README.
- Cross-family calls and the `family_unconfirmed` flag (#188, #189; curator decisions of 2026-10-29; `src/MATPredict/detect/family_status.py`, `--cross-family-policy {off,on}` and `--preset v0.8` = `--preset v0.7` plus the policy): a locus called by two families is resolved after searching (the scoped family wins; otherwise the family is decided once per genome by the higher mean best core-gene identity, with no margin, and each such call carries `family_ambiguous` with both identities, #219; an out-of-scope MTL call with fewer than two flank genes cannot win against a family whose own flank gene, SLA2 or DIC1, is in the call: it is set aside as `cross_call` with `flank_rule`, left out of the genome mean, and `family_ambiguous.reason` is `cross_call_flank_rule`), an MTL call outside Serinales is `family_unconfirmed` unless two of PAP1, OBP1, PIK1 lie within 15 kb, and the Ogataea MATa and MATalpha cassettes are indicated (`cassette_undetermined`). Default off; the flag is separate from `confidence`. Evidence and replays: `results/2026-10-29_yeast_routing/`.
- MATog scope widened to the genus *Ogataea* (461281) with `cassette_report`: 22 genomes change confidence in the replay and one MTL-only call (GCA_019232045.1) is no longer searched.
- Family `MATpk` for *Pichia kudriavzevii* (scope: the species, 4909; heterothallic) with the CBS 573 MATa record and the MATalpha record MH260578 (#189, #216): 42 -> 54 of the 56 *P. kudriavzevii* genomes called; other Pichia species unchanged (`results/2026-10-29_pichia_family/`). Label code `PicKud` (frozen).
- Labels: `registry_version` (`<registry_date>+<8-hex hash>`, stamped by `scripts/idiomorph_labels.py --stamp`) on every labelled output, frozen lineage and system codes with aliases, and a `label_status` column and explorer marker with the six values (including `family_unconfirmed`); the earlier mapping column is `label_mapping` (#210, #211, #212, #217, #218). System codes frozen by the curator: `PicKud`, `SchPom`, `SchCom` (families `mat3`, `Abeta`) and, on 2026-10-30, `SacSch`, `PetFab`, `BarCal`.
- Families `MATss` (*Saccharomycopsis schoenii*, scope: the species 44080; homothallic, NCBI genetic code 12), `MATpf` (*Petasospora fabianii*, species 36022; heterothallic) and `MATbz` (*Barnettozyma*, genus 599802), all scope-only, with five tier-3 candidate records (`db/candidates`, `needs_review`; our own gene models, because the *S. schoenii* authors' GTF is not public and no RNA-seq exists for the other two genera; #220, #189). Whether a *Saccharomycopsis* species other than *S. schoenii* is heterothallic (expected for *S. lassenensis*) is not encoded.
- Supporting profiles for MATa1, MATalpha2 and MATa2 (`results/2026-10-30_mat_hmms/`; 53 of 54 held-out sequences classified; a score difference, not referenced by any rule) that corrected the first *S. schoenii* gene names; evidence for the SLA2-type genera in `results/2026-10-29_near_tie/` and `results/2026-10-29_ss_records/`.
- Browser tool page: a "Try an example" list that fetches two *Neurospora crassa* MAT-locus regions from NCBI, with download links (#215).
- `scripts/dedupe_family_calls.py` opt-in overlap mode for family calls with different boundaries (#209).
### Changed
- The site pages count visits with [GoatCounter](https://www.goatcounter.com) (#225): anonymous page views only, no cookies. A footer note says so. The per-genome HTML reports written by `matpredict` carry no counter.
- Site wording: "Gating" becomes "Checks and confidence" (nav item "Confidence", table headers Check / What it tests, pipeline step "Verify"); the old `#gates` anchor still resolves. Wording only: no data or code change.
- Locus maps (explorer) mark the outcome of a resolved idiomorph pair: the gene that won the call keeps its solid fill and the gene set aside (`lost_to` in the report) is drawn pale with a solid outline and labelled "(set aside)"; the tooltip gives the basis and scores (HMM bitscore or alignment identity, for example "won over sexM (HMM score 364 vs 71)"), the locus line names the winner, and the legend has a new entry. The locus shards carry the verdict as a seventh item of a gene entry (`[W|L, other gene, basis, scores]`); a gene outside a resolved pair keeps the six-item form, so readers of the old shards are unaffected. Shards rebuilt for v0.7.0 data.
- Mucoromycotina tables: a blank order or family is filled from the genus (`scripts/mucoromycotina_genus_lineage.tsv`, NCBI taxonomy of the 51 genera, option `--genus-lineage`); a given order or family is never overwritten. The ZyGoLife LCG and Jena genomes had no order or family, and none of the Mucoromycotina inputs carried a family, so they sat under "(none)" in the explorer. Only the `order` and `family` columns of `genomes.tsv` and `loci.tsv` change (685 and 974 rows); calls, loci and `by_class.tsv` are identical. `build_site_data.py` now carries the family into the taxa rows and the genome tables (it was hard-coded blank), so the explorer groups Mucoromycotina by order, family and genus.


### Fixed
- A `high` call is no longer lowered to `medium` because of a stand-in margin (#137). When the idiomorph vote produced no ranking (a single candidate), the reported margin fell back to the identity gap between two different genes on one HMG box, and the tier cap for a thin margin then applied to it. That gap is reported but no longer lowers the tier. On the 381 `margin_fallback` genomes of the v0.7.0 run (228 *Macrophomina*, 98 *Pyricularia*, others), 380 loci return from `medium` to `high` and nothing else changes. In the full run the same change raises 645 loci in all (639 Ascomycota, 6 Basidiomycota), every one with an identity-gap margin below 5.
- Gene model selection is deterministic when one reference record aligns equally well at two places in a polish window (a gene that flanks a MAT locus on both sides). The tie went to whichever alignment the tool listed first, so the `sla2` model of a withheld MATsc locus in GCA_030583325.1 landed on the left copy in 2 of 10 repeat runs and on the right copy in 8, moving the `core_span` and `beyond_core_bp` of that locus. After score, identity and reference record id, the lower start and then the strand decide. Calls, confidence and the GFF3 did not change between the repeat runs.

### The definitive run
The v0.8.0 data come from commit `a2933f7` (the v0.8.0 release-prep branch: main at `6959657` plus the version bump; `src`, `db` and `scripts` are identical to main), run with `--preset v0.8` on database digest `70afa14e884604f06fa359b083a43317d83458f60b1b0c8254be561d7a131676` (184 records; candidate records excluded). The lists are those of the v0.7.0 run: 19,415 Ascomycota, 3,275 Basidiomycota and 974 Mucoromycotina genomes. Folders: `results/2026-10-31_v080_{ascomycota,basidiomycota,mucoromycotina}`.

| Campaign | Genomes listed | Called | Uncalled | No report (suppressed) | Loci | Loci high / medium / low |
|---|---|---|---|---|---|---|
| Ascomycota | 19,415 | 17,607 | 1,736 | 72 | 20,792 | 14,559 / 6,126 / 107 |
| Basidiomycota | 3,275 | 2,977 | 277 | 21 | 5,964 | 3,805 / 2,149 / 10 |
| Mucoromycotina | 974 | 868 | 105 | 1 | 914 | 689 / 212 / 13 |

Against v0.7.0, by locus: 475 calls lost (453 Ascomycota, 22 Basidiomycota) and 260 new (all Ascomycota), 645 confidence up (639 Ascomycota, 6 Basidiomycota) and none down, no idiomorph flips; Mucoromycotina is unchanged. 401 of the lost calls are the cross-family policy reporting one call where v0.7.0 reported two families on the same span (381 Ascomycota, 20 Basidiomycota), and the 260 new calls are the `MATog` (136), `MATpk` (88), `MATss` (15), `MATbz` (12) and `MATpf` (9) halves of such pairs; 69 are the *S. pombe* `mat2` copies of #227, one more is the `PR`-over-`MAT` swap of #227, and 3 are the *Ogataea* routing loss of #228. All 645 confidence rises are the #137 fix: each is a single-candidate locus whose stand-in identity-gap margin was below 5 and that was `medium` in v0.7.0 (380 of them carry the `margin_fallback` note, 265 had no note). 171 loci (67 Basidiomycota, 104 Ascomycota) have the same length, identity, span and call but a different exon placement, as the #224 tie-break intends. One low-confidence Basidiomycota partial `HD` locus (GCF_022376455.1) is called in v0.7.0 and not in v0.8.0; the same code calls it in one of three repeat runs, so it is the run-to-run variance of #161 that #224 removes. Nothing is unexplained. The determinism check (24 Ascomycota and 24 Basidiomycota genomes drawn with seed 20261102, run twice) gave no differences between the runs or against the campaign.

### Known issues
- Tuning of confidence on weak or short hits is held for v0.9.0: the single-candidate margin (#145), the vote on short fragments and tiny margins (#138), and the unpolished-gene noise floor (#149).
- The *Hyphopichia burtonii* *MTLA1* candidate record (#148) and the *Pichia* reference record (#189) are held for v0.8.1; until then 94 of the 158 *Pichia* genomes in the full run still call `MATsc` (54 call `MATpk`, 9 are uncalled and 1 has no report).
- In *S. pombe*, the silent `mat2`-P cassette can be reported as `mat1` (#227, held for v0.8.1). When two families are both in scope and tie on every evidence term, the cross-family rule falls back to the family name and `mat1` wins over `mat2`. In the v0.8.0 run this replaced 69 `mat2` calls in 54 genomes (taxid 4896) that v0.7.0 reported; in 54 of the 69 `mat1` says `M` where `mat2` said `P`. The cassette is still recorded as `cross_call` on the kept call.
- The same tie rule can swap a MAT call for a pheromone-receptor call: GCA_017654745.1 loses its `MAT` alpha call to an undetermined `PR` call on the same span (#227). With `--cross-family-policy off` both calls return.
- Genomes now routed by lineage to `MATog` (*Ogataea*) are no longer searched with `MATsc` and `MTL`. Where the `MATog` hits are withheld, a call v0.7.0 made is lost (#228, held for v0.8.1): *O. degrootiae* GCA_019232045.1 (a high-confidence `MTL` alpha call), GCA_030557025.1 and GCA_003707685.2 (three loci in all).
- The known issues of the v0.7.0 section below still apply unchanged: run-to-run variance on equal-scoring copies (#161), the rust genomes that are not searched (#165; the same 21 Basidiomycota genomes have no report), the b-locus intron cap (#153, #154), the defective-gene flag, and the rollout path without policy options.

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
