# Localize engines: what depends on the hit scores

Date: 2026-10-12. Branch: `verticillium-reference-records` at `39a65a2279`.
Type: design-implications review. Read-only on `src/`, `db/`, `results/`.

> Note (2026-10-13): the measured numbers in sections 2 and 4 (13 genome pairs,
> 85% retention) were taken while the benchmark was still running. Final numbers
> (15 genomes, 86.6% retention, regression-panel results) are in
> `2026-10-12_tblastn-fast-engine.md`. This report is model output from a
> subagent. Its code line references were not all re-verified.

## Question

The localize step (`search_localize`, `src/MATPredict/detect/search.py:430`) runs one
genome-wide `tblastn -seg no` with the BLAST default e-value of 10. No `-evalue`,
no `-num_threads` and no `-task` are passed (`search.py:455-458`). The plan is to make
this step pluggable: `tblastn` (default), `tblastn -task tblastn-fast`, and maybe
MMseqs2 translated search. This note lists every downstream use of the localize
scores, says which uses break with a second engine, compares ways to track scores
from more than one engine, and proposes a validation plan on final calls.

Every claim below cites a file. Anything not checked is marked **not verified**.

## Method

- Grep of `src/MATPredict/detect/*.py` for `identity`, `coverage`, `bitscore`,
  `evalue`, `align_length_aa`, `reference_length_aa`, and method names; each hit read
  in context.
- Read of analysis and results notes on distant lineages, the bar, holdouts and the
  flank bitscore floor.
- One small measurement on the benchmark files already on disk
  (`analysis/tblastn_fast/out/*.tsv.zst`, 13 genome pairs complete at the time of
  reading; the run is still in progress). Script:
  scratchpad `samehsp.py` (not committed). It joins HSPs by
  `(qseqid, sseqid, sstart, send)` and compares scores. It is preliminary.

## 1. Inventory: where localize scores are consumed

Kind: **A** = absolute threshold; **R** = relative comparison inside one run's own
hits; **X** = comparison across sources (methods, runs, or calibration sets).

`SearchHit` fields are defined at `search.py:94-133`. The docstring at
`search.py:111-117` says nothing ranks on bitscore/evalue in `_rank`; that is still
true for `_rank`, but other code now reads bitscore (rows 9-13).

| # | Site | Field | Kind | Use |
|---|---|---|---|---|
| 1 | `search.py:136`, `search.py:155-164` `drop_low_quality_hits` | genomic span `end-start+1` | A | Drop HSPs shorter than 90 bp before clustering. Not a score, but the HSP extent is engine-dependent. |
| 2 | `search.py:455-458` | e-value | A (implicit) | No e-value cutoff beyond BLAST's default 10. Another engine's default cutoff changes the hit set. |
| 3 | `clustering.py:17-32` `cluster_hits` | coordinates | A | Single-linkage by `max_gap`. No scores. HSP extent and count still matter. |
| 4 | `pipeline.py:675`, `:950` evidence floor `min_identity` | identity | A | Default `None` (off). CLI flag `--min-identity` (`cli.py:391`). Admission itself is gene count + core role (`pipeline.py:667-676`). |
| 5 | `pipeline.py:792-810` `_polish_rank` | distinct genes, max identity, hit count | R (but X across methods) | Order clusters for the per-family polish cap of 6 (`pipeline.py:707`). Max identity is taken over all the family's hits, any method. Hit count depends on how many HSPs an engine reports. |
| 6 | `pipeline.py:725`, `:877-907` strong-identity tier | max identity | A | Every cluster with best identity >= 50.0% is polished past the cap. Calibrated on 2,722 Dothideomycete genomes with tblastn (comment `pipeline.py:709-724`). |
| 7 | `pipeline.py:813-861` `strong_fragment_scores` | HSP span (translated) | A | Translates tblastn HSP spans (`method == "tblastn_genome"`, `:846`) and scores them with the family HMM against `gate_min_score` (default 100.0, `mat_gene_gate.py:68`; Mucoromycota manifest 97.8, `db/Mucoromycota/classifiers/MAT/manifest.yaml:1192`). Ranks clusters for the cap. HMM bits are engine-independent; the HSP extent is not. |
| 8 | `pipeline.py:2516-2533` classifier on `hsp_fragment` | HSP span (translated) | A (`min_margin`) | When no core gene was modelled, the idiomorph classifier scores tblastn HSP translations. Filter is the literal `"tblastn_genome"` method name. |
| 9 | `pipeline.py:976-993` `sub_floor_genes`, used at `:2657-2660` | bitscore | A | A gene whose every hit is below the family flank floor (default 39.0 bits) does not cap the tier at medium. |
| 10 | `flank_carried.py:76-127` `apply_flank_carried_rule`; floor `family_registry.py:86`, per-family override `family_registry.py:403` | bitscore (best raw, via `GeneEvidence.bitscore`) | A | A flank-carried call (no modelled core gene) is kept only if its strongest core hit is >= 39.0 bits and near the flank span; else withheld. |
| 11 | `idiomorph.py:102-140` `idiomorph_candidates`, used by `assign_idiomorph` (`idiomorph.py:168-190`) | bitscore, else identity | R (and X) | Idiomorph vote: the idiomorph with the highest best-hit bitscore wins. This is the call path for every family without a classifier. Only `db/Mucoromycota/classifiers/MAT` exists today. Falls back to identity when bitscore is None. |
| 12 | `split_locus.py:93-94` `_strength`; `split_locus.py:116-146` | bitscore, identity | R and X | Paralog guard (strongest hit for the core gene must be this locus) and flank choice. Flank raw identity is compared with modelled core identity: `best.identity < core.identity - 10` (`:141`). Core must be modelled at >= 95% (`split_locus.py:52`). |
| 13 | `pipeline.py:1772-1778` in `_gene_evidence` | evalue, bitscore | X | Per gene, the min e-value and max bitscore over ALL live hits of that gene in the cluster, any method. These go to the report and to rows 9-10. |
| 14 | `pipeline.py:1675-1682`, `:1786-1793` `_RAW_METHOD_PREFERENCE` | method, then identity | R within one method | Choose the raw hit shown as a gene's evidence. Method order first; identity breaks ties only inside one method. |
| 15 | `idiomorph.py:277-301` `_rank`; `resolve_idiomorph_overlaps` `idiomorph.py:584-650` | method, then identity | R within one method | Pre-polish sexM/sexP overlap resolution. Proteome method beats tblastn; identity decides inside one method. Exact tie left unresolved. |
| 16 | `idiomorph.py:662-720` `resolve_idiomorph_by_models` | polish model score | R within one tool | Post-polish resolution. "One tool's score is never compared with the other's" (docstring, `idiomorph.py:675-685`). Not a localize score. |
| 17 | `tiering.py:18`, `:25`, `:28-97` `allele_absent_genes_to_ignore`; identities built at `pipeline.py:2699-2704` | identity | A and X | Ignore a weak other-allele gene if its identity < 50.0 and >= 10.0 below the called allele's best modelled core. Identities come from `GeneEvidence`: polished identity for modelled genes, raw hit identity for unpolished genes. |
| 18 | `assembly_gap.py:69`, `:126-135` | identity | A | Flank hit >= 50.0% anchors an assembly gap. Report only (`pipeline.py:469`). |
| 19 | `receptor_arrays.py:51`, `:119-148`, `:167` | `align_length_aa / reference_length_aa`, else coverage | A | Receptor locus must cover >= 50% of one reference. Report-only flag (module docstring). Sums HSP lengths, so HSP fragmentation matters. |
| 20 | `mat_gene_gate.py:69`, `:81-88`, `:141-150` | identity of modelled flanks | A | Flank support >= 40.0% identity. Reads modelled evidence (polish or diamond), so not localize scores, except through which genes got modelled. |
| 21 | `pipeline.py:3062-3069` below-fraction-floor record | identity | report | `best_identity` in `withheld_detail` and diagnostics. |
| 22 | `report.py:396-414` JSON `gene_evidence`; `report.py:201-206` GFF3 attrs; `report/html.py:128`, `:345-354` | identity, coverage, evalue, bitscore | report/persisted | Written to `detection_report.yaml`. No engine name is written. `provenance.py:75-97` records `matpredict_version` and db digest, not BLAST version or engine. |
| 23 | `pipeline.py:1272-1305` evidence diagnostics JSONL | gene_count, hit_count, best_identity | persisted, X across runs | The replay dataset for cap decisions (comment `pipeline.py:807-809`; replays in `results/2026-10-05_dothideomycetes_full/cap_*_replay.py`). |
| 24 | `pipeline.py:1308-1345` idiomorph diagnostics | winner/loser identity, coverage, model scores | persisted, X across runs | Named as the calibration dataset for `DEFAULT_MIN_OVERLAP_FRACTION` (`idiomorph.py:12`, `:48`). Overlap fraction is computed from HSP coordinates. |
| 25 | `regression.py:131-139`, `:183-212` | report `identity`, classifier scores | X across runs | Regression diff prints core-model identity changes and classifier margin shifts >= 5.0 bits. Does not read bitscore/evalue. |

`classifier_build.py` does not read localize hits. It trains HMMs from curated
record proteins plus `training_extra.faa` (`classifier_build.py:8-17`). No change
there from an engine switch.

Raw localize hits are not persisted. The only persisted forms are the per-gene best
values in the report (row 22) and the per-cluster summaries in the diagnostics
(rows 23-24). An engine comparison cannot be replayed from old runs; it needs reruns.

## 2. What breaks or shifts with a second engine

### Measured: same alignment, same scale (tblastn vs tblastn-fast only)

On the 13 genome pairs on disk (preliminary, run in progress):

- 12,067 fast HSPs have a default HSP with the same `(query, contig, sstart, send)`.
- 11,379 of those rows are identical in pident, length, e-value and bitscore.
- 218 differ in bitscore; max |delta| = 4.6 bits; 98 lower, 120 higher in fast.
  Default bitscores of these rows span 27.7-386.0 bits.
- 94 of 4,002 matched rows with default bitscore 30-50 differ in bitscore.
- 680 differ in e-value; 470 of those have the same bitscore. The cause is
  **not verified**. A likely cause is tblastn sum statistics, where an HSP's e-value
  depends on the other HSPs it links with.

Reading: for tblastn-fast the bitscore scale is nearly the same. The main effect is
which HSPs exist, and their extent, not the scale. E-values are less stable than
bitscores. For MMseqs2, the score scale is **not verified**; it uses its own
scoring and filters, so assume a different scale until measured.

Benchmark caveats (from `analysis/tblastn_fast/run_bench.sh`): it uses
`-num_threads 4`, while `search_localize` passes no thread count; it does not pass
`-db_gencode`, which `search_localize` passes for non-standard codes (`search.py:459`).

### (a) Absolute cutoffs that need per-engine calibration

| Cutoff | Value | Where | Risk with a new engine |
|---|---|---|---|
| Flank-carried floor | 39.0 bits | `family_registry.py:86`; `flank_carried.py:120-122` | High. Tuned on 8 real loci (33.5-68.2 bits) and 52 noise loci (<= 30.4 bits). gzUmbVina2 passes at 39.7 bits, a 0.7-bit margin (`results/2026-09-27_flank_bitscore_floor/NOTE.md`). The measured bit deltas (up to 4.6) exceed that margin. |
| Sub-floor noise genes | same 39.0 bits | `pipeline.py:976`, `:2657` | Medium. Changes tier (medium vs high), not presence. |
| Strong-identity polish tier | 50.0% | `pipeline.py:725` | Medium. Fewer or shorter HSPs can drop a cluster's best identity under 50%, so the true locus falls back to the cap rank. The cap lost 17 true loci without this tier (comment `pipeline.py:719-721`). |
| Minimum HSP length | 90 bp | `search.py:136` | Medium for pheromone precursors and short HSPs; the floor was set on tblastn HSPs. |
| Implicit e-value cutoff | 10 | `search.py:455-458` | High for MMseqs2 (different default, **not verified** here); none for tblastn-fast if left at default. |
| HSP-fragment classifier gate | 97.8 / 100 HMM bits, `min_margin` | `pipeline.py:813-861`, `:2516-2533` | Medium. HMM bits do not depend on the engine, but the translated span does. |
| Allele-absent guard | < 50.0% and gap >= 10.0 | `tiering.py:18`, `:25` | Low-medium. Only raw (unpolished) identities come from localize. |
| Split-locus flank gap | 10 points | `split_locus.py:55`, `:141` | Low-medium. Flank raw identity vs core modelled identity. |
| Gap anchor / receptor coverage | 50% / 0.5 | `assembly_gap.py:69`; `receptor_arrays.py:51` | Low. Report-only. |

### (b) Relative uses that are scale-safe

These compare scores from one engine's output in one run. A monotone rescaling does
not change them. A change in which HSPs exist still can.

- `_RAW_METHOD_PREFERENCE` tie-break by identity inside one method (`pipeline.py:1786-1793`).
- `_rank` inside one method (`idiomorph.py:277-301`).
- `_polish_rank` ordering (`pipeline.py:792-810`), if all of a genome's hits come
  from one engine.
- `idiomorph_candidates` bitscore vote (`idiomorph.py:133-136`), if both idiomorphs'
  hits come from one engine and all carry a bitscore.
- `split_locus._strength` paralog guard (`split_locus.py:93`, `:126-130`), if all
  hits come from one engine.

Condition for all of these: **one engine per run**. If a run mixes engines (for
example tblastn-fast first, then default tblastn for a rescue), they become X uses.

### (c) Places that already mix scores from different sources

These show how the code handles more than one scoring source today.

1. **Named method preference, not score comparison.** `_RAW_METHOD_PREFERENCE`
   (`pipeline.py:1670-1682`) says: "This is a named method preference, NOT a score
   comparison: their identity numbers come from different searches and are not on a
   common scale." `_rank` (`idiomorph.py:277-301`) does the same for idiomorph
   overlaps. `resolve_idiomorph_by_models` never compares miniprot and exonerate
   scores (`idiomorph.py:675-698`).
2. **Mixing without a guard.** These sites mix sources today:
   - `_gene_evidence` takes min e-value and max bitscore over all methods in a
     cluster (`pipeline.py:1772-1778`). A diamond bitscore and a tblastn bitscore
     go into one max. The result feeds the 39-bit floor (row 10).
   - `idiomorph_candidates` uses bitscore when present and identity when not
     (`idiomorph.py:133`). A polished rescue hit is built with no bitscore
     (`pipeline.py:1640-1649`), so its identity (0-100 %) can be compared with
     another idiomorph's bitscore. Diamond and tblastn bitscores are also compared
     here.
   - `split_locus._strength` gives -1.0 to a hit with no bitscore
     (`split_locus.py:94`), so a polished rescue hit always loses to any raw hit.
     `split_locus.py:141` compares raw flank identity with modelled core identity.
   - `tiering.allele_absent_genes_to_ignore` compares polished and raw identities
     in one dict (`pipeline.py:2699-2704`).
   - Filters on the literal method `"tblastn_genome"` (`pipeline.py:846`, `:2524`)
     would silently exclude hits from a new engine if it gets a new method name.
     `_RAW_METHOD_PREFERENCE` would rank a new method name last
     (`pipeline.py:1679-1682`).

Design consequence: if a new engine gets a new `method` string, rows 7, 8, 14 and 15
change behaviour by name, not by score. If it keeps `"tblastn_genome"`, the engine
is invisible in every report.

## 3. Options for tracking scores from more than one engine

All options start with the same step: add `engine` (name + version + task/options) to
`SearchHit` and to the report run block. Without it no later analysis can separate
engines.

| Option | Cost | Must measure | Failure modes |
|---|---|---|---|
| **O1. Raw scores + engine tag; per-engine table of the absolute cutoffs in §2(a)** | Low. One field; a small cutoff table keyed by engine (default = current values). | Paired runs per engine on labelled sets (flank floor set, Dothideomycete cap set, regression panel). | Cutoffs multiply per engine. A cutoff not in the table silently uses tblastn values. |
| **O2. Coverage-weighted identity** (`identity * coverage/100`) | Very low. Fields exist. | Distribution per engine on true vs noise loci. | Still engine-dependent through HSP extent. `qcovhsp` is per HSP (`search.py:356-361`), so multi-exon genes score low. Identity is a weak signal: minus-strain identities run 25.9-43.5% (`search.py:148-151`). |
| **O3. Bits per aligned residue** (`bitscore / align_length_aa`) | Very low. | Same as O2. | Matrix- and engine-specific. Noisy for short HSPs, which are the pheromone-precursor case. |
| **O4. Blast score ratio** (`bitscore / self-bitscore` of the reference, per engine) | Low-medium. One self-search of `reference.faa` per engine per db build; store a table. Named as the candidate in `search.py:111-117`. | Self-scores per engine; replay of the idiomorph vote and the 39-bit floor against Zygo 23 and the flank-floor labels. | Still depends on HSP extent. A split gene gives low BSR per HSP. The search.py comment notes idiomorph margins as small as 0.59 points must be rechecked before BSR replaces any ranking. |
| **O5. Rank within query** | Very low. | Nothing new. | Scale-free but loses absolute strength. Cannot express a noise floor, which rows 9-10 need. |
| **O6. Per-engine empirical e-value from decoys** (shuffled or reversed genome) | High. An extra search per genome or per calibration genome. | Decoy score distribution per engine and genome size. | Decoys model random noise, not the real noise. The known noise here is paralogs and artefacts: Debaryomyces-type Serinales calls, the Mucor P1 HMG paralog (`classifier.py:36-47`). The code already moved off e-values because they grow with genome size (`family_registry.py:89-93`). |

### Recommendation

Use **O1** now, and store **O4 (BSR)** as an extra derived field without using it in
any decision.

Reasons:

1. The absolute cutoffs that read localize scores are few (§2(a): about nine, two of
   them report-only). Each one was tuned on a labelled set that is on disk. A
   per-engine replay of those sets is cheaper than a new normalised scale.
2. For tblastn-fast the measured bitscore scale is nearly the same for the same
   alignment (§2). The bigger risk is lost or shorter HSPs, which no normalisation
   repairs. Only final-call comparison can measure that.
3. The code's existing rule for mixed sources is "named preference, never cross-scale
   comparison" (§2(c)1). O1 keeps that rule. A normalised score would invite
   cross-engine comparisons the curated rules were never tested on.
4. BSR is the project's own named candidate (`search.py:111-117`). Recording it costs
   little and gives data for a later decision. Using it to rank needs the Zygo 23
   check first, as that comment says.

Decision that I cannot support with data: whether MMseqs2 is worth adding at all.
There is no MMseqs2 measurement in the repo.

Also fix, independent of the engine choice, so mixing is visible: give
`_gene_evidence` best-score fields per method, or record which method gave the best
bitscore (`pipeline.py:1772-1778`). This is a suggestion only; it is not decided.

## 4. "Maybe the distant hits are ambiguous anyway?"

Short answer: partly true for presence calls in distant lineages; false for some
real calls that already rest on weak hits.

### Evidence that weak/distant localize hits rarely become calls on their own

- Admission needs >= 2 distinct genes with a core gene (`pipeline.py:640-676`), and a
  reported locus needs >= 2 modelled genes (`MIN_POLISHED_GENES`,
  `pipeline.py:679-686`). The cost note there: on 2026-09-22 Pezizomycotina panels the
  bar cut 46,647 loci to 5,648 and kept all 2,957 high-confidence calls.
- Gap-lineage pilots and holdout rerun (`docs/notes/2026-09-24_gap-lineage-pilots-and-holdout-rerun.md:79-98`;
  `results/2026-09-24_pilots/summary.txt`): "Every record not found was withheld by the
  bar at the right place." "The withheld evidence is often thin -- a single unmodelled
  HSP." "In these lineages the bar, not the search, is what stands between a genome
  and a call." Example: K. lactis MATa, one 159 bp MATA1 HSP, 0 modelled.
- Synteny test (`docs/notes/2026-09-24_bar-synteny-test.md:45-59`): in Pezizales (69/81),
  Dipodascales (35/46), Saccharomycodales (15/18) the withheld locus is next to SLA2,
  so localize found the right place and the bar withheld it. In Orbiliales (5/40) and
  Serinales (3/73) the locus was not localised at all (reference gap).
- Once a cluster is admitted, polish rescues missing core genes in the cluster window
  (`pipeline.py:2343-2397`). The window comes from cluster coordinates and reference
  protein length, not from scores (`pipeline.py:1578-1606`). So a weak hit only has to
  seed a cluster; polish rebuilds the gene. How often polish recovers a gene whose
  localize hit was lost: **no data**.
- With a proteome, localize runs only as a rescue for missing core genes
  (`pipeline.py:2115-2131`). Engine choice matters less for annotated genomes.
  Size of that effect: **no data**.

These notes are from code of 2026-09-24. The bar is unchanged in current code
(`pipeline.py:686`). I found no later re-measurement and no bar-only-tier decision in
`analysis/decisions.md` or `analysis/open-questions.md`.

### Evidence that the weak tail does change final calls

- Flank-carried calls rest on the strongest raw core hit's bitscore. Real loci score
  33.5-68.2 bits; the floor is 39 bits; one real locus passes by 0.7 bits
  (`results/2026-09-27_flank_bitscore_floor/NOTE.md`). These are Umbelopsis-type HMG
  fragments of 33-44 bits (`family_registry.py:89-93`). A small shift in that band
  changes kept vs withheld.
- HSP-fragment classifier calls: when no core gene is modelled, the idiomorph label
  comes from tblastn HSP translations (`pipeline.py:2516-2533`; `classifier.py:13-19`,
  Umbelopsis sp. M5902 and Mucor hiemalis gzMucHiem1).
- Idiomorph vote for families without a classifier compares best bitscores between
  idiomorphs (`idiomorph.py:102-140`). The Actinomucor example decided on a 4.6-bit
  gap (`idiomorph.py:151-153`), the same size as the largest bit delta measured in §2.
- Polish cap rank uses hit counts and max identity (`pipeline.py:792-810`). The
  Dothideomycetes run showed true loci ranked 7th-12th and never polished in 15
  genomes before the 50% tier (`pipeline.py:709-724`). Fewer HSPs per cluster change
  the rank.
- Withheld loci are written to `suppressed_loci` and diagnostics. If the curator
  later adds a bar-only tier, the weak-hit population becomes call-relevant.

Conclusion: there is no data that the weak tail is irrelevant to final calls. There
is data for the opposite in four mechanisms above. The 85% retention of strong best
hits from the partial benchmark does not answer this; only a final-call comparison
does (§5).

## 5. Calibration and validation plan (final calls, existing assets)

Principle: compare the two engines at the **same commit and same db**. Old cohort runs
used older code, so a diff against them mixes code changes with engine changes.

### Step 0. Prerequisite code change (not done; for the implementer)

Engine tag on `SearchHit` and in `provenance.run_block`; one engine per run; decide
the method string (§2(c)).

### Step 1. Hit level, truth assemblies (in progress)

- Assets: `results/truth_assemblies.json` (18 records, 16 assemblies, verified),
  `testset/coordinate_truth/ascomycota.yaml`, `analysis/tblastn_fast/`.
- Add to `analyze.py`: per truth gene, best bitscore in each engine; count of HSPs
  lost in the 30-50 bit band; bit delta on matched HSPs (as in §2).
- Rerun with `-db_gencode` as detect does, and with the thread count detect uses, if
  timing is to be quoted.

### Step 2. Final calls, fixed regression panel

- Assets: `testset/regression_panel.tsv` (30 ascomycota, 33 basidiomycota, 80
  mucoromycota, 23 record_sources rows, verified by count) and Zygo 23
  (`scripts/run_zygo_regression.slurm`, 46 reports in prior runs).
- Run `scripts/run_regression_panel.sh` twice at one commit: base = default tblastn,
  candidate = tblastn-fast. Diff with
  `scripts/regression_check.py diff --pair ...`. It lists `call_lost`, `call_gained`,
  `idiomorph_changed`, `confidence_changed`, `gene_set_changed`, `core_model_changed`,
  `classifier_input_changed`, `classifier_shift` (`regression.py:183-212`).
- Prior bases exist, verified by listing, each with 30/33/80/23/46 reports:
  `results/2026-10-11_regression_vd11/base/{ascomycota,basidiomycota,mucoromycota,record_sources,zygo23}`
  (diff summary names code `e9060daf32`). Also
  `results/2026-10-12_regression_vocab/base_ascomycota` and `cand/ascomycota` (30).
  Use them only if the candidate runs at the same code and db; otherwise run a fresh
  base.

### Step 3. Weak-tail sets

- Flank-carried floor set: `results/2026-09-27_flank_bitscore_floor/` (`loci.tsv`,
  `genomes.tsv`, `hit_one.sh`, `hits/` with 318 entries). Rerun `hit_one.sh` with
  `-task tblastn-fast`. Check the 8 real loci still score >= 39 bits and the 52 noise
  loci stay below.
- Holdout radii: `scripts/run_holdout_benchmark.py`, prior results
  `results/2026-09-23_holdout/` and `results/2026-09-24_holdout/`
  (`record/species/genus/family/order.yaml`). Distant radii are where weak hits
  matter. Score with no_reference cases removed, as ruled on 2026-09-26.

### Step 4. Large cohort (only if Steps 2-3 show no or explained changes)

Candidate frames, verified by listing:

| Directory | Content | Code | Use |
|---|---|---|---|
| `results/2026-10-05_dothideomycetes_full/` | `reports_all.tar.zst`, 2,731 listed, 2,728 reports (`compare_vs_v060_summary.md`) | frozen `run-84361bc` (`analysis/2026-10-05_dothideomycetes-full-run.md:10`) | Best fit: the 50% tier and cap replays were measured on it (`cap_*_replay.py`). |
| `results/2026-10-03_ascomycota_v060/` | `reports_all.tar.zst`, `genomes.tsv` (19,415 rows), `loci.tsv` | v0.6.0, `run-7c7ed99` (README) | Sampling frame by order; too old to be a direct baseline. |
| `results/2026-10-03_basidiomycota_v060/` | `genomes.tsv` (3,270 rows), `loci.tsv`, `by_order.tsv` | v0.6.0, `run-7c7ed99` (README) | Sampling frame for Basidiomycota. |

Plan: draw a stratified sample by order from one frame, run both engines at the
current commit, diff with `regression_check.py`. Sample size and acceptance limits
are curator decisions (below). Runtime for the cohort: **not measured**; production
runs single-threaded tblastn, the benchmark used 4 threads.

### What to report

Per engine pair: counts of each `change_types` value; for every `call_lost` and
`idiomorph_changed`, the locus, its best raw bitscore under each engine, and which
mechanism in §4 moved it; runtime per genome on the same thread count.

## Open questions for curator

1. What change in final calls is acceptable for a faster engine: zero lost calls on
   the regression panel + Zygo 23, or a stated tolerance?
2. Should an engine change keep the method string `tblastn_genome`, or get its own?
   (§2(c): the name drives rows 7, 8, 14, 15.)
3. Should the 39-bit flank floor be re-tuned per engine, or should a faster engine be
   barred from runs where flank-carried calls matter (Mucoromycota, Serinales)?
4. Is a mixed mode wanted (fast engine first, default tblastn for families or
   genomes with no admitted cluster)? This turns scale-safe uses into mixed ones.
5. Should withheld (bar-only) loci count in the engine comparison, given the open
   bar-only-tier question from 2026-09-24?
6. Is MMseqs2 in scope, given no measurement exists in the repo?
7. Large-cohort step: which frame (Dothideomycetes full, Ascomycota v0.6.0,
   Basidiomycota v0.6.0) and what sample size?
