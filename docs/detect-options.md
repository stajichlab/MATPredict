# `matpredict detect`: options, policies, presets and environment variables

This page lists every setting that changes what `matpredict detect` does, so a run can be reproduced and a difference between two runs can
be traced. The tables under "Options" are generated from the command's own `--help` (`matpredict detect --help` is the authority; a test
checks that this page names every policy, preset and environment variable the code defines). Written for v0.8.1.

## The short version

- **Code defaults reproduce the oldest behaviour.** Every policy is `off` unless asked for, so a bare `matpredict detect` run is not the run
  behind the published tables.
- **A preset turns on the policies of one release.** `--preset v0.8.1` is what the v0.8.1 definitive run used. An option given on the command
  line beats the preset, even when it names the default value. The report records the preset (and each policy) it ran with.
- **Policies change calls; options tune thresholds; environment variables point at files or change speed.** Only the policies and the
  thresholds can change a call.

## Presets

A preset fills the policy options that were not passed explicitly.

| Preset | partial_core | homothallic | confidence_cap | strong_single_gene | scoped_roster | cross_family | fallback_after_withheld |
|---|---|---|---|---|---|---|---|
| (none) | off | off | off | off | off | off | off |
| `v0.7` | coverage_profile | both_genes | off | off | on | off | off |
| `v0.8` | coverage_profile | both_genes | off | off | on | on | off |
| `v0.8.1` | coverage_profile | both_genes | off | off | on | on | on |
| `v0.9` | coverage_profile | both_genes | off | off | on | keep_labelled | on |

Each preset adds to the one before it: `v0.8` adds the cross-family policy to `v0.7`, `v0.8.1` adds the fallback after withheld hits to
`v0.8`, and `v0.9` adds the three **report-only** rules of the v0.9.0 cycle to `v0.8.1` (curator's ruling 2026-11-07): `cross_family_policy=keep_labelled`
(both calls kept where the kept one has no idiomorph label and the other has one, #255), `idiomorph_margin_policy=separate` (a single-candidate
identity gap reported as `idiomorph_identity_gap`, #145) and `out_of_scope_label_policy=on` (a Basidiomycota MAT call outside its scope on STE3, PAN6
and RPL39 alone marked `family_unconfirmed`, #170; the command-line flag is `--scope-label-policy`). `v0.9` turns on no rule that changes a tier or an idiomorph label: the weak-fragment,
weak-flank, short-fragment-vote, single-region and weak-coverage rules stay off in every preset until the curator approves their genome lists. Two policies stay off in every preset: `confidence_cap_policy` and `strong_single_gene_policy` (the latter waits for a 12-locus
sequence check, #261).

A full run passes the preset through the batch scripts, for example `DETECT_ARGS="--preset v0.8.1"` with `scripts/run_clade_panel.slurm`.
Mucoromycotina runs use `--phylum Mucoromycota` instead of a taxid.

## Options

### Policies

| Option | Values | Default | What it does |
|---|---|---|---|
| `--partial-core-policy` | off,coverage_profile | `off` | off: the best-bitscore idiomorph vote stands (default). coverage_profile: a core hit that covers < 80% of its reference does not beat a full-length (>= 90%) core gene of the other idiomorph when the alpha-box/HMG-box profiles agree; where they disagree or are unavailable the call is `undetermined` (see `idiomorph_vote_note`). A partial hit near a contig end is treated as truncation and the vote is kept. Independent of --idiomorph-policy. |
| `--homothallic-policy` | off,both_genes | `off` | off: a locus with an alpha-box and an HMG-box gene is reported as one idiomorph unless the reference-coverage rule of `classify_locus` passes (default). both_genes: a locus whose modelled core genes include a COMPLETE alpha-box and a COMPLETE HMG-box protein, each confirmed by the gene-identity profiles (not an HMG paralog or MAT1-1-3), within the family's separation limit, is classed `homothallic_candidate` (see `homothallic_note`). A candidate, never an assertion; the idiomorph is reported `undetermined` (the vote's winner is kept in the note as `vote_idiomorph`) and the confidence is unchanged. |
| `--confidence-cap-policy` | off,partial_weak,weak_only | `off` | off: no extra cap (default). partial_weak: a `high` call drops to `medium` when its best core hit is a contig-edge partial (< 80% of its reference, within 500 bp of a contig end, no full-length core gene) or a weak hit (< 45% identity) with no flank gene. weak_only: the second criterion alone. Changes confidence only; the reason is written as `confidence_cap`. |
| `--cross-family-policy` | off,on,keep_labelled | `off` | off: every family's calls are reported as they are (default). on: two families' calls at the same locus (>= 90% of the shorter) are one locus, kept by the family whose taxonomic scope contains the genome (else the better score), the other recorded as `cross_call`; an MTL call in a genome outside Serinales is flagged `family_unconfirmed`, or `family_confirmed` when two of PAP1, OBP1, PIK1 lie within 15 kb of its core genes. Confidence is never changed. keep_labelled: as on, and where the kept call has no idiomorph label while another family's call at the same locus has one, both calls are kept (#255). |
| `--fallback-after-withheld-policy` | off,on | `off` | off: a genome routed by taxid or lineage searches only the families scoped to it (default). on: when that search reports no locus because the modelled-gene bar withheld every hit, the phylum-fallback families are searched as well (an Ogataea genome whose MATog hits are withheld is then searched with MATsc and MTL); the report's `fallback_search` says so. Needs --taxid. |
| `--weak-fragment-cap-policy` | off,on | `off` | off: no cap (default). on: an Ascomycota:MAT call whose strongest core gene (MAT1-*) is not a model on which both polish tools agree and scores under 50 bits is capped at tier `low`; the idiomorph, the locus and the class are never changed (#275). The report's `weak_fragment_cap` says which gene decided. |
| `--weak-coverage-noise-policy` | off,on | `off` | off: only a gene scoring under the family's noise floor (39 bits) is ignored when an unpolished gene would cap the tier (default). on: a gene whose every hit covers < 30% of its reference is ignored too, whatever its bitscore, so a fragment crossing 39 bits no longer lowers a call to `medium` (#149). Tier only (medium to high). |
| `--scope-label-policy` | off,on | `off` | off: no label (default). on: a Basidiomycota:MAT call in a genome outside the family's taxonomic scope whose genes are only STE3, PAN6 and RPL39 is marked `family_unconfirmed` (reason `mat_outside_scope_weak_genes_only`) in `family_status`; the call, idiomorph and tier are never changed (#170). An unknown lineage is never flagged. |
| `--short-fragment-vote-policy` | off,undetermined | `off` | off: the best-bitscore idiomorph vote stands (default). undetermined: a vote won by less than the family's minimum margin (5 bits) while the best core gene covers < 60% of its reference is reported `undetermined`; the vote's winner, margin and coverage go into `idiomorph_vote_note` (#138). The locus, coordinates and tier are not changed. |
| `--weak-flank-cap-policy` | off,on | `off` | off: no cap (default). on: a `high` Ascomycota:MAT call whose every flank gene is weak (unpolished, and under 50% identity or 30% coverage) is lowered to `medium`; the idiomorph, the locus and the class are never changed (#129). The report's `weak_flank_cap` names the weak genes. |
| `--withheld-detail-policy` | off,on | `off` | off: a withheld locus lists its genes by name only (default). on: each entry of `suppressed_loci` also lists `genes` (name, role, polish status, identity, coverage, bitscore), which shows what a locus under the polished-gene bar consisted of. Report only; no call changes. |
| `--single-region-cap-policy` | off,on | `off` | off: no cap (default). on: a `high` Basidiomycota call whose live core genes fall in fewer than two distinct regions (models overlapping by half or more are one region) and that has no polished flank gene is lowered to `medium`; the idiomorph, the locus and the class are never changed (#167). The report's `single_region_cap` says what was found. |
| `--idiomorph-margin-policy` | off,separate | `off` | off: `idiomorph_margin` is reported as it is (default). separate: for a locus with a single idiomorph candidate the identity-point gap the margin falls back to is reported as `idiomorph_identity_gap` and `idiomorph_margin` is null (#145). Report only: no call, coordinate or tier changes. |
| `--strong-single-gene-policy` | off,on | `off` | off: a lone core gene with no flank and no second modelled gene is withheld (default). on: in a family with no other report, one modelled core gene at >= 85% identity covering >= 90% of its reference is reported as `idiomorph_gene_only`, confidence `low`, flagged `strong_single_gene`; two such loci of opposite idiomorph are both reported and the genome gets a `two_idiomorphs` statement. Independent of --idiomorph-policy and --partial-core-policy. |
| `--scoped-roster-policy` | off,on | `off` | off: a roster gene that declares its own `taxonomic_scope` in order.yml is not searched (default; the run is the one it was before such genes were added). on: it is searched, and counts, for a genome whose lineage reaches its scope, and nowhere else (prototype: GAP1/BGL1 as optional bonus flanks of Ascomycota MAT in Pleosporaceae). |

### Routing and taxonomy

| Option | Values | Default | What it does |
|---|---|---|---|
| `--exclude-records` |  | `` | comma-separated curated record_ids to withhold from this run, for leave-one-out recall. They leave the query set AND every db-derived quantity, so the genome a withheld record came from becomes a genuine held-out test rather than a self-consistency check. See MATPredict.detect.holdout for choosing a set by taxonomic radius. |
| `--taxid` |  | `None` |  |
| `--phylum` | Ascomycota,Basidiomycota,Mucoromycota,auto,auto-weak | `None` | Restrict detection to one phylum's curated families outright, skipping taxid-based routing entirely. Use when the genome's phylum is known but its taxid routes badly -- without it, a taxid no family's taxonomic_scope covers falls back to searching every family in every phylum. The choices are read from the database root at startup, so they always match what is actually curated. `auto` asks GmG (https://github.com/stajichlab/GmG) which phylum the genome is in and searches that phylum; a weak or non-fungal call, or a phylum with no curated family, is NOT searched and the report says why. `auto-weak` also accepts a weak call. An explicit --taxid wins over both. |
| `--gmg-index` |  | `None` | With --phylum auto/auto-weak: the GmG reference index. A path, a release file name (downloaded once to ~/.cache/matpredict/gmg and checked), `latest`, or a URL. Default: $MATPREDICT_GMG_INDEX, else the index MATPredict is validated with; a warning says when a newer one exists. |
| `--genetic-code` |  | `None` | NCBI translation table for this genome (e.g. 12 for the CUG-Ser1 clade). Omit to derive it from --taxid, which reads the same cached taxonomy document routing already fetches; falls back to 1 when that lookup cannot answer. An explicit value always wins. |
| `--genetic-code-map` |  | `None` | TSV/CSV file mapping a genome id to its translation table, for batch runs where each genome differs. Two columns, id and code, with an optional header; the id is matched against the genome FASTA's basename with extensions stripped. Used only when --genetic-code is not given. The BFD sample sheet's ASMID/TRANSL_TABLE columns are exactly this shape. |
| `--exhaustive` |  | `False` | Search every family in every phylum when no curated family covers the genome (no taxid, an unresolvable lineage, or an uncurated phylum). Default: such a genome is NOT searched and its report says so (curator's ruling 2026-09-26; the chytrid control ran ~31 min per genome under exhaustive and called nothing). Never widens a route that matched. |

### Search and evidence

| Option | Values | Default | What it does |
|---|---|---|---|
| `--evidence-diagnostics` |  | `None` |  |
| `--min-hits` |  | `2` | Minimum DISTINCT genes a family needs in a cluster to reach Stage 2 polishing (default: 2) |
| `--min-identity` |  | `None` | Minimum best-hit percent identity for a family to reach Stage 2 polishing (default: no identity cutoff) |
| `--require-core-role, --no-require-core-role` |  | `True` | Require at least one core_MAT hit for a family to reach Stage 2 polishing (default: True) |
| `--supported-min-bitscore` |  | `33.0` | Bitscore floor for the report-only `supported_span`: a locus's own hits at or above it, plus its modelled genes, define the supported extent; weaker hits do not. Does not change calls or the cluster span. Default: 33. |

### Polishing

| Option | Values | Default | What it does |
|---|---|---|---|
| `--annotation` |  | `None` | GenBank annotation of the genome, shown as a gene-model track in the report. A GenBank file given as --genome is used the same way, and its CDS proteins also serve as the gene calls. |
| `--helixerlite` |  | `None` | helixer-lite executable (default: $MATPREDICT_HELIXERLITE). When set, it runs on the locus windows only and its models are a report track; calls are unchanged. |
| `--helixer-lineage` |  | `fungi` | helixer-lite lineage (default: fungi) |
| `--helixer-cpus` |  | `4` | CPUs for helixer-lite (default: 4) |
| `--max-polished-clusters-per-family` |  | `6` | Polish at most N admitted clusters per family, ranked before polishing by distinct genes, then best identity, then hit count. Default: 6 (curator's ruling 2026-09-26; measured -45% compute, 4 of 613 calls lost). 0 = no cap. |
| `--polish-strong-identity` |  | `50.0` | Also polish every admitted cluster whose best identity is at or above this percentage, even past the cap; the rest of the cap is filled from the usual rank. Stops a true locus that hits few genes at high identity being outranked by noise that hits more genes at low identity. Default: 50. 0 = the plain cap. |
| `--polish-core-pair, --no-polish-core-pair` |  | `True` | Also polish every admitted cluster with two distinct live core genes of one idiomorph (e.g. MTLA1 + MTLA2) whatever its rank, on top of the cap (it takes no slot). The property is the cluster's own, so adding families that admit more clusters cannot push such a cluster out of the cap. Default: on. --no-polish-core-pair restores the earlier behaviour. |

### Inputs and outputs

| Option | Values | Default | What it does |
|---|---|---|---|
| `--genome` |  | `None` |  |
| `--proteins` |  | `None` |  |
| `--html, --no-html` |  | `None` | Write report.html in --out-dir (default: on; $MATPREDICT_HTML=0 turns the default off) |
| `--pdf` |  | `False` | Also write report.pdf (needs WeasyPrint or Chrome/Chromium; off by default) |
| `--sample` |  | `None` | Sample name for the report (default: the genome file name up to the first dot) |
| `--organism` |  | `None` | Organism name for the report, e.g. 'Fusarium oxysporum f. sp. lactucae'. Shown as given; routing still uses --taxid |
| `--out-dir` |  | `None` |  |
| `--preset` | v0.7,v0.8,v0.8.1,v0.9 | `None` | Set the five policy options below to a named release setting: v0.7: partial-core-policy=coverage_profile, homothallic-policy=both_genes, confidence-cap-policy=off, strong-single-gene-policy=off, scoped-roster-policy=on; v0.8: partial-core-policy=coverage_profile, homothallic-policy=both_genes, confidence-cap-policy=off, strong-single-gene-policy=off, scoped-roster-policy=on, cross-family-policy=on; v0.8.1: partial-core-policy=coverage_profile, homothallic-policy=both_genes, confidence-cap-policy=off, strong-single-gene-policy=off, scoped-roster-policy=on, cross-family-policy=on, fallback-after-withheld-policy=on; v0.9: partial-core-policy=coverage_profile, homothallic-policy=both_genes, confidence-cap-policy=off, strong-single-gene-policy=off, scoped-roster-policy=on, cross-family-policy=keep_labelled, fallback-after-withheld-policy=on, idiomorph-margin-policy=separate, out-of-scope-label-policy=on. A policy flag given explicitly overrides the preset. Without --preset every policy keeps its own default (off). |
| `--emit-cds-fasta` |  | `False` | Also emit real CDS features with translation= attributes in the GFF3, and write a companion detected_loci.fasta next to it. The companion FASTA contains the FULL sequence of every contig the detection results reference, not just the locus spans, so it can be large (a whole chromosome-scale contig per referenced contig). Off by default. |

## Environment variables

| Variable | Effect |
|---|---|
| `MATPREDICT_DB_ROOT` | the curated database folder (default `db/` of the checkout). A run must use the database of the same tree as its code. |
| `MATPREDICT_CACHE_DIR` | cache for taxonomy lookups and the GmG index (default `.matpredict_cache`). |
| `MATPREDICT_CONFIG` | path of the config file (else `$XDG_CONFIG_HOME/matpredict/config.toml`). |
| `MATPREDICT_NCBI_EMAIL`, `MATPREDICT_NCBI_API_KEY` | NCBI E-utilities identity for taxonomy lookups (else the `[ncbi]` table of the config file). |
| `MATPREDICT_TAXONOMY` | a local taxonomy table; used first, with no network. |
| `MATPREDICT_OFFLINE` | `1`: a taxid the local table does not hold is an error instead of an NCBI lookup (then give `--phylum` and `--genetic-code`). |
| `MATPREDICT_GMG_INDEX`, `MATPREDICT_GMG_CACHE` | the GmG index used by `--phylum auto`, and where it is kept. |
| `MATPREDICT_HELIXERLITE` | the helixer-lite executable (default for `--helixerlite`). |
| `MATPREDICT_MAX_INTRON_BP` | overrides the intron cap of the polish calls for a run (a family may set its own cap in `order.yml`). |
| `MATPREDICT_POLISH_WORKERS` | number of parallel polish calls per genome (default the smaller of 4 and the CPU count; `1` runs them inline in the serial order). Does not change calls. |
| `MATPREDICT_TBLASTN_SEG` | `yes` runs the genome-wide `tblastn` as two calls: `-seg no` for most queries and `-seg yes -soft_masking true` for the selected ones (default `no`: one `-seg no` call, as in v0.8.1). In `main`, not in a release yet (#276). |
| `MATPREDICT_TBLASTN_SEG_GENES` | comma-separated reference gene names that get the filtered search (recommended `SLA2,APN2,MAT1-1-1`); replaces the length rule. |
| `MATPREDICT_TBLASTN_SEG_MIN_AA` | with `MATPREDICT_TBLASTN_SEG=yes` and no gene list: queries of at least this many residues get the filtered search (default 250). |
| `MATPREDICT_HTML` | `0`, `false`, `no` or `off` turns `report.html` off unless `--html` is given. |
| `MATPREDICT_BFD_SUPPRESS` | overrides the path of the shared BFD suppress list (`detect suppress-filter`). |

## In `main`, not in a release yet

These are merged but unreleased; none changes a default, and none is in a preset (the preset change belongs in the release-prep pull request).

| Item | What it is |
|---|---|
| `--weak-fragment-cap-policy off\|on` (#275, PR #288) | an `Ascomycota:MAT` call at medium or high whose strongest core gene is not a both-tools-agree model and scores under 50 bits is capped at tier `low`; the idiomorph label is kept. On the v0.8.1 run: 108 loci match, 68 would change tier (23 medium, 45 high). `docs/notes/2026-11-06_weak-fragment-cap.md`. |
| `--weak-coverage-noise-policy off\|on` (#149) | an unpolished gene whose every hit covers under 30% of its reference is noise whatever its bitscore, so it no longer caps the call at `medium`. Tier only. |
| `--scope-label-policy off\|on` (#170) | a `Basidiomycota:MAT` call outside the family's scope that rests on STE3, PAN6 and RPL39 alone is marked `family_unconfirmed` in `family_status`. Report only. |
| `--short-fragment-vote-policy off\|undetermined` (#138) | a close idiomorph vote (< 5 bits) whose best core gene covers < 60% of its reference is reported `undetermined`, the vote kept in `idiomorph_vote_note`. |
| `--weak-flank-cap-policy off\|on` (#129) | a `high` `Ascomycota:MAT` call whose every flank gene is weak (unpolished, under 50% identity or 30% coverage) is lowered to `medium`. Changes the tier only. |
| `--withheld-detail-policy off\|on` | each withheld locus also lists its genes with role, polish status, identity, coverage and bitscore. Report only. |
| `--single-region-cap-policy off\|on` (#167) | a `high` Basidiomycota call on fewer than two distinct core regions and no polished flank gene is lowered to `medium`. Changes the tier only. |
| `--idiomorph-margin-policy off\|separate` (#145) | a locus with a single idiomorph candidate reports the identity-point gap its margin falls back to as `idiomorph_identity_gap` and `idiomorph_margin` as null. Report only: no call, coordinate or tier changes. |
| `MATPREDICT_TBLASTN_SEG`, `MATPREDICT_TBLASTN_SEG_GENES`, `MATPREDICT_TBLASTN_SEG_MIN_AA` (#276, PR #284) | opt-in `tblastn` filtering for the slowest genomes. With the gene list `SLA2,APN2,MAT1-1-1` it takes two outlier genomes from 3,063 s and 1,685 s to 468 s and 288 s with identical loci; it gives no speed-up on typical genomes. `docs/notes/2026-11-06_tblastn-seg-policy.md`. |

## Tried and not built

| Item | Outcome |
|---|---|
| `--idiomorph-rank-basis identity\|bitscore` (draft PR #285, closed unmerged) | compared bitscores instead of identity in the overlap resolver. It changed 13 of 300 Ascomycota calls, dropped silent cassette loci in Saccharomyces-type yeasts and moved the tier at 73 of 347 loci in both directions. The code stays on the `resolver-robust` branch. `docs/notes/2026-11-06_idiomorph-rank-basis.md` (on that branch). |
| `--polished-vote-policy` (#275; branch `polished-vote`, no pull request) | a polished idiomorph beating an unpolished one in the vote. It flipped a correct call (XylFL0594, MAT1-2 confirmed by reads) and changed none of the 47 pilot assemblies, including the *T. viride* case it was meant for. `docs/notes/2026-11-07_polished-vote-not-built.md`. |

## Adding a policy

1. Add the option with `action=_PolicyFlag`, `choices` and a default of `off` in `src/MATPredict/detect/cli.py`, and pass it into
   `run_pipeline` through `pipeline_kwargs` (the `parameters` block of the report is only a record, it does not reach the pipeline).
2. Add a test that the option reaches `run_pipeline` (see `tests/test_cli_smoke.py`).
3. Put it in a preset only in the release-prep pull request, with the measured effect in the CHANGELOG.
4. Update this page; `tests/test_docs_detect_options.py` fails until it names the new option.
