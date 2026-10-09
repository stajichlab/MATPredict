# Changelog

All notable changes to MATPredict. Versions follow [semantic versioning](https://semver.org/);
release tags may carry a clade suffix naming the lineage whose infrastructure that
release completed.

## [Unreleased]

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
- When two curated reference records tie exactly in alignment score for one gene, the polisher now names the one with the lower
  record id instead of whichever alignment the tool listed first (the tool's order among equal scores is not guaranteed). Found on
  Leppa1, where one gene model was attributed to the A1163 record in most runs and to the Af293 record in some. Affects only the
  record named in a model and in `reference_records`; calls, spans and gene sets are unchanged.
- PDF report: gene-evidence tables no longer run past their card. In print, header, number and status cells may wrap, and
  table cells break long names (for example `fungal_mating_type_pheromone`). Found on real campaign reports (Serpula
  lacrymans, 78 px past the card). New layout test `tests/report/test_report_layout.py` checks every fixture, including one
  real campaign report.
- `detect --genome` on a gzip or zstd file now stops at once with a message that says to decompress it (it used to fail
  minutes later with a BLAST or index error). Recognised by the first bytes, not the file name.

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
