# MAT Locus Visualization — Follow-ups Handoff (2026-09-19)

Read this first in a new session before touching this area. It orients you
to exactly what's left on the visualization work, not a full project
status — for that, see `docs/HANDOFF.md`.

## What's already done and merged (do not redo)

A full plan (`docs/superpowers/plans/2026-09-18-mat-locus-visualization.md`,
now fully executed — commits `8bb132e..1c43f2c` plus a final-review fix
wave `c1e7769`/`574c06a`) built:

- `src/MATPredict/db/gff_export.py::write_genbank` — now writes real
  segment sequence (was an all-N placeholder) and real `CDS`/`translation`
  features (was `gene`-only), including a real `CompoundLocation` for
  multi-exon genes built from `gene["exons"]` (transcript-order convention:
  ascending for `+` strand, descending for `-` strand — already stored that
  way in the schema, do not reverse it again).
- `src/MATPredict/detect/report.py::write_detection_gff3` — gained an
  optional `genome_fasta` parameter that, when given, writes a companion
  FASTA + `CDS`/`translation` features for detection results (reuses
  `benchmark.py`'s already-tested `_extract_translated_gene`/
  `_splice_transcript`).
- `matpredict curate-db draw-locus --phylum <p> --order-or-family <o>
  --record-id <id> --out <path>` — single-locus gene-structure diagram via
  pyGenomeViz (`src/MATPredict/db/draw.py`). Single-segment records only;
  fragmented (multi-segment) records raise `NotImplementedError` on
  purpose (a real pyGenomeViz limitation found by hands-on trial, not an
  oversight — see the research doc below for why).
- `matpredict curate-db draw-synteny --record-ids <id1> <id2> ... --out
  <path>` — multi-locus comparison diagram via clinker
  (`src/MATPredict/db/synteny.py`). Handles fragmented records fine
  (clinker's own capability, confirmed live — asymmetric with
  `draw-locus` on purpose, already documented in `synteny.py`'s
  docstring). Auto-generates clinker's `--gene_functions`/`--colour_map`
  CSVs from curated `role`/`gene_class` schema fields. Warns (doesn't
  raise) when a resolved record has zero `CDS` features.
- `pygenomeviz`/`clinker` are lazy-imported (inside `_cmd_draw_locus`/
  `_cmd_draw_synteny` in `src/MATPredict/db/cli.py`, not module scope) and
  declared as a `[project.optional-dependencies] viz` extra in
  `pyproject.toml` — this was a real, measured ~50% CLI startup-time fix,
  don't move the imports back to module scope.
- `draw_locus` walks GenBank features positionally (not keyed by gene
  name) specifically so records with duplicate gene names (real
  Basidiomycota B-locus multi-copy pheromone/receptor cases — normal
  biology, see the domain-knowledge memory) draw every gene, not just one
  per name. Don't reintroduce a name-keyed dict here.

**Full research context**: `docs/notes/2026-09-18_mat-locus-visualization-research.md`
(includes a real, evidence-based §7 appendix from hands-on tool trials —
exact clinker/pyGenomeViz source line citations, real error messages,
real package names). Read this before touching either tool's integration
again; it already answers most "does X work with Y" questions empirically.

**Tests**: 277 passing as of this handoff. Run with `pixi run -e test
pytest -v` — **not** bare `pixi run pytest`, which can resolve to a stray
system Python lacking `pygenomeviz`/`clinker` (a real, confirmed
environment gotcha, cost real debugging time twice already).

## What's left — 7 real, named follow-ups, none blocking, none started

Ordered roughly by value/effort, not strict priority — use judgment.

### 1. DB-wide `backfill-gff` staleness sweep (the big one)

**56 of 59 accepted curated records still have pre-fix `locus.gbk` files**
— all-N placeholder sequence, zero `CDS` features, generated before the
`write_genbank` fix above existed. This is the dominant state of the DB
today, confirmed by a full scan during the visualization plan's final
review. Two records were caught and fixed one-off during that plan's own
work (`199306_rmscc1040_MAT_MAT1-1`, `2903220_liq80xsp_MAT_combined`) —
the other 54 haven't been touched.

**Fix**: run `matpredict curate-db backfill-gff` — wait, check first
whether that command's own `find_records_missing_proteins_faa` check
(which only looks for a missing/empty `proteins.faa`, not a stale
`locus.gbk`) will actually re-process these 56 records or skip them
(most of them likely already have a real, non-empty `proteins.faa` from
the earlier `proteins.faa`-backfill plan, which would make
`backfill-gff` skip them even though their `locus.gbk` is still stale).
If it skips them, you'll need a small new script/command that regenerates
`locus.gbk` specifically for every record with real `exons` data in
`metadata.yaml` but no `join(` in its current `locus.gbk` — reuse
`build_gff_for_record` (`src/MATPredict/db/cli.py`) directly, don't write
a second regeneration path.

**This needs real NCBI/UniProt network access** (fetching each record's
protein sequences) — budget real time, and expect some records to fail
individually (isolate failures per-record, matching this project's
established convention, don't let one bad record abort the whole sweep).

### 2. Wire `genome_fasta` into `detect`'s real CLI/batch-runner path

`write_detection_gff3`'s `genome_fasta` parameter (built and tested) is
currently **unreachable in production** — neither `src/MATPredict/detect/
cli.py`'s `_cmd_detect` nor `src/MATPredict/detect/batch_runner.py`'s
`run_batch` passes it, and no CLI flag exposes it. This means the whole
"real CDS/translation for detection results" deliverable has never
actually been produced by a real run. Needs: a new `--emit-cds-fasta`-style
flag (name it whatever fits) on `matpredict detect`, threaded through to
`run_pipeline`'s `write_detection_gff3` call, plus the analogous wiring in
`batch_runner.run_batch` (which already writes each genome's own
`genome_fasta` path — reuse it, don't refetch).

**Bundle with this** (both currently unreachable/latent until #2 lands, so
fix together): `_extract_translated_gene` re-parses the whole genome FASTA
once per gene plus once more for the companion FASTA (N+1 linear scans —
fine for a handful of genes, real cost at scale); the companion FASTA
writes each contig as one giant unwrapped line (wrap at 60-80 cols,
matching standard FASTA convention and this project's line-length norms
elsewhere).

### 3. Per-record clinker cluster naming (cosmetic)

Every `locus.gbk` file is literally named `locus.gbk`, so a multi-record
`draw-synteny` comparison shows every cluster captioned "locus" in
clinker's output (only the inner contig/accession label distinguishes
them). Fix: in `src/MATPredict/db/synteny.py`, stage each resolved
`locus.gbk` into a temp dir as `<record_id>.gbk` before invoking clinker,
so the cluster names are meaningful.

### 4. Missing `draw-synteny` CLI smoke test

`tests/test_cli_smoke.py` has smoke tests for `draw-locus`'s argparse
wiring but not `draw-synteny`'s (`--record-ids` nargs parsing,
`_cmd_draw_synteny`). Add one matching the existing `draw-locus` smoke
tests' style.

### 5. CSV-clobber risk in `draw_synteny`

`draw_synteny` writes its generated `--gene_functions`/`--colour_map` CSVs
beside `out_path` using `out_path.stem` — two concurrent runs sharing an
output directory/stem would overwrite each other's CSVs mid-run. Low risk
in practice (this is a manual/on-demand command, not a batch process
today), but a real bug if it's ever automated. Fix: write the generated
CSVs to a `tempfile.mkdtemp()`-style scratch dir instead.

### 6. Duplicated `_FALLBACK_COLOR` constant

`src/MATPredict/db/draw.py:43` and `src/MATPredict/db/synteny.py:81` both
define the same fallback color independently. `synteny.py` already
imports `ROLE_COLORS` from `draw.py` — import `_FALLBACK_COLOR` the same
way and delete the duplicate.

### 7. (Optional, lowest priority) `matpredict draw-locus` for fragmented records

Currently raises `NotImplementedError` on purpose. pyGenomeViz's real
GFF3 parser has 2 confirmed bugs with this project's fragmented-locus
shape (documented in the research doc's §7: silent multi-seqid dropping,
and a `FeatureRangeError` from absolute-vs-relative pragma coordinates).
If this is ever wanted, the real fix is either (a) rewrite the
`##sequence-region` pragma to `1 <length>` and rebase feature coordinates
to match before handing the GFF3 to pyGenomeViz, or (b) bypass `Gff`
entirely and build `FeatureTrack`/`FeatureSegment` objects directly from
parsed GFF3 data. Not started; genuinely more work than it looks, budget
accordingly.

## Relevant persistent memory (auto-loaded in Claude Code, but listed here
in case you're not using memory or want to read it directly)

- `~/.claude/projects/-bigdata-stajichlab-jstajich-projects-MATPredict/memory/`
  — `matpredict_pr_locus_domain_knowledge.md` (why duplicate gene names are
  normal, relevant to item 1's fix above and to why `draw_locus`'s
  positional-not-name-keyed fix matters), `matpredict_esm2_future_direction.md`
  and `matpredict_pfam_domain_future_direction.md` (unrelated future
  directions, not blocking, just context).

## Workflow convention for this project (keep following it)

Every non-trivial code/data change goes through the full
subagent-driven-development loop (`superpowers:subagent-driven-development`
skill): fresh implementer subagent → task review → fix round(s) if needed
→ final whole-branch review on the most capable model → one fix wave → one
scoped re-review. Every review should independently re-verify claims
against live data/real files rather than trusting the implementer's report
— this caught real bugs every single time it was done seriously across
this whole project's history so far.
