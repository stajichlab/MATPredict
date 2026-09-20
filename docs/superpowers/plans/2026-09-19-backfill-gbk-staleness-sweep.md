# Plan: DB-wide `locus.gbk` staleness sweep (visualization follow-up #1)

## Spec

`docs/HANDOFF-visualization-followups.md` §"1. DB-wide `backfill-gff` staleness
sweep". That section is the binding authority for this plan.

## Context (verified against live data on 2026-09-19, not assumed)

A full scan of `db/*/*/*/metadata.yaml` (59 accepted, non-candidate records):

- **56 of 59** records have a `locus.gbk` with **zero `CDS` features** and an
  all-`N` `ORIGIN` block — written before `gff_export.write_genbank` gained real
  sequence + real `CDS`/`translation` output. The 3 already-good records are
  `199306_rmscc1040_MAT_MAT1-1` (4 CDS), `2903220_liq80xsp_MAT_combined` (4 CDS),
  and `199306_silveira_MAT_MAT1-2` (1 CDS).
- **All 59** records have a non-empty `proteins.faa`. Therefore
  `find_records_missing_proteins_faa` returns `[]` and `curate-db backfill-gff`
  as it exists today would skip **every one** of the 56 stale records. The
  handoff's suspicion is confirmed: a new selector is required.
- Segment `sequence_source.type` across the 61 segments: **56 `insdc_nucleotide`,
  5 `assembly`**. `write_genbank` only fetches real nucleotide sequence for
  `insdc_nucleotide`, so the 5 `assembly`-sourced segments would stay all-`N`
  even after a re-run. Those 5 segments each carry a real, individually
  fetchable `seq_region` contig accession — verified live against NCBI efetch:
  `NW_026089539.1` and `JAAGWA010000001.1` both return real sequence.
  Affected records: `5334_h4-8_Aalpha_4`, `5334_h4-8_Balpha_3`,
  `5334_h4-8_Bbeta_2`, `5346_a43-b43-okayama-7_HD_A43`,
  `5346_a43-b43-okayama-7_PR_B43`.
- Network and the warm `.matpredict_cache` (834 entries, 43 MB) are both live.

## Global Constraints

- **Do not write a second locus-regeneration path.** Reuse
  `build_gff_for_record` (`src/MATPredict/db/cli.py`) — it is already the single
  place this logic lives, and the single-record `build-gff` command and the
  batch backfill command both call it.
- **Per-record failure isolation.** One record's live-fetch failure must never
  abort the sweep, matching `backfill_missing_proteins_faa`'s existing
  discipline. Failures are reported, never silently swallowed.
- **Never fabricate sequence.** An unfetchable segment falls back to the all-`N`
  placeholder for that segment only, exactly as `write_genbank` already does.
- Tests run with `pixi run -e test pytest -v` — **not** bare `pixi run pytest`.
- Staleness is defined by **`CDS` feature count == 0**, not by absence of
  `join(`. Many curated genes are single-exon and legitimately never produce a
  `join(` even when freshly regenerated (verified: 39 of the 56 stale records
  have zero multi-exon genes). A `join(`-based test would mis-classify them as
  permanently stale.

## Task 1 — staleness selector, `assembly` seq_region fallback, CLI flag

Files: `src/MATPredict/db/cli.py`, `src/MATPredict/db/gff_export.py`,
`tests/db/test_backfill_gff.py`, `tests/test_cli_smoke.py`.

Write tests first (project uses TDD).

1. **New selector `find_records_with_stale_locus_gbk(db_root) -> list[tuple[str, str, str]]`**
   in `src/MATPredict/db/cli.py`, beside `find_records_missing_proteins_faa` and
   sharing its shape and conventions (skip the `candidates/` tree; sorted
   `db_root.glob("*/*/*/metadata.yaml")`; return `(phylum, order_or_family,
   record_id)` triples).

   A record is **stale** when it has at least one `present` gene AND its
   `locus.gbk` is missing, OR parses to zero `CDS` features across all its
   `SeqRecord`s. Parse with `Bio.SeqIO.parse(path, "genbank")`, not a text
   search — a substring check on `"     CDS  "` is column-fragile. A `locus.gbk`
   that fails to parse counts as stale (regenerating it is the fix).

   Required tests:
   - a record whose `locus.gbk` has real `CDS` features is NOT returned;
   - a record whose `locus.gbk` has only `gene` features IS returned;
   - a record with no `locus.gbk` at all IS returned;
   - a record under `candidates/` is never returned;
   - a record whose only genes are `present: false` is NOT returned (nothing to
     regenerate).

2. **`write_genbank` `assembly` fallback** in `src/MATPredict/db/gff_export.py`.
   Today the fetch is gated on `source.get("type") == "insdc_nucleotide"`. Widen
   it: when `type` is `"assembly"` (an assembly-level `GCA_`/`GCF_` accession
   `NcbiClient` cannot resolve), fetch using `source["seq_region"]` instead —
   the contig/scaffold accession the coordinates are actually relative to.
   For `insdc_nucleotide`, keep using `source["accession"]` exactly as now.
   If neither yields a usable accession, fall back to the all-`N` placeholder
   for that segment, unchanged.

   Document in the docstring WHY this is correct: the record's `start`/`end`
   coordinates are relative to `seq_region`, not to the assembly accession, so
   `seq_region` is the right thing to subrange-fetch in both branches — the
   assembly accession was never a usable efetch target. Both real values were
   verified live (see Context).

   Required tests (mock `ncbi`, no network in tests):
   - an `assembly`-typed segment fetches with its `seq_region` value and the
     returned sequence lands in the written GenBank;
   - an `insdc_nucleotide` segment still fetches with its `accession` value
     (regression guard);
   - a segment whose fetch raises still writes the all-`N` placeholder.

3. **CLI wiring.** Give `curate-db backfill-gff` a `--stale-gbk` flag that
   selects records via `find_records_with_stale_locus_gbk` instead of
   `find_records_missing_proteins_faa`. Default (flag absent) behavior must be
   byte-identical to today's — existing tests must keep passing unchanged.

   Refactor `backfill_missing_proteins_faa` minimally so the two selectors share
   one failure-isolating driver loop rather than duplicating it. Keep the
   existing public function name working (other tests import it).

   Add a `tests/test_cli_smoke.py` smoke test for `--stale-gbk` argparse wiring,
   matching the existing `draw-locus` smoke tests' style.

4. Also print, per succeeded record, a warning when the regenerated `locus.gbk`
   still contains a placeholder-only segment, so the sweep's output names
   exactly which records remain without real nucleotide sequence.

## Task 2 — run the sweep over the real database

Not a code task. Runs live network.

1. Record the pre-state: the 56 record ids and their `CDS` counts (0).
2. Run `matpredict curate-db backfill-gff --stale-gbk` against `db/`. Expect
   real NCBI/UniProt traffic and real wall-clock. Expect some individual records
   to fail; they are reported, not fatal.
3. Re-scan `db/` and report, factually: how many records now have ≥1 `CDS`
   feature, how many still have zero and why, how many still carry an all-`N`
   segment, and any per-record failures with their messages.
4. Spot-verify at least 3 regenerated records by eye against their
   `metadata.yaml`: CDS coordinates match the curated gene/exon coordinates, the
   `translation` qualifier matches the corresponding `proteins.faa` entry, and a
   multi-exon minus-strand gene renders as `complement(join(...))` in
   descending genomic order (the transcript-order convention — do NOT reverse
   it).
5. Confirm `matpredict curate-db draw-locus` now produces a gene-structure
   diagram for a record that previously had no `CDS` features.
6. Commit the regenerated `db/**/locus.gbk` files (and any `locus.gff3`/
   `proteins.faa` that legitimately changed) as a data commit, separate from
   Task 1's code commit.
