# Plan: make detection CDS/translation output reachable (handoff item 2)

## Spec / authority

`docs/HANDOFF-visualization-followups.md` §2, "Wire `genome_fasta` into
`detect`'s real CLI/batch-runner path". That section is binding, with the
corrections and additions recorded under Context below.

## Context (read from the real code on 2026-09-19, not assumed)

- `report.write_detection_gff3(outcome, out_path, genome_fasta=None)` already
  implements everything: a companion FASTA at `out_path.with_suffix(".fasta")`
  and a sibling `CDS` feature per gene carrying `translation=`. It is tested
  (`tests/detect/test_report.py:301,339,356`).
- **It is unreachable in production.** Both production call sites omit
  `genome_fasta`: `detect/cli.py:41` and `detect/batch_runner.py:226`. No CLI
  flag exposes it. The deliverable has never been produced by a real run.
- **Correction to the handoff**: it says to thread the parameter "through to
  `run_pipeline`'s `write_detection_gff3` call". `run_pipeline` does NOT call
  `write_detection_gff3` — the CLI and the batch runner each call it directly
  after `run_pipeline` returns. So no `run_pipeline` signature change is
  needed; wiring is at the two call sites only.
- `batch_runner` already has the genome on local disk as `decompressed_path`
  (`batch_runner.py:219`), still valid at the write site (the `finally` that
  removes it runs after). Reuse it; do not refetch or re-decompress.
- **Storage, not in the handoff**: the companion FASTA writes each referenced
  contig's FULL sequence (`report.py`, `contig_sequences[record.id] =
  str(record.seq)`), not the locus slice. This is deliberate — the GFF3 carries
  absolute contig coordinates, which only line up against a full contig. The
  consequence is that a locus on a large chromosome yields a FASTA of that
  chromosome's whole length, per genome. This is why the feature must be
  opt-in at BOTH call sites rather than always-on.

## Global Constraints

- **Default behavior must not change.** With the new flag absent, `detect` and
  `run_batch` must produce byte-identical output to today. Existing rollout
  artifacts must stay reproducible.
- **Never fabricate sequence.** The existing graceful-degradation behavior is
  binding: a contig the GFF3 references but the genome FASTA lacks is logged
  and omitted; a gene whose sequence cannot be extracted gets no `CDS` feature.
  Neither is an error and neither aborts the write.
- **One bad genome must not sink a batch.** `run_batch`'s existing per-genome
  try/except discipline is binding; emitting sequence must stay inside it.
- Tests run with `pixi run -e test pytest -v` -- never bare `pixi run pytest`.
- `GeneEvidence.exons` in the DETECTION schema is always ASCENDING genomic
  order regardless of strand, and `_splice_transcript` performs the explicit
  minus-strand reversal. This is the OPPOSITE convention to sub-project 1's
  curated `gene["exons"]` (transcript order). Do not conflate them; do not
  "fix" either.

## Task 1 -- wire the flag, fix the N+1 scan, wrap the FASTA

Files: `src/MATPredict/detect/cli.py`, `src/MATPredict/detect/batch_runner.py`,
`src/MATPredict/detect/report.py`, and their tests.

Tests first.

1. **CLI flag.** Add `--emit-cds-fasta` (a `store_true`) to the `detect`
   subparser. When set, `_cmd_detect` passes `genome_fasta=Path(args.genome)`
   to `write_detection_gff3`. When absent, it passes nothing, exactly as today.
   Document in the flag's help that it additionally writes a companion FASTA
   containing the full sequence of every contig the results reference, and that
   this can be large.

2. **Batch runner.** Give `run_batch` an `emit_cds_fasta: bool = False`
   keyword-only parameter, defaulting False so every existing caller is
   unchanged. When True, pass the already-on-disk `decompressed_path` as
   `genome_fasta`. Do not re-decompress, do not refetch, and keep the call
   inside the existing per-genome try/except.

3. **Fix the N+1 genome scan.** `write_detection_gff3` currently calls
   `_extract_translated_gene(genome_fasta, ...)` once per gene -- each call
   opens and parses the entire genome FASTA -- and then parses it once more for
   the companion FASTA. Parse the genome ONCE per call and reuse it.

   Do this without duplicating `benchmark.py`'s extraction/splicing logic,
   which is already tested (`_extract_translated_gene`/`_splice_transcript`,
   including the minus-strand exon-order handling). The preferred shape is to
   read the needed contigs once into a `dict[str, str]` (the companion-FASTA
   path already builds exactly this) and have the per-gene extraction work from
   that in-memory mapping. If that requires a small, additive refactor of
   `benchmark.py` -- e.g. splitting the sequence-fetch from the
   splice-and-translate so both callers share the latter -- do that rather than
   copying the logic. `benchmark.py`'s own callers and tests must keep passing
   unchanged.

   Required test: a genome FASTA is opened/parsed exactly ONCE for an outcome
   with 2+ genes. Assert on a real counter (e.g. monkeypatch the open/parse
   helper and count calls), not on timing.

4. **Wrap the companion FASTA.** Sequence lines are currently written
   unwrapped, one contig per line. Wrap at 60 columns (standard FASTA, and what
   the rest of this project emits). Required test: no line in the companion
   FASTA exceeds 60 characters except the `>` headers.

5. **Measure, do not guess.** In the report, state the real companion-FASTA
   size produced for at least one real genome, and the real wall-clock
   difference the N+1 fix makes for a multi-gene outcome. If no real genome is
   available in the environment, say so plainly rather than estimating.

## Task 2 -- prove it end to end on a real genome

Not a code task.

1. Find a real genome already on disk in this repo/environment (check
   `testset/`, and any rollout output directory). Do NOT download one if a
   suitable genome already exists; say which you used.
2. Run `matpredict detect --genome <g> --out-dir <tmp> --emit-cds-fasta` and
   confirm: `detected_loci.gff3` gains `CDS` lines with `translation=`
   attributes, the companion `detected_loci.fasta` exists, its sequence lines
   are <= 60 cols, and its contig names match the GFF3's seqids exactly.
3. Verify at least one emitted `translation=` is a real protein: no internal
   stop codons, and it corresponds to the gene's reported coordinates.
4. Run the SAME command WITHOUT the flag and confirm the GFF3 is byte-identical
   to what the current `main` produces -- this is the no-regression check that
   matters most, since existing rollout artifacts must stay reproducible.
5. Report the measured companion-FASTA size and the genome's size, so the
   storage cost of a rollout-scale run is a real number rather than a guess.
6. Do NOT commit any genome, FASTA, or detection output into the repo. Work in
   the scratchpad directory.
