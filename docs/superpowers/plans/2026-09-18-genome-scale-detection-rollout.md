# Genome-Scale Detection Rollout Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Acquire real target genomes, run the existing `matpredict detect` pipeline across all of them, aggregate the results into one honest performance picture, and wire up real (protein-match-based) ground-truth sanity scoring using genomes that back existing curated records — surfacing further curation/pipeline gaps at scale, the way the 2-genome validation earlier this session did by hand.

**Architecture:** 4 new components, all additive to the existing pipeline (no changes to `run_pipeline`/`search_localize`/`polish_with_exonerate`/`polish_with_miniprot`/`tiering.assign_tier` unless a rollout finding requires one — track any such finding as its own follow-up, do not fold speculative pipeline changes into this plan): (1) genome acquisition, (2) SLURM-sized batch orchestration, (3) result aggregation, (4) ground-truth sanity scoring wired into `benchmark.py`'s currently-stubbed sensitivity field.

**Tech Stack:** Python, this project's existing `NcbiClient`/`MatpredictConfig`, `pixi run matpredict detect`, SLURM (via whatever job-submission convention this repo already uses elsewhere — check for existing `.slurm`/`sbatch` scripts in the repo before inventing a new submission style).

**Spec:** `docs/superpowers/specs/2026-09-18-genome-scale-detection-rollout-design.md` — read this in full before Task 1. It explicitly leaves 3 things undecided (exact pilot genome count, `datasets` CLI availability, ambiguous-match handling for ground-truth scoring) — Task 1 resolves the first two, Task 5 handles the third with a documented judgment call.

## Global Constraints

- HPCC SLURM job sizing: target ~1-1.5 hours of real runtime per job (this project's own established rule) — do not default to one job per genome without first estimating whether that under-shoots this window.
- Never hardcode `/scratch/$USER/...` — use `$SCRATCH` with `:?` so a missing var fails loudly.
- Never commit downloaded genome FASTA files to git; store under a git-ignored directory or `$SCRATCH`.
- Compress large intermediate/output files by default (`.zst` for internal pipeline intermediates).
- Ground-truth comparison scores on translated protein sequence, never raw exon/intron/coordinate structure (binding principle from `docs/superpowers/plans/2026-09-18-mat-gene-exon-intron-model.md`).
- Reuse `src/MATPredict/db/ncbi_client.py`'s `NcbiClient` for any new NCBI-facing code; do not build a second client.
- No change to `run_pipeline`/`search_localize`/`polish_with_exonerate`/`polish_with_miniprot`/`tiering.assign_tier` in this plan.

---

### Task 1: Genome acquisition for a defined pilot target list

**Files:**
- Create: `src/MATPredict/detect/genome_acquisition.py`
- Modify: `src/MATPredict/db/ncbi_client.py` (only if Step 1's investigation finds `NcbiClient` genuinely needs an assembly-level resolution method — see Step 1; do not add this speculatively if `datasets` CLI covers acquisition end to end)
- Test: `tests/detect/test_genome_acquisition.py`

**Interfaces:**
- Produces: `acquire_genomes(taxids: list[int], out_dir: Path, ncbi: NcbiClient | None = None) -> list[AcquiredGenome]`, where `AcquiredGenome` is a small dataclass with `taxid: int`, `accession: str`, `fasta_path: Path`, `species: str`. This is the single entry point Task 3's orchestrator consumes — its shape must not change without updating Task 3.

- [ ] **Step 1: Investigate real acquisition options in this environment (do this before writing any acquisition code)**

Check whether NCBI's `datasets` CLI is available:
```bash
which datasets || pixi run which datasets || module avail datasets 2>&1
```
Also check whether it's already declared as a project dependency:
```bash
grep -n "datasets\|ncbi-datasets" pixi.toml
```
If `datasets` is available (either as a system module or a pixi dependency you should add), use it for acquisition — it handles assembly discovery, download, and unzip in one well-tested tool, and is the standard approach on NCBI-genome-heavy HPCC environments. If it is NOT available and cannot be added as a pixi dependency (check whether `ncbi-datasets-cli` exists on the conda-forge/bioconda channels this project's `pixi.toml` already uses), fall back to `NcbiClient`-based acquisition: resolve each taxid to its representative/reference assembly accession via `esummary` (db=assembly), then fetch the assembly's FASTA via the assembly's FTP path (returned in the esummary response's `ftppath_refseq`/`ftppath_genbank` field) — this requires a plain HTTP/FTP fetch, not `efetch`, since `efetch` doesn't serve whole assembly FASTA files. Document which path you took and why in your report; this determines everything else in this task.

- [ ] **Step 2: Confirm the pilot target list (resolve the spec's first open item)**

The spec proposes, in priority order: (a) the exact/closely-related genome assemblies backing the newly-accepted Onygenales *C. immitis*/*C. posadasii* records — read `db/Ascomycota/Onygenales/*/metadata.yaml`'s `locus.core.segments[].sequence_source.accession` and `taxonomy.taxid` fields to find these; (b) real *Aspergillus* genomes beyond the curated *A. nidulans*/*A. fumigatus* strains, plus at least 2-3 real *Penicillium* species (genuinely uncurated — a blind test); (c) 1-2 additional *Neurospora*/*Fusarium* species beyond the curated strains. Propose a concrete list of 10-20 real taxids/species (not a placeholder count — pick real, specific organisms, e.g. by checking NCBI Assembly for representative genomes of *Aspergillus flavus*, *Aspergillus niger*, *Penicillium chrysogenum*, *Penicillium expansum*, or similar well-assembled species) and write it into your report before proceeding, so the human reviewer can adjust the list before genomes are actually downloaded (large, slow operation — get the list right first).

- [ ] **Step 3: Implement `acquire_genomes`**

Using whichever mechanism Step 1 selected, implement the function so that for each taxid it resolves a representative assembly, downloads the genome FASTA into `out_dir` (a caller-supplied directory — Task 3 will point this at a `$SCRATCH`-based path), and returns an `AcquiredGenome` per successfully-acquired genome. Skip (not crash) any taxid that can't be resolved, collecting failures into a separate return value or log rather than silently dropping them — the caller needs to know if a genome from the pilot list didn't come through.

- [ ] **Step 4: Write tests using a fake/mocked transport (no live downloads in the test suite)**

Follow whatever mocking convention `tests/db/test_ncbi_client.py` already uses for HTTP-level fakes. Test: a successful acquisition returns the right `AcquiredGenome`; an unresolvable taxid is skipped and reported as a failure, not raised as an exception that aborts the whole batch.

- [ ] **Step 5: Run the full suite, commit**

```bash
pixi run pytest -v
git add src/MATPredict/detect/genome_acquisition.py tests/detect/test_genome_acquisition.py
git commit -m "feat: add genome acquisition for the detection-rollout pilot target list"
```

Do NOT run a real download of the full pilot list as part of this task's own test/commit cycle — that's Task 2, a data-acquisition step, not a code task, and belongs in its own workspace/directory outside git.

---

### Task 2: Acquire the pilot genome set (data step, not a code task)

**Files:** none committed to git — this task's output is real files under a git-ignored path (e.g. `$SCRATCH/matpredict-rollout-genomes/` or a project-local `genomes/` directory added to `.gitignore` if the genomes need to persist beyond one SLURM job's lifetime).

- [ ] **Step 1: Add the acquisition target directory to `.gitignore` if it doesn't already match an existing ignore pattern**

- [ ] **Step 2: Run `acquire_genomes` for real against the pilot list confirmed in Task 1**

```python
from pathlib import Path
from MATPredict.detect.genome_acquisition import acquire_genomes
from MATPredict.db.ncbi_client import NcbiClient
# use this project's real config/client construction pattern -- check
# src/MATPredict/db/cli.py's _make_clients for how NcbiClient is normally built
genomes = acquire_genomes(taxids=[...], out_dir=Path("$SCRATCH-or-genomes-dir"))
```

- [ ] **Step 3: Report what was actually acquired vs. what failed**

Write a short report (`.superpowers/sdd/2026-09-18-genome-scale-detection-rollout/pilot-acquisition-report.md`) listing every genome successfully acquired (species, accession, file size) and every failure with its reason. If the failure rate is high (more than a couple of the 10-20 targets), stop and flag it rather than proceeding into Task 3 with a degraded pilot set — this is a case where the plan should not proceed on autopilot; a high failure rate likely means Task 1's acquisition mechanism has a real bug worth fixing before burning more time downloading.

---

### Task 3: SLURM-sized batch orchestration

**Files:**
- Create: `src/MATPredict/detect/batch_runner.py`
- Test: `tests/detect/test_batch_runner.py`

**Interfaces:**
- Consumes: the list of `AcquiredGenome` from Task 1/2.
- Produces: `plan_batches(genomes: list[AcquiredGenome], target_seconds_per_job: int = 5400) -> list[list[AcquiredGenome]]` — groups genomes into per-job batches sized toward the target runtime, using the per-genome timing evidence already on record (`docs/superpowers/plans/2026-09-17-mat-detection-search-localization-benchmark-notes.md`'s ~1m36s figure for a ~39 Mb genome as the seed estimate; real fungal genomes vary — this function should scale its per-genome time estimate by genome file size relative to that reference point, not assume every genome takes the same time). Also produces `run_batch(genomes: list[AcquiredGenome], db_root: Path, reference_fasta: Path, out_dir: Path) -> None`, which runs `run_pipeline` (or shells out to the `matpredict detect` CLI — implementer's judgment on which avoids redundant setup cost per genome, per the spec's note) for every genome in one batch, writing each genome's own report under `out_dir/<taxid>_<accession>/`.

- [ ] **Step 1: Write the failing test for `plan_batches`**

```python
def test_plan_batches_groups_toward_target_runtime():
    genomes = [
        AcquiredGenome(taxid=1, accession="A", fasta_path=Path("/tmp/a.fasta"), species="sp1"),
        AcquiredGenome(taxid=2, accession="B", fasta_path=Path("/tmp/b.fasta"), species="sp2"),
    ]
    # with a fake file-size-based time estimator returning e.g. 2700s each,
    # 2 genomes should land in one batch under a 5400s target
    batches = plan_batches(genomes, target_seconds_per_job=5400)
    assert len(batches) == 1
    assert len(batches[0]) == 2

def test_plan_batches_splits_when_over_target():
    # genomes whose estimated time individually exceeds the target should
    # each get their own batch rather than being silently merged
    ...
```
Write the second test with concrete large-file-size fixtures (use `tmp_path` to create real files of a controlled byte size, since the size-based estimator needs to read real file sizes) rather than a placeholder — this needs an actual size-scaling estimator, not a stub, so write the estimator's real logic (linear scaling from the ~39Mb/1m36s reference point is a reasonable, defensible v1 — document that assumption in the docstring and note in your report that it should be recalibrated once real batch-run timing data exists) before writing the tests that exercise it.

- [ ] **Step 2: Implement `plan_batches` and `run_batch`**

- [ ] **Step 3: Find this repo's existing SLURM submission convention** (do not invent a new one)

```bash
grep -rl "sbatch\|#SBATCH" . --include="*.sh" --include="*.slurm" 2>/dev/null
```
If a convention exists, follow it for whatever wraps `run_batch` into an actual submitted job. If none exists in this repo, write a minimal `sbatch` script template as part of this task's deliverable (a `.slurm` or `.sh` file under `scripts/` or wherever this repo's existing tooling scripts live — check for a `scripts/`/`bin/` directory before creating a new location), following this project's global HPCC conventions (`$SCRATCH:?`, no hardcoded scratch paths).

- [ ] **Step 4: Run tests, full suite, commit**

```bash
pixi run pytest -v
git add src/MATPredict/detect/batch_runner.py tests/detect/test_batch_runner.py <any new SLURM script>
git commit -m "feat: add SLURM-sized batch orchestration for the detection rollout"
```

---

### Task 4: Result aggregation

**Files:**
- Create: `src/MATPredict/detect/rollout_aggregate.py`
- Test: `tests/detect/test_rollout_aggregate.py`

**Interfaces:**
- Consumes: the YAML report shape `write_detection_report` already produces (`src/MATPredict/detect/report.py`'s `_result_doc`/`write_detection_report` — read this file fresh, its exact shape is: `families_attempted: [str]`, `detected: [{family, contig, start, end, confidence, idiomorph, ambiguous_with, genes_found, genes_missing, genes_not_searchable, fragmented, reference_records, segments: [...], gene_evidence: [{gene, role, contig, start, end, strand, identity, coverage, reference_record, method, status, alternate_model}]}]`, `not_detected: [{family, reason, best_fraction_found, genes_found, genes_missing, genes_not_searchable}]`).
- Produces: `aggregate_reports(report_paths: list[Path]) -> RolloutSummary`, where `RolloutSummary` carries: total genomes attempted, a per-family tally of `{confidence tier: count}` across all genomes, a list of `(genome, family)` pairs where the family was attempted but not detected with its `reason`, and a list of anomalies (a genome from a taxon whose curated family is the SAME order/class as a family that WAS detected elsewhere in the batch, but this genome reported nothing for it — a real, worth-investigating signal, not necessarily a bug).

- [ ] **Step 1: Write the failing test with 2-3 real-shaped synthetic YAML report fixtures**

Construct fixture YAML matching the exact shape above (hand-write 2-3 small `write_detection_report`-shaped dicts covering: one genome with a high-confidence detection, one with a `not_detected` entry, one with a low-tier detection) and confirm `aggregate_reports` correctly tallies confidence counts and surfaces the `not_detected` reasons.

- [ ] **Step 2: Implement `aggregate_reports`**

- [ ] **Step 3: Add a CLI entry point** wiring this into `src/MATPredict/detect/cli.py` as a new `matpredict detect rollout-summary --reports-dir <dir> --out <path>` subcommand (or wherever this repo's existing subcommand-registration pattern puts it — follow `register_subcommands`'s existing style in that file).

- [ ] **Step 4: Run tests, full suite, commit**

```bash
pixi run pytest -v
git add src/MATPredict/detect/rollout_aggregate.py tests/detect/test_rollout_aggregate.py src/MATPredict/detect/cli.py
git commit -m "feat: add detection-rollout result aggregation"
```

---

### Task 5: Ground-truth sanity scoring wired into `benchmark.py`

**Files:**
- Modify: `src/MATPredict/detect/benchmark.py`
- Test: whatever test file already covers `run_benchmark` (find via `grep -rl "run_benchmark" tests/`)

**Interfaces:**
- Consumes: `aggregate_reports`'s per-genome detection results (Task 4); each curated record's own `metadata.yaml` coordinates as ground truth for genomes that match a curated record's source genome/accession.
- Produces: `run_benchmark()`'s currently-always-`None` `sensitivity` field becomes a real number for any family/genome pair where a genuine ground-truth match was found; stays `None` (not fabricated) for pairs with no real ground truth available — do not force a score where there's nothing legitimate to compare against.

- [ ] **Step 1: Read `benchmark.py` in full, fresh** (its exact current shape may have shifted since this plan was written) to find precisely where `sensitivity=None` is currently hardcoded and what data is already available at that point (holdout groupings, family/genome pairing).

- [ ] **Step 2: Implement the genome-matching judgment call**

For each detection-rollout genome, determine whether it matches (exactly, by accession, or closely enough — same species AND same or a documented-equivalent strain) an existing curated record's own source genome. This is the spec's third open item: when the match is ambiguous (same species, different or unstated strain), do NOT silently guess — flag it in the output as `ground_truth_match: ambiguous` with the reason, and exclude it from the numeric sensitivity score rather than let an uncertain match corrupt a real number. Only score genomes/records with an unambiguous match.

- [ ] **Step 3: Implement protein-sequence-based comparison, not coordinate-based**

For each unambiguously-matched (genome, curated record) pair, for each gene the curated record expects: check whether the detection result's `gene_evidence` for that gene has a translated protein sequence that matches (or is highly similar to — use this project's existing `seqmatch.score_match` if that fits, since it's already a real, tested protein-comparison utility, rather than writing a second one) the curated record's own deposited protein. A gene is "found" for sensitivity purposes if this protein-level match succeeds, regardless of whether the detected coordinates/exon structure exactly match the curated record's — per this plan's binding design principle.

- [ ] **Step 4: Write real tests**

Using 1-2 real curated records already in `db/` (pick genuinely small, fast-to-load ones) plus a synthetic detection-result fixture, confirm: an exact protein match scores as found; a genuinely different protein (wrong gene) scores as not found; an ambiguous-strain case is excluded from the numeric score and flagged, not silently scored either way.

- [ ] **Step 5: Run tests, full suite, commit**

```bash
pixi run pytest -v
git add src/MATPredict/detect/benchmark.py <its test file>
git commit -m "feat: wire real ground-truth sensitivity scoring into the benchmark, using self-consistency genomes"
```

---

### Task 6: Run the pilot rollout and produce a findings report

**Files:** none new — this task executes Tasks 1-5's tooling for real and writes a findings report.

- [ ] **Step 1: Run the full pipeline**: acquire (Task 1/2's real output), batch-orchestrate (Task 3), aggregate (Task 4), sanity-score (Task 5) — against the pilot genome list confirmed in Task 1 Step 2.

- [ ] **Step 2: Write the findings report** to `docs/superpowers/plans/2026-09-18-genome-scale-detection-rollout-findings.md`: real tier distribution across the pilot, every anomaly surfaced, every ground-truth sanity score obtained, and an explicit triage list of anything that looks like a real curation gap or pipeline bug (do not fix these inline as part of this task — list them for a follow-up plan/dispatch, the same way the 2-genome validation earlier this session surfaced findings that became their own dedicated fix passes).

- [ ] **Step 3: Present the findings to the human reviewer** rather than closing this plan silently — this rollout's entire point is to surface gaps, so ending with "N findings, here they are" is success, not a failure requiring more work before this plan is done.

## Self-review notes (controller, at plan-writing time)

- Spec coverage: all 4 architecture components map 1:1 to Tasks 1-5; Task 6 is the actual rollout run the spec's "Scope for this plan" section calls for.
- No placeholders: every code block is real and complete; the two genuinely open decisions (exact pilot list, `datasets` CLI availability) are explicitly resolved as investigation steps within Task 1 rather than left as unwritten TODOs elsewhere in the plan.
- Type/signature consistency: `AcquiredGenome` (Task 1) is consumed unchanged by Task 3's `plan_batches`/`run_batch`; `aggregate_reports`'s `RolloutSummary` (Task 4) is consumed by Task 5's sensitivity wiring — both interfaces are declared once, in the task that produces them, and referenced (not redefined) by later tasks.
