# CLAUDE.md: working rules for MATPredict

Rules for any agent (or person) working in this repository. They are short on purpose; the reasons are in the notes under `docs/notes/`. Add to this file when a rule has cost a run.

## Where work happens

1. **All analysis and data processing runs on the HPCC. The laptop is an orchestrator only.** It edits code, opens PRs, reads results and submits jobs over SSH. It does not run panels, searches or downloads.
2. **Explorations and one-off work stay inside `~/projects/MATPredict`** (on the HPCC: `/bigdata/stajichlab/jstajich/projects/MATPredict`), in a named folder that says what it is and when it can be purged: `results/<YYYY-MM-DD>_<topic>/` for results, `.claude/worktrees/<name>/` for code. A folder that is a scratch test gets a name that says so (`*_test`, `*_probe`) so it can be deleted later.
3. **Temporary files go to `$SCRATCH`** (fast node-local SSD, not billed, cleared when the job ends). That covers downloads (for example `.sra` from NCBI before `parallel-fastq-dump`), uncompressed genomes and intermediate files. Keep only the results worth keeping under `results/`. `/scratch` is node-local: never put a worktree or anything a later job needs there.
4. **Scripts live in `scripts/`, library code in `src/MATPredict/`, tests in `tests/`.** Any compiled helper (for example Rust) goes in its own clearly named place, not mixed into the Python package. Finished code is documented (docstring or `docs/`) and has a test. New options go on the options page `docs/detect-options.md`; a test fails when the page and the CLI disagree.
5. **Dataset metadata comes from the PopSeq repository** (`/bigdata/stajichlab/shared/projects/Fungi/PopSeq`). MATPredict is a consumer of that metadata and of SRA accessions. If the metadata needs a fix, or another NCBI/SRA query is needed, consider fixing it in PopSeq rather than here.

## Git and release

- Work in a git worktree on a branch, never on the main checkout. Check `git branch --show-current` before editing.
- **Open a PR; the curator merges.** Never merge a PR, push to `main` or create a tag or release unless asked in that message.
- Do not edit a tree a Slurm job is running from. A job executes `$SRC` directly, and a half-saved file breaks every genome still to run. Use a worktree and point the job at it.
- Freeze each full run: a worktree at a commit (`run-<sha>`), a frozen copy of `db/`, and the result folders recorded in memory or a note.
- Never use bare `git stash`; the stash is shared across worktrees.

## HPCC practice

- Over SSH use a login shell (`bash -l`), or `sbatch`/`squeue` are not on the PATH. Prefer short commands; long ones go to the background.
- Compute nodes cannot reach GitHub over SSH. Fetch on the login node (`git fetch`, then `git rev-parse FETCH_HEAD`, then `git worktree add --detach`).
- Use `--constraint=milan` (EPYC 7713) for timing-sensitive jobs. The `c01`-`c03` Opteron nodes are about 2.15x slower, so never compare timings across node types.
- Do not run more than one copy of a submit command. A hung SSH call has once duplicated a job into the same output folder.
- Another session may be running reads work from `run-wave3` and `run-expand300` (jobs named `wave-*`, `genus-*`). Do not move those worktrees.

## Detector rules

- **A behaviour change is opt-in.** New logic goes behind a policy flag that defaults to `off`, is switched on by a named preset only after it passes, and is described in `docs/detect-options.md`.
- **Plumb every option.** An option must reach `pipeline_kwargs` in `detect/cli.py`. The `parameters` block of the report is only a record. A CLI passthrough test guards this.
- **Acceptance rule for a speed or refactor change: no change at all** in family, idiomorph, tier or coordinates on the 510 panel genomes (300 Ascomycota, 120 Basidiomycota, 90 Mucoromycotina) plus the 313 FungiDB genomes. Compare `genes_found` as a set (its order can vary). Any difference is examined and approved genome by genome by the curator.
- **Never present a speed number across different node types or different batching.** Compare like with like; state the node type.
- The FungiDB annotation is an independent check, not the truth. Report agreement and disagreement both.
- No quiet fixes: a call that changes is listed with its genome, old value and new value.

## Tests

Run from a worktree:

```
PYTHONPATH=src <env>/bin/python -m pytest tests -q
```

The suite must pass before a PR is opened. Report failures with their output; do not hide a skipped step.

## Reporting

- End a long task with what was done, where it lives (path, branch, PR) and the next command.
- Say plainly when a result is not verified. Retract a wrong claim where it was made.
- Dates are absolute (`2026-11-07`), not relative.
