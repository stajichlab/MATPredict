# Uncommitted edits found in HPCC run worktrees (saved 2026-10-05)

These two frozen run worktrees on HPCC had edits that were never committed. They were saved here so they do not depend on that disk. Nothing here is on a merge path.
Apply with `git checkout <base>` then `git apply <diff>`.

| File | Worktree | Base commit | What it is |
|---|---|---|---|
| `run-capprot_pipeline.diff` | `run-capprot` | see `run-capprot_base.txt` (`a2fe1b4`) | Cap-protection experiment in `detect/pipeline.py`: an `MATP_CAP` environment override for the per-family polish cap, plus the gate minimum score passed through. 62 lines added, 6 removed. |
| `run-f017ecb-nocaax_order.diff` | `run-f017ecb-nocaax` | see `run-f017ecb-nocaax_base.txt` (`f017ecb`) | "Finder-off arm": removes the `pheromone_precursor_scan` block from `db/Basidiomycota/order.yml` to isolate the CAAX scan's effect on calls. |

The frozen worktrees `run-541a658-dothideo` and `run-7277e10-dothideo` showed 40 staged additions each: the ten Dothideomycete records, which are now on `main` (PR #36). Nothing else was saved from them.

Also archived on 2026-10-05, as branches `archive/*` on the remote: `aligner-candidate-measure`, `backup-curation-umbelopsis-19e1d12`, `backup-curation-umbelopsis-pre-076afe4`, `merge-measure`, `merge-measure2`, `nf-measure`,
`review-umb-noguard`, `umb-noguard`. Each held commits that existed only on HPCC.
