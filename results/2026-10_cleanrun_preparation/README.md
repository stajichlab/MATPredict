# Clean-run preparation: judged findings (2026-10-22/23)

Written by subagents (model output, reviewed by the curator in the session) and kept here as the record behind section 10 of `docs/notes/2026-10-22_clean-run-plan.md`.

- `regression_judgement.md`: the regression gate (v0.6.0 code `7c7ed99` against the candidate, 163-genome panel plus Zygo 23, defaults and clean-run flags). Every call-touching change traced to a commit by bisecting; finds `0d216c4` (fixed by #140).
- `pilot_judgement.md`: the pilot (565 Ascomycota/Basidiomycota genomes plus 25 *Mucoromycotina*) against v0.6.0, with causes established by re-running 32 genomes with flags toggled, and the supported / tie-break split of `high` calls.

Paths in the reports (`scratchpad/`, `results/2026-10_cleanrun_*`) point at HPCC runs under `/bigdata/stajichlab/jstajich/projects/MATPredict/results/`; the raw per-genome reports are not committed. Some statements in them are marked as not verified by their authors; read those markers. The determinism check (0 differences in 198 genomes, three runs) is summarised in the plan.
