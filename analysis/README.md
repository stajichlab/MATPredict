# analysis/ — study reports

This folder holds one report per study, an index, the curator's decisions and
the open questions. It is the place to start when investigating a result.

- `INDEX.md` — every study, grouped by theme, with status and key numbers.
- `decisions.md` — curator rulings (J. Stajich) in date order, with evidence.
- `open-questions.md` — pending items and what each waits on.
- `<date>_<topic>.md` — one report per study (related runs are grouped).

## Paths

`results/<folder>/` paths refer to the main checkout
(`/bigdata/stajichlab/jstajich/projects/MATPredict/results/`). Summaries and
small tables of most folders are also committed on the PR #9 branch
(`polish-scope-cuts`) or on a curation branch; a report names the branch when
a folder is committed only there. Per-genome TSVs open in any spreadsheet.

## Report format (use for every new study)

```
# <Title>
Status: decided | open | superseded | running

## Question
One or two sentences.

## Data and code version
- Code: commit <sha> (frozen worktree .claude/worktrees/run-<sha>)
- Database: same tree as the code unless stated
- Inputs: genome list path, panel size

## Method
Short numbered steps. Name the rule, threshold and tool.

## Results
Numbers and small tables. Every number cites its source file.

## What changed in detection
Commits that implement a change, if any.

## Limits
What the data cannot show; sample sizes; circularity.

## Curator decisions
Made (date, ruling) / open (what is needed).

## Files
Paths to NOTE.md, TSVs, PDFs.
```

## Rules

- Plain short sentences. No inflated claims. Say what was verified and what
  was inferred. Mark a number "unverified" if no source file shows it.
- Mark superseded conclusions explicitly, with what replaced them.
- Unpublished manuscript content (`resource/MBE_202608/`, git-excluded) must
  never be copied here. Cite the preprint DOI 10.1101/2025.09.11.675505 only
  for facts it states, and only aggregate comparisons made by this project.

## Pre-sign-off regression check (curator ruling 2026-09-30)

Every new or changed curated record, classifier rebuild, paralog class,
roster/scope change or detection-rule change must attach a regression check
before sign-off. Run `scripts/run_regression_panel.sh` on baseline and
candidate frozen worktrees (or `scripts/regression_check.py diff` on existing
runs), and link `diff/summary.md` plus each panel's `regression_summary.md` in
the study report. The curator reviews every call_lost, call_gained,
idiomorph_changed, confidence_changed and core_model_changed line. The panel is
`testset/regression_panel.tsv` (163 genomes + Zygo 23 both inputs). See
docs/regression-check.md.
