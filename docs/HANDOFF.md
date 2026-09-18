# MATPredict Handoff — 2026-09-18

Read this first in a new session before touching anything. It orients you to
where the project actually is, not where any single spec says it should be.

## What MATPredict is

A general-purpose MAT (mating-type) locus annotator/classifier for fungi,
decomposed into 6 sub-projects. Only sub-projects 1 and 2 have real code so
far.

- **Sub-project 1 — Reference database** (`db/`, `src/MATPredict/db/`):
  curated, evidence-tiered MAT-locus records used as search references and
  as evaluation ground truth.
- **Sub-project 2 — Detection pipeline** (`src/MATPredict/detect/`):
  `matpredict detect --genome <fasta> [--proteins <fasta>] [--taxid <id>]
  --out-dir <dir>` runs a localize-then-polish search against one genome
  and writes a GFF3 + YAML report.
- Sub-projects 3-6 (ab initio prediction hook, ML/homology corpus,
  cross-taxon evaluation harness, whatever comes after) do not exist yet as
  code. Do not assume any scaffolding for them.

## Sub-project 1 — current state

**Coverage** (`db/<Phylum>/<Order-or-Family>/<record_id>/metadata.yaml`,
one directory per accepted record): 3 phyla (Ascomycota, Basidiomycota,
Mucoromycota), 17+ families as of this handoff, including — added this
session — Eurotiales (*A. nidulans*, *A. fumigatus*), Onygenales
(*Coccidioides immitis*, *C. posadasii*), Hypocreales (*Fusarium
graminearum*, *Gibberella fujikuroi*), Sordariales (*Neurospora crassa*),
Helotiales (*Sclerotinia sclerotiorum*, *S. trifoliorum*, *Botrytis
cinerea*), Teloschistales + Lecanorales (6 lichen genera), Pezizales
(*Tuber* — first Pezizomycetes family), plus the pre-existing
Saccharomycetales/Schizosaccharomycetales/Pneumocystidales/Taphrinales/
Ophiostomatales/Diaporthales (Ascomycota) and Agaricales/Ustilaginales
(Basidiomycota) families.

**Held, not accepted:** `db/candidates/Ascomycota/1301609_fen-198_MAT_combined`
(*Xanthocarpia feracissima*) — no NCBI assembly exists for the underlying
reads; revisit if one is deposited.

**Curation workflow:** `matpredict curate-db propose|validate|accept|reject
|build-gff|build-duckdb`. Never self-accept a proposed record — always a
human decision. `validate` runs live NCBI/UniProt checks plus a real
independent-translation check (`sequence_match`, see below) and a
gene-vocabulary check against `db/<Phylum>/order.yml`.

**Schema highlights:** 1-based fully-closed coordinates throughout. Tier-1
(real experimental evidence) vs. tier-2 (genome-annotation/bioinformatic
only) evidence bar, tracked per `locus_existence`/`boundaries`/
`idiomorph_assignment`. Genes can now carry a real exon/intron/frame model
(`exons`, `codon_start`, `transl_table` — added this session, see
`docs/superpowers/plans/2026-09-18-mat-gene-exon-intron-model.md`) instead
of a single naive span; this is what `sequence_match` actually verifies
against a live-fetched deposited protein.

**Known, disclosed residual debt in sub-project 1 (nothing hidden, all
flagged in-record or in this file):**
- `seqmatch.score_match`'s pass/warn/fail status gates on `percent_identity`
  alone; `coverage` is never thresholded. A short-amplicon record can show
  `pass` at very low coverage — read `pass` as "no frameshift detected in
  the covered fragment", not "full locus confirmed." (`src/MATPredict/db/
  seqmatch.py`)
- `gff_export.write_gff3` still emits one flat `gene` span per curated gene
  even where real `exons` now exist — multi-exon curated records don't get
  multi-feature GFF3 output yet. Not wired to anything downstream currently,
  so low urgency, but a real gap if anyone starts consuming curated GFF3
  for exon-level work.
- Lineage-aware `taxonomic_scope` routing (added this session,
  `detect/family_registry.py`) was live-verified on 8 of the ~50 previously
  mis-routed taxids, not all of them. Should hold for the rest (same
  mechanism, no per-species special-casing), but hasn't been swept.
- `gene_class` (functional cross-family gene identity, orthogonal to the
  structural `role` field) is populated in `db/Ascomycota/order.yml` and
  `db/Basidiomycota/order.yml` only where a curated record's own citations
  support the call. Left unclassified: `SXI1`/`SXI2`, `bE`/`bW`, and the
  `MATsc`/`MATyl`/`MTL`/`PM`/`mat1`-family gene lists.
- A **paralog-identity/reporting-model limitation** is diagnosed but not
  fixed: when a family has multiple real gene copies (e.g. duplicated
  pheromone-receptor alleles), detection currently collapses them to one
  reported evidence entry. Flagged as a larger architectural change,
  deliberately deferred.
- Design principle, binding on all future work here: **protein-sequence
  match is ground truth; exon/intron structural divergence between two
  records, or between a curated record and a detection result, is expected
  and never itself an error.** See `matpredict_esm2_future_direction.md` in
  the memory system for the full statement and its origin.

## Sub-project 2 — current state

**Architecture** (`docs/superpowers/specs/2026-09-17-mat-detection-search-
localization-design.md` is the binding spec): Stage 0 routes a genome to a
fast path (diamond, proteome-only, cheap) or a full localize-then-polish
path; Stage 1 is batched `tblastn` genome-wide localization; Stage 2 polishes
each candidate window with BOTH `exonerate --refine region` and `miniprot`
(agreement is reported per gene but does NOT feed confidence tiering by
design); Stage 3 classifies each gene's status
(`polished_agree`/`polished_disagree`/`polished_single`/`unpolished`/
`STATUS_NOT_POLISH_CANDIDATE`); Stage 4 emits a GFF3 + YAML report
(`src/MATPredict/detect/report.py`).

**Validated empirically** (not just unit-tested) against 2 real genomes
this session: a ~39 Mb *Coprinopsis cinerea*-class genome (1m36s vs. a
>300s timeout pre-revision) and a second real genome for correctness
(a medium→high tier improvement from fixing a real multi-gene-window
polishing bug). See `docs/superpowers/plans/2026-09-17-mat-detection-
search-localization-post-fix-validation-notes.md` and the sibling
`-benchmark-notes.md`.

**Known, disclosed residual debt in sub-project 2:**
- `src/MATPredict/detect/benchmark.py`'s `run_benchmark()` computes
  species/genus holdout groupings but **always reports `sensitivity=None`**
  — there is no wired-up recall/precision scoring against curated ground
  truth yet. This is the single most relevant gap for whatever you build
  next (see the companion spec/plan this handoff was written alongside).
- Fast-path zero-hit rescue, cluster-aware rescue-eligibility, and the
  fragmented-locus reporting path all have real, committed fixes and unit
  tests, but are under-exercised by real curated data — the investigation
  this session found one real usable example for the cheap fast-path case
  and confirmed no real candidate yet exists in the curated DB for the
  zero-hit-rescue or fragmented-locus scenarios specifically.
- No genome-acquisition or batch-orchestration tooling exists yet. `detect`
  runs against exactly one genome per invocation; there is no script that
  fetches real target genomes or runs `detect` across many of them.

## Where to look for detail, not summary

- Ledgers (per-plan execution history, every ruling, every review finding):
  `.superpowers/sdd/<plan-basename>/progress.md` while a plan is in
  progress; deleted once a plan's final review is clean (git history is the
  record after that — `git log` for anything you don't find in a ledger).
- Specs: `docs/superpowers/specs/`. Plans: `docs/superpowers/plans/`.
- Memory (persistent across sessions, auto-loaded):
  `~/.claude/projects/-bigdata-stajichlab-jstajich-projects-MATPredict/memory/MEMORY.md`
  and its linked files — covers user domain-expertise notes (e.g. PR-locus
  gene-duplication is normal, don't default to low confidence for it), the
  annotation-gap lesson (prefer classic locus-specific GenBank deposits
  over modern whole-genome annotation for tiny MAT genes), and the
  ESM2/structural-divergence future direction.

## Workflow conventions established this session (keep following them)

- Work in place on `main`, no worktree, for this project.
- Every non-trivial code/data change goes through the full
  subagent-driven-development loop: fresh implementer → task review → fix
  round(s) if needed → final whole-branch review on the most capable model
  → one fix wave → one scoped re-review. Every review in this session
  independently re-verified claims against live NCBI/PubMed data rather
  than trusting the implementer's report — keep doing that; it caught real
  bugs every single time it was done seriously.
- Curated-record proposals are never self-accepted. A human (you) makes the
  accept/reject/hold call every time, informed by the tier table and a
  named reason for anything held or rejected.
- Attribution for commits/PRs from this session: see the standing
  system-reminder in-session for the exact `Co-Authored-By`/`Claude-Session`
  lines — do not guess, re-check the current reminder each session since
  this has changed format before.
