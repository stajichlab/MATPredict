# Detect targeting + Mucoromycota validation — handoff (2026-09-20)

Read this before touching detection scope, the evidence floor, or the Zygo test set.
Branch `detect-targeting`, **13 commits ahead of `main`, 377 tests passing, NOT MERGED**.

## TL;DR

The targeting work is **done and reviewed**. The validation run is **done**. The fixes
that validation revealed are **NOT implemented** — nothing has been re-run with them.
That is the next step.

## What shipped on this branch (all reviewed clean)

| # | Change | Effect |
|---|---|---|
| 1 | Phylum routing + routing-aware reference FASTA | query set **181 -> 19** proteins for Mucoromycota; `--phylum` override; `routing_mode` recorded in every report |
| 2 | `max_cluster_gap_bp` per locus in `order.yml` | Mucoromycota `MAT` = 50000, others default 25000; run takes the MAX over routed families |
| 3 | `EvidenceFloor` defaults `min_hits=2, require_core_role=True` | polish admission = >=2 distinct genes with >=1 core_MAT |
| 4 | `run_batch` + `run_detection_batch.py` + `.slurm` wrapper | a batch pays the same narrowed cost; `evidence_floor` opt-out; wrapper tail-forwards flags |

Folded in from the abandoned `detect-emit-cds` branch: `--emit-cds-fasta`, an N+1
genome-parse fix, 60-column FASTA wrapping. **First proven end to end in this run.**

## Measured results (real runs, not estimates)

**23 ground-truth Mucoromycota genomes** (`testset/Zygo/{plus,minus}.abspres.csv`):

- 23/23 loci detected **on the correct ground-truth scaffold**
- per-gene recall **100%**: `tptA` 23/23, `rnhA` 23/23, `sexP` 16/16, `sexM` 7/7
- mean **41.8 s/genome**, max 55.1 s, 0 errors
- for contrast: the same pipeline pre-fix was killed after **>24 min on ONE small yeast
  genome without completing**

**44 previously-unrepresented genera** (1 representative each, 43 completed):
6 loci found at 66-89% identity — *Ellisomyces*, *Rhizomucor*, *Thamnidium*,
*Blakeslea*, *Gongronella*, *Choanephora*.

Full analysis: `docs/notes/2026-09-20_zygo-mucoromycota-validation.md`.

## The four open defects, in priority order

### 1. Idiomorph calling fails in 23/23 — ONE root cause, fix is measured

`sexM` and `sexP` share an HMG domain, so both curated references hit the SAME locus gene
(92-100% coordinate overlap, all 23 genomes). Consequences, both from this one defect:

* `idiomorph=undetermined` in 23/23 — the tool cannot call Plus vs Minus.
* Because both are "found", two idiomorphs are always named, so
  `expected_genes_for_idiomorph` (`family_registry.py:164-181`) never narrows the expected
  roster — its documented conservative fallback — inflating `fraction_found`'s denominator
  from ~5-6 to 7 (`scoring.py:35-40`).

**Measured fix candidate: take the higher-identity member of an overlapping, mutually
exclusive idiomorph pair. Correct 23/23 (100%) on the ground-truth set.**

Margins are wide for Plus (sexP 45-100% vs sexM 27-36%) and **thin for Minus
(2.3-13.0 points)** — narrowest `Cunninghamella_bertholletiae_NRRL_1376`, sexP 28.3 vs
sexM 30.6. Do not treat the rule as safe without also fixing #4.

This is an ALGORITHMIC/biological decision, deliberately left to the curator.

### 2. The ambiguity floor rejects fragmented loci

All 37 "no locus" genera in the discovery sweep were rejected by the **ambiguity floor
(0.50)** — not the evidence floor, not for lack of hits. **15 were at exactly 0.429
(3 of 7), one gene short**; 12 of those found `rnhA` + `sexM` + `sexP`.

Clearest case: *Rhizopus acetoinus*, a cluster at **92.6% identity, 3 genes, admitted to
polishing**, discarded at 0.429, locus at `scaffold_516:467-7397` — 467 bp from the contig
edge, so its flanking genes are on adjacent contigs.

**Fixing #1 does NOT rescue this one.** Recomputed resolved-to-Plus: expected 6, found
`sexP, btbA` = 0.33, still under 0.50. Genuinely fragmented loci need the **relaxed second
pass the curator proposed on 2026-09-19 and which was never settled** — specifically
whether it runs always, or only when the strict pass finds nothing genome-wide. That
question is still open and blocks this.

### 3. One `exonerate` SIGSEGV costs an entire genome

*Gilbertella persicaria*: `exonerate --model protein2genome ... --refine region ... exited
-11` on window `scaffold_71:24582-37340`. MATPredict propagated it as fatal, losing every
locus already found in that genome. Observed rate **1/44 (2.3%)**. A bad polish window
should downgrade that gene to `unpolished`, not abort the run. `run_batch`'s per-genome
try/except contains this to one genome in a rollout, but nothing degrades within a genome.

### 4. Reference coverage is why Minus margins are thin

The curated DB holds 3 Minus + 3 Plus Mucoromycota records (*Mucor circinelloides*,
*Phycomyces blakesleeanus*, *Rhizopus arrhizus*) and **none close to *Cunninghamella*,
*Absidia*, *Chaetocladium* or *Actinomucor* — 5 of the 7 genera tested**. All 7 Minus
genomes came out `confidence: medium`; 14 of 16 Plus came out `high`.

The 23 validated loci (plus the 6 new genera) carry real CDS features with translations
and are **candidates for promotion into the curated DB**, which addresses this directly.

### 5. Minor: evidence-diagnostics rows are written twice

60 rows = 30 unique `(cluster, family)` pairs, inflating the calibration dataset 2x.
Not investigated.

## Deliberately NOT done, and why

* **`min_identity` is still `None`.** On one genome the separation looked clean (true
  locus 75.9% vs 20 spurious admissions at 30.6-50.0%, a 25.9-point margin). But the
  Minus-strain identities across the 23 are **25.9-43.5%**, so a global cutoff near 55%
  would destroy Minus detection entirely. Any identity floor must be role- or gene-aware.
  Setting one from n=1 is exactly the mistake `EvidenceFloor`'s original docstring warned
  against.
* **The relaxed second pass** — the curator answered the flank question ("another gene,
  not necessarily a flank") but never the always-vs-only-on-empty question.
* **Promoting any detected locus into `db/`** — curation is the curator's call.

## Recommended next steps, in order

1. **Curator decision on #1** (resolve overlapping idiomorph calls by identity?) and on
   #2's relaxed-second-pass policy. Both are biology decisions, not engineering.
2. **Implement #1**, then **re-run the 23** and confirm idiomorph accuracy goes 0/23 ->
   23/23 and that `fraction_found` rises for every genome. This is the re-run that has NOT
   happened.
3. **Implement #3** (isolate polish failures per gene) — cheap, and it is pure loss today.
4. **Re-run the 44-genus sweep** after #1 and #2 and see how many of the 15 near misses
   are recovered.
5. **Only then** consider a `min_identity`, calibrated on the full diagnostics corpus.
6. Merge `detect-targeting`. It is green and reviewed; it is unmerged only because the
   validation run raised these questions.

## Reproducing the runs

Inputs are in `$SCRATCH` (node-local, NOT backed up, and gone when the job ends):

* `$SCRATCH/zygo/fasta/*.{fna,faa}` — 67 converted genomes. **Proteome deflines MUST carry
  `contig:start-end:strand`** or the diamond fast path rejects them; that cost a failed run.
  Converter: `$SCRATCH/zygo/convert2.py`.
* `$SCRATCH/zygo/runs/<organism>/` — per-genome GFF3, companion FASTA, report, diagnostics.
* `$SCRATCH/zygo/reps.txt` — the 44 genus representatives.

Command shape:
```
matpredict detect --genome <g>.fna --proteins <g>.faa --phylum Mucoromycota \
  --out-dir <out> --evidence-diagnostics <out>/evidence_diagnostics.jsonl --emit-cds-fasta
```

This SLURM allocation had **4 CPUs / 16 GB**; the sweep ran at concurrency 3. Disk stayed
at 1% of a 3.2 TB `/scratch`. Nothing was written to `/bigdata` outside the repo.

## Process note worth keeping

A real `detect` run was launched against the working tree while a reviewer was authorised
to temporarily rename a CLI flag in it. The run died with
`'Namespace' object has no attribute 'min_hits'` and was nearly reported as a Critical bug.
It was not. **Never run a real experiment against a tree an agent may mutate** — use a git
worktree, or serialise the two.
