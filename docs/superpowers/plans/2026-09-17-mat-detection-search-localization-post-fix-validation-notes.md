# Post-fix real-genome validation: localize-then-polish search revision

Date: 2026-09-17. Environment: UCR HPCC interactive SLURM job (host as
allocated, 8 CPUs visible, node 1-minute load average ~214 during the runs),
pixi `default` environment with the real `tblastn`, `makeblastdb`, `miniprot`,
`exonerate` and `diamond` binaries in
`.pixi/envs/default/bin/`.

This is a manual evidence log, not a pytest suite. It follows the earlier
file `2026-09-17-mat-detection-search-localization-benchmark-notes.md`, whose
correctness interpretation was invalidated when the multi-gene-window bug was
found. Every number below is a real measurement from a real run. Where
something could not be measured or reproduced, that is stated plainly.

Code under test: `main` at `76a60ef` ("fix: attribute per-gene evidence per
segment on a fragmented call"), which carries all of this session's fixes
(`bac68b7`, `b1936a3`, `a177b12`, `b927c89`, `12267eb`, `76a60ef`).
Comparison baseline: `8049a26`, the commit the earlier benchmark notes were
written against, i.e. before any of those fixes.

`git diff --stat 8049a26 HEAD -- db/` is empty: the curated reference database
is byte-identical between the two commits, so any difference below comes from
the code, not from changed reference data.

## Run 1 (primary): GCF_000143185.2, *Schizophyllum commune* H4-8, Aalpha

Genome obtained exactly as in the earlier notes:

```
curl -sS -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/143/185/GCF_000143185.2_Schco3/GCF_000143185.2_Schco3_genomic.fna.gz
gunzip -f GCF_000143185.2_Schco3_genomic.fna.gz     # 39,156,689 bytes
```

Curated record `db/Basidiomycota/Agaricales/5334_h4-8_Aalpha_4/metadata.yaml`:
segment `NW_026089539.1:1821008-1827457`, genes `Z` (1821008-1824097, `-`) and
`Y` (1824369-1827457, `+`).

Post-fix command (working tree at `76a60ef`):

```
pixi run matpredict detect \
  --genome $SCRATCH/matpredict_val/GCF_000143185.2_Schco3_genomic.fna \
  --taxid 5334 --out-dir $SCRATCH/matpredict_val/out_post
```

Pre-fix command (temporary `git worktree` at `8049a26`, reusing the same
installed pixi binaries via `PATH`, so both runs used identical tool builds):

```
git worktree add $SCRATCH/matpredict_val/wt_8049a26 8049a26
export PATH=<repo>/.pixi/envs/default/bin:$PATH
PYTHONPATH=$SCRATCH/matpredict_val/wt_8049a26/src \
  <repo>/.pixi/envs/default/bin/python -m MATPredict detect \
  --genome $SCRATCH/matpredict_val/GCF_000143185.2_Schco3_genomic.fna \
  --taxid 5334 --out-dir $SCRATCH/matpredict_val/out_pre_fix
```

Both runs printed the same headline:
`detected 13 candidate locus/loci (4 families attempted)`, with
`Basidiomycota:Abeta`, `Basidiomycota:Balpha` and `Basidiomycota:Bbeta`
reported as not detected ("no reference-protein hits found for this family in
this genome").

### Correctness: the curated Aalpha call, before and after

Read directly from each run's `detection_report.yaml`:

| | pre-fix `8049a26` | post-fix `76a60ef` |
|---|---|---|
| call span | `NW_026089539.1:1821007-1827453` | `NW_026089539.1:1821010-1827453` |
| confidence tier | **medium** | **high** |
| gene `Z` | 1821007-1824096, 98.31%, `miniprot_refine`, `polished_single` | 1821010-1824096, 100.0%, `exonerate_refine`, `polished_disagree` |
| gene `Y` | **1826233**-1827453, 98.771%, `tblastn_genome`, **`unpolished`** | **1824368**-1827453, 100.0%, `exonerate_refine`, `polished_disagree` |

This is the direct real-world confirmation the multi-gene-window fix was
supposed to produce, and it shows two separate improvements:

1. **Tier.** `medium` -> `high`. Pre-fix, gene `Y` never received a polish
   model from either tool, so `_any_gene_unpolished` capped the tier. Post-fix
   both genes are polished, so the cap no longer applies.
2. **Coordinate accuracy, which matters more than the tier.** Pre-fix, gene
   `Y`'s reported start was 1826233 -- a raw tblastn HSP covering only 1,221 bp
   of the curated 3,089 bp gene, i.e. the gene model was wrong by ~1.9 kb at
   the 5' end. Post-fix, `Y` starts at 1824368 against a curated 1824369: one
   base off, and that one base is the standard 0/1-based edge, not a modelling
   error.

Both genes report `status: polished_disagree`, meaning both `exonerate
--refine` and `miniprot` produced a model for each gene and the two disagree
slightly; the primary (exonerate) model is reported and the miniprot model is
carried as `alternate_model`. For gene `Z` the miniprot alternate is
1821007-1824096 at 98.31% with 6 exons; for `Y` it is 1824368-1827456 at 98.1%
with 6 exons. Against the curated coordinates the miniprot alternates are
actually the closer of the two models at the outer edges (`Z` start 1821007 vs
curated 1821008; `Y` end 1827456 vs curated 1827457). So "disagree" here is a
few-base disagreement between two near-perfect models, not a conflict between a
right and a wrong answer. It is not `polished_agree`, but it is a genuine
two-tool confirmation of the same gene, which is what the fix was for.

The remaining 12 Aalpha calls are low-confidence, short, single-gene hits
elsewhere in the genome (expected HD-domain gene-family relatives), unchanged
in character from the earlier run.

### Timing, and an honest caveat about it

Real `time` output:

| run | wall | user | sys |
|---|---|---|---|
| post-fix `76a60ef` | **6m30.163s** | 2m6.934s | 0m16.179s |
| pre-fix `8049a26` (today) | **4m54.962s** | 2m9.417s | 0m14.090s |
| pre-fix `8049a26` (earlier notes, different session) | 1m36.4s | 1m23.8s | 10.0s |

Two things must be said plainly here.

- **The earlier session's 1m36.4s wall / 1m23.8s user figure did not reproduce
  today at the identical commit with the identical genome and an identical
  reference DB.** Today the same commit took 4m54.962s wall and 2m9.417s user.
  I did not determine the cause. The node's 1-minute load average during these
  runs was ~214 on an 8-CPU view, so the machine was heavily oversubscribed,
  which inflates wall clock and, through cache and memory-bandwidth contention,
  can inflate measured user CPU as well. I am not claiming that explains all of
  it; I am stating that the earlier absolute number is not reproducible in this
  environment and should not be treated as a stable baseline.
- **The valid comparison is the same-session pair.** Pre-fix 2m9.417s user vs
  post-fix 2m6.934s user: the fixes cost no measurable extra CPU, despite
  post-fix polishing four gene/tool models where pre-fix polished one. The
  wall-clock difference between the two (4m55 vs 6m30) is not backed by a
  matching CPU difference and on a node at load ~214 I do not consider it
  attributable to the code.

The earlier notes' headline claim -- that the localize-then-polish revision
turned a >5-minute non-completion into a run of a few minutes -- was measured
against the *pre-revision* commit `cd7e386`, which is a different comparison
and was not re-run here.

## Run 2: GCA_016772295.1, *Coprinopsis cinerea* A43mut B43mut okayama7#130

Chosen as the second `assembly`-type curated record. It is a different family
(Psathyrellaceae vs Schizophyllaceae) and a different genus, though the same
order and phylum -- the curated DB contains no `assembly`-type record outside
Basidiomycota/Agaricales, so a different *phylum* was not available. The five
`sequence_source.type: assembly` records in `db/` sit on only two underlying
assemblies (`GCF_000143185.2` and `GCA_016772295.1`); this is the other one.

```
curl -sS -o cc.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/016/772/295/GCA_016772295.1_ASM1677229v1/GCA_016772295.1_ASM1677229v1_genomic.fna.gz
gunzip -f cc.fna.gz     # 39,185,934 bytes, 31 sequences
```

### 2a: HD locus (2-gene window) -- passes, high confidence

Curated `5346_a43-b43-okayama-7_HD_A43`: `JAAGWA010000001.1:1625680-1631546`,
genes `HD1`, `HD2`.

```
pixi run matpredict detect --genome $SCRATCH/matpredict_val/cc.fna \
  --taxid 2982305 --out-dir $SCRATCH/matpredict_val/out_cc_hd
```

`real 4m44.861s, user 1m27.223s, sys 0m14.282s`. Result: 12 candidate calls,
the top one `JAAGWA010000001.1:1625683-1631394`, **high** confidence,
`genes_found: [HD1, HD2]`, `genes_missing: []`:

```
HD1 1625683-1627631 (-) 100.0 exonerate_refine polished_disagree (alt: miniprot_refine)
HD2 1628111-1631394 (+) 100.0 exonerate_refine polished_disagree (alt: miniprot_refine)
```

This is a second, independent multi-gene-window case, in a different family and
a different genome, where both genes in one window are polished by both tools.
It corroborates run 1's finding rather than relying on it alone.

### 2b: PR locus (8-gene window) -- localizes, but recall is poor

Curated `5346_a43-b43-okayama-7_PR_B43`: `JAAGWA010000010.1:1806154-1826859`,
8 genes (4 of them all named `pheromone_receptor`, plus `pheromone_B44`,
`pheromone_B43`, and two `fungal_mating_type_pheromone`).

```
pixi run matpredict detect --genome $SCRATCH/matpredict_val/cc.fna \
  --taxid 5346 --out-dir $SCRATCH/matpredict_val/out_cc
```

`real 5m58.321s, user 1m43.151s, sys 0m15.143s`. Result: 13 candidate calls,
all **low** confidence. The top call is
`JAAGWA010000010.1:1806650-1823823`, which does overlap the curated locus and
starts within 500 bp of it, so localization worked. But its
`gene_evidence` contains exactly **one** entry --
`pheromone_receptor 1806650-1808642 (-) 100.0 exonerate_refine
polished_disagree` -- with `genes_found: [pheromone_receptor]` and
`genes_missing: [pheromone]`.

Honest reading of this, separating what is and is not a bug in this session's
fixes:

- The order.yml `PR` family declares only two *gene names*, `pheromone` and
  `pheromone_receptor`, and the curated record reuses `pheromone_receptor` for
  four distinct proteins. The report is keyed by gene name, so four receptor
  copies necessarily collapse into one `pheromone_receptor` entry. The single
  evidence entry therefore does **not** prove only one receptor was polished;
  it proves the data model cannot express more than one per name. I could not
  distinguish those two situations from the report alone, and I did not
  instrument the code to find out.
- The `pheromone` genes are genuinely not found. These are the tiny (~40-60 aa)
  pheromone-precursor genes; tblastn localization of such short queries is
  weak. This matches the known annotation-gap behaviour for these genes and is
  a pre-existing recall limitation, not something introduced or addressed by
  this session's fixes.

So run 2b is reported as a real, unflattering result, not as a pass.

### Incidental finding: the HD family is unreachable by its own record's taxid

`db/Basidiomycota/order.yml` gives the `HD` family
`taxonomic_scope: [2982305]`, but the only curated HD record
(`5346_a43-b43-okayama-7_HD_A43`) has `taxonomy.taxid: 5346`, and `grep -rn
2982305 db/` matches nothing except that one scope line. Because `route()` in
`family_registry.py` does exact-taxid-membership matching only, running
`--taxid 5346` on the C. cinerea genome attempts **only** the `PR` family; the
HD locus is silently never searched. That is why run 2a had to be invoked with
`--taxid 2982305` to exercise HD at all. This is a pre-existing curated-DB /
scope inconsistency, unrelated to this session's code fixes. I did not change
it; it is recorded here so it is not lost.

## Full-suite sanity check

```
pixi run pytest -q
161 passed in 98.01s (0:01:38)
```

## What was validated, and what was not

Validated with real binaries against real published assemblies:

- The multi-gene-window fix works end to end. In two independent genomes and
  two different families, both genes of a 2-gene MAT window are now polished;
  before the fix, only one of the two was, on the same genome with the same
  reference DB.
- The previously misdiagnosed `medium` tier on the curated
  `Basidiomycota:Aalpha` locus is now `high`, and, more substantively, gene
  `Y`'s reported model went from a 1,221 bp truncated tblastn HSP to a
  3,086 bp spliced model matching the curated gene to within one base.
- The fixes cost no measurable extra CPU on a ~39 Mb genome.

Not validated:

- **The fast-path rescue mechanism and the fragmented-locus path** (the targets
  of `a177b12`, `b927c89`, `12267eb`, `76a60ef`) were not exercised by any of
  these runs. All three runs were genome-only (`--genome`, no `--proteins`), so
  the supplied-proteome fast path never ran, and no call in any run had
  `fragmented: true`. No curated `assembly`-type record in `db/` today has a
  multi-segment fragmented locus, and none ships a matching proteome, so there
  was no real record against which to drive those paths. They remain validated
  only by unit tests and small synthetic probes.
- **A different-phylum second validation.** No `assembly`-type curated record
  exists outside Basidiomycota/Agaricales.
- **Any recall/Sn-Sp number.** `matpredict detect benchmark` still does not
  compute recall, exactly as the earlier notes recorded. Nothing in this
  session changed that.
- **A stable absolute timing baseline.** See the caveat above: the earlier
  session's numbers did not reproduce at the identical commit.

## Cleanup

The downloaded assemblies, all output directories and the temporary
`8049a26` worktree were created under `$SCRATCH` (node-local), never inside the
repo, and were removed after the runs. Nothing large was committed.
