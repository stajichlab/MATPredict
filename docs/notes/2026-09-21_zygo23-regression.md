# The 23-genome Zygo regression, after the 2026-09-21 changes

First genome run against any of the `detect-scope-and-gaps` work. Purpose: prove
the nine unvalidated changes did not break the validated Mucoromycota result.

Reports archived at `results/2026-09-21_zygo23/reports.tar.zst` (38 KB). Genome
FASTAs were node-local scratch and are gone with the job; they are reproducible
from `testset/Zygo/query/*.gbk` with the converter described below.

## Result — no regression

| | before (2026-09-20) | after |
|---|---|---|
| locus on truth scaffold | 23/23 | **23/23** |
| idiomorph correct | 23/23 | **23/23** |
| total loci reported | 31 | **31** |
| `locus_class` | 31 `mat_locus` | **23 `mat_locus` + 8 `partial_locus`** |
| confidence | 18 high / 13 medium | **18 high / 13 medium** |
| `detection_pass` | all strict | all strict |
| mean wall time | 41.8 s | **18.5 s** (min 13.5, max 23.9) |

23 genomes, 0 errors, 425 s of CPU total at concurrency 3 on 4 cores.

## `partial_locus` landed on exactly the intended calls

All 8 `partial_locus` calls are the same element, and nothing else was
relabelled:

```
gene set : sexP|sexM|glrA          x8  (no other combination)
span     : 6,461 - 6,487 bp
conf     : medium                  x8
idiomorph: undetermined            x8
```

Actinomucor A-23671, Circinomucor circinelloides NRRL 22899, Mucor alternans
A-15142, and Mucor sp. A-14906 / A-21230 / A-21232 / A-21236 / A-25970 — one per
genome, none on a truth scaffold.

This is the conserved HMG-paralog region that
`2026-09-20_mucoromycota-scale-testing.md` (Addendum 4) identified as the
precision cost of taking `btbA` out of the denominator: 3 found of 6 expected =
exactly 0.500, surviving because the floor test is `<` and not `<=`. The
curator's ruling of 2026-09-21 kept `<` and labelled the tie instead. The label
now fires on precisely those 8 and on nothing else.

**The 23 true loci were untouched** — still `mat_locus`, still 18 high / 5
medium, still the correct idiomorph in every case.

## What this does NOT establish

* **It is a regression guard, not a new capability test.** Every change under
  test here is Mucoromycota-neutral by construction: the run uses
  `--phylum Mucoromycota`, which routes one family, so the order-level
  Basidiomycota scope, the per-locus Basidiomycota gaps and the routing-aware
  gap floor are all inert on this set. They remain unvalidated.
* **The speedup is not attributable.** 41.8 s -> 18.5 s is a real measured
  difference but the two runs used different hardware and this one ran three
  genomes at a time on four cores. The two caches are the likely cause and the
  component measurements support that (19.4x on window extraction, 201x on the
  curated-DB parse), but this run does not isolate them. A before/after on
  identical hardware was not done.
* **Nothing about `partial_locus` outside Mucoromycota.** 4 of the 5 families
  that can be `mat_locus` at exactly 0.500 are Ascomycota, and no Ascomycota
  genome has been run.
* **Nothing about the relaxed pass.** All 31 calls came from the strict pass, so
  the "relaxed never emits `mat_locus`" rule was never exercised here.

## Reproducing

`testset/Zygo/query/<organism>.gbk` symlinks into the ZyGoLife LCG annotation
tree; all 23 ground-truth organisms resolve. Convert each to `.fna` + `.faa`,
with proteome deflines carrying `contig:start-end:strand` — without that the
diamond fast path raises `ProteomeDeflineError`, which the previous handoff
records as having cost a whole run. Validate the deflines against
`search._parse_proteome_location` before converting the whole set.

```
matpredict detect --genome <g>.fna --proteins <g>.faa --phylum Mucoromycota \
  --out-dir <out> --evidence-diagnostics <out>/evidence_diagnostics.jsonl
```
