# Gap-lineage pilots and the holdout re-run, 2026-09-24

Runs: pilots from frozen worktree `run-a90b97d` (Pezizomycetes and some
Taphrinomycotina / Dothideomycetes / uncurated-Pezizomycotina / Dipodascomycetes
reports from `run-6d15944`, kept only where routing was `lineage` or `direct`);
holdout from `run-6d15944`. Raw results:
`results/2026-09-24_pilots/<lineage>/runs/`, `results/2026-09-24_pilots/summary.txt`,
`results/2026-09-24_holdout/`. Lists: `results/2026-09-24_pilots/lists/`,
stratified by family with >=1 genome per family, seed 20260924.

## Pilot call rates

"Called" = at least one reported locus. "Bar only" = nothing reported, but the
modelled-gene bar withheld a locus (now written to `suppressed_loci`).
Wall time is measured per genome (`wall_seconds`), one CPU slot each.

| lineage | genomes | called | bar only | neither | loci/genome | routing | median s | p90 s |
|---|---:|---:|---:|---:|---:|---|---:|---:|
| Dothideomycetes (0 records) | 100 | **90** | 10 | 0 | 0.91 | lineage 100 | 622 | 948 |
| uncurated Pezizomycotina orders | 100 | **82** | 18 | 0 | 0.88 | lineage 100 | 620 | 856 |
| Taphrinomycotina | 103 | 80 | 23 | 0 | 1.74 | lineage/direct | 13 | 69 |
| Saccharomycetes, not Saccharomycetaceae | 100 | 56 | 44 | 0 | 0.80 | phylum_fallback 100 | 516 | 900 |
| Dipodascomycetes + small classes | 100 | 48 | 51 | 1 | 0.63 | mostly phylum_fallback | 549 | 1029 |
| Orbiliomycetes | 67 | 27 | 40 | 0 | 0.42 | phylum_fallback 67 | 336 | 582 |
| Pichiomycetes | 99 | 18 | 81 | 0 | 0.37 | lineage 75 / fallback 24 | 4 | 663 |
| Pezizomycetes | 172 | 27 | 81 | 64 | 0.16 | lineage | 22 | 41 |

By order, where it matters:

* **Dothideomycetes generalise with no in-class record.** Pleosporales 35/37,
  Mycosphaerellales 11/12, Botryosphaeriales 8/8, Dothideales 8/8,
  Venturiales 5/5. Idiomorph split 46 MAT1-1 / 45 MAT1-2. Not verified against
  truth coordinates: this is a call rate, not a recall.
* **Serinales 2/75**, with 73 bar-only. Nearly every withheld locus is one
  100-260 bp MTLA1 HSP with 0 modelled genes. `MTL` has ONE record,
  D. hansenii MTLa; there is no alpha reference. Curation gap.
* **Taphrinales 5/22**, Saitoella 0/2, Neolecta 0/1, all the rest bar-only.
  Withheld clusters carry matMc with matPi or matMi in one window -- the
  expected primary-homothallic layout -- but 0-1 genes could be modelled from
  the single T. deformans reference.
* **Saccharomycodales 3/21** (Hanseniaspora; the literature reports possible
  MAT gene loss -- unverified here), Ascoideales 11/26.
* **Pezizales outside Tuberaceae** 15/160: withheld clusters hold ONE MAT
  core gene plus a Tuber flanking ORF.
* Schizosaccharomycetales 66/67, Pneumocystales 9/11, Pichiales 16/24.

Families outside their home clade are being called under `phylum_fallback`:
the Pezizomycotina `MAT` family in 37 Dipodascales and 26 Orbiliales loci,
`MATsc` in Phaffomycetales. HMG and alpha-box homology makes this plausible;
whether those calls are right is unmeasured.

## The taxonomy defect the first pilot launch exposed

The first launch (8 panels x 16 concurrent runs) routed one taxid three ways
and sent Leotiomyceta genomes to `exhaustive`. Cause and fix in `a90b97d`:
non-atomic cache writes, a 12 s retry window against NCBI's keyless 3 req/s,
and silent degradation. After the fix and a serial warm-up, all 521 pilot
reports written by `a90b97d` have `routing_error: null` and
`genetic_code_error: null`. The other 320 were kept from the first launch
because they routed `lineage` (252) or `direct` (68), so they did not degrade;
they predate the field. The one warm-up failure (taxid 45151) was an HTTP 400
from NCBI, and its genome's own lookup then succeeded. Taphrinomycotina's slowest
genome went from 817 s to 94 s.

## Holdout, re-run on current code and re-scored

`results/2026-09-24_holdout/`, scored by `holdout.score_holdout`. Reports from
`6d15944`; re-scored offline with `680435b`'s scorer
(`rescored_680435b.txt`), valid because no run degraded (78 lineage, 5 direct).

| radius | found | right idiomorph | bar only (right / wrong idiomorph) | no_reference_family |
|---|---:|---:|---:|---:|
| record | 13/16 | 12 | 2 / 1 | 1 |
| species | 12/16 | 12 | 3 / 1 | 1 |
| genus | 12/16 | 12 | 3 / 1 | 1 |
| family | 9/11 | 9 | 1 / 1 | 6 |
| order | 9/11 | 9 | 1 / 1 | 6 |

**There is no genuine miss and no `no_reference_idiomorph` at any radius.**
Every record not found was withheld by the bar at the right place. The
2026-09-23 table's conclusions change in three places:

* *S. pombe* mat1-M is not "structurally unfindable". With routing fixed it is
  found at every radius, idiomorph M, by the PM family's matMc -- and withheld,
  0 modelled genes.
* The "only two genuine misses" (*K. lactis* MATa, S288C HMRa) are bar losses.
  K. lactis: one 159 bp MATA1 HSP, 0 modelled. S288C HMRa at species/genus:
  MATA1 + MATALPHA2 over 1 kb, 1 modelled, labelled MATa.
* S288C HMRa at record radius is found but labelled **MATalpha**, a
  `wrong_idiomorph`. *S. octosporus* mat1-P is withheld and labelled M.

The withheld evidence is often thin -- a single unmodelled HSP -- so "bar
loss" means "the right place held a signal the bar judged too weak", not "a
good locus was thrown away".

## What this says about the bar

In these lineages the bar, not the search, is what stands between a genome
and a call. It was set on Pezizomycotina panels (46,647 loci; one-modelled-gene
loci produced 0 high-confidence calls), where references are close. Where the
nearest reference is distant or single -- Serinales, Taphrinales, Pezizales
outside Tuberaceae, Orbiliales -- genes localise but cannot be modelled.

Two routes, not decided here:

1. **Curation.** Add references for the lineages above; the verified
   candidates are in `2026-09-24_mat-reference-gap-literature.md`. Expected to
   raise modelled-gene counts directly. Unmeasured until added.
2. **A tier for bar-only loci.** Report them as a distinct, capped class so a
   sweep does not read them as absence. This is a curator decision: it
   reverses ruling (3) of 2026-09-22 for a subset.
