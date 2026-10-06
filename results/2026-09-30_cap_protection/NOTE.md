# Polish-cap protection for classifier-strong clusters (2026-09-30)

## Question
Should the per-family polish cap (top 6 clusters) exempt a cluster whose
unmodelled tblastn HSP fragment already scores at or above the family's gate
threshold on the idiomorph classifier? Motivating case: Mucor (Rhizomucor)
pusillus NRRL A-13674 lost a Plus call after the Umbelopsis records added
competing clusters (results/2026-09-29_umbelopsis_merge/NOTE.md).

## Data and code version
- Code/db: frozen worktree `.claude/worktrees/run-capprot` at a2fe1b4 plus an
  UNCOMMITTED patch to `_polish_rank` use in `pipeline.py` (env `MATP_CAP_EXEMPT`).
- Baseline: the a2fe1b4-equivalent runs in results/2026-09-29_umbelopsis_merge/
  (Mucoromycota_f59353c, lcg_runs, jena/runs_scaffolds, zygo23_f59353c).
- Panels: Mucoromycota 293, LCG Mucoromycotina 621, Jena 64, Zygo 23 (both
  inputs); controls: 50 Ascomycota (cap panel) + 50 Basidiomycota (<750 Mb),
  baseline vs V1. SLURM short, jobs 29300917-29300935 (jobs.txt).

## Method
A cluster is "strong" when the classifier, scoring its HSP translations, gives a
best MAT score >= the manifest gate threshold and no paralog class.
- V1: polish every strong cluster beyond the cap.
- V2: polish at most one extra strong cluster per family per genome.
- V3: rank strong clusters first within the cap (no extra polishes).

## Results
| Set | Variant | Called (base -> var) | Lost / gained / changed |
|---|---|---|---|
| Mucoromycota 293 | V1, V2, V3 | 253 -> 253 | 0 / 0 / 0 |
| LCG 621 | V1, V2, V3 | 535 -> 536 | 0 / 1 / 0 |
| Jena 64 | V1, V2, V3 | 61 -> 61 | 0 / 0 / 0 |
| Asco 50 / Basidio 50 (control) | V1 | unchanged | 0 / 0 / 0 |

- The one gain in every variant is M. pusillus NRRL A-13674: Plus, medium,
  model-typed, margin 90.5, scaffold_161:16,276-25,616 (fragment 224.8 bits).
- Strong clusters: 965 across all genomes; 954 were already inside the top 6.
  V1/V2 newly polished 11 clusters (exempt_V1.log): M. pusillus A-13674,
  Lichtheimia hyalospora NRRL 1305, Benjaminiella poitrasii, Mycotypha africana,
  Jena CBS210_80's BFD copy, and six Syncephalastrum loci at 102-104 bits. Only
  M. pusillus produced a reported call; the others stay withheld by later rules.
- Zygo 23: 23/23 locus and idiomorph, scaffolds and contigs, all variants.
- No paralog-class or weak-margin call appeared. Clean strain labels unchanged
  (no labelled genome changed).
- Controls: no change; families without a classifier are unaffected, as designed.

## Limits
- Runtime cannot be measured from these runs: the variant jobs ran on less
  loaded nodes and finished faster than the baseline, although V1/V2 only add
  work. The added work is 11 extra polishes over ~980 genomes.
- One recovered call; the evidence for benefit rests on a single genome.

## Extension (curator-approved, 2026-09-30): stress test, the 11 clusters, early-diverging

Code: same run-capprot copy; added env MATP_CAP (cap override, 0 = off) to the
uncommitted patch. 84 jobs (jobs2.txt), all COMPLETED; every genome has a report
(Mucoromycota 293, LCG 621, Jena 64 in all 9 configs; ED 290 in all 4); no FAIL lines.
Script: compare_stress.py -> compare_stress_output.txt.

### Stress test against a cap-off reference
| Set, cap | Cap-off calls lost by plain cap | V1 recovered | V2 recovered | V3 recovered |
|---|---|---|---|---|
| Mucoromycota, 3 | 1 | 1 | 1 | 1 |
| Mucoromycota, 2 | 4 | 4 | 4 | 4 |
| LCG, 3 | 3 | 3 | 3 | 3 |
| LCG, 2 | 11 | 10 | 9 | 10 |
| Jena, 3 | 0 | - | - | - |
| Jena, 2 | 1 | 1 | 1 | 1 |
| Total | 20 | 19 (95%) | 18 (90%) | 19 (95%) |

- False promotions (a call not in cap-off, a different label, or a paralog class):
  V1 0, V3 0, V2 1 -- Mucor indicus NRRL 13082 at cap 2: V2 polished the best-scoring
  extra cluster (scaffold_101) and reported Plus/medium partial_locus (margin 247.8),
  while cap-off reports Minus/high at scaffold_1916. With the true locus unpolished,
  the relaxed path reported the other cluster. V3 gave the cap-off answer.
- V3 cost: at cap 2 only, it displaced two calls that plain cap 2 kept -- Mucor
  genevensis NRRL 1821 scaffold_188 and NRRL 1412 scaffold_980, both Minus, low,
  partial_locus, margins 31.4-31.7 (the weak ~30-bit second Minus calls flagged as
  possible paralogs in results/2026-09-29_two_idiomorphs/). Still lost under V3:
  M. genevensis NRRL 1756 scaffold_619 (same kind). No V3 drop at cap 3 or 6.
- Zygo 23: 23/23 locus and idiomorph, scaffolds and contigs, at every cap and variant
  (cap-off; caps 3 and 2 with plain, V1, V2, V3).
- Newly polished clusters beyond the cap: V1 21 (cap 3) / 36 (cap 2); V2 20 / 29;
  V3 0 by design.

### The 11 extra-polished clusters (first replay, V1 at cap 6)
- Called: M. pusillus A-13674 (the gain); Mycotypha africana NW_026515730.1 (already
  called identically in the baseline).
- below_fraction_floor (0.4): Benjaminiella poitrasii (undetermined, margin 0.1);
  the BFD copy of Jena CBS210_80 (Plus 245.4; the genome keeps its Minus/high call);
  S. racemosum CBS 440.59 (Minus 34.2).
- modelled_gene_bar: Lichtheimia hyalospora NRRL 1305 (Plus 118.6; Lichtheimiaceae
  gap); Syncephalastrum Synrac1 (flank genes only, no core); four Syncephalastrum
  genomes (Minus, margin 30.4, sexP+sexM).
- None was withheld by the paralog class or the gate. The later rules still stop
  them; but the cap-2 M. indicus case shows V1/V2 can surface a wrong cluster when
  the true locus is also displaced. V3 cannot add clusters, only reorder them.

### Mortierellomycota + Kickxellomycota (290, --phylum Mucoromycota)
- 0 calls at cap 6 and cap 3, plain and V3 alike (current code withholds all of
  them). 35 strong-fragment clusters exist; V3 promoted none to a call. No false
  promotions. (Discovery-only lineages; any call would be unverified.)

### Runtime
Wall-clock totals are not comparable across these jobs: node load varied (for
example Mucoromycota cap 3: plain 49,096 s, V1 43,193 s, V2 71,560 s), and V1/V2
only add work. The one like-for-like pair, ED cap 6 Mortierellomycota, was 10,801 s
plain vs 10,932 s V3 (+1%). V3 adds the fragment scoring only.

## Recommendation
Adopt V3 (rank classifier-strong clusters first, inside the cap). It recovered 19 of
20 displaced cap-off calls with 0 false promotions, adds no polishing, kept Zygo
23/23, and changed nothing in the controls or the early-diverging lineages. Its one
cost appears only at cap 2: two weak low-confidence partial Minus calls of the
likely-paralog kind. Do not adopt V2 (1 false promotion). Implement test-first in
PR #9 if the curator agrees.

## Files
compare.py, compare_output.txt, summary.tsv, exempt_V{1,2,3}.log, jobs.txt,
lists/, V1/ V2/ V3/ ctrl_base/ ctrl_V1/ run directories; extension:
compare_stress.py, compare_stress_output.txt, jobs2.txt, stress/ (9 configs,
exempt_*.log), ed/ (4 configs).
