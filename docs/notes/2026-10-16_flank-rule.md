# Flank-only helixer rule (2026-10-16)

Code: `src/MATPredict/detect/helixer_rule.py`; evaluation: `scripts/panel_helixer_metric.py score` (prints a
`flank rule:` line per genome). The rule is a function that returns a refined model beside the homology
model; nothing calls it from `detect` yet, so no call or report changes.

Rule: a `flanking_*` gene takes helixer-lite's ends when a helixer model matches it (same strand, 50%
reciprocal overlap) and is not an explicitly incomplete ORF. Everything else, including every `core_MAT`
gene and every gene with no helixer model, keeps its homology model.

## Panel result (eight genomes, same helixer tables as the 2026-10-15 scoring; descriptive)

| genome | flank genes | homology within 30 bp | rule within 30 bp | gained | lost |
|---|---|---|---|---|---|
| phycomyces | 2 | 0 | 0 | 0 | 0 |
| rhizopus | 2 | 0 | 2 | 2 | 0 |
| mucor_sp1006 | 4 | 2 | 3 | 2 | 1 |
| fusarium | 3 | 0 | 3 | 3 | 0 |
| aspergillus | 3 | 3 | 1 | 0 | 2 |
| ustilago | 1 | 1 | 1 | 0 | 0 |
| mucor_mucedo | 3 | 0 | 0 | 0 | 0 |
| saccharomyces_r64 | 0 | - | - | - | - |
| **total** | **18** | **6** | **10** | **7** | **3** |

Core MAT agreement is unchanged by construction (15/17 stays 15/17). The adoption test from the metric
section, beating the baseline on flank without losing core MAT agreement, is met.

Caveats: the rule was chosen after seeing these counts, so this is an in-sample check, not a held-out one;
n is 18 flank genes; NCBI is a prediction; the helixer tables predate the `complete_orf` column, so the
ORF condition was not exercised (unknown passes). Aspergillus is the counter-example (2 of 3 flank genes
lose agreement) and is worth reading before merging. A held-out panel is the right next test.
