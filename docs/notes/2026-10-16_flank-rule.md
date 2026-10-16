# Flank-only helixer rule (2026-10-16)

Code: `src/MATPredict/detect/helixer_rule.py`; evaluation: `scripts/panel_helixer_metric.py score` (prints a
`flank rule:` line per genome). The rule is a function that returns a refined model beside the homology
model; nothing calls it from `detect` yet, so no call or report changes.

Rule: a `flanking_*` gene is **extended** to helixer-lite's ends when a helixer model matches it (same strand,
50% reciprocal overlap) and is not an explicitly incomplete ORF. Each end moves on its own and only outward;
helixer never shortens a homology model. Extend-only is the curator's choice until a larger gene set
(Aspergillus, Neurospora) is tested; a first version that also shortened scored 10/18 and lost 2 of 3 on
Aspergillus. Everything else, including every `core_MAT`
gene and every gene with no helixer model, keeps its homology model.

## Panel result (eight genomes, same helixer tables as the 2026-10-15 scoring; descriptive)

| genome | flank genes | homology within 30 bp | rule within 30 bp | gained | lost |
|---|---|---|---|---|---|
| phycomyces | 2 | 0 | 0 | 0 | 0 |
| rhizopus | 2 | 0 | 2 | 2 | 0 |
| mucor_sp1006 | 4 | 2 | 3 | 2 | 1 |
| fusarium | 3 | 0 | 3 | 3 | 0 |
| aspergillus | 3 | 3 | 2 | 0 | 1 |
| ustilago | 1 | 1 | 1 | 0 | 0 |
| mucor_mucedo | 3 | 0 | 0 | 0 | 0 |
| saccharomyces_r64 | 0 | - | - | - | - |
| **total** | **18** | **6** | **11** | **7** | **2** |

Core MAT agreement is unchanged by construction (15/17 stays 15/17). The adoption test from the metric
section, beating the baseline on flank without losing core MAT agreement, is met.

Caveats: the rule was chosen after seeing these counts, so this is an in-sample check, not a held-out one;
n is 18 flank genes; NCBI is a prediction; the helixer tables predate the `complete_orf` column, so the
ORF condition was not exercised (unknown passes). Aspergillus is still the counter-example (1 of 3 flank genes
loses agreement) and Mucor sp. 1006Ph loses 1; both are worth reading before merging. A held-out panel is the right next test.

## Query check (added after the two extend-only losses)

The two genes extend-only lost (Aspergillus APN2, Mucor sp. 1006Ph rnhA) were both cases where helixer's
extended 5′ end replaces the N-terminus with a different start (a frame-shifted start codon plus an extra
intron in APN2; an 11 bp micro-exon 147 bp upstream in rnhA) while the rest of the protein is identical to
NCBI's. The rule now also requires helixer's protein to align the curated query protein better than the
homology protein does (BLOSUM62 local score gain of at least 10, `MIN_SCORE_GAIN`, fixed before this run).

`scripts/panel_rule_check.py`, same eight genomes, flank genes matching NCBI, ends within 30 bp:

| | within 30 bp | changed | gained | lost |
|---|---|---|---|---|
| homology only | 6/18 | – | – | – |
| extend-only, unverified | 11/18 | 17 | 7 | 2 |
| extend-only, query-checked | 10/18 | 5 | 4 | 0 |

The query check removes both losses and keeps 4 of 7 gains. It rejects 12 extensions; most of those
(Rhizopus rnhA and glrA, Fusarium COX13, Mucor mucedo) were helixer ends within 30 bp of NCBI that the query
cannot confirm because they only add the stop codon or UTR-side residues, so the query check is conservative
about extensions that are probably right. Aspergillus's query is its own genome's curated record, so its
agreement is partly self-consistent.
