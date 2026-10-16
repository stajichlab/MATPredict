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

## Pair rule decision (2026-10-16)

The profile-bitscore rule (pyhmmer, each gene against its own idiomorph's profile) stays as the sexM/sexP pair
rule for families with a classifier, replacing the 2026-09-26 model-score rule. Basis: the old and new rules
agree on 6 of 8 pair decisions in five Mucoromycota genomes; the 2 that differ (Mucor sp. 1006Ph and
*M. circinelloides*, identical numbers: profile 76 against 57 for sexM, model score 129 against 101 for sexP) are
low-identity pairs. Families without a classifier keep the model-score rule.

## Held-out test (pre-registered 2026-10-16, before any of these genomes is run)

Genomes (RefSeq, NCBI-annotated, none in the 8-genome panel): *Neurospora crassa* OR74A `GCF_000182925.2_NC12`,
*Aspergillus nidulans* FGSC A4 `GCF_000011425.1_ASM1142v1`, *Aspergillus niger* CBS 513.88 `GCF_000002855.4_ASM285v2`.
The curated database holds records from the same strains (`5141_74-ors-a`, `162425_fgsc-a4`); each is excluded
from the search with `--exclude-records` so the query protein for every gene comes from another strain or species.

Pipeline unchanged: detect (phylum Ascomycota), helixer-lite on the locus ±5 kb windows (fungi lineage), the
rule with the query check, scored against the NCBI annotation with the panel metric (50% reciprocal match, 30 bp
tolerance, strata core MAT / flank).

Pass criteria, fixed now: (1) core MAT agreement is unchanged by the rule (true by construction; reported);
(2) on flank genes the rule gains more genes than it loses (gained minus lost is above zero) and loses none that
homology had right and helixer's extension moves more than 30 bp from NCBI without query support (the query check
should prevent this; any loss is read gene by gene); (3) no more than one extension in three is rejected by the query
check on genes where the extension ends within 30 bp of NCBI. With three genomes the counts are descriptive; no
significance is claimed. Performance is recorded for each stage: detect wall time and peak memory, helixer-lite
window wall time and peak memory, and the rule check wall time.

## Held-out result (2026-10-16; criteria above were fixed before the run)

*N. crassa* OR74A, *A. nidulans* FGSC A4 and *A. niger* CBS 513.88, with the strain's own curated records excluded from
the search. `scripts/panel_rule_check.py`, scores in `scripts/panel_helixer_metric.py`.

**Metric correction made after the first held-out look.** In these annotations the NCBI gene feature includes
UTRs (N. crassa: 256 to 568 bp 5′ and up to 429 bp 3′ beyond the CDS), while the pipeline and helixer models are
coding-only, so the gene-span comparison penalised correct models. Both scripts now take `--ncbi-extent cds`
(the span of the gene's CDS features). It applies equally to homology, helixer and the rule. The panel was rescored
the same way. Both versions are shown; the CDS extent is the like-for-like comparison and is the one to read.

Flank genes within 30 bp of NCBI (CDS extent):

| set | n | homology | rule, unverified | rule, query-checked | gained | lost |
|---|---|---|---|---|---|---|
| panel (8 genomes) | 18 | 8 | 15 | 13 | 5 | 0 |
| held-out | 9 | 2 | 6 | 6 | 4 | 0 |
| both | 27 | 10 | 21 | 19 | 9 | 0 |

(Unverified, both sets: 13 gained, 2 lost; the two losses are on the panel and are removed by the query check.)
With the gene extent the held-out set is 1/9 → 5/9 (query-checked: 4 gained, 0 lost) and the panel 6/18 → 10/18.

Pre-registered criteria: (1) core MAT unchanged: yes, by construction. (2) gained minus lost above zero and nothing
lost: 4 gained, 0 lost on the held-out genes. (3) at most one in three within-30-bp extensions rejected by the
query check: 0 of 4 rejected. All three met.

**Core MAT is a different story with the self-records out.** Held-out core genes (CDS extent): homology within 30 bp
1 of 3 matched (nidulans 1/2, niger 0/1; N. crassa MAT1-1-1, -1-2 and -1-3 were only short tblastn fragments of
105, 321 and 117 bp, none matching an NCBI gene), helixer matched 3 of 3 and within 30 bp 2 of 3 (nidulans 2/2).
On the panel homology was better on core (15 of 17 against 13 of 15). The panel includes genomes whose own record is
in the database (A1163), which flatters homology. The rule still leaves core on homology (curator's decision), but
this is evidence to revisit with a held-out core test, n = 3 here.

**Performance** (one genome each, wall time, peak memory):

| genome | detect | helixer-lite on the locus windows | rule check |
|---|---|---|---|
| N. crassa | 496 s, 0.31 GB | 25 s, 1.0 GB | 16 s |
| A. nidulans | 316 s, 0.31 GB | 75 s, 0.96 GB | 23 s |
| A. niger | 600 s, 0.31 GB | 139 s, 0.83 GB | 30 s |

Detect is single-threaded (CPU time equals wall time). helixer-lite on the windows adds 5 to 25% to detect (CPU time
18 to 26 s; most of the wall time is model loading), against 34 min and 3.9 GB for the whole Phycomyces genome. The
rule check includes a second helixer run (no cache) and takes 16 to 30 s.
