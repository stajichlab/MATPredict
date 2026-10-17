# Remnant status for genes: decision paper (2026-10-17)

Branch `remnant-status` from `origin/main` 5b9a9e0 (v0.6.3). Analysis and notes only: nothing in `src/`, `db/` or the schema was changed.
Context: the curator accepted the Botrytis SAS405 record (`40559_sas405_MAT_MAT1-2`, PR #95) with its truncated MAT1-1-1 remnant stored as
gene 2, `present: false`, no coordinates, no protein, and decided to revisit a schema `remnant` status after the coverage rule merged.

## 1. Inventory

The database has 168 records (`db/*/*/*/metadata.yaml`). A grep for remnant|relic|pseudogene|degenerate|truncated|partial hits 90+ files, almost all
because of "partial" (completeness) and "truncated" (annotation or sequence-boundary truncation). The nine records the curator listed, plus the
`present: false` genes and a second grep for truncat|fragment|frameshift|stop, give this classification.

### 1a. True remnants, relics or pseudogenes of an idiomorph's gene (9 records)

| record | gene | stored as | what it is | protein in `proteins.faa` (a search reference) |
|---|---|---|---|---|
| 40559_sas405_MAT_MAT1-2 | MAT1-1-1, 5'-truncated, ~nt 3267-3960, residues 140-353 of 353 (61%), alpha box missing | `present: false`, no coords | cross-idiomorph remnant in a MAT1-2 locus | none |
| 28548_06cwm-g22_MAT_combined (S. trifoliorum, fused L allele) | MAT1-1-1, split into 5' and 3' fragments by MAT1-2-1/-4 (Xu 2016) | `present: false` | gene exists only as two fragments | none |
| 1708542_exf-10342_wallMAT_v2 (Wallemia) | HMG, 43-aa box fragment, unannotated | `present: false` | truncated HMG in the inverted locus version (Gostincar 2019) | none |
| 746128_a1163_MAT_MAT1-1 | MAT1-2-1, degenerate 3' remnant (HMG-box 5' part absent or frameshifted) | `present: true`, role core_MAT, EDP53119.1 | cross-idiomorph relic in a MAT1-1 locus, deliberately kept "so it is not lost" | 250 aa (Af293's full MAT1-2-1 is 322) |
| 230073_uamh-11059_MAT_MAT1-2 | MAT1-1-1, alpha-box domain deleted (Tsui 2013) | `present: true`, AGH03258.1 | cross-idiomorph remnant, "almost certainly non-functional" | 224 aa (sas56 full MAT1-1-1 is 353) |
| 578113_sxlc146_MAT_MAT1-2 | MAT1-1-2, "co-located remnant from the opposite idiomorph" | `present: true`, AQT26493.1 | cross-idiomorph remnant per the curator's note, but the stored protein is 391 aa; whether it is truncated was NOT checked | 391 aa |
| 107463_unknown-1_MAT_MAT1-2 (Xanthoria) | MAT1-2-1, CDS "truncated protein, interrupted by stop codon at beginning of HMG box" | `present: true`, CAI59774.2 | same-idiomorph disabled gene (stop-truncated), not a cross-idiomorph remnant | 104 aa |
| 101119_atcc11807_MAT_Plus and _MAT_Minus (Mucor, homothallic) | rnhA (Plus) and glrA (Minus) annotated as pseudogenes, no translation | no gene entry at all, prose only | true pseudogenes of FLANKING (variable/conserved) genes, not core | none |

### 1b. Mentioned in prose, not recorded (1 record)

- 27339_unknown-1_MAT_MAT1-2 (Lecanorales): a deposit `misc_feature complement(880..1405) "mat1-1-1 like region"` inside the MAT1-2 idiomorph, deliberately
  not recorded as a gene (no CDS, no protein). A possible cross-idiomorph remnant that nobody has verified.

### 1c. Other uses of the words (3 records, no action)

- 1873960_unknown-1 and -2 (Mycosphaerellales): "degenerate PCR" is a method.
- 5270_521_aLocus_a1 (Ustilago maydis): "remnant pheromone-gene sequence in a2" is in the evidence text for the a1 vs a2 comparison. There is no a2 record in
  the database, so nothing is recorded; it is a possible FUTURE remnant record (a2 carries a remnant of the a1 pheromone gene).
- Not hits but related: 5180_44ba1 describes the Inv+ alleles' truncated MAT1-1-1 (not recorded, the record is the intact Inv- allele); 498019_b11243 and
  4952_w29 are annotation or sequence-boundary truncations (completeness `partial`), not biological remnants.

### 1d. All `present: false` genes (10 records, 20 gene entries): three different meanings

`present: false` is used for (i) true absence: S. trifoliorum 06CWM-G27 MAT1-2-1/-10 (a demonstrated deletion), Mucorales tptA/algA/glrA/btbA (4 records),
Umbelopsidales rnhA (2 records), Teloschistales 115814 SLA2/APN2; (ii) a remnant or fragment: sas405, 06cwm-g22, Wallemia HMG (3 core genes); and
(iii) an unannotated, never-found gene (Wallemia SXI1: absent per Gostincar 2019). Core genes: 6 entries, 4 records (3 are remnants, 1 is a true deletion
plus 2 genes of the same deletion). Flanking genes: 14 entries. Only the three remnant entries share a meaning the curator wants to distinguish.

Counts across the database (genes by role and present): core_MAT present 334, absent 6; flanking_conserved present 97, absent 7; flanking_variable present 44,
absent 7.

### 1e. Other data sets

| source | what | numbers | verified here? |
|---|---|---|---|
| Botrytis genomes (`docs/notes/2026-10-17_botrytis-curation.md`, `results/2026-10-17_partial_coverage_vote/botrytis_*`) | MAT1-1-1 remnant in MAT1-2 genomes, 61-62% coverage, scores 410-421 vs MAT1-2-1 222-331 bits | 34 genomes in the evaluation (17 stated MAT1-2 + T4; 15 + B05.10 MAT1-1); 17 flip MAT1-1 to MAT1-2 with the rule on; 31 Botrytis notes in the 243-genome campaign | the note and per-genome tables, read; not re-run |
| Botrytis MAT1-2-1 relic in MAT1-1 genomes | note cites Sl8 (GCA_015148665.1) at 16% of reference, 87 bits | brief says ~11%; **I could not find an ~11% figure**. `results/2026-10-17_botrytis_mislabels` and `docs/notes/2026-10-17_botrytis-mislabels.md` are not in the repo or on `main` | NOT verified |
| A. fumigatus (`analysis/2026-10-07_afum-reads-and-assembly.md`, `analysis/figures/afum_hybrid_loci.tsv`) | degenerate MAT1-2-1 remnant at the same place, 117 bp from MAT1-1-1, in all MAT1-1 assemblies; only a 36% 3' stretch in MAT1-2 assemblies | 160 of 160 MAT1-1 assemblies carry it (coverage 0.964-1.0, median identity 99.9%); 128 MAT1-2 assemblies carry the 36% stretch; 168 MAT1-1 loci were labelled `homothallic_candidate` on 2026-10-07 because the record stores the remnant as a core gene; hybrid loci TSV has 20 rows (IFM_59359, IFM_61407: MAT1-1-1 + remnant on a 2.4-kb contig, no flanks) | read; the 168 figure predates the full-length-model condition (see 2) and was NOT re-run on current main |
| Morchella note (2026-09-25) | the homothallic rule simulated without the full-length condition would label 113 loci homothallic in heterothallic-rich panels (MAT1-1-3 vs MAT1-2-1 cross-matches; truncated MAT1-1-1 remnants in MAT1-2 idiomorphs); with it, 0 | 113 loci | read |
| Tremellales note (2026-09-21) | "Degenerate remnant, or two paralogs?" a 38% and 27% pair 70 kb apart (Figure C) and two ~32% HMG fragments (Figure G): these are detection-side partials, not curated cases | 2 figures | read; no curated remnant in Tremellales |
| Campaign scan (`results/2026-10-17_partial_coverage_vote/scan/`) | 2,015 contested loci (1,962 genomes) from 19,400 v0.6.0 Ascomycota reports; labels MAT1-1 1,681, MAT1-2 334; classes mat_locus 1,737, homothallic_candidate 201, idiomorph_gene_only 77 | 54 loci at PARTIAL < 80 and FULL >= 90 (53 are MAT1-1 calls with a partial MAT1-1-1 and full MAT1-2-1; 1 is MAT1-2); partial coverage 60-70% in 41 of the 54, below 50% in 8; both genes full in 678; the rule's flips: 46 in the campaign (Botrytis 31, Ophiostoma 5, Sclerotinia 4, 6 other genera), 2 `undetermined`, 1 edge-guarded true MAT1-1-1 (A. fumigatus) | tables read and recounted from the TSV; correctness of the 44 non-Botrytis flips is unverified (no truth), as the vote note says |

Reading: cross-idiomorph remnants are common in Sclerotiniaceae and Ophiostomatales (Abate 2018 dMAT1-1-1; Tsui 2013), and the 54 triggered loci are the
population a remnant status would describe in genomes. Curated remnants are 9 records.

## 2. What `present`, `role` and the new field touch today

Verified by grep over `src/`:

- `genes[].present` is read ONLY in `src/MATPredict/db/` and `detect/benchmark.py`:
  - `gff_export.py` (lines 26, 47, 219, 301): `present: false` genes get no GFF3 gene, no GenBank feature, no `proteins.faa` entry.
  - `cli.py` (170, 250): the backfill staleness test and protein fetch skip them.
  - `schema.py:80` (`validate_gene_vocabulary`) and `validate.py:161`, `synteny.py:149`: skipped.
  - `build_duckdb.py:140`: stored as the `present` column.
  - `detect/benchmark.py:585`: the self-consistency benchmark's expected set is the `present` genes; a `present: true` remnant that detect does not report counts as a real "not found".
- `detect/search.py`, `family_registry.py`, `reference_fasta.py`, `scoring.py`, `idiomorph.py` and `report.py` do NOT read the record's `present`. They read
  `proteins.faa` and `order.yml`. Consequence: **a `present: false` gene has no effect on detection**, only that it contributes no sequence. The report's
  `present=false` GFF3 attribute is generated from the family roster's `genes_missing`, not from curated records.
- `role` (`core_MAT|flanking_conserved|flanking_variable`) is read by `reference_fasta`, the `proteins.faa` header and the classifier; a remnant keeps its
  gene identity's role (sas405 gene 2 is `core_MAT`).
- The coverage rule (`detect/partial_core.py`, `pipeline.model_reference_coverage`, `full_length_models`) measures a model against the length of the
  reference RECORD it was built from (`record_lengths[(record_id, gene)]`), not the longest. **A remnant stored as `present: true` with a protein is a
  reference of its own length, so a model built from it reads ~100% "full length".** This affects three records: a1163's 250-aa MAT1-2-1 (78% of Af293's 322),
  uamh-11059's 224-aa MAT1-1-1 (63% of 353), and 107463's 104-aa MAT1-2-1. Whether a model is ever built from them (best reference) in real genomes was
  NOT measured; the Afum 168 `homothallic_candidate` loci suggest it happens for A1163. The sas405 remnant, being `present: false`, cannot have this effect.
- Held-out tooling (`detect/holdout.py`): reads `mating_type.idiomorphs` and the locus spans, never the genes. Withholding a remnant record changes nothing
  about how remnants are scored; `benchmark.py` is the only consumer that needs a change (below).
- `additionalProperties` is not set on gene items, so a new gene key already validates today; the gap is vocabulary and consumers, not validation.

## 3. Options

### Option A. Keep `present: false` (status quo)

Schema diff: none. Migration: none. Records touched: 0 (sas405 stays as merged). Tests: none.
For the coverage rule: nothing changes; the rule works from sequence evidence and never read the record's gene status. For held-out tooling: nothing.
Costs: the three remnant `present: false` entries are indistinguishable from true absences (deletion, missing flank) except in prose; the four `present: true` remnants
(a1163, uamh-11059, sxlc146, 107463) cannot be identified at all and a1163 and uamh-11059 remain references of their own truncated length; `benchmark.py`
counts them as expected genes; the Afum `homothallic_candidate` open question stays open. Loci affected: 0.

### Option B. Gene-level `status` (recommended)

```yaml
# db/_schema/metadata.schema.yaml, genes[].items.properties (additive; absent = intact)
        present: {type: boolean}
        # status: what the recorded gene IS. Absent means intact (a gene model that
        # encodes a functional protein). Orthogonal to `present`: `present` says whether
        # the record carries a gene model with sequence; `status` says what that gene is.
        #   remnant     a truncated or degenerate copy of an idiomorph's gene, typically
        #               the cross-idiomorph copy left in the other idiomorph
        #               (SAS405 MAT1-1-1, A1163 MAT1-2-1). Not a functional gene.
        #   pseudogene  a gene with its own coding sequence disabled (stop, frameshift,
        #               annotated pseudogene): atcc11807 glrA/rnhA, 107463 MAT1-2-1.
        # Not for annotation or assembly truncation (use locus.core.completeness).
        status: {enum: [remnant, pseudogene]}
        # remnant_note: optional, one sentence of evidence (span, % of reference).
        remnant_note: {type: string}
```

Semantic rules (checked by a new `validate_gene_status`): `status` set and `present: false` means no model recorded (as sas405 today); `status` set and `present: true` means
a model is recorded, and its protein is NOT used as a search reference (B2); `status` on `role: core_MAT` or flanking genes is allowed.

Migration (records touched: 7, plus 2 optional):

| record | edit |
|---|---|
| 40559_sas405_MAT_MAT1-2 | gene 2: add `status: remnant` (stays `present: false`); drop the "present: false so that it contributes no sequence" wording from the note; record_version +1 |
| 28548_06cwm-g22_MAT_combined | gene 5 MAT1-1-1: `status: remnant`, record_version +1 |
| 1708542_exf-10342_wallMAT_v2 | gene 3 HMG: `status: remnant` (SXI1 stays absent without status), record_version +1 |
| 746128_a1163_MAT_MAT1-1 | gene 2 MAT1-2-1: `status: remnant`; changes the note's "recorded as present core_MAT so the remnant is not lost" |
| 230073_uamh-11059_MAT_MAT1-2 | gene 1 MAT1-1-1: `status: remnant` |
| 578113_sxlc146_MAT_MAT1-2 | gene 4 MAT1-1-2: `status: remnant` ONLY after the protein is checked (391 aa; not verified as truncated) |
| 107463_unknown-1_MAT_MAT1-2 | gene 1 MAT1-2-1: `status: pseudogene` |
| 101119_atcc11807 Plus and Minus (optional) | add rnhA / glrA as `present: false, status: pseudogene`, `flanking_*`, or leave as prose |

Code changes, in two steps:

- B1 (metadata only, no change in any detection output): `schema` enum; `validate_gene_status`; `build_duckdb.py` + `duckdb_schema.sql` add a `status` column;
  `benchmark.py:585` expected set = `present` and no `status`; `gff_export.py` writes `/pseudo` on a `status` gene's feature (GenBank convention) and still skips
  `present: false`; docs line in `docs/`. Loci affected in genomes: 0.
- B2 (behavioural, separate decision after measurement): `gff_export.write_proteins_fasta` omits `status` genes, so `reference_fasta` and
  `_curated_protein_lengths` stop treating them as references (the remnant stops reading as "100% full length"; a1163's, uamh-11059's and 107463's proteins
  leave the search set). Expected effects to measure: Afum MAT1-1 loci (160 to 169; the `homothallic_candidate` class, 168 on 2026-10-07 on an older rule), Ophiostomatales
  MAT1-2 loci, Xanthoria loci. B2 could remove the only MAT1-1-1 reference for a clade (uamh-11059's paired MAT1-1 record is 230073_uamh-1363, whose MAT1-1-1 was not checked; 107463 has
  two other Teloschistales MAT1-2 records (115814, 1301501_voucher-267; lengths not checked)), so the held-out run at SPECIES/GENUS radius is the check.

Tests needed: schema accepts `status: remnant|pseudogene` and rejects other values and `partial` (`tests/db/test_schema.py`); `validate_gene_status` (status
without a reason note allowed, status on an intact-looking gene with protein flagged as a warning only); `test_build_duckdb` column; `test_gff_export` (`present:
false` skipped, remnant `/pseudo`, `proteins.faa` omission in B2); `test_benchmark` expected set; a database-wide test that the migrated records carry the flag (a
list in the test, like `test_every_accepted_record_gene_name_is_declared`); B2 adds a regression: Afum MAT1-1 locus class before and after, Botrytis 34 genomes
unchanged (they use sas405 `present: false`, no effect), the 2,015 contested loci unchanged except Afum/Ophiostoma.
Coverage rule: unchanged in B1; in B2 the denominator for models built on remnant references changes (more honest). Held-out tooling: unchanged except
the `benchmark.py` expected set. Cost: about 1 day for B1 including migration of 7 records and tests.

### Option C. Locus-level field

```yaml
# locus.core
          cross_idiomorph_remnants:
            type: array
            items: {type: object, required: [gene, of_idiomorph],
                    properties: {gene: {type: string}, of_idiomorph: {type: string},
                                 span: {type: string}, reference_coverage: {type: number}, note: {type: string}}}
```

Fits the cross-idiomorph pattern only: sas405, uamh-11059, a1163, sxlc146, 06cwm-g22 (5 records). Does not fit Wallemia HMG (the idiomorph is
unknown), 107463 or the atcc11807 pseudogenes (same-gene disabling, not cross-idiomorph), so those 4 records need another mechanism. It also has no link to
`genes[]`, so a1163's and uamh-11059's proteins still enter the references, and the `/pseudo` export and `benchmark.py` need a lookup by name. Detection and the coverage rule: no change.
Its one advantage: it can describe a remnant that has no gene entry (27339's `misc_feature`, the Ustilago a2 pheromone remnant), with a span and coverage.
Records touched: 5 to 6.

### Comparison

| | A keep | B gene `status` | C locus field |
|---|---|---|---|
| Records changed (curated) | 0 | 7 (9 with the optional pseudogenes) | 5 to 6 |
| Remnant cases it can express | 3 of 9 (present false only) | 9 of 9 | 5 of 9 |
| Distinguishes remnant from true absence | prose only | yes | partly (only listed) |
| Marks stored remnant proteins as non-reference | no | yes (B2) | no |
| Effect on detect in B1 | none | none | none |
| Genome loci touched by behaviour | 0 | B1: 0; B2: Afum (~160 to 169), Ophiostomatales, Xanthoria, not measured | 0 |
| Coverage rule | unchanged | unchanged (B1); remnant denominators fixed (B2) | unchanged |
| Held-out tooling | unchanged | `benchmark.py` expected set | needs name lookup |

## 4. Recommendation

Option B, B1 now and B2 only after a measured run.

1. It is the only option that expresses all nine true cases, including a remnant that has a protein, which is where the real hazard is
   (a1163's 250-aa and uamh-11059's 224-aa remnant proteins are references that make a model look full length to the coverage rule).
2. It does not change detection: no `src/MATPredict/detect` code reads the record's `present` outside `benchmark.py`, so B1 is metadata, exports and one benchmark line.
3. It is additive with no default change; existing records stay valid and the 3 `present: false` remnants get the label they lacked.
4. `status` is separate from `present` and from `completeness`, so "no sequence recorded", "what it is" and "assembly truncation" stay three different things.
5. It does not make the coverage rule depend on the DB: the rule keeps working from alignments in genomes it has no record for (the 31 Botrytis flips and 54 triggered loci).
   A remnant status is documentation and a reference filter, not a replacement for the rule.

Do not add `partial` to the enum. In the database, "partial" already means assembly or annotation truncation (`completeness`; 498019_b11243, 4952_w29, sas405's locus), and a
gene-level `partial` would mix that with a biological remnant.

## 5. What I could not verify

- The ~11% MAT1-2-1 relic coverage and the "~34 genomes / B05.10" breakdown of that relic: the `botrytis_mislabels` data is not in the repo. The vote note cites 16% (Sl8).
- Whether the 168 Afum `homothallic_candidate` calls still hold on current main (older rule; the full-length condition may already have cleared them), and whether any
  model in a real genome is built from the a1163, uamh-11059 or 107463 proteins as its best reference: not run (no HPCC job; reports are on /bigdata).
- That sxlc146's MAT1-1-2 (391 aa) is truncated or a remnant at the protein level; the curator's note says remnant.
- That 06cwm-g22's MAT1-1-1 and the Wallemia HMG fragment are the same kind of event as sas405's (read from the notes, not from sequence).
- The correctness of the 44 non-Botrytis campaign flips (no truth), and the edge case Rf2 (remnant on a separate contig, which the rule cannot see; a remnant sequence record with
  coordinates might catch it, untested).
- The duckdb DDL (`db/_schema/duckdb_schema.sql`) column list was not opened; the B1 `status` column is an assumption of one added column.
- Counts of "90+ files" in the first grep were not itemised; only the keyword-specific greps were classified.

## 6. Questions for the curator

1. Is `status` (orthogonal to `present`) acceptable, or should a remnant be `present: false` by definition?
2. Should a remnant carry coordinates and a protein when they exist (A1163 and UAMH do), with `status` stopping it being a reference (B2), or should those two records drop the protein?
3. Enum: `remnant | pseudogene` only, or also `partial` (my advice: no)?
4. Should B2 (remnant proteins leave the search references) wait for a held-out and Afum rerun? Is the Afum `homothallic_candidate` label the intended consequence to clear?
5. sxlc146 MAT1-1-2 (391 aa) and 107463 MAT1-2-1 (stop-truncated): remnant, pseudogene or intact?
6. Do you want atcc11807's glrA/rnhA pseudogenes and the 27339 `misc_feature` recorded (new curation), or left in prose?
7. Should the SAS405 remnant get coordinates (nt 3267-3960 of KF944386.1) so a future decoy or locus-split check could use the sequence?
8. Do you want a separate Ustilago a2 record to capture the remnant a1 pheromone gene?
