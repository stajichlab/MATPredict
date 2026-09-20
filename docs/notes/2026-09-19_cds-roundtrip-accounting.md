# CDS round-trip accounting — all 181 curated CDS features (2026-09-19)

First generated after merging `backfill-gbk-staleness-sweep` (`b5e5b91`).
Re-measured on branch `gene-class-join` after all 59 accepted records were
regenerated so `locus.gbk` carries the new `/gene_class` CDS qualifier. Every
round-trip number below was re-derived from the regenerated files and is
unchanged; only the `gene_class` section at the end is new. The per-gene table
`docs/notes/2026-09-19_cds-roundtrip-accounting.tsv` gained a populated
`gene_class` column (it replaces the always-empty `gclass` column); all its
other columns are byte-identical to the `b5e5b91` measurement.

## Method

For every `CDS` feature in every accepted record's `locus.gbk`: extract the
feature's nucleotide sequence from the record's own written sequence, apply the
feature's `/codon_start` and `/transl_table`, translate, and compare to that
same feature's `/translation` qualifier. PASS means the file is internally
consistent — its coordinates reproduce its own protein.

`excess_codons` = (CDS nt / 3) − protein length. For a gene whose introns are
not modelled, this approximates the unmodelled intron burden in codons.

## Result

| | Count |
|---|---|
| Total CDS features | 181 |
| **PASS** | **103** |
| **FAIL** | **78** |

### Every failure has one cause

| Failure class | Count |
|---|---|
| Genes with modelled exons that fail | **0** |
| Genes with `exons` unmodelled (`exons: None`) that fail | **78** |
| …of those, showing in-frame internal stop codons | **78 (all)** |
| …of those, failing only on length with no internal stop | 0 |

All 78 failures are the same thing: the gene has introns, the curated record
does not model them, so the CDS location spans intron sequence and translating
it hits stop codons. This is a **curation-completeness gap in `metadata.yaml`,
not a code defect** — the `/translation` qualifier remains authoritative and
correct in every case (all 181 match their `proteins.faa` entry).

Of the 130 genes with no modelled exons, **52 genuinely have no introns** and
round-trip cleanly. The remaining 78 need exon boundaries.

### Where the 78 sit

| Cut | Count |
|---|---|
| `core_MAT` | 53 |
| `flanking_conserved` | 25 |
| Implied burden <20 codons (~one short intron) | 9 |
| 20–60 codons | 24 |
| 60–150 codons | 29 |
| >150 codons (large, or a mis-specified gene span) | 16 |

Median implied burden 74.5 codons; range 1.7 to 1266.

### Worth looking at first

Cheapest (likely a single short intron), all `core_MAT` unless noted:

- `230073_uamh-11059_MAT_MAT1-2` SLA2 — 1.7 excess codons, 2 stops (*flanking*)
- `55307_bsh1-11_MATtub_MAT1-1`, `57749_ym1-1_MATtub_MAT1-1`,
  `752769_k467_MATtub_MAT1-1` MAT1-1-1 — 17.0 each
- `746128_a1163_MAT_MAT1-1` MAT1-1-1 — 18.0
- `5346_a43-b43-okayama-7_PR_B43` fungal_mating_type_pheromone — 18.3

Largest, which may indicate a mis-specified span rather than introns:

- `2903220_liq80xsp_MAT_combined` / `2903222_liq146xsp_MAT_combined` SLA2 —
  1266 / 1265 excess codons. These two are the MAG-derived records with no
  protein accessions; their SLA2 and APN2 (544.7 each) dominate the tail.
- `5346_a43-b43-okayama-7_PR_B43` pheromone_receptor ×4 — 148 to 295 each.

## `gene_class` coverage after the order.yml join

Superseded finding. This note originally recorded `gene_class` as unset for all
181 CDS. That is no longer true. `write_genbank` now resolves each present
gene's `gene_class` from its phylum vocabulary `db/<phylum>/order.yml`, keyed on
`(mating_type.locus_name, gene.name)`, and emits it as a `/gene_class` CDS
qualifier. All 59 accepted records were regenerated through
`matpredict curate-db build-gff`; the only change to any `locus.gbk` was 112
added `/gene_class` lines across 45 records (0 deletions), and `locus.gff3` and
`proteins.faa` were byte-identical.

| | Count |
|---|---|
| CDS features | 181 |
| **CDS carrying `/gene_class`** | **112** |
| CDS without it | 69 |

| `gene_class` | CDS |
|---|---|
| HMG_box | 21 |
| pheromone_precursor | 20 |
| alpha_box | 19 |
| sla2_homolog | 17 |
| apn2_homolog | 16 |
| pheromone_receptor | 9 |
| HD1 | 5 |
| HD2 | 5 |
| **Total** | **112** |

### Why the other genes carry no qualifier

Counted over the 182 present genes (181 have a CDS feature; one does not).

| Reason | Genes |
|---|---|
| Resolves to an `order.yml` entry that declares no `gene_class` | 70 |
| Gene name not declared in `order.yml` at all (lookup miss) | 0 |
| Has a `gene_class` but no CDS feature to hang it on | 0 |

So the qualifier is limited only by vocabulary coverage, never by a failed join.
The single present gene with no CDS feature (`4754_b80_PM_combined`) is in the
70 unclassified, so no classified gene lost its qualifier.

Of the 78 round-trip failures above, 57 now carry a `gene_class`.
`db/synteny.py`'s clinker `--gene_functions`/`--colour_map` generation can now
colour by 8 `gene_class` values on 112 of 181 CDS, instead of by `role` alone
(2 values). A visualization or Pfam/domain plan that needs per-`gene_class`
granularity has data for roughly 62% of curated CDS, not 0%.
