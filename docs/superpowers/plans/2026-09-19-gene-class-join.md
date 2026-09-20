# Plan: populate `gene_class` on exported features (join from order.yml)

## Spec / authority

User decisions taken on 2026-09-19 (binding, recorded verbatim):

1. **HD mapping**: `bW` -> `HD1`, `bE` -> `HD2`, `SXI1` -> `HD1`, `SXI2` -> `HD2`.
   The `HD1`/`HD2` enum values are kept; no generic `homeodomain` value is added.
2. **Flanking classes**: `APN2`/`apn2` -> `apn2_homolog`; `SLA2`/`sla2` ->
   `sla2_homolog`. Named after the canonical gene, not the function. This is a
   deliberate departure from the schema comment's "functional cross-family
   identity" framing for these two entries; the user chose it for extensibility
   to `COX13`, `tptA`, `rnhA` and the other flanking genes later.
3. **Source of truth**: `gene_class` stays ONLY in `db/<phylum>/order.yml`. It is
   joined at export time on `(phylum, mating_type.locus_name, gene.name)`.
   `metadata.yaml` and `metadata.schema.yaml` are NOT given a `gene_class` field.

## Context (measured on 2026-09-19, not assumed)

- `gene_class` is defined in `db/_schema/order.schema.yaml` as an optional
  per-gene enum: `[alpha_box, HMG_box, pheromone_receptor, pheromone_precursor,
  HD1, HD2]`. It is orthogonal to `role` (structural position).
- It is populated for **38 of 91** `order.yml` gene entries across 3 phyla /
  19 loci. It is absent from **every** `metadata.yaml` and from
  `metadata.schema.yaml`.
- `gff_export.write_genbank` writes a `/gene_class` qualifier conditionally from
  `gene.get("gene_class")` -- a key that by schema can never exist on a record
  gene. Result: **0 of 181 CDS features in the whole DB carry `/gene_class`.**
- The join was measured over all 59 accepted records: **0 present genes fail to
  resolve** to an `order.yml` entry. 75 present genes already resolve to a set
  `gene_class`; 107 resolve to an entry with none set.
- The only downstream consumer that reads the qualifier is `db/draw.py`, which
  builds its gene label from the `gene_class` qualifier. It gains the label once
  the qualifier is written. `db/synteny.py` does **not** read `locus.gbk`
  annotation at all: `_gene_function_rows` (`synteny.py:140-154`) reads
  `(name, role)` from `metadata.yaml`, and `synteny.py:57-63` documents labeling
  clinker's `--gene_functions`/`--colour_map` by `role` ALONE as a deliberate
  anti-fragmentation choice. **This branch therefore improves `draw.py` labels
  only and has zero effect on clinker/synteny output.**

## Global Constraints

- **Do not add `gene_class` to `metadata.schema.yaml` or to any `metadata.yaml`.**
  Decision 3 above is binding.
- **Do not invent `gene_class` values.** Only the 9 `order.yml` entries named in
  Task 1 gain one. Every other unset entry stays unset.
- The join key is `(phylum, mating_type.locus_name, gene.name)` and is
  **case-sensitive** -- `APN2` and `apn2` are distinct entries in distinct loci
  and both must be curated.
- A record whose gene does not resolve, or resolves to an entry with no
  `gene_class`, must write NO `/gene_class` qualifier -- absent, never empty.
- Tests run with `pixi run -e test pytest -v` -- never bare `pixi run pytest`.
- Duplicate gene names within a record are normal biology; the lookup is by name
  and legitimately returns the same class for every copy. Do not deduplicate.

## Task 1 -- schema, curation, join code, tests

Files: `db/_schema/order.schema.yaml`, `db/Ascomycota/order.yml`,
`db/Basidiomycota/order.yml`, `src/MATPredict/db/gff_export.py`,
`src/MATPredict/db/cli.py`, tests.

Tests first.

1. **Extend the enum** in `db/_schema/order.schema.yaml` to
   `[alpha_box, HMG_box, pheromone_receptor, pheromone_precursor, HD1, HD2,
   apn2_homolog, sla2_homolog]`. Update the explanatory comment above it to say
   that the two `*_homolog` values are named after a canonical gene rather than a
   function, and why (user decision, extensibility to other flanking genes).

2. **Curate exactly these 9 `order.yml` entries** (verified present, currently
   unset; usage counts are curated present genes affected):

   | phylum | locus | gene | gene_class | uses |
   |---|---|---|---|---|
   | Ascomycota | MAT | APN2 | apn2_homolog | 15 |
   | Ascomycota | MAT | SLA2 | sla2_homolog | 15 |
   | Ascomycota | MATyl | apn2 | apn2_homolog | 1 |
   | Ascomycota | MATsc | sla2 | sla2_homolog | 1 |
   | Ascomycota | MATyl | sla2 | sla2_homolog | 1 |
   | Basidiomycota | MAT | SXI1 | HD1 | 1 |
   | Basidiomycota | MAT | SXI2 | HD2 | 1 |
   | Basidiomycota | bLocus | bW | HD1 | 1 |
   | Basidiomycota | bLocus | bE | HD2 | 1 |

   Total newly classified: 37 curated present genes.

3. **Implement the join.** Add a resolver to `src/MATPredict/db/gff_export.py`
   (or a small module beside it) that loads a phylum's `order.yml` once and
   answers `(locus_name, gene_name) -> gene_class | None`. Then thread the
   resolved classes into `write_genbank` as a `gene_classes: dict[int, str]`
   parameter keyed by `gene_index` -- deliberately mirroring the existing
   `sequences: dict[int, str]` parameter's shape, so the function keeps taking
   plain data and stays testable without touching the filesystem.

   `write_genbank` must prefer the passed-in mapping and must no longer rely on
   `gene.get("gene_class")`. Keep the parameter optional with a default so every
   existing caller and test keeps working unchanged.

   `build_gff_for_record` (`src/MATPredict/db/cli.py`) is the single place that
   loads `order.yml` and builds the mapping, since it already knows the phylum
   and already owns the one regeneration path.

4. **Required tests** (no network; mock `ncbi`):
   - a gene resolving to a class gets the `/gene_class` qualifier;
   - a gene resolving to an entry with no class gets NO qualifier (not empty);
   - a gene not present in `order.yml` at all gets NO qualifier and does not raise;
   - two genes sharing a name in one record both get the qualifier (the real
     Basidiomycota multi-copy case);
   - a resolver test proving the lookup is case-sensitive (`APN2` vs `apn2`);
   - a schema test that every `gene_class` value used across all `order.yml`
     files is a member of the enum -- this is the permanent net against a typo'd
     class silently disappearing from figures.

## Task 2 -- regenerate and verify

1. Regenerate every one of the 59 accepted records via the existing per-record
   path (`curate-db build-gff`), so `/gene_class` lands in `locus.gbk`. Do not
   write a new regeneration path.
2. Report, measured: how many of the 181 CDS now carry `/gene_class`, broken
   down by value; how many do not and why (resolves-but-unset vs unresolvable).
   Expected from the measured join: 75 already-classified + 37 newly classified
   = **112 of 182 present genes**, with the CDS count depending on how many of
   those genes have a translation.
3. Confirm no `locus.gbk` lost a CDS, a translation, a `/codon_start` or a
   `/transl_table` in the regeneration -- diff the feature counts against the
   pre-change state and report any record whose CDS count changed.
4. Re-run the CDS round-trip accounting
   (`docs/notes/2026-09-19_cds-roundtrip-accounting.tsv` was generated at
   `b5e5b91`: 181 CDS, 103 PASS, 78 FAIL, all 78 `exons: None`). Confirm the
   PASS/FAIL split is unchanged -- this task must not alter translation at all.
5. Commit code+tests separately from regenerated data.
