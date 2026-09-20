# Changelog

All notable changes to MATPredict. Versions follow [semantic versioning](https://semver.org/);
release tags may carry a clade suffix naming the lineage whose infrastructure that
release completed.

## [0.5.0] — 2026-09-20 — `v0.5.0-mucoromycota`

Mucoromycota detection becomes usable end to end: idiomorph calling works, the
scoring denominator is honest, and the reference set covers both idiomorphs.

**Measured on 23 ground-truth Mucoromycota genomes** (`testset/Zygo`, curated
absence/presence tables):

| | before | after |
|---|---|---|
| loci on the correct scaffold | 23/23 | 23/23 |
| **idiomorph called correctly** | **0/23** (all `undetermined`) | **23/23** |
| Plus `fraction_found` | 0.571 | **0.833** |
| Minus `fraction_found` | 0.571 | 0.800 |
| roster genes with no reference | `algA`, `glrA` | **none** |
| confidence | — | 18 high / 5 medium |

**Measured on a 44-genus discovery sweep**: strict locus calls rose from
**6 genera to 14**, with `algA`/`glrA` — previously unfindable — present in 12.

### Added

- **Idiomorph resolution.** `sexM` and `sexP` share an HMG box, so one locus
  gene draws both curated references; this happened in 23/23 genomes and made
  the idiomorph uncallable. Overlapping pairs are now collapsed to one gene,
  with the loser annotated `superseded_by` rather than deleted so the ambiguity
  stays in the report.
- **Search-path tiebreak.** The winner is decided by *which search found the
  gene* — the annotated proteome outranks a genome-wide tblastn rescue —
  falling back to identity when both come from the same search. Identity alone
  scored 19/23 once new references were added; this scores 22/23, and 16/16 for
  both on the smaller reference set.
- **`locus_class`**: `mat_locus`, `homothallic_candidate`, `idiomorph_gene_only`.
  Says *what* was found, orthogonally to `detection_pass` (*how* it was
  admitted). A lone `sexM` or `sexP` with no flank is **kept deliberately** as
  training material for a future per-idiomorph HMM.
- **Homothallic detection.** Both idiomorphs in one locus is real biology —
  *Syzygites megalocarpus* encodes both HMG transcription factors, each flanked
  by its own intact gene with the other pseudogenised (Idnurm 2011). Requires
  both genes proteome-supported and within `max_homothallic_separation_bp`.
- **Relaxed second pass**, run only when the strict pass finds nothing
  genome-wide, admitting on the existing `EvidenceFloor` bar and capped at
  medium confidence.
- **Per-locus curation thresholds** in `order.yml`: `min_idiomorph_margin`
  (5.0) and `max_homothallic_separation_bp` (20 kb), alongside the existing
  `max_cluster_gap_bp`.
- **Report fields**: `idiomorph_margin`, `idiomorph_resolutions` (both members'
  identity, coverage and overlap), `locus_class`, `detection_pass`.
- **Diagnostics corpus**: every row now carries `run_id` and `genome_id`, and
  each idiomorph resolution is logged with the measurements a recalibration
  needs. Without `genome_id` no per-genome statistic could be computed from a
  batch at all.

### Fixed

- **`fraction_found` counted genes that could never be found.** `algA` and
  `glrA` were in the Mucoromycota roster with no reference protein anywhere in
  `db/`, so every genome was scored against an unreachable ceiling. The
  denominator now counts only genes present in the reference set the run
  actually searched with. Without this, resolving the sexM/sexP overlap alone
  would have moved Plus to exactly 0.500 — the rejection boundary.
- **A polish crash cost an entire genome.** `exonerate --refine` segfaults
  deterministically on some windows (verified: `--refine full` also crashes,
  no-refine succeeds, miniprot succeeds). Signal deaths now retry without
  `--refine`, then fall back to miniprot; non-zero exits stay fatal so a
  missing binary is still loud. Observed rate 1 in 44 genomes.
- **Relaxed calls carried no gene evidence** — no coordinates, identities or
  references — making them impossible to verify.
- **Resolution ran only before polishing**, but polishing moves coordinates, so
  pairs that overlap only afterwards escaped collapse. It now runs again after.
- **Alignments under 90 bp are dropped.** A 27 bp "gene" (nine codons) was
  being reported as a locus. The bar is on the alignment, well below the 60 aa
  short-ORF floor, so small MAT genes still pass. It is not an identity bar:
  Minus identities run 25.9–43.5%.
- **Reported evidence could be a superseded hit** while the gene counted as
  found via a different live one.
- **`score_match` passed on identity alone**, ignoring the coverage it
  computed: a record validated as `identity=100.0, coverage=0.86, pass`. Now
  warns below 80% coverage — a warning, not a failure, because 7 accepted
  records legitimately sit at 3.0–86.7%.
- **Evidence floor counted one gene as two** when a cluster's only hits were an
  unresolved `sexM`/`sexP` pair on the same protein.

### Database

Mucoromycota grew from 6 to 15 accepted records (9 Minus / 6 Plus, 2
homothallic); every gene validates at 100% identity **and** 100% coverage.
All are tier-1 published deposits.

| record | accession | contribution |
|---|---|---|
| *Mooraboolomyces wintlei* | `OR965930.1` | the only source of `algA` and `glrA` proteins |
| *Absidia urquhartii* ×2 | `PP971768/9` | first Cunninghamellaceae; first Minus outside *Mucor*/*Phycomyces*/*Rhizopus* |
| *Mucor mucedo* | `JN587498.1` | the alginate lyase the `algL`/`algA` question concerns |
| *Rhizopus azygosporus* | `MG967659.1` | second `glrA` |
| *Blakeslea trispora* | `HG939558.1` | `sexM` only — its `tptA` (38 aa) and `rnhA` (82 aa) are fragments |
| *Syzygites megalocarpus* ×2 | `JN112239/40` | **first homothallic reference**; third `glrA` |
| *Parasitella parasitica* | `KY081664.1` | Parasitellaceae Minus |

Assembly-derived loci for *Cunninghamella* and *Chaetocladium* were considered
and **rejected**: they are tier-2 homology-inferred, which the database design
spec excludes in this phase, and a literature search confirmed no MAT-locus
deposit exists for those genera (nor for *Actinomucor* or 11 others checked).

**Reference balance is not neutral.** Adding two Plus references flipped four
*Cunninghamella* Minus genomes to Plus; check per-idiomorph counts before any
ingest.

### Documentation

- `docs/notes/2026-09-20_algA-glrA-mucorales-literature.md` — establishes that
  `algA` and `algL` are one gene (62.8% identity between the two deposited
  proteins; the same author uses both spellings in one publication), that the
  2017 review carries **no citation** for either gene, and that `glrA` traces
  to Idnurm 2011.
- `docs/HANDOFF-detect-scoring-idiomorph.md` — every threshold with the
  measurement behind it, and the open questions.

### Notes

Every threshold is provisional and documented at its definition with its
evidence. `min_identity` remains `None`: Minus identities run 25.9–43.5%, so
any global cutoff would destroy Minus detection.

---

## [0.4.0] — 2026-09-20 — detection targeting (branch `detect-targeting`)

- Phylum routing: the query set narrows from 181 to 19 proteins for a
  Mucoromycota run; `--phylum` override; `routing_mode` in every report.
- `max_cluster_gap_bp` became per-locus curation data in `order.yml`
  (Mucoromycota `MAT` = 50 kb); a run takes the maximum over routed families.
- `EvidenceFloor` defaults to the curator's ruling: ≥2 distinct genes including
  ≥1 `core_MAT`.
- `run_batch` and a SLURM wrapper, so a batch pays the same narrowed cost.
- `--emit-cds-fasta`, an N+1 genome-parse fix and 60-column FASTA wrapping.

Validated on 23 ground-truth genomes at 100% per-gene recall, mean 41.8 s per
genome; the same pipeline pre-fix had been killed after >24 minutes on one
small yeast genome without completing.

## [0.3.0] — 2026-09 — localize-then-polish search

- Genome-wide `tblastn` localization followed by per-gene refinement with both
  `exonerate --refine region` and `miniprot`, classified into four explicit
  polish statuses; only genuinely unconfirmed genes cap confidence at medium.
- Fast-path rescue for core genes a supplied annotation does not contain — the
  blind spot the pipeline exists to close, since small pheromone-precursor
  genes are routinely missing from whole-genome annotation.
- Fragmented-locus handling: multi-segment calls with per-segment
  `contig_edge_distance`, scoped per contig for valid GFF3 parentage.
- Genome acquisition, SLURM-sized batch orchestration and rollout aggregation.

## [0.2.0] — 2026-09 — detection pipeline

- `matpredict detect`: family routing by taxid, diamond fast path against a
  supplied proteome, gap-based clustering, fractional per-family scoring with
  cross-family ambiguity detection, confidence tiering, idiomorph assignment.
- GFF3 and YAML report writers; families that fall short are reported as "not
  detected" with a reason rather than silently omitted.
- Short-ORF genes distinguished from real absences via `genes_not_searchable`.
- Leave-one-out benchmark suite.

## [0.1.0] — 2026-09 — curated reference database

- `matpredict curate-db`: propose → validate → accept/reject lifecycle, with
  GFF3, GenBank and protein FASTA export and a DuckDB query cache.
- Record schema with per-claim evidence tiers and citations, multi-exon
  `exons`/`codon_start`/`transl_table`, and independent translation
  verification — validation re-derives each protein from the record's own
  claimed coordinates rather than comparing a fetched protein to itself.
- Cached NCBI E-utilities and UniProt clients with retry on 429/5xx.
- `order.yml` controlled vocabulary per phylum, with `taxonomic_scope` routing
  and an optional `gene_class` for cross-family gene grouping.
- Locus and synteny diagram commands (`draw-locus`, `draw-synteny`).

---

## Project state at 0.5.0

- **68 accepted records** — 43 Ascomycota, 10 Basidiomycota, 15 Mucoromycota
- **19 loci** defined across three phyla
- **463 tests** passing, 46 test modules, 35 source modules
