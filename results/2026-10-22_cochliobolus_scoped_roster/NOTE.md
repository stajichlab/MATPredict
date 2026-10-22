# Pleosporaceae-scoped GAP1/BGL1 roster flanks: prototype and measurement (2026-10-22)

Branch `cochliobolus-scoped-roster` (commits dcbdd3a, 7828dff, 24d01fa on origin/main 2d14712). Follows
`analysis/2026-10-18_cochliobolus-mat-recall.md` (branch `cochliobolus-mat-recall`), whose global prototype recovered 90 of 101
uncalled Bipolaris/Curvularia/Exserohilum/Ophiosphaerella genomes but lowered 24 of 94 Dothideomycetes controls from high to medium.
"Measured" = a number from the runs below. "Judgement" is labelled.

## 1. Design

No existing mechanism scopes a roster gene to a taxon: `taxonomic_scope`, `fallback_subphylum` and `fallback_searchable` act on
whole families; `db/taxon_overrides.tsv` is about mislabelled genomes. A second Pleosporaceae family would move the B. maydis records
out of `Ascomycota:MAT` and rename every Pleosporaceae call. So the prototype adds, in the style of the other policies:

* **Per-gene `taxonomic_scope`** in `order.yml` (`db/_schema/order.schema.yaml`). `family_registry.route()` removes a scoped gene from
  the routed family unless `--scoped-roster-policy on` and the genome's lineage (taxid + NCBI ancestors, the cached document routing
  already read) reaches the scope (`scope_roster_genes`, `_apply_gene_scopes`). Operator routes (`--phylum`, `--exhaustive`), no taxid
  or a failed lookup remove it. `reference_fasta.build_reference_fasta(families=...)` leaves a scoped gene's proteins out of the
  query set unless the routed roster still carries it (callers that pass no `families` never query them). **Default off**: with the
  policy off, or outside the scope, the roster and the query set are byte-for-byte those of main.
* **GAP1 and BGL1** on the `Ascomycota:MAT` roster: `flanking_conserved`, `taxonomic_scope: [28556]` (Pleosporaceae; NCBI esummary:
  family, Pleosporaceae), `optional: true` and each in a **`slot` of its own**. Optional keeps an absent GAP1/BGL1 out of
  `fraction_found` (what cost the global prototype its 24 controls: the denominator grew from 4 to 6); the one-member slot adds +1/+1
  when found, so a lone MAT gene with both flanks scores (1+2)/(4+2) = 0.5, the same as the non-optional global prototype, and a
  call that lacks them keeps its score. Reuses the slot mechanism approved for SXI1/SXI2 with one member per slot: needs a ruling.
* **Per-gene `min_identity`** (`search.drop_below_gene_min_identity`, before clustering), 50% for GAP1/BGL1. Added after the first run
  (v1): BGL1 is a beta-glucosidase family, and 27.7%-identity BGL1 fragments (13.8% coverage) beside an HMG paralog (54%, 21% coverage)
  and a COX13 fragment (32%) on another contig made **8 spurious `high` second MAT calls in Alternaria** (`changes_v1_vs_v2.tsv`).
  Every GAP1/BGL1 hit beside a real MAT gene is >= 70% (`flank_ev_proto2.tsv`: 571 rows, GAP1 70-100%, BGL1 80-100%, all >= 50%
  coverage, all within 10 kb of the core gene; v1 had 8 more BGL1 rows at 20-30%, < 50% coverage, far from core). Only scoped genes
  declare a floor, so no other run changes. The 50 is the gap's midpoint, chosen on these data (no held-out set).

ORF1 (AAB84003.1) was not added.

## 2. References (verified)

C. heterostrophus (= B. maydis) deposits, Turgeon et al. 1993 (PMID 8479433): C4 MAT-2 AF027687.1, GAP1 **AAB84002.1** (714 aa,
complement join 633..3006, 5 exons), BGL1 **AAB84005.1** (870 aa, 8898..11742, 4 exons); C5 MAT-1 AF029913.1, GAP1 **AAB82943.1**
(714 aa), BGL1 **AAB82946.1** (870 aa, 9022..11866). Each CDS join translated from the deposit (`verify_cds.py`, Biopython) equals the
NCBI protein record (no internal stop). C4 and C5 differ at 1 (GAP1) and 3 (BGL1) residues. The C4 pair is the same sequence the
earlier prototype used (`mat_adjacent_refs.faa`). Added to `5016_c4_MAT_MAT1-2` and `5016_c5_MAT_MAT1-1` as genes 1/2 (segment
extended to 633..end), locus files regenerated with `curate-db build-gff` (`edit_records.py`).

## 3. Runs

Every Dothideomycetes genome of `results/2026-10-18_dothideomycetes_audit/genome_list.tsv` (2,731; `lists/all.tsv`), HPCC `short`,
`worker.slurm`, `--partial-core-policy coverage_profile --confidence-cap-policy off --strong-single-gene-policy off`, no
`--homothallic-policy`, lineage routing from the taxid, `timeout 3600`.

| variant | code | flags | genomes |
|---|---|---|---|
| base | origin/main 2d14712 | - | 2,731 (2 time out in every variant: GCA_982436715.1, GCA_986281625.1) |
| proto (v1) | dcbdd3a | `--scoped-roster-policy on` | 2,731 |
| **proto2 (v2, final)** | 24d01fa | `--scoped-roster-policy on` | 2,731 |
| off | branch (7828dff/24d01fa) | policy off | 297 (`lists/sub.txt`: clade, 94 controls, 100 random) |
| baserep | main again | - | the same 297 |
| ssg | main | `--strong-single-gene-policy on` | 707 (`lists/ssg.txt`: clade, controls, Pleosporaceae, audit no-calls) |

Scope membership: 477 genomes have 28556 in their NCBI lineage (`lineage_scope.py`, MATPredict's own resolver; all lineages resolved).
Comparison (`collect.py`, `analyze.py`): `report_sha` = SHA-256 of the whole `detection_report.yaml` minus its `run` block; locus key =
family, contig, start, end, idiomorph, confidence, locus_class, genes_found.

## 4. Results (proto2)

**Outside the scope: 2,252 of 2,252 genomes with both runs have byte-identical reports** (everything except the `run` block), so
idiomorph, confidence, genes_found, class and every suppressed/not-detected entry are unchanged. The policy-off branch equals main on
297 of 297 and main equals itself on 297 of 297 (no run-to-run noise in this setting).

**The 101 uncalled clade genomes** (`clade_recovery.tsv`): 90 called (Bipolaris 57/57, Curvularia 23/23, Exserohilum 10/10,
Ophiosphaerella 0/11, Phaeosphaeriaceae, out of scope). 79 `medium/partial_locus`, 11 `high/mat_locus` (one genome, B. gigantea
GCA_035079505.1, gets a MAT1-1 and a MAT1-2 locus). Against the pipeline-independent box-ORF scan: 86 agree, 4 agree with one of two
scanned idiomorphs (B. gigantea with both ORFs on different contigs), **0 disagree**. The four C4/C5 strains (strain identity per
Turgeon 1993; the records come from these strains, so not independent) are all called correctly (MAT1-1 for C5, MAT1-2 for C4,
medium). No clade BioSample states a mating type. 7 of the high calls list the opposite core gene (a superseded cross-hit) in
genes_found, as in the earlier prototype.

**Inside the scope, outside the clade** (385 genomes, `changes_inside_nonclade.tsv`): 379 change, none downwards.

| change | n | detail |
|---|---|---|
| GAP1 added to genes_found only (idiomorph, confidence, class unchanged) | 308 | |
| medium -> high | 50 | 44 Alternaria + 2 Pyrenophora `partial_locus` MAT1-x + APN2 -> `mat_locus` + GAP1 (2/4 = tie with the floor -> 3/5); 4 Stemphylium `idiomorph_gene_only` -> `mat_locus` + GAP1/BGL1 |
| class only | 4 | 3 Pyrenophora, 1 Stemphylium `idiomorph_gene_only` -> `partial_locus` (medium kept) |
| new call (medium, partial_locus) | 17 | 9 Alternaria, 7 Stemphylium, 1 Paradendryphiella; **17 of 17 agree with the audit's genome-wide box scan** |
| idiomorph changed / call lost / high -> lower | 0 | (v1 had 8 spurious second calls; removed by the 50% floor) |

Of the earlier prototype's 94 controls: 84 out of scope, identical; 10 in scope: 4 medium->high, 1 class change, 5 GAP1 added. The
global prototype changed 26 of them (24 high->medium, 2 idiomorph sets).

**Fallback, strong-single-gene rule** (main, `--strong-single-gene-policy on`, 707 genomes): 46 of the 101 called (all Bipolaris,
`low/idiomorph_gene_only`; Curvularia/Exserohilum/Ophiosphaerella 0), 46 agree with the scan, the 4 C4/C5 correct; plus 3 new low
Parastagonospora calls; 0 changes to any existing call (656 of 656 other reports identical; 94/94 controls identical).

| | scoped roster (v2, policy on) | single-gene rule |
|---|---|---|
| clade genomes recovered (of 101) | **90** (79 medium, 11 high) | 46 (all low) |
| disagree with the scan | 0 | 0 |
| changes outside Pleosporaceae (2,252) | **0** | 0 on its 707 (by design it touches only families with no call) |
| changes to existing calls inside Pleosporaceae | 362 (all additive/upward: 308 genes, 50 medium->high, 4 class) | 0 |
| new calls outside the clade | 17 (all scan-consistent) | 3 Parastagonospora |

## 5. Tests

`tests/db` 215 passed. `tests/detect` (final code): 1,073 passed, 6 skipped, 1 failed:
`test_pipeline.py::test_fast_path_missing_gene_is_polished_in_the_existing_clusters_window`, the known intermittent; it also fails
in a full `tests/detect` run of unmodified main (1,063 passed, same failure), passed in the first branch run, and passes 5/5 alone on
both. New: `tests/detect/test_scoped_roster.py` (10 tests). `test_not_searched_routing.py`'s route stub gained the new keyword.

## 6. Limits

* Truth is the box-ORF scan and 4 strain identities; no BioSample mating type exists in the clade. The 101-of-101 locus presence and
  the scan's idiomorph depend on its 100/80-bit thresholds.
* The 50% floor and the 0.5-tie behaviour were set on these genomes (no held-out set). The 50 medium->high Alternaria/Pyrenophora
  promotions rest on GAP1 at 70-90% next to the core gene; they are changes the curator may or may not want.
* Only Dothideomycetes were run; the scope cannot reach other classes by construction (outside Pleosporaceae the roster and query set
  are main's), but other phyla were not run.
* `batch_runner` without `--phylum` builds one shared query set and passes no `families`, so scoped genes are never searched there.
* `tables_v1/` is the v1 comparison made while `base` was still running (2,716 of 2,731); superseded by the top-level tables.
* Wall times are cluster-noisy (the same queries ran 25% slower in one array); no runtime claim is made.
