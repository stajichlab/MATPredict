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
* **Per-gene `min_identity`** (`search.drop_below_gene_min_identity`, before clustering), measured at 50% for GAP1/BGL1; **BGL1 raised
  to 60% by curator ruling 2026-10-22 (GAP1 stays 50); see section 4b for the re-run**. Added after the first run
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
Ophiosphaerella 0/11, Phaeosphaeriaceae, out of scope). As reported: 79 `medium/partial_locus`, 11 `high/mat_locus` (one genome,
B. gigantea GCA_035079505.1, gets a MAT1-1 and a MAT1-2 locus). **Correction (review 2026-10-22, corrected again 2026-10-22 with `high_support.py`): none of the 11 recovered `high` calls is high
on real flank evidence.** The review's "1 of 11" (and this note's previous version) counted Curvularia kusanoi GCA_011058905.1, but
that genome was already called `high` on main (APN2 64.8%/60.6% and COX13 54.2%/91.2%, both modelled; it only gained GAP1) and is
not a recovery. **All 11 recovered `high` calls reach `high` only because a noise third gene breaks the 0.5 tie**: a lone core gene +
GAP1 + BGL1 scores 3/6 = 0.5 = the floor and is capped `partial_locus`/medium; an unmodelled APN2 HSP at 5-6% coverage (5
Exserohilum) or an unsuperseded opposite-idiomorph core cross-hit at 10-21% coverage (6 Curvularia, which then list MAT1-1-1 and
MAT1-2-1 together in genes_found) adds a required gene, the score becomes 3/5 or 4/7, and the found GAP1/BGL1 flank lets the tier
reach `high`. Split both ways: **as called, 79 medium + 11 high; on real flank evidence, 90 medium-equivalent + 0 high.** The slot
is what moves these clusters over the floor (on main they scored 1/4 or 2/5), so the curator should decide whether such `high`
calls are acceptable before the policy goes on in a clean run. Against the pipeline-independent box-ORF scan: 86 agree, 4 agree with one of two
scanned idiomorphs (B. gigantea with both ORFs on different contigs), **0 disagree**. The four C4/C5 strains (strain identity per
Turgeon 1993; the records come from these strains, so not independent) are all called correctly (MAT1-1 for C5, MAT1-2 for C4,
medium). No clade BioSample states a mating type. 7 of the high calls list the opposite core gene (a superseded cross-hit) in
genes_found, as in the earlier prototype.

**Inside the scope, outside the clade** (385 genomes, `changes_inside_nonclade.tsv`): 379 change, none downwards.

| change | n | detail |
|---|---|---|
| GAP1 added to genes_found only (idiomorph, confidence, class unchanged) | 308 | |
| medium -> high | 50 | 44 Alternaria + 2 Pyrenophora `partial_locus` MAT1-x + APN2 -> `mat_locus` + GAP1 (2/4 = tie with the floor -> 3/5); 4 Stemphylium `idiomorph_gene_only` -> `mat_locus` + GAP1/BGL1. The 44 Alternaria are near-identical strains (same APN2/GAP1 identities), one measurement rather than 44. 2 Pyrenophora (GCA_012365135.1, GCA_025504585.1) promote on an APN2 covering 14.7% of its reference plus GAP1: the same shape as the noise-tie `high` calls above |
| class only | 4 | 3 Pyrenophora, 1 Stemphylium `idiomorph_gene_only` -> `partial_locus` (medium kept) |
| new call (medium, partial_locus) | 17 | 9 Alternaria, 7 Stemphylium, 1 Paradendryphiella; **17 of 17 agree with the audit's genome-wide box scan** |
| idiomorph changed / call lost / high -> lower | 0 | (v1 had 8 spurious second calls; removed by the 50% floor, and still absent at 60%) |

Of the earlier prototype's 94 controls: 84 out of scope, identical; 10 in scope: 4 medium->high, 1 class change, 5 GAP1 added. The
global prototype changed 26 of them (24 high->medium, 2 idiomorph sets).

**Fallback, strong-single-gene rule** (main, `--strong-single-gene-policy on`, 707 genomes): 46 of the 101 called (all Bipolaris,
`low/idiomorph_gene_only`; Curvularia/Exserohilum/Ophiosphaerella 0), 46 agree with the scan, the 4 C4/C5 correct; plus 3 new low
Parastagonospora calls; 0 changes to any existing call (656 of 656 other reports identical; 94/94 controls identical).

| | scoped roster (v2, policy on) | single-gene rule |
|---|---|---|
| clade genomes recovered (of 101) | **90** (79 medium, 11 high as called; 1 high on real flank evidence, 10 high only via a noise third gene) | 46 (all low) |
| disagree with the scan | 0 | 0 |
| changes outside Pleosporaceae (2,252) | **0** | 0 on its 707 (by design it touches only families with no call) |
| changes to existing calls inside Pleosporaceae | 362 (all additive/upward: 308 genes, 50 medium->high, 4 class) | 0 |
| new calls outside the clade | 17 (all scan-consistent) | 3 Parastagonospora |

Not observed is not impossible: no in-scope existing call moved down, but GAP1/BGL1 can still change an existing call
non-additively through the per-family polish cap (`_polish_rank` counts distinct genes, so a cluster gaining GAP1/BGL1 can outrank
and displace another when > 6 clusters are admitted), by bridging two clusters within the 25 kb gap, or by admission (a core
fragment + one flank now clears `min_hits=2`, which is how the 17 new calls and the 8 v1 spurious calls arose). None fired in
these 477 genomes.

**BGL1 floor margin (as measured at 50%; now 60%, see 4b).** The 50% floor sits only 3-6 points above a full-length BGL1 paralog (44.4-46.7% identity over 802-890 aa,
on another contig) that the review found in every one of 15 sampled recovered genomes (GAP1's nearest paralog: 36.9% over 103 aa).
A genus whose BGL1 paralog drifts above 50% would rebuild the v1 spurious second calls. The value is not changed here (every real
flank hit is >= 70%, so a floor near 60 would keep recall with a wider margin; curator's call); it must be re-validated for every
new genus brought into the scope.

## 4b. BGL1 floor 60% (curator ruling 2026-10-22), re-run

Code: origin/main after PR #132 (2cb453f: BGL1 `min_identity: 60`, GAP1 50). Variant `proto60`, same clean-run settings and
`--scoped-roster-policy on`, on 722 genomes (`lists/f60.txt`): the 103 clade genomes, the 94 controls, all 477 in-scope genomes
(including the 385 in-scope non-clade genomes and the 8 Alternaria with the v1 spurious second call) and 150 random out-of-scope
genomes (245 out-of-scope genomes in all). Tables: `floor60/`.

* **Recoveries: 90 of 101, unchanged** (Bipolaris 57/57, Curvularia 23/23, Exserohilum 10/10, Ophiosphaerella 0/11); none lost;
  same 79 medium + 11 high; 0 disagree with the scan; the four C4/C5 strains still correct.
* **No call changed against the 50% run.** 623 of 722 reports are byte-identical to proto2 (minus `run`); the other 99 have the same
  calls (family, contig, start, end, idiomorph, confidence, class, genes_found as a set) and differ only in (a) the order of
  `genes_found` in 48 detected calls, (b) withheld loci (`suppressed_loci`: 50 genomes, fewer or smaller withheld clusters, because
  BGL1-paralog hits at 50-60% no longer seed them) and `suppressed_unpolished` counts (4), and (c) one `not_detected` entry. All 571
  GAP1/BGL1 rows in detected calls are the same as at 50% (lowest BGL1 87.4%).
* **Inside the scope, non-clade (385): the same counts as at 50%**: 308 GAP1 added, 50 medium->high, 4 class only, 17 new calls
  (17/17 scan-consistent), 0 downwards. The 8 v1 spurious second calls stay absent.
* **High split** (`high_support.py`: "supported" = a SLA2/APN2/COX13 flank modelled over >= 50% of its reference): recovered clade
  highs 0 supported / 11 tie-break; of the 50 promotions, 44 Alternaria supported (APN2 51-53%, 95-98% coverage; near-identical
  strains), 2 Pyrenophora tie-break (APN2 14.7% coverage) and 4 Stemphylium tie-break (no SLA2/APN2/COX13; high through GAP1+BGL1
  plus an opposite-idiomorph cross-hit).
* **Outside the scope: 245 of 245 reports identical to main** (minus `run`); controls 84/84 identical.
* Tests on a compute node at 2cb453f: `tests/detect` 1,076 passed, 6 skipped, 1 failed (the known intermittent
  `test_fast_path_missing_gene_is_polished_in_the_existing_clusters_window`, which also fails in full runs of main and passes alone).

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
* Reachability: `run_pipeline` and `batch_runner` have no `scoped_roster_policy` parameter (`run_pipeline` self-routes with the
  default off; `batch_runner` routes by `--phylum`, an operator route, and builds the query set without `families`), so the policy is
  reachable only through per-genome `matpredict detect`. The clean run uses per-genome `detect` (scripts/submit_binned_panel.py ->
  scripts/run_clade_panel.slurm, flags via `DETECT_ARGS`), so `--scoped-roster-policy on` can be passed there.
* Scoping outcome is recorded: the run block's `parameters.scoped_roster` says whether the scoped genes were searched; a failed
  lineage lookup is never silent (`SKIPPED ...` there and appended to `routing_error`). `matpredict detect audit-scope` also lists any
  record that carries a scoped gene from outside that gene's scope (none now).
* Records `5016_c4_MAT_MAT1-2` and `5016_c5_MAT_MAT1-1` are now `record_version: 2` (material change: two genes added, segment
  moved); **the curator must re-review and sign off version 2**; the existing sign-off is for version 1.
* `tables_v1/` is the v1 comparison made while `base` was still running (2,716 of 2,731); superseded by the top-level tables.
* Wall times are cluster-noisy (the same queries ran 25% slower in one array); no runtime claim is made.
