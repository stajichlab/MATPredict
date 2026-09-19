# Genome-scale detection rollout — findings (Task 6, 2026-09-18)

## Status

The 13-genome pilot rollout was started for real (Tasks 1-5's real tooling,
no mocks) and was **stopped early** after a severe, root-caused performance
regression made completing all 13 genomes impractical (projected ~1 day of
wall-clock for a batch estimated at ~19 minutes). Only genome 1 of 13
(*Coccidioides immitis*, taxid 5501, `GCA_004115165.2`) ran to completion.
Genome 2 (*C. posadasii*, taxid 199306) was ~1.5h into its own run, with no
report written, when the process was killed. The remaining 11 genomes were
never started.

This is reported as a genuine, reproducible pipeline finding — not a
network/hardware flake, and not something waiting longer would have fixed
(see the root cause below, which affects most of the pilot list, not a
transient slowdown on one genome).

## What actually ran

- **Task 1/2 (acquisition)**: `acquire_genomes` re-run for real against all 13
  pilot taxids. **All 13/13 resolved via the local BFD genome library**, zero
  network calls, zero acquisition failures.
- **Task 3 (batch planning)**: `plan_batches(genomes, target_seconds_per_job=5400)`
  produced **1 batch of 13 genomes**, estimated total ~1147s (~19 minutes)
  by the gzip-size-linear estimator. Run directly in-session (not via sbatch)
  per the brief's guidance, since 1 batch well under the 5400s target.
- **Actual run**: genome 1 (5501) took **~2h20m wall clock** (15:23→17:43),
  roughly **87x** the estimator's own 71.5s prediction and roughly **145x**
  the ~96s single-genome reference benchmark this project has on record.
  Genome 2 was still running after ~1.5h with no output. The batch was
  killed at this point rather than waited out.

## Headline finding: `taxonomic_scope` routing is broken for most of the Ascomycota `MAT` family, causing near-total loss of taxonomic routing across the pilot

`family_registry.route()` is supposed to narrow each genome's search to only
the families whose `taxonomic_scope` (declared in `db/<Phylum>/order.yml`)
covers the genome's taxid, directly or via NCBI lineage ancestry. When no
family's scope matches, `route()` deliberately falls back to searching
**every** curated family (documented, intentional fallback behavior in
`route()` itself — see its docstring and final line). That fallback is what
is firing, pervasively, across this pilot.

**Live-reproduced, with real code, real data, no speculation:**

```
$ route(taxid, load_all_families(Path("db")))   # 19 total curated families
5501   (C. immitis)        -> 19/19 families  <-- EXHAUSTIVE FALLBACK
199306 (C. posadasii)      -> 19/19 families  <-- EXHAUSTIVE FALLBACK
162425 (A. nidulans)       -> 19/19 families  <-- EXHAUSTIVE FALLBACK
746128 (A. fumigatus)      -> 19/19 families  <-- EXHAUSTIVE FALLBACK
5061   (A. niger)          -> 19/19 families  <-- EXHAUSTIVE FALLBACK
5059   (A. flavus)         -> 19/19 families  <-- EXHAUSTIVE FALLBACK
5076   (P. chrysogenum)    -> 19/19 families  <-- EXHAUSTIVE FALLBACK
27334  (P. expansum)       -> 19/19 families  <-- EXHAUSTIVE FALLBACK
36651  (P. digitatum)      -> 19/19 families  <-- EXHAUSTIVE FALLBACK
5141   (N. crassa)         -> 1/19 families   (correctly routed: Ascomycota:MAT)
5518   (F. graminearum)    -> 19/19 families  <-- EXHAUSTIVE FALLBACK
5507   (F. oxysporum)      -> 19/19 families  <-- EXHAUSTIVE FALLBACK
510951 (N. tetrasperma)    -> 1/19 families   (correctly routed: Ascomycota:MAT)
```

**11 of 13 pilot genomes (85%) fall into the exhaustive fallback** — searched
against all 19 curated families instead of the 1 actually relevant to them.
Only the two *Neurospora* genomes route correctly.

**Root cause, verified independently (own live checks, not inherited from
the coordinator's initial report):** `db/Ascomycota/order.yml`'s `MAT` locus
declares:

```yaml
taxonomic_scope: [222544, 5180, 28548, 40559]
```

All four of these taxids resolve (confirmed via `default_lineage_taxids`) to
species in **Helotiales** (the same order as this project's `28548`/`40559`/
`5180` curated Helotiales `MAT` records). Meanwhile, the `Ascomycota:MAT`
family's own curated records (`db.glob("*/*/*/metadata.yaml")`, excluding
`candidates/`) span **at least 6 different orders**: Helotiales, Onygenales
(*Coccidioides*), Eurotiales (*Aspergillus*), Hypocreales (*Fusarium*,
*Metarhizium*), Sordariales (*Neurospora*), and others. A systematic audit
(scripted, checking every family's own curated records' taxids against that
family's declared `taxonomic_scope`, by both direct membership and live
NCBI-lineage-ancestor overlap) found:

- **17 of 22 curated `Ascomycota:MAT` records (77%) are NOT covered by their
  own family's `taxonomic_scope`** — including every ground-truth-tier
  record (5501 x2, 199306 x2, 162425, 746128 x2) and several blind-tier-
  relevant ones (5518, 5127 x2, 27339, 107463, 115814, 1301501 x2, 2903220,
  2903222).
- Every **other** family in all three `order.yml` files (`PM`, `mat1`,
  `mat2`, `mat3`, `MATsc`, `MATyl`, `MTL`, `MATtub` in Ascomycota; all 9
  Basidiomycota loci; Mucoromycota's `MAT`) checked out fine — every one of
  their own curated records' taxids is covered, directly or via lineage, by
  that family's declared scope. **This is specific to the `Ascomycota:MAT`
  locus's `taxonomic_scope` value, not a systemic problem with `route()`
  itself or with every family's curation.**

This is the exact, previously-flagged, unswept risk from `docs/HANDOFF.md`:
> "Lineage-aware `taxonomic_scope` routing... live-verified on 8 of the ~50
> previously mis-routed taxids, not all of them... hasn't been swept."

This rollout is that sweep, and it found a large, reproducible instance.

### Runtime consequence

Until this is fixed, per-genome detection cost for any taxid not well-covered
by `order.yml`'s scope data is dominated by the exhaustive-fallback effect —
searching all 19 families (each against every gene in every curated record
for that family) instead of the 1-4 actually relevant. This makes the
project's existing ~96s/39Mb reference timing figure a **massive
underestimate** for 11 of 13 pilot genomes, and by extension for most
Eurotiomycetes/Hypocreales/Onygenales taxa in general — not a one-off. The
batch-runner's own `estimate_seconds_from_gzipped_size` (Task 3) has no way
to account for this, since it estimates from file size only, not family-count
routing outcome — it will keep underestimating every affected genome's real
runtime until the scope data is fixed.

**Recommendation: re-run the full 13-genome rollout as a fast-follow only
after `taxonomic_scope` data across all three `order.yml` files is corrected
and re-swept with the same live-lineage audit used here** — not attempted
again as-is; a second attempt without the fix will reproduce the same
multi-hour-per-genome cost.

## Genome 1's real, complete result (5501, *C. immitis*, `GCA_004115165.2`)

- `families_attempted`: 19 (all — the exhaustive-fallback effect above)
- `detected`: 67 candidate-locus entries
- Confidence tally across those 67: **57 low, 10 medium, 0 high**
- Per-family hit counts (all spurious except `Ascomycota:MAT`, the genome's
  real family): `Basidiomycota:bLocus` 14, `Ascomycota:mat2` 13,
  `Basidiomycota:HD` 11, `Basidiomycota:Aalpha` 5, `Ascomycota:PM` 4,
  `Ascomycota:MATyl` 4, `Basidiomycota:PR` 3, `Basidiomycota:MAT` 3,
  `Ascomycota:mat1` 3, `Ascomycota:MATsc` 3, `Mucoromycota:MAT` 2,
  `Basidiomycota:Bbeta` 1, **`Ascomycota:MAT` 1** (the real family).
- `not_detected` (6 entries, all real "no hits"/"below ambiguity floor"
  outcomes, not errors): `Ascomycota:mat3`, `Ascomycota:MTL`,
  `Ascomycota:MATtub`, `Basidiomycota:Abeta` (no reference-protein hits);
  `Basidiomycota:Balpha` (best cluster 0.25 of expected genes, below 0.50
  floor); `Basidiomycota:aLocus` (0.33, below floor).

**The real (`Ascomycota:MAT`) family's own result is weak**: only 1 gene
(`APN2`, a flanking-conserved gene, not a core idiomorph gene) was found, at
low confidence, matched against reference record `578113_sxlc146_MAT_MAT1-2`
— a **Diaporthales** record, not one of this genome's own curated Onygenales
records (`5501_h538-4_MAT_MAT1-1`, `5501_rs_MAT_MAT1-2`). None of the core
idiomorph genes (`MAT1-1-1/2/3/5`, `MAT1-2-1/10/4`, `SLA2`, `COX13`,
`CIMG_00407`, `matA-1/2/3`, `mt a-1`) were found. This is a second, separate
signal (not just the scope-routing noise) that real detection sensitivity
for this genome's true locus is currently poor — worth its own follow-up
investigation, independent of the routing bug (see triage list).

The other 66 detections are, with essentially total confidence, cross-phylum
noise: e.g. Basidiomycota `bLocus`/`HD`/`PR`/`Aalpha` genes (real only in
Basidiomycota) and yeast-specific `MATsc` genes (`MATALPHA1/2`, `alpha1/2`)
matching an Ascomycete filamentous fungus genome at 22-51% identity with no
biological basis — a direct downstream symptom of the routing bug, not
independent evidence of anything real.

## Ground-truth sanity scoring (Task 5 wiring), real result

`match_ground_truth("5501_GCA_004115165.2", db_root)` returned **2 matches,
both `status="ambiguous"`, 0 `"exact"`**:

- vs `5501_h538-4_MAT_MAT1-1`: same taxid, but rollout genome accession
  (`GCA_004115165.2`) ≠ curated source accession (`EF472259.1`), and rollout
  strain (`WA_211`, from the BFD manifest) ≠ curated strain (`H538.4`).
- vs `5501_rs_MAT_MAT1-2`: same taxid, accession ≠ `NW_004504310.1`, strain
  `WA_211` ≠ `RS`.

Since `score_self_consistency` only scores `"exact"` matches (by design —
see `benchmark.py`'s own documented judgment call), and there are zero exact
matches for this genome, **the real, honest numeric sensitivity result is:
0 `FamilyBenchmark` entries scored.** This is the exact, already-known
outcome from Task 5's development (different isolate/strain than every
curated Coccidioides record) — re-confirmed live here, not a new problem and
not evidence the pipeline is broken. No other ground-truth-tier genome
(199306, 162425, 746128) completed, so no additional ground-truth scores
exist from this run.

## Aggregation (Task 4), real result

`aggregate_reports` / `matpredict detect rollout-summary` run over the
`pilot-rollout-out/` directory (2 genome directories that exist: 1 complete,
1 empty from the killed run):

- `total_genomes`: 2 (both attempted genome directories, per Task 4's
  documented "count every attempted genome" contract)
- `genome_errors`: 1 — `199306_GCA_020976775.1` (`detection_report.yaml`
  does not exist — the process was killed before this genome finished; this
  is Task 3's known, documented "empty directory on failure/interruption"
  gap working exactly as designed, not a new bug)
- `anomalies`: 0 (anomaly detection needs ≥2 genomes with completed reports
  to compare across a shared taxonomic group; only 1 genome completed, so
  there is nothing to compare against — an honest "not enough data," not a
  clean bill of health)
- `not_detected`: 6 entries (all from genome 5501, listed above)

## Triage list (real findings for follow-up — none fixed here, per task constraints)

1. **[HEADLINE, high priority] `Ascomycota:MAT`'s `taxonomic_scope` in
   `db/Ascomycota/order.yml` is stale/incomplete**: `[222544, 5180, 28548,
   40559]` covers only Helotiales-adjacent taxa and misses Onygenales,
   Eurotiales, Hypocreales, and others that this same family's own 22
   curated records actually span (17/22, 77%, fall outside it). Real,
   reproducible, live-verified. **Fix should replace this with a scope that
   actually covers Pezizomycotina (or whatever the correct common ancestor
   of this family's real curated diversity is)**, then be re-verified with
   the same live-lineage-overlap audit used here, extended to **every**
   curated record in **every** family across all three `order.yml` files
   (this pilot's audit script already does this check cheaply — no genome
   download or detection run required, just `default_lineage_taxids` calls
   against existing metadata — and found no other family with this problem,
   but a full DB-wide sweep beyond just the pilot's 13 taxids has not been
   done and should be, since new curated records get added over time and
   could silently reintroduce this).
2. **[Medium priority] Real detection sensitivity for genome 5501's actual
   MAT locus looks weak, independent of the routing bug**: only 1 of ~10+
   core idiomorph genes was found, matched against the wrong-order reference
   record (Diaporthales, not this genome's own curated Onygenales records),
   at low confidence. Worth a dedicated investigation once the routing bug
   is fixed and the noise from irrelevant families is gone — is this a real
   sensitivity problem for Onygenales-family detection, or an artifact of
   this specific highly-diverged genome/strain (`WA_211`, not a strain any
   curated record represents)?
3. **[Informational] The batch-runner's runtime estimator
   (`estimate_seconds_from_gzipped_size`) has no way to account for
   scope-routing outcome**, only file size. Once the scope bug above is
   fixed, this estimator should be revisited with real, multi-genome timing
   data from a corrected run (as the estimator's own docstring already
   anticipates) — the current single 96s/39Mb data point plus this pilot's
   badly-inflated numbers are not a usable calibration set together.
4. **[Informational, not actionable without more genomes] Zero anomalies
   detected** is an artifact of only 1 genome completing, not a real "clean"
   result — re-run once more genomes complete.

## What was not attempted

- Genomes 2-13's detection runs (162425, 746128, 5061, 5059, 5076, 27334,
  36651, 5141, 5518, 5507, 510951) — not run, per the decision to stop early
  once the root cause was identified rather than wait out a ~1-day batch.
- Ground-truth sanity scoring for 199306, 162425, 746128 — blocked on their
  detection runs not having completed.
- No source code was modified as part of this task, per the task's explicit
  constraint — the `taxonomic_scope` fix and any detection-sensitivity
  investigation are follow-up work.

## Artifacts

- Real per-genome detection report (only completed genome): repo-relative
  `.superpowers/sdd/2026-09-18-genome-scale-detection-rollout/pilot-rollout-out/5501_GCA_004115165.2/detection_report.yaml`
  (gitignored working directory, not committed)
- Rollout summary YAML:
  `.superpowers/sdd/2026-09-18-genome-scale-detection-rollout/rollout_summary.yaml`
  (gitignored, not committed)
- This findings document (committed): `docs/superpowers/plans/2026-09-18-genome-scale-detection-rollout-findings.md`
