# Plan: make `detect` targeted — phylum routing, per-locus cluster gap, real polish gate

## Spec / authority

Curator rulings by J. Stajich on 2026-09-19/20, recorded verbatim:

1. **Phylum routing.** "it is reasonable to do a test by a phyla - so the run
   will know to search one of Ascomycota, Basidiomycota, Mucoromycota".
   When no family's `taxonomic_scope` matches, **fall back to the query taxon's
   phylum**, and record in the report that the phylum fallback was used.
2. **Zygo scope.** "for the current zygo testing - let's just focus on this
   running the Mucoromycota db against the candidate genomes."
3. **Cluster distance.** "for ascomycetes the distance is fine, for
   mucoromyocta lets start with 50kb and we will learn empirically after we
   have built up a few more examples."
4. **Polish trigger.** "polish trigger could be another gene not necessarily
   flank - at least at start until we have a bigger collection of examples of
   fragments." i.e. >=2 DISTINCT genes with >=1 `core_MAT`; a flanking partner
   is NOT required.

## Context (read from the real code on 2026-09-19, not assumed)

- `route()` (`detect/family_registry.py:113`) tries direct `taxonomic_scope`
  membership, then NCBI lineage ancestors, then `return list(families)` --
  an exhaustive fallback over EVERY family in EVERY phylum.
- This fires on a common case, not an exotic one. Measured live: taxid 294748
  resolves its lineage fine (`... 4890 Ascomycota, 716545 saccharomyceta,
  147537 Saccharomycotina ...`) but **no family scope intersects it**
  (Ascomycota `MAT` is scoped to 716546 leotiomyceta, `MATsc` to 4893, `MTL`
  to 4959), so it routed to **19 of 19 families across all 3 phyla, 91 genes**.
- `build_reference_fasta(db_root, out_path)` (`detect/reference_fasta.py:14`)
  takes NO family argument -- it globs every accepted `proteins.faa`. The
  live run's query set was **181 proteins, the whole DB**, regardless of
  routing. Routing today affects only attribution and polishing, never the
  tblastn query set.
- Per-phylum scale, measured: Ascomycota 9 families / 47 genes / 122 proteins;
  Basidiomycota 9 / 37 / 40; Mucoromycota 1 / 7 / 19. Total 19 / 91 / 181.
  So a Mucoromycota-only run is a **9.5x smaller query set**.
- `cluster_hits` (`detect/clustering.py:17`) groups **by contig only**, not by
  family, with one global `max_gap=25_000`. Because routing will be
  phylum-restricted, one run sees one phylum, so a per-run gap suffices --
  `cluster_hits`'s shape does not need to change.
- `EvidenceFloor` (`detect/pipeline.py:285`) already gates polishing per
  (cluster, family) BEFORE Stage 2, and `min_hits` already counts DISTINCT
  genes by `gene_name`, never raw HSPs. Defaults are `min_hits=1`,
  `require_core_role=False`, `min_identity=None` -- i.e. one HSP of any gene,
  any role, any identity admits a family to the two-tool polish loop. The
  docstring states these defaults are deliberate and deferred "until
  db/Ascomycota/order.yml's taxonomic_scope fix ... lets most genomes run with
  correctly-narrowed routing instead of the exhaustive fallback". That fix is
  this plan.
- Real cost of the status quo: a 14.7 MB single-genome run was killed after
  ~24 minutes without completing, serially polishing spurious cross-phylum
  candidates.

## Global Constraints

- **Never silently widen scope.** Any fallback (phylum, or exhaustive) MUST be
  recorded in the detection report so a coverage gap surfaces in output rather
  than being absorbed into noisy results.
- **Never fabricate a phylum.** If the lineage cannot be resolved, do not guess
  -- fall through to today's exhaustive behaviour and say so in the report.
- **One bad genome must not sink a batch** -- `run_batch`'s per-genome
  try/except discipline is binding.
- Tests run with `pixi run -e test pytest -v`, never bare `pixi run pytest`.
- Do NOT change detection semantics beyond scope/gating: no change to
  clustering's contig-only grouping, to `_splice_transcript`'s ASCENDING
  detection-schema exon convention, or to polish tool invocation.
- Duplicate gene names within a family are normal biology.

## Task 1 -- phylum routing and a routing-aware reference FASTA

1. **`--phylum` override** on `matpredict detect` (choices limited to the
   phyla present under `db_root`, discovered at runtime, not hardcoded). When
   given, it restricts routed families to that phylum outright and skips taxid
   routing. This is what the Zygo run will use.
2. **Phylum fallback in `route()`.** When neither direct nor lineage matching
   yields a family, resolve the query taxon's phylum from the lineage already
   fetched and return that phylum's families. Only if the phylum cannot be
   determined does it fall through to today's `list(families)`.
   `route()` must report WHICH path it took -- return it alongside the
   families, or expose a small result object; do not make callers re-derive it.
3. **Narrow the query set.** `build_reference_fasta` gains an optional family
   (or phylum) filter so the tblastn query set contains only the routed
   families' curated proteins. Default (no filter) behaviour unchanged.
   Verify with a real count that a Mucoromycota-restricted run queries 19
   proteins, not 181.
4. **Record the routing decision in the report** (`write_detection_report`):
   which families were searched, and whether routing was direct, lineage,
   phylum-fallback, or exhaustive.

## Task 2 -- per-locus cluster gap

1. Add an optional `max_cluster_gap_bp` per locus to
   `db/_schema/order.schema.yaml`, defaulting to 25000 when absent.
2. Set it to **50000** on the Mucoromycota `MAT` locus, with a comment
   recording the curator ruling and that it is a starting value to be revised
   empirically.
3. `run_pipeline` derives the run's `max_gap` as the **MAXIMUM** over routed
   families' values. Rationale, and it must be in the code comment:
   under-splitting is recoverable because the evidence floor and polishing
   still discriminate within a cluster, whereas over-splitting silently
   destroys a real locus by cutting it in two.

## Task 3 -- flip the evidence-floor defaults

1. `EvidenceFloor` defaults become `min_hits=2`, `require_core_role=True`.
   `min_identity` stays `None`.
2. Update the docstring: the deferral it describes is now resolved by this
   plan's routing fix, and the new defaults are the curator's ruling -- >=2
   distinct genes with >=1 core_MAT, a flanking partner deliberately NOT
   required until more fragmented-locus examples exist.
3. The CLI flags `--min-hits`/`--require-core-role` keep working and still
   override.
4. Report the real effect on at least one genome: how many (cluster, family)
   pairs are admitted under the old vs new defaults.

## Task 4 -- make the batch path pay the same narrowed cost

Found by Task 1's implementer, and it blocks the stated goal: Tasks 1-3 make the
SINGLE-GENOME CLI path targeted, but the 813-genome Zygo rollout goes through
`run_batch`, which **builds ONE unrestricted reference FASTA for the whole batch
and self-routes per genome**. So a rollout still queries all 181 curated proteins
per genome and none of Task 1's narrowing applies to it.

This conflicts with `run_batch`'s deliberate build-once design, so it needs a real
decision, not a patch.

### Context (measured)

- `src/MATPredict/detect/batch_runner.py` builds the reference FASTA once before
  the per-genome loop, by design, to avoid rebuilding it 813 times.
- `scripts/run_detection_batch.py` is the SLURM entry point and exposes no way to
  restrict scope.
- Task 3's implementer additionally found `run_batch` has **no `evidence_floor`
  parameter at all**, so every batch rollout silently inherited the new stricter
  defaults with no opt-out.
- Per-phylum query sets, measured: Ascomycota 122 proteins, Basidiomycota 40,
  Mucoromycota 19, unrestricted 181.

### Requirements

1. **`run_batch` gains an explicit scope.** A `phylum: str | None = None`
   parameter that, when set, restricts routing AND the reference FASTA for every
   genome in the batch -- built ONCE for that phylum, preserving the build-once
   design rather than rebuilding per genome. Default `None` keeps today's
   behaviour.
2. **`run_batch` gains `evidence_floor: EvidenceFloor | None = None`**, threaded
   to `run_pipeline`, so a batch can opt out of or tighten the new defaults.
   Default `None` means "use `run_pipeline`'s default", NOT a hardcoded floor.
3. **`scripts/run_detection_batch.py` exposes both** as CLI flags, named to match
   `matpredict detect`'s (`--phylum`, `--min-hits`, `--require-core-role` /
   `--no-require-core-role`).
4. **A mixed-phylum batch must still work.** If `phylum` is None, per-genome
   self-routing is unchanged. Do not make the batch path require a phylum.
5. Report the measured query-protein count for a `--phylum Mucoromycota` batch
   (expect 19) versus unrestricted (expect 181).

### Constraints

- `run_batch`'s per-genome try/except isolation is binding -- one bad genome must
  never sink the batch.
- The reference FASTA must still be built ONCE per batch, not per genome.
- Do not change routing (Task 1), gap derivation (Task 2), or `EvidenceFloor`'s
  defaults (Task 3).
