# Why three Parastagonospora genomes have no call (2026-10-17)

Analysis only. `src/` and `db/` are untouched; the proposed diff in section 6 is a sketch and has not been run.
Code: origin/main 624048d, run on the HPCC (pixi env mat-prof, checkout /bigdata/stajichlab/jstajich/projects/parmiss, jobs `parmiss-*`).
Genomes: GCA_003497985.1 and GCA_003503065.1 (P. nodorum, taxid 13684), GCA_022559985.1 (P. avenae, taxid 2814595).
Baseline `detect --taxid T --evidence-diagnostics ...` reproduces the campaign: 0 loci, the same `not_detected` text, for all three.

## 1. Short answer

The evidence floor is a cause, but not the only one. The true gene is a lone, near-complete core gene that is
either on a short contig or far from any strong flank. It is dropped at three successive gates, and removing
only the first one does not recover a call:

| gate | where | what it does to the true cluster |
|---|---|---|
| 1. evidence floor | `_families_meeting_evidence_floor` (pipeline.py ~945) | needs 2 distinct non-superseded genes; the diagnostics row says `gene_count: 2` because MAT1-1-3 also hit, but that hit is `superseded_by` the winner and is not counted, so 1 gene: `admitted: False`, never polished |
| 2. modelled-gene bar | `MIN_POLISHED_GENES = 2` (`_withhold`) | with the floor relaxed the cluster is polished (1 modelled gene) and withheld as `modelled_gene_bar` |
| 3. fraction floor | `ambiguity_floor = 0.5` | after the bar is also relaxed, 1 of 4 roster genes = 0.25 < 0.5: `below_fraction_floor` (the relaxed pass rescued it only in the genome that had no strict call) |

Not the cause: tblastn localization (all three clusters exist), clustering window (the cluster is the whole gene),
the polish cap (`--max-polished-clusters-per-family 0` changes nothing: the true cluster was never admitted, and
once admitted it is polished anyway by the strong-identity tier), the MAT-gene gate (`suppressed_mat_gene_gate: 0`
in every run), gene-model building (miniprot/exonerate built the model, `polished_genes: 1`).

## 2. Per genome

### GCA_003497985.1 (P. nodorum, expected MAT1-2)
- True hit: MAT1-2-1, 99.39% over 165 aa (own-taxon record SN436GA98 is in the reference set), contig QVWR01000152.1
  (60,679 bp), cluster 57,395-58,482 (1,087 bp, the whole gene), 2.2 kb from the contig end.
- Flank genes are on other contigs and weak: APN2 and COX13 at about 74% on QVWR01000146.1 (start of contig), SLA2 78.6% on QVWR01000006.1.
  So the split-locus rule (core >=95%, <=500 bp from the end, 2 flanks within 10 points of the core identity) fails on edge distance (2.2 kb) and flank identity.
- Baseline: 243 clusters, 36 admitted; this one `admitted: False` (1 gene after supersession).
- `--min-hits 1`: admitted, polished (1 gene), withheld `modelled_gene_bar` and `below_fraction_floor` (0.25 < 0.5). Still 0 calls.
- `--min-hits 1` + bar 1 only: reports QVWR01000056.1 MAT1-2-1 + APN2 at 55.6% (a paralog) as `partial_locus`, medium. Right idiomorph by luck, wrong locus.
- `--min-hits 1` + bar 1 + fraction floor 0.2: reports the true locus (MAT1-2, `idiomorph_gene_only`, medium, 99.4%) plus the 55.6% paralog as `mat_locus`, **high** (a false high-confidence call).

### GCA_003503065.1 (P. nodorum, expected MAT1-2)
- Same hit, 98.9-99.4%, contig QVWJ01004591.1 (3,621 bp), cluster 318-1,405, 318 bp from the contig start. APN2/COX13 at about 74% on QVWJ01000028.1, SLA2 74% on QVWJ01000221.1.
- Edge distance passes (318 < 500) but flank identity (74% vs the 89.4% needed) does not, so no split-locus call.
- Gates 1, 2 and 3 as above. `--min-hits 1` + bar 1 recovers it through the relaxed pass (MAT1-2, medium) with a second junk locus (QVWJ01001800.1, 53% MAT1-2-1, 361 bp).

### GCA_022559985.1 (P. avenae, expected MAT1-1)
- True hit: MAT1-1-1 85.6% over 271 aa (SN435PL98 is a P. nodorum record), JAFJTR010000080.1 (135 kb), cluster 3,562-4,681, 3.5 kb from the contig start. The cluster also carries a superseded MAT1-2-1 hit.
- Flanks: APN2 70% and COX13 73% on JAFJTR010000071.1, SLA2 70.6% on JAFJTR010000083.1; weaker than in nodorum.
- Baseline `admitted: False`. 85.6% is below the split-locus 95%, so even that rule would not fire.
- Only `--min-hits 1` + bar 1 + fraction floor 0.2 recovers it (MAT1-1, medium, `idiomorph_gene_only`). **The same run also reports a wrong-idiomorph call**: JAFJTR010000003.1 MAT1-2-1 + APN2 at 55.6%, `mat_locus`, high.

### GCA_003503165.1 (P. avenae, N50 5 kb): assembly gap
Not re-run. The campaign report has 5 suppressed loci, all on 0.3-3 kb contigs with 0 modelled genes. Assembly gap, nothing to fix.

## 3. Relaxations tried (3 genomes x each; all 12 + 9 jobs under 1.1 min)

| setting | GCA_003497985.1 | GCA_003503065.1 | GCA_022559985.1 |
|---|---|---|---|
| baseline | none | none | none |
| `--min-hits 1` | none | none | none |
| `--max-polished-clusters-per-family 0` | none | none | none |
| both of the above | none | none | none |
| `--min-hits 1` + bar 1 (`min_polished_genes=1`) | wrong locus (55%), MAT1-2 | none | wrong locus (55%), MAT1-2 (wrong idiomorph) |
| bar 1 only | wrong locus (55%), MAT1-2 | none | wrong locus (55%), MAT1-2 (wrong idiomorph) |
| `--min-hits 1` + bar 1 + fraction floor 0.2 | true locus MAT1-2 + junk high | true locus MAT1-2 + junk | true locus MAT1-1 + junk high MAT1-2 |

`min_polished_genes` and `ambiguity_floor` are not CLI options; they were set with a wrapper
(`scratchpad/mp_relax.py`, monkeypatching `run_pipeline`; not committed). Admitting with `--min-hits 1` raises admitted clusters
from 36/13/37 to 172/162/224 per genome, but the polish cap still holds polishing to about 6 + strong-identity clusters,
so the compute cost on these genomes is about 1 minute either way (wall time only; CPU not measured).

The relaxations alone are not an acceptable fix: with the bars off, 55% HMG paralogs become calls, one of them high confidence with the wrong idiomorph.
What separates the true gene is identity: 99.4, 98.9, 85.6% against 40-56% for paralogs.

## 4. Is this specific to three genomes? No.

Records-run Dothideomycetes reports (reports_all.tar.zst, 2,728 reports), Pleosporales (order from results/2026-10-03_ascomycota_v060/genomes.tsv):
1,012 genomes, 875 called, **137 uncalled**. In 131 of the 137 at least one cluster was admitted and polished (reason text: 32 "0 modelled genes",
50 "1 modelled gene", 49 flank-carried withheld), 6 had no suppressed loci. So the evidence floor is not what leaves these genomes uncalled in general; the bar and flank-carried rule are.

Panel (seeded random, `scratchpad/panel.tsv`): 24 uncalled + 24 called Pleosporales, run with the permissive setting (min-hits 1, bar 1, fraction floor 0.2),
then filtered offline by the proposed rule (a single-gene locus is kept only if the cluster's best identity >= T; multi-gene loci kept as they are):

| T | uncalled genomes recovered (of 24) | called genomes with an extra locus (of 24) |
|---|---|---|
| 95 | 3 | 0 |
| 90 | 8 | 0 |
| 80 | 14 | 0 |
| 70 | n/a | 0 |
| 50 | n/a | 8 genomes, 9 loci (mostly MAT1-2 where the call is MAT1-1) |
| 0 (nothing but the permissive run) | n/a | 10 genomes, 11 loci |

At T=80 one recovered genome (GCA_035078975.1) gets both MAT1-2 (84.6%) and MAT1-1 (94.5%) single-gene loci, which needs curation; at T=90 it gets only the 94.5% MAT1-1 locus.
Of the 24 uncalled, 10 have no single core cluster >=80% at all (diverged or absent), so about 40% are unrecoverable by identity.
Two of the "recovered" genomes (GCA_059628625.1, GCA_052323635.1) are multi-gene (87.7 and 88.6%): they come from the lowered fraction floor / min-hits 1, not from the single-gene rule.
Extrapolating the T=90 rate (8/24) to 137 uncalled gives roughly 45 genomes; this is a projection from a 24-genome sample, not a measurement.

The earlier study docs/notes/2026-09-24_saccharomyces-cassette-diagnosis.md reached the same shape for yeast: 90 single-gene >=90% clusters
fail `min_hits=2`, "at most 11" more fail the bar; its fixes 2 and 4. The 2026-09-20 Absidia note is the reason for the floor
(HMG paralogs at 30-50%), which an identity threshold respects. Studies read: analysis/2026-09-26_polish-cap.md (cap is not involved here),
2026-09-27_confidence-rules.md (fallback identity does not separate good from bad fallback calls, so the confidence must not rest on identity alone; here it is a gate, not a tier),
2026-09-28_mat-gene-gate-and-validation.md (the gate types a gene; it never ran for these single-gene clusters). I did not re-run their replays.

## 5. What would be checked before adopting

- Sequence check that the 99.4% gene is a full-length MAT1-2-1 (the tBLASTn HSP is 165 aa and the cluster 1,087 bp; not translated or inspected).
- Whether the weak flank (about 74%) on the neighbouring contigs is the true flank of the locus: not verified. If so, adding P. nodorum APN2/COX13/SLA2 flank records to db would let the split-locus rule fire without code change (a `db/` curation option).
- The sample is Pleosporales only; no other order or phylum was run. N50s in the brief were not re-checked.
- The brief's check note docs/notes/2026-10-17_parastagonospora-check.md is not on origin/parastagonospora-check; the numbers above were reproduced from my own runs and results/2026-10-05_dothideomycetes_full/ reports.

## 6. Proposed diff (sketch, untested)

Add a rule modelled on `split_locus`, applied only to a family that is still unreported after the strict, relaxed and split-locus passes, so a genome that already works is not changed. Parameters as an `EvidenceFloor`-style frozen dataclass with `None` meaning off, default off until the curator rules.

```diff
--- a/src/MATPredict/detect/pipeline.py
+++ b/src/MATPredict/detect/pipeline.py
@@ class EvidenceFloor:
     min_hits: int = 2
     min_identity: float | None = None
     require_core_role: bool = True
+    #: A cluster with ONE distinct core gene still reaches polish when that gene's best
+    #: identity is at least this (percent). None = off (current behaviour).
+    strong_single_identity: float | None = None
@@ def _families_meeting_evidence_floor(...):
         distinct_genes = {h.gene_name for h in own_hits}
-        if len(distinct_genes) < floor.min_hits:
-            continue
+        if len(distinct_genes) < floor.min_hits:
+            strong = (floor.strong_single_identity is not None and own_hits
+                      and any(h.role == "core_MAT" and h.identity is not None
+                              and h.identity >= floor.strong_single_identity for h in own_hits))
+            if not strong:
+                continue
@@ after the split-locus block (results = results + split_calls ...)
+    # Strong single core gene: one MODELLED core gene at >= strong_single_identity, the genome's
+    # best hit for that gene, in a family with no report. Reported as idiomorph_gene_only,
+    # confidence low, flagged `strong_single_core`; skips the 2-gene bar and the fraction floor.
+    strong_calls = _strong_single_core_calls(results, clusters, families, _build,
+                                             identity=evidence_floor.strong_single_identity)
+    strong_calls, strong_paralogs = apply_mat_gene_gate(strong_calls, classifier_specs)
```

`_strong_single_core_calls` would mirror `_split_locus_calls`: for each cluster of an unreported family, require `polished_genes >= 1` with the core gene modelled,
the cluster's best core identity >= the parameter, and the gene the genome's best bitscore hit for that gene; drop a call that overlaps a reported locus;
if two single-gene loci of different idiomorphs qualify in one genome, report neither as a call and list both (the GCA_035078975.1 case).
CLI: `--strong-single-identity` (default None). Suggested value for review: 90 (8/24 recovered, 0/24 collateral, no sign of the 84.6% paralog-like second idiomorph).
At 85 the third genome (GCA_022559985.1, 85.6%) would also be recovered; at 90 it is not.

## 7. Questions for the curator

1. Is a lone core gene at >=90% identity, with no flank and no second modelled gene, a reportable locus (class `idiomorph_gene_only`, confidence `low`), or should it stay a diagnostic only?
2. Which threshold: 95 (the split-locus ruling, recovers 3/24), 90 (8/24), or 85 (adds GCA_022559985.1)? The 85.6% P. avenae hit is a different species from its nearest record.
3. Should a genome with two qualifying single-gene loci of opposite idiomorph (GCA_035078975.1) be reported as two calls, as undetermined, or not at all?
4. Prefer curation to code for the Parastagonospora cases: add P. nodorum APN2/COX13/SLA2 to the roster so the split-locus rule (95% / 500 bp / flanks within 10 points) can work, and widen `max_edge_bp` (2.2 kb here) or flank gap? Note the flanks here are about 74%, not within 10 points of 99%.
5. Do the 49 Pleosporales uncalled genomes withheld by the flank-carried rule, and the 32 with 0 modelled genes, deserve the same investigation? I did not look at them.
