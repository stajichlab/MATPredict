# Why the C. cinerea B43 record changes 164 withheld loci and a called HD span
Status: open (cause found; the code question needs a curator ruling)

## Question
Branch `curate-cinerea-b43` (record `5346_a43-b43-okayama-7_PR_B43` v2, commit 6983433) left calls unchanged but its regression (`results/2026-10-06_regression_cinerea_b43/`) changed 164 withheld loci in the Basidiomycota panel and the span of the called HD locus in T48-F (GCA_982397435.1). Why, and is it a defect?

## Data and code version
- Baseline: frozen worktree `run-2d1860a` (code and database). Candidate: `run-6983433` (adds only the record change, the `order.yml` aliases and a 5-line `gff_export.py` edit).
- Three genomes chosen from the most-changed in the regression: Mixia osmundae (33 changed loci), Dacryopinax primogenitus (18), Leucosporidium creatinivorum (18). 850 loci in the baseline reports. Jobs 29617229-31, 29618151-52, 29618482 (exfab), each through `scripts/run_clade_panel.slurm`.
- Variants (new frozen-style worktrees off 6983433, `proteins.faa` of the record edited, nothing else): `noNew` (4 new peptides removed), `oldSeq6` (gene 6 given back its old sequence), `both` (both).

## Method
1. Read the diff of the record. Sequence comparison old versus new: the 5 receptor proteins are unchanged; gene 6 changed from 123 aa (`MDASVSTPIYP...SSIPSL`) to a different 72-aa protein (`MSDLFASLD...CTIA`, named `phb3.1`; the old entry was wrong according to the curation note); four peptides are new (`phb1.1` 51 aa, `phb2.3` 53 aa, `phb2.1` 61 aa, `phb3.3` 85 aa). The other edits are names and aliases.
2. Controls and variants run on the same three genomes; `cmp.py` compares every locus (called and withheld) by family and contig: span, withheld reason, genes found.
3. For two changed loci, BLAST (tblastn, `-seg no`, e-value 10, the settings `detect` uses) of old and new queries against the genome around the changed boundary.

## Results
Loci that differ from the baseline (`comparisons.txt`; 850 loci in the baseline):
| Run versus baseline | identical | span changed | gene set changed | reason changed | loci only in one run |
|---|---|---|---|---|---|
| baseline repeated (same code, same database) | 850 | 0 | 0 | 0 | 0 |
| candidate | 789 | 56 | 11 | 1 | 3 |
| candidate without the 4 new peptides | 814 | 32 | 5 | 1 | 1 |
| candidate with the old gene-6 sequence | 822 | 27 | 7 | 0 | 2 |
| candidate with both reverted (aliases and names kept) | 850 | 0 | 0 | 0 | 0 |
- The change is not noise: two baseline runs agree on all 850 loci.
- Both sequence changes contribute (removing either leaves part of the effect); reverting both gives the baseline exactly. So the alias and rename handling changes nothing.
- Mechanism, traced for two loci. `detect` runs tblastn at BLAST's default e-value (10) and clusters hits of all routed families together (`cluster_hits` over every hit; the code notes that "a cluster routinely mixes families"). A very weak hit from a short precursor query inside the 25-kb gap of another family's locus joins the cluster and moves the reported span.
  - Mixia, HD locus on NW_014575563.1: baseline 1,108,035-1,121,335, candidate 1,096,033-1,121,335. The new 72-aa `phb3.1` has a hit at 1,096,033-1,096,128 (32 aa, 40.6% identity, e-value 3.0, 24.6 bits), the new left boundary. With the old gene-6 sequence the span returns to the baseline value.
  - T48-F, called HD locus on CEVXIV010000004.1: baseline 118,043-130,301, candidate 118,043-153,196 (same four genes, same confidence). The new `phb3.3` has a hit at 153,008-153,196 (63 aa, 28.6% identity, e-value 4.1, 25.8 bits), the new right boundary. That is a called locus whose span grows 23 kb from one noise-level hit.
- The other 164 changed loci in the regression are of the same kind in the three genomes I decomposed; I did not trace each of them to a hit. The second Mixia HD span change (2,169,865 to 2,234,630 start) is reproduced by the gene-6 swap but I did not find the individual hit.
- Calls: in the regression, 1 call touched (the C. cinerea CC3 genome, the curated one). No call changed in T48-F (HD high, PR medium before and after); only the HD span.

## What changed in detection
Nothing. This note adds the experiment; no source change.

## Limits
- Three genomes of 33 were run for the decomposition; two loci were traced to individual hits. The mechanism is shown, not counted over all 164.
- The hit coordinates come from my own tblastn runs with detect's settings, not from detect's hit table; the match with the cluster boundaries is exact in both cases.
- Whether the spans of withheld loci matter to anyone is not established; they only appear in reports.

## Update: `core_span` (curator ruled option 3, report only) and what it shows
Built on branch `core-span` (code, tests and results there): each called locus gets `core_span` (extent of its own gene models on its contig) and `beyond_core_bp` (part of the cluster span outside them), in `detection_report.yaml`, on the GFF3 `MAT_locus` line and in the HTML report. Calls and clustering are unchanged.

Results (`results/2026-10-08_core_span_panel/` and `results/2026-10-08_core_span_trace3/` on `core-span`):
- **T48-F**, the case that motivated it: HD span 118,043-153,196, own genes 118,043-130,301, 22,895 bp beyond. This matches the baseline run (span 118,043-130,301) exactly, so for this locus the field recovers the pre-edit extent.
- **The 33-genome Basidiomycota regression panel, old and new database** (new code both times): 33 genomes each, 47 called loci each. 44 of the 45 loci that match one-to-one have identical span and core under both databases. The exception is the curated C. cinerea genome (CC3), where span and core both moved because its own gene models changed. Two loci could not be matched one-to-one (several loci on one contig and family). So this panel contains no case of "span moved, core unchanged"; T48-F is the only one, and it is outside the panel. The earlier 164 changed loci are withheld loci, which carry no gene coordinates in the report, so `core_span` cannot be computed for them.
- **How far spans extend past their genes** (new database, 47 called loci): 37 loci (79%) have `beyond_core_bp` > 0; the median is 948 bp (13% of the span); 23 loci are 1 kb or more beyond, 11 are 10 kb or more (largest 43.5 kb, a PR locus on NCVV01000007.1); 10 loci have at least half of their span outside their own genes. Why each extends was traced for only two loci (below), so these numbers show the extent, not that the extra span is noise: some may be real unmodelled genes.
- **Leucr1 bLocus on MCGR01000028.1** (span 145,579-160,941; own genes 145,636-147,187; 13.8 kb beyond). tblastn of the run's whole reference set over that contig: the stretch 147,188-161,500 holds 52 hits; 50 have e-value 0.088 to 10 and come from a dozen unrelated queries (bW, pheromone, receptor, bE, mfa1, HD2, PAN6 and others). The only two better hits (e-value 5e-5, 29 bits, 30% identity over 80 aa) come from one bE query at 160,622-160,941, and the span ends at 160,941, the end of that hit. No modelled gene lies in the stretch. Whether that bE hit is real is not established.

Open: (a) `core_span` for withheld loci needs gene coordinates added to `suppressed_loci` in the report (a pipeline-output change, not done). (b) Which of the 10 loci with half their span outside their genes are noise and which hold real unmodelled genes.

### Hit tally in the stretches beyond the core (11 loci with 10 kb or more beyond; curator asked, 2026-10-08)
Method (`scripts` in `results/2026-10-08_core_span_panel/tally.py`, job 29637408; files on branch `core-span`): for each locus, the run's own `_reference.faa` was run through tblastn (detect's settings, the genome's genetic code) over the contig window; hits lying wholly inside a stretch outside the core are counted as strong (e-value < 1e-3) or weak (>= 1e-3), and the best hit is listed. Hits that straddle the core boundary belong to the core genes and are not counted. (My first version counted overlapping hits and overstated the strong hits next to the core; that file is kept as `tally_newdb_overlap_criterion.tsv`.) 14 stretches of 1 kb or more.

| Group | Loci | Beyond the core | What lies in the extra stretch |
|---|---|---|---|
| Noise only | GCA_000715385.1 (Rhizoctonia solani) PR (left 16.1 kb), GCA_001542265.1 redHD LNKU01000002.1 (left 14.2 kb, right 3.4 kb), GCA_023273805.1 MAT CP096880.1 (right 16.9 kb), Stehi1 HD (left 11.7 kb) | 4 loci, 62.3 kb | 0 strong hits in all 5 stretches (22 to 421 weak hits each; best e-value 1e-3 to 7e-3). The span edge is a weak hit (e-value 1e-3 to 7e-3). |
| Nearly noise | Leucr1 bLocus (right 13.8 kb) | 13.8 kb | 2 strong hits, both from one bE query at the span end (e-value 4e-6, 30% identity over 80 aa); 419 weak hits. |
| Real hits, same family | PR loci: ASM209295v1 NCVV01000007.1 (right 43.5 kb), Gabo G3 CM035310.1 (left 13.7 kb, right 20.4 kb), Trametes versicolor NW_007360328.1 (right 40.4 kb), Heterobasidion irregulare NW_009258203.1 (left 10.6 kb) | 4 loci, 128.6 kb | 152, 7 and 135, 111 and 178 strong hits; the best are PR receptor queries (e-values 1.5e-55, 9.9e-95, 5e-142, 0; Heterobasidion 74% identity over 291 aa). These look like further receptor genes in a receptor array. |
| Real hits, other family | Balpha JAGVSI010000976.1 (left 17.4 kb, right 20.8 kb); HD NW_006763290.1 (right 23.5 kb) | 2 loci, 61.7 kb | 38, 69 and 35 strong hits, best e-values 1e-24, 1e-21, 5e-15, all from PR (receptor) queries: genuine hits of another family inside this family's span. |

What it means:
- The curator's lean towards trimming is supported for 5 of 11 loci (4 noise-only plus Leucr1), 76 kb of extra span in total, where trimming to the core would remove only weak hits. With the two earlier traced loci (T48-F, Mixia) that is 7 loci with noise-defined spans.
- Trimming everything to the core would be wrong for the other 6 loci: in 4 PR loci the extra stretch holds strong hits of the locus's own family (receptor genes that the gene list does not include), and in 2 loci (Balpha, HD) it holds strong hits of another family. In the 4 PR loci, trimming would cut real receptors out of a receptor array; whether a PR "locus" should include them is a curation question (the 2026-10-06 handoff notes that B-locus receptors sit in one array and that array span does not identify a B locus).
- So the data point to a rule on hit quality, not on the core: build the reported span from hits with e-value below a floor (strong hits), so a stretch carried only by weak hits is dropped and a stretch with real hits is kept. With a floor of 1e-3 the 4 noise-only loci would be trimmed, and Leucr1 would not (its edge hit is 4e-6); a floor of 1e-5 would trim Leucr1 too. I have not tested any floor on the regression panel, so the number is open.

Limits: 11 loci from one 33-genome panel; "strong" is my cut-off (e < 1e-3), and a strong hit shows similarity to a curated protein, not that a gene is real (non-mating STE3 receptors also hit); tblastn at e-value 10 on the run's query set approximates detect's hit table, whose boundary hits matched in every locus where I checked the edge (in all 16 stretches, including two under 1 kb that are not in the table, a tblastn hit lay exactly at the span edge); I did not tally the 36 loci with less than 10 kb beyond.

## Update: `supported_span` prototype, and which review claims were checked (curator asked, 2026-10-08)
An independent review (Opus, read-only) advised against trimming to the core and for a second, report-only span built from the modelled genes plus own-family hits above a bitscore floor. I built that as `supported_span` on branch `core-span` (code, tests, results). Calls, clustering and the cluster span are unchanged. Floor option: `--supported-min-bitscore` (default 39, the flank-carried floor, kept for now).

What was checked, and how:
| Claim from the review | Status |
|---|---|
| `locus_merge` compares each result's own `start`/`end`, so trimming the span before the merge would split merged A/B calls | **Verified by running it.** Two B-locus calls with identical spans merge (1 call); the same calls with disjoint supported spans do not (2 calls), and a 2 kb overlap of 15 kb does not either. `supported_span` is a separate field and leaves `start`/`end` alone, so it cannot do this. |
| Report-only changes cannot change calls | **Verified by 0 differences**: 33 genomes, five floors against the earlier run of the same database (family, contig, span, confidence, idiomorph, class, genes, merged_from, withheld loci), and 0 differences between floors. |
| Reported genes must lie inside the span | **Checked: 160 of 160 inside at every floor.** The first run found 5 outside: a real bug of mine (the supported span was stale after `locus_merge` added genes). Fixed test-first by taking the union of the members' supported spans; re-run. |
| A weak hit from another family can bridge a family's own hits across more than the cluster gap (a membership effect that report-only fields do not touch) | **Measured as a proxy, not proven.** Over the three v0.6.0 campaigns, 290 of 27,822 called loci (1.0%; not merged calls) have their own genes further apart than the family's cluster gap: Ascomycota MAT 90, Basidiomycota PR 83, Ascomycota MTL 40. This cannot say whether the bridge is another family's hit or an own-family hit that was not reported. |
| The polish window, rescue acceptance and the gene-count gate depend on the cluster span or membership | **Not tested.** Only the polish window use (`_padded_window` uses `cluster.start/end`) was confirmed by reading the code. These matter only if clustering or the cluster span is changed; `supported_span` changes neither. |
| An e-value floor depends on genome size, a bitscore floor is the existing precedent | **Precedent confirmed in `flank_carried.py` (curator ruling 2026-09-27); the genome-size dependence was not re-measured here.** |

Results (34 genomes: the 33-genome Basidiomycota panel plus T48-F; 49 called loci; floors 30, 33, 36, 39 and 50 bits; `results/2026-10-08_supported_span/`):
- **Sweep.** Loci with any span beyond the supported extent: 17, 20, 21, 23, 25 (floors 30, 33, 36, 39, 50); 7 loci with 10 kb or more beyond up to floor 36, 9 at 39 and 50. Total span beyond the supported extent: 143, 146, 146, 183 and 186 kb, against 349 kb beyond the core. So the supported span keeps about half of the extra span that the core span would cut, and it is nearly flat from 30 to 36 bits and steps up at 39.
- **The 14 hand-tallied stretches** (from the previous tally; this check uses the same 11 loci the idea came from, so it is not independent). Noise-only stretches: 5 of 5 dropped at every floor. Strong own-family hits kept: 7 of 8 at floors 30 to 36, 6 of 8 at 39 and 50. Lost at every floor: the Leucr1 `bE` hit (29 bits, 30% identity over 80 aa; doubtful). Lost only from 39: a Gabo G3 PR pheromone hit (46 aa, 48% identity, 38 bits). Short real genes cannot reach high bitscores, so a floor of 39 can drop them; 30 to 36 keeps it. The noise stretches' best hits were 25 to 31 bits.
- **T48-F** (the case that started this): cluster 118,043-153,196; core and supported 118,043-130,301 at every floor, so the 22.9 kb carried by one weak hit is dropped. **Leucr1 bLocus:** cluster 145,579-160,941; supported 145,579-147,340 / 147,223 (floor 30 / 39): the 13.7 kb stretch is dropped.
- All other-family strong hits (receptor hits inside HD and Balpha spans) are excluded from the supported span by design.

Limits: 49 called loci from one panel; the 14-stretch check shares its loci with the idea; the floor is chosen on very few short genes (one pheromone); withheld loci have no `supported_span` (no gene coordinates in the report); no check on Ascomycota or Mucoromycota genomes; the bridging number is a proxy.

Recommendation (for the curator): keep `supported_span` report-only; consider lowering the default floor from 39 to about 33 (plateau 30 to 36), which keeps the short pheromone and still drops all five noise stretches. The Ascomycota and Mucoromycota panels were then tested (next update): the floor makes no difference there. Ruled and applied: default 33 (see the last update). Changing clustering or the reported `start`/`end` is not supported by this evidence and could change calls (merge, polish window).

## Update: other groups, and withheld loci (curator asked, 2026-10-08)
Code: branch `core-span` (draft PR #59): `supported_span` and `core_span` now also on withheld loci. Full suite: 1,170 passed.

**Other groups** (Ascomycota 30 + Mucoromycota 80 regression-panel genomes; baseline = current `main`; supported_span code at floors 33 and 39; `results/2026-10-08_supported_span_other_groups/`):
- 0 call differences on all 110 genomes at both floors (family, contig, cluster span, confidence, idiomorph, class, genes, merged_from, withheld loci). Every reported gene inside its supported span; the core always inside the supported span.
- Ascomycota (26 called loci): extra span beyond the core 65.2 kb (20 loci), beyond the supported span 62.9 kb: the supported span changes little, because few loci have noise stretches. Mucoromycota (52 called loci): 147.1 kb beyond the core in 39 loci, 95.3 kb beyond the supported span in 9: the supported span keeps the extra stretches that hold hits and drops the rest.
- The 8 called loci with 10 kb or more beyond the supported span (3 Ascomycota, 5 Mucoromycota; 10 to 41 kb): tblastn of the run's query set over the dropped stretches finds 0 strong hits (e < 1e-3) in all 8, with 37 to 241 weak hits each (best e-value 2e-3 to 2e-2). So the part dropped is noise-level hits only.
- Floor 33 and floor 39 give the same totals on called loci in these groups (extra span 62.9 and 95.3 kb at both); only withheld-locus counts differ slightly (beyond-supported >= 10 kb: Ascomycota 304 at 33 and 305 at 39, Mucoromycota 59 and 64). Floor 39 also has 0 strong hits in the dropped stretches. The floor choice therefore matters in the Basidiomycota set (the short pheromone), not here.

**Withheld loci, old versus new database** (the question this note started from; Basidiomycota panel + T48-F, 34 genomes, new code at floor 33, old = the database before the C. cinerea record, new = current; `results/2026-10-08_supported_span_withheld/`):
- 4,021 withheld loci under the new database and 3,995 under the old. Matched one-to-one: 1,625; 1,590 identical in span, supported and core.
- At the genome/family/contig level (2,205 keys in both databases, 2,185 with the same number of loci): **89 keys had a changed cluster span, and in 86 (97%) the supported span and the core span did not change.** So the span changes of the 164-locus regression are almost entirely noise stretching the cluster span. The 3 exceptions are real changes (a strong own-family hit now extends the locus; example: Fomme1 PR NW_006760402.1, 1,111,569-1,121,523 to 1,111,569-1,141,382, with the supported span moving the same way). 3 more keys changed their supported span although the cluster span did not (gene-set changes). 20 keys have a different number of loci and 24 keys exist in only one database (22 loci only in the new, 4 only in the old); these were not compared.
- The curated record change therefore does what it should (the gene set changes a few loci) and the cluster-span noise it exposes is separable with `supported_span`.

Limits: 34 and 110 genomes from the fixed regression panels, not a random sample; the old database is the one before the record change plus whatever else changed between those trees (db diff not isolated); the hit tallies use tblastn on the run's query set as an approximation of detect's hit table; the 164 and 89 are counted differently (different matching), so they are not the same number.

## Update: default floor set to 33 and re-tested (curator ruled 2026-10-08)
`--supported-min-bitscore` now defaults to 33 (was 39; test first, then the one-line change; branch `core-span`, draft PR #59). Re-run with no flag from a fresh frozen tree (`results/2026-10-08_supported_span_default33/`):
- **Calls unchanged:** 0 differences against current `main` on the 110 Ascomycota and Mucoromycota genomes, and 0 differences against the earlier new-code run on the 34 Basidiomycota genomes (family, contig, cluster span, confidence, idiomorph, class, genes, merged_from, withheld loci). The reports record a floor of 33 in every locus.
- **Spans identical to the explicit 33-bit runs** in every field except one gene's `reference_record` in Leppa1 (below): 34 of 34 Basidiomycota reports and 109 of 110 others. (A first comparison against a floor-33 run made before withheld loci had the new fields differed only in those new fields.)
- **Containment:** 160 and 348 reported genes, 0 outside the supported span.
- **Full suite:** 1,170 passed, 0 failed.
- The effect of 33 against 39 is as measured before: same totals on Ascomycota and Mucoromycota called loci; on the Basidiomycota panel 33 keeps the 46-aa Gabo G3 pheromone hit (38 bits) that 39 drops, and still drops the five noise-only stretches.

**A separate finding: an unstable `reference_record` on an exact tie.** In Leppa1 (GCA_001692735.1) one gene model (58,355-60,583, `exonerate_refine`, 59.6% identity) is attributed to the A1163 or to the Af293 curated record in different runs of the same code: current `main` 0 of 17 runs chose Af293, the new code 2 of 17 (1 of 16 clean repeats). Same coordinates, same identity: two curated references tie exactly (the A1163 remnant and the Af293 MAT1-2-1 share their C-terminal region). Nothing in a call, span or gene set changed; only the record named for that gene and `reference_records`. The numbers do not show that my change is involved (clean repeats 1/16 against 0/16, Fisher p = 1.00; 2/17 against 0/17, p = 0.48), and it does not touch hit ordering, but I cannot rule it out. A deterministic tie-break (for example, the lower record id) would remove it; that changes reported records for tied hits, so it needs a ruling and is not done. `results/2026-10-08_leppa1_tiebreak/`.

## Curator decisions
Ruled 2026-10-08 (J. Stajich): option 3, a report-only `core span` (built, see the update above). Options considered:
1. Keep the curation (the old gene-6 entry was wrong and the four peptides are real) and accept span noise from short queries as a known property. Cost: reported spans of called loci can include a 20-kb or larger stretch for no biological reason (T48-F).
2. Change the code so that a hit joins a cluster only above a floor (a minimum bit-score or e-value, or a per-family cluster span built from that family's own hits). This changes detection for every family and needs the regression panel and a ruling on the threshold.
3. Report a "core span" next to the cluster span (report-only), so nothing in calls changes.
Either way the `curate-cinerea-b43` branch itself is clean: its changes are the two sequence edits, and aliases are harmless.

## Files
- `results/2026-10-08_cinerea_trace/` (`cmp.py`, `comparisons.txt`, `jobs.tsv`, `lists/three.tsv`, `gene6_old_vs_new.faa`).
- Variant worktrees kept on HPCC until a ruling: `.claude/worktrees/trace-noNew`, `trace-oldSeq6`, `trace-both`. Run outputs: `results/2026-10-08_cinerea_trace/{baseA,baseB,cand,noNew,oldSeq6,both}/` (main checkout, not tracked).
