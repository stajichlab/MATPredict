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

## Curator decisions
Ruled 2026-10-08 (J. Stajich): option 3, a report-only `core span` (built, see the update above). Options considered:
1. Keep the curation (the old gene-6 entry was wrong and the four peptides are real) and accept span noise from short queries as a known property. Cost: reported spans of called loci can include a 20-kb or larger stretch for no biological reason (T48-F).
2. Change the code so that a hit joins a cluster only above a floor (a minimum bit-score or e-value, or a per-family cluster span built from that family's own hits). This changes detection for every family and needs the regression panel and a ruling on the threshold.
3. Report a "core span" next to the cluster span (report-only), so nothing in calls changes.
Either way the `curate-cinerea-b43` branch itself is clean: its changes are the two sequence edits, and aliases are harmless.

## Files
- `results/2026-10-08_cinerea_trace/` (`cmp.py`, `comparisons.txt`, `jobs.tsv`, `lists/three.tsv`, `gene6_old_vs_new.faa`).
- Variant worktrees kept on HPCC until a ruling: `.claude/worktrees/trace-noNew`, `trace-oldSeq6`, `trace-both`. Run outputs: `results/2026-10-08_cinerea_trace/{baseA,baseB,cand,noNew,oldSeq6,both}/` (main checkout, not tracked).
