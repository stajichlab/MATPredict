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
- The other 164 changed loci in the regression are of the same kind in the three genomes I decomposed; I did not trace each of them to a hit. The second Mixia HD span change (2,169,865 to 2,234630 start) is reproduced by the gene-6 swap but I did not find the individual hit.
- Calls: in the regression, 1 call touched (the C. cinerea CC3 genome, the curated one). No call changed in T48-F (HD high, PR medium before and after); only the HD span.

## What changed in detection
Nothing. This note adds the experiment; no source change.

## Limits
- Three genomes of 33 were run for the decomposition; two loci were traced to individual hits. The mechanism is shown, not counted over all 164.
- The hit coordinates come from my own tblastn runs with detect's settings, not from detect's hit table; the match with the cluster boundaries is exact in both cases.
- Whether the spans of withheld loci matter to anyone is not established; they only appear in reports.

## Curator decisions
Open, with the options:
1. Keep the curation (the old gene-6 entry was wrong and the four peptides are real) and accept span noise from short queries as a known property. Cost: reported spans of called loci can include a 20-kb or larger stretch for no biological reason (T48-F).
2. Change the code so that a hit joins a cluster only above a floor (a minimum bit-score or e-value, or a per-family cluster span built from that family's own hits). This changes detection for every family and needs the regression panel and a ruling on the threshold.
3. Report a "core span" next to the cluster span (report-only), so nothing in calls changes.
Either way the `curate-cinerea-b43` branch itself is clean: its changes are the two sequence edits, and aliases are harmless.

## Files
- `results/2026-10-08_cinerea_trace/` (`cmp.py`, `comparisons.txt`, `jobs.tsv`, `lists/three.tsv`, `gene6_old_vs_new.faa`).
- Variant worktrees kept on HPCC until a ruling: `.claude/worktrees/trace-noNew`, `trace-oldSeq6`, `trace-both`. Run outputs: `results/2026-10-08_cinerea_trace/{baseA,baseB,cand,noNew,oldSeq6,both}/` (main checkout, not tracked).
