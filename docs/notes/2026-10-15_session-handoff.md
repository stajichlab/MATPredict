# Handoff, 2026-10-15

## Where things stand

- **main**: v0.6.2 (tag `v0.6.2`, GitHub release published). Merged since: #64 (browser page, profile-bitscore pair rule, winners-only verdict), #65 (version bump), #66 (docker publishes on release only), #67 (locus downloads: FASTA, GenBank, proteins).
- **Live Pages site**: deployed at v0.6.2 (04:35 UTC). It does **not** yet include #67's downloads. Redeploy: manual run of "Deploy browser page" from `main`, or tag `v0.6.3`.
- **Open PR #68**: runs the independent exonerate and miniprot calls concurrently (serial and 4-worker reports and GFF identical on Phycomyces; 133 s to 100 s). Ready for review.
- **Branch `annotation-tracks`** (`25ba354`, pushed, no PR): `detect/tracks.py`, the data model for annotation tracks (helixer-lite GFF3, GenBank input, extension summary) with six passing tests. Not yet wired into calls or the report.

## Measurements

- **Saccharomyces, one detect run** (cProfile, serial): 6 min 2 s, 329 MB peak. Exonerate polishing is about 84% of the time (302 s of 359 s profiled). Database loading and reference build are about 4% (16 s). Pruning `db/` would save little and would change the query set, so it is not worth doing for speed.
- **helixer-lite, fungi lineage, 8 CPUs** (HPCC, job 29977074):
  - Saccharomyces (GCA_000151485.1): 8 min 49 s wall, 3.86 GB peak memory, 3.9 MB GFF.
  - Phycomyces (GCF_001638985.1): still running at the time of writing; see the job log in `/bigdata/stajichlab/jstajich/projects/profile/helixerlite/phycomyces/`.
- **Phycomyces locus window** (helixer-lite, laptop): all four polished MAT genes get complete ATG-to-stop ORFs. Helixer and the NCBI annotation (PHYBLDRAFT) both extend the polished models at the ends, in the same direction. 10 of 12 helixer genes in the window match an NCBI gene.

## Decisions taken

- Refined gene ends from helixer-lite are **reported alongside** the polished models, not replacing them.
- The helixer-lite runs are on HPCC (the lab venv at `/bigdata/stajichlab/jstajich/projects/helixer/venv` already has helixer-lite and the fungi weights).
- GenBank input is a planned **second route for candidate genes** and an **annotation track**. The report will show three tracks: the homology and exonerate evidence track, the helixer-lite track, and the GenBank annotation track when one is supplied, with a summary of how the tracks agree.
- The browser build will eventually include TensorFlow and the helixer-lite weights. Peak memory (3.9 GB on Saccharomyces) matters for that.

## Next steps

1. Read the Phycomyces helixer-lite result and record its timing alongside the Saccharomyces numbers.
2. Wire `tracks.py` into the report: a track per source, plus the summary table.
3. Add GenBank input: parse it, use its CDS proteins as candidate queries, and show it as the annotation track.
4. Decide whether to redeploy Pages now, or wait for a `v0.6.3` tag.
5. Merge PR #68 after review.

## Known issues

- `tests/db` fails at collection (unrelated to the detection changes; not investigated).
- Twelve tests are skipped on macOS because they need Linux (mafft from pixi, `/usr/bin/bash`). They run on HPCC.
- The NCBI e-mail warning on HPCC (`MATPREDICT_NCBI_EMAIL` is not set). Runs still complete.
- The Ascomycota panel representatives (*Fusarium*, *Aspergillus*) are not confirmed MAT1-1/MAT1-2 pairs.
- The pair rule (profile bitscore) overrides the 2026-09-26 curator ruling; the curator has not signed off.

## helixer-lite on locus windows (2026-10-15, later)

- Goal changed: helixer-lite runs only on windows around the loci the search found, to propose better gene
  models, not to annotate the genome. `detect/helixer_windows.py` widens each locus by ±5 kb, merges overlapping
  windows per contig, runs helixer-lite once on the windows, and maps coordinates back. Predictions are a track;
  they do not change calls.
- Phycomyces, whole genome (HPCC): 34 min 38 s, 3.9 GB peak, 17,037 gene lines. The window run takes 9.3 s.
- Phycomyces locus window (laptop, ±5 kb): 7 models. Polished vs helixer-lite 3′/5′ (bp, on the gene's strand):
  - algA: 3′ +189, complete ORF
  - tptA: 3′ +648, 5′ +12, complete ORF
  - sexM: 5′ 0, 3′ +3, complete ORF
  - rnhA: 3′ +2107, complete ORF (NCBI PHYBLDRAFT_80075 ends at 3991627; helixer 3991273)
- **Undercalling:** set-aside sexP (3985731–3985973) has no helixer model; it lies inside the sexM model.
- **Window sensitivity:** the ±5 kb window and the 33 kb window disagree on the ends (algA 3′, rnhA 5′). The window
  policy must be fixed and recorded before any model is accepted or compared across runs.

## Decision question (not yet decided)

Which model to keep when homology (polished) and helixer-lite disagree. Candidate rules, to be tested on the
panel rather than picked by eye:

1. Keep homology boundaries where exonerate and miniprot agree; report helixer's extension separately.
2. Extend a homology end to helixer's ATG or stop only when the extension is a complete ORF on the same strand
   and the other tool (miniprot) does not contradict it.
3. Never let helixer remove or replace a homology gene; it can only add or extend, and each change is recorded.

Evidence needed: for MAT genes (where helixer undercalls), homology should win; for flanking genes, helixer's
ends may be better. Check against NCBI annotation on the panel genomes.

## Decision rule: refinements from the quiz (2026-10-15)

Settled before any panel run, so the rule is not tuned against the panel:

- **Silence is not contradiction.** A tool not reaching a model's end (for example miniprot stopping at the
  polished boundary) is no evidence against the extension. Contradiction means an explicit conflict: an
  intron, a stop codon, or a frame change inside the proposed extension.
- **Missing helixer models are not evidence against homology.** Helixer undercalls MAT genes (sexP is absorbed
  into sexM). Homology keeps a gene when helixer has no model for it.
- **One window policy, fixed in advance.** The window is the locus ±5 kb, merged per contig, for every locus.
  No per-gene or per-run choice of window (in particular not by ORF length, which biases results upward).
- **The metric is fixed before scoring** (below).

## Metric (pre-registered; to be applied to the panel, not yet run)

Reference: the NCBI annotation of the same assembly (PHYBLDRAFT and equivalents). NCBI is a prediction, so
agreement is evidence, not ground truth. Each result is reported with that caveat.

1. **Gene match.** A pipeline or helixer model matches an NCBI gene when both are on the same strand and
   the overlap is at least 50% of the longer of the two.
2. **Boundary agreement.** For a matched pair, the 5′ and 3′ differences are compared with NCBI, in bp, strand-aware.
   A boundary agrees when the difference is within 30 bp (ten codons; a tolerance chosen before the run).
3. **Strata.** Every count is reported separately for core MAT genes (`role core_MAT`) and flanking genes
   (`flanking_conserved`, `flanking_variable`). The rules are judged per stratum.
4. **Per-rule score.** For each rule, count the genes where the rule's chosen boundary agrees with NCBI, and the
   genes where it changes a boundary away from NCBI. Compare against the baseline (keep the homology boundary).
5. **Undercalling (recall).** Fraction of NCBI MAT genes with any helixer model overlapping them.
6. **Report, do not hide.** Counts (n) with every number. With nine genomes the MAT strata are small (one to
   three genes per genome), so results are descriptive; no significance claim is made from them.

A rule is adopted only if it beats the baseline on the flanking stratum without losing MAT agreement. The
MAT stratum is the guard against helixer's undercalling.
