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
