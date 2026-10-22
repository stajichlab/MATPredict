Pilot (2026-10-22): core-only curated records from GenBank MAT1-1-1/MAT1-2-1 deposits, and their effect on detection.
Full write-up: the agent report handed to the curator (not committed). Nothing here is a curator decision.

generator_dryrun_106.tsv   scripts/deposit_to_record.py on all 106 selected deposits (no records written): verdict, reasons
generator_pilot.tsv        the 14 pilot deposits (13 records written, 1 NEEDS_HUMAN)
generator_dryrun_106_v2.tsv  generator v2 (curator rulings of 2026-10-22) on the same 106 deposits, offline from the
                           pilot NCBI cache, against the pre-pilot db (2d14712, the db the v1 dry run saw; the v1 run
                           reproduces byte-identically from it). No records written. Columns: v1 + species, order.
                           After the second review: 9 CONFIRMED / 56 CAVEAT / 41 NEEDS_HUMAN (v1 11/73/22).
                           Candidates proposed earlier in the run are kept apart from db/ in messages ("in-batch
                           candidate (not written)", reason names batch_*, `source` column in conflicts.tsv), but the
                           unknown-N slug counts db/ AND batch record_ids (the db rule; v1 did the same): two strain-less
                           isolates of one taxid get unknown-1 / unknown-2, so a strain-less record_id DEPENDS ON BATCH
                           ORDER and on the verdicts of earlier deposits in the run.
write_run_v1_vs_v2.tsv     real write runs (--out to a scratch directory, not db/) of v1 (c9f69ce) and v2 on the 106:
                           verdict and record_id per deposit. Every record_id is identical; v2 writes 65 records, a
                           subset of v1's 84 (e.g. Berkeleyomyces basicola 124036_unknown-1_MAT_MAT1-2 + unknown-2_MAT1-1).
generator_dryrun_106_v2_changes.tsv  per deposit: v1 and v2 verdict, moved to NEEDS_HUMAN, and the ruling(s) that
                           cause each NEEDS_HUMAN (scripts/attribute_rulings.py)
conflicts.tsv              ruling 8 conflicts table of that run, with a `source` column (db = curated record, batch =
                           in-batch candidate); header only: no deposit disagrees with a curated record
generator_dryrun_106_v2_branchdb.tsv  the same run against this branch's db (13 pilot records present): the 13 pilot
                           deposits come back NEEDS_HUMAN as identical-model duplicates of curated records (no conflict)
replay_genomes.tsv         283 BFD genomes (fixed seed 20261022): pilot species, 60 Hypocreales, 60 Helotiales,
                           16 other pilot orders, 86 Ascomycota controls, 12 Basidiomycota, 8 Mucoromycota
replay_changes.tsv         every difference between `detect` on origin/main (2d14712) and this branch, per genome
                           (detected and suppressed loci; base/pilot JSON), from scripts/compare_runs.py
replay_pilot_species_calls.tsv  per pilot-species genome: base vs pilot call and core models
scripts/                   selection, slurm worker (HPCC, --partial-core-policy coverage_profile
                           --confidence-cap-policy off --strong-single-gene-policy off), comparison, curate flow,
                           post-checks, truncation alignments. Run outputs: HPCC
                           /bigdata/stajichlab/jstajich/projects/MATPredict/curate_pilot_replay/runs/{base,pilot}
