# Handoff: Verticillium reference records and what still needs validation

Date: 2026-10-11. Commit `8e34533fac` (main, not pushed) holds the three accepted records and the analysis note.
Read first: `docs/notes/2026-10-11_verticillium-reference-records.md` (withheld-locus analysis, both regression panels).
Consumer project: `/bigdata/stajichlab/shared/projects/KassonLab/Verticillium/nf/pipeline/mating_type/README.md`.

## State

| Item | State |
|---|---|
| `db/Ascomycota/Glomerellales/27337_tv103_MAT_MAT1-2` (AB505214.1) | accepted, committed, signed off by J. Stajich |
| `.../27337_nbrc6119_MAT_MAT1-1` (AB505215.1) | accepted, committed, signed off |
| `.../27337_vdls-17_MAT_MAT1-2` (RefSeq GCF_000150675.1, genome-derived) | accepted, committed, signed off |
| 8 records in `db/candidates/Ascomycota/` (F. fulva 2, P. chrysogenum 2, E. crustaceum 2, F. sacchari 2) | validated, status `needs_review`, NOT accepted, NOT committed |
| `db/candidates/Ascomycota/1301609_fen-198_MAT_combined` | pre-existing, untouched |

Effect on the four Verticillium genomes: MAT calls unchanged (MD4 MAT1-1, WVU/S13/S32 MAT1-2). Core hit identity rose from 34-50% to 94.8-100%. Withheld loci rose by 3 to 5 per genome.

## To verify (none of these is done)

1. **Five call changes in the second panel.** Not judged as improvement or error. Genomes: `GCA_001557505.1` (Reticulascus), `GCA_025802115.2` (Bertia moriformis, withheld -> MAT1-2), `GCA_004918105.1` (Aureobasidium, confidence medium -> high), `GCA_026023515.1` (Claussenomyces, MAT1-1-2 model added), `GCA_026023555.1` (Chrysothrix, model changed). Need: for each, find which reference hit changed and check the gene models by hand. All 11 records were added together, so re-run with records added one at a time (or by group) to assign each change.
2. **Withheld-locus rise, mechanism.** I did not extract the exact hit set that `detect` clustered. "A new weak hit completes a two-gene cluster" is the likely mechanism, not shown. Check on one locus (MD4 `MTK05_134`) with the detect debug output.
3. **273 withheld-only changes in the Ascomycota panel.** Not examined. The four Verticillium genomes only were checked.
4. **Stricter tblastn e-value.** Not tested. Need to decide whether to change the default (now 10) and measure the effect on withheld counts and on calls across all panels. This is a detect-behaviour change, so it needs the full regression panel.
5. **Tests and DB checks on the live repo.** I do not have a record of running the pytest suite, `curate-db validate` over the full live db, or `build-gff` after the three records were added. Run them on `main` at `8e34533fac`.
6. **Regression coverage of the live db.** Panel 1 covered only the two Usami records. Panel 2 used a worktree with all 11 records, not the committed state (3 records). No panel has run on the committed state.
7. **Self-call check.** VdLs.17 is not in the panel, so the check reports `no_report`. A direct detect run called MAT1-2 (high) at its own genome. TV103 and NBRC 6119 have no source assembly (`not_possible`).
8. **Eight proposed records.** Decision pending (curator only). Open question: F. sacchari MAT1-2 omits MAT1-2-3 because it is not in the `order.yml` MAT vocabulary. Decide whether to add MAT1-2-3 to the vocabulary or accept without it.
9. **Short et al. 2014 (PMID 25383550) follow-up.** My notes show eight records from other species were curated. I did not record whether every other locus accession in that paper was checked. Re-check the paper's accession list against the db before saying this is complete.
10. **Optional.** V. dahliae JR2 (`CP009079.1:3830361-3865909`) has the same arrangement and a MAT1-2 call. It is not a record. No MAT1-1 genome-derived record exists for V. dahliae (NBRC 6119 is the only MAT1-1 record, MAT1-1-1 only).

## Housekeeping

- Frozen worktrees in `.claude/worktrees/`: `run-e9060daf32-base`, `run-e9060daf32-vdahliae` (old `unknown-1` candidate name, now moot), `run-e9060daf32-vd11`. Remove when done.
- Result dirs: `results/2026-10-11_regression_vdahliae`, `results/2026-10-11_regression_vd11` (base symlinked), `results/2026-10-11_vdahliae_genomes`, `results/2026-10-11_vdahliae_genomes_final`.
- Scratch helpers (make_records.py, gen.py, gen_vdls17.py) are in node-local `$SCRATCH` and may be gone. The records are in git, so the generators are not needed to reproduce them, but their method is not recorded elsewhere.
- Rule: only the curator (J. Stajich) accepts records.
