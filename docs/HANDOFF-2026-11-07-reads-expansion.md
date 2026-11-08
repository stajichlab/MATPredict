# Handoff 2026-11-07: raw-reads MAT typing, expansion to 300 per species and waves 3-4

Written 2026-11-07 ~10:00 PST. Read with `results/2026-11-06_reads_expand300/NOTE.md` (the decisions and the first numbers) and `datasets/reads/availability/README.md` (the PopSeq join).

## The project

MATPredict types mating-type idiomorphs from short reads (`wave_driver.py` + `genotype_mat_one_accession.sh`: 2M spots per run, k-mer panels cut from BFD assemblies, cores and flanks per lineage). Waves 1-2 typed up to 50 isolates per species. On 2026-11-06 the lab's read catalog **PopSeq** (`/bigdata/stajichlab/shared/projects/Fungi/PopSeq`) was joined to what MATPredict had typed (`datasets/reads/availability/`, PR #282) and showed ~45,000 untyped eligible isolates of species that already pass, plus species never tried. The goal: estimate idiomorph frequencies from many more isolates, and type new species.

## Curator decisions (2026-11-06/07)

1. **Cap 300 isolates per species** (typed + new), first with at most 25 new runs per study, then also the 2,782 runs that per-study cap left out (keep the run lists: the curator wants the same runs used later in a **PopSketch** analysis).
2. **Gate rulings** (`wave_driver.py`): `ZYGOSITY` = *Candida albicans*, *Pichia kudriavzevii* (heterozygous a/alpha diploids: `both` is a call and agrees with a collapsed assembly; words `PASS_ZYGOSITY`/`PROVISIONAL_ZYGOSITY`). `RELEASED` = *S. paradoxus, S. eubayanus, S. pastorianus, S. mikatae, Ogataea haglerorum* (HML/HMR cassettes confuse the panel; run flagged `UNVALIDATED:released_cassettes`). The 50% rule (UNVALIDATED gets the full run) was already in `decide()`.
3. **Wave 3** (new species with records): *Cryptococcus* (6), *Rhodotorula* (2), Mucorales (3), *Pichia kudriavzevii*, plus *Pyricularia oryzae* re-planned with lineage panels.
4. **Species groups** (`GROUPS`): *Cryptococcus gattii* complex (gattii, deuterogattii, bacillisporus), *Blumeria graminis* complex (graminis, hordei), *Sporothrix schenckii* complex (schenckii, brasiliensis).
5. **Wave 4**: all 28 `needs_reference` species that have loci of both idiomorphs (own or genus) in v0.8.1 (`results/2026-11-06_reads_expand300/needs_reference_fit.tsv`).
6. Partitions: `preempt` (`-p preempt -A preempt`) and `short` (QOS `short`, 2 h max) may be used. Curator: this would have been better as a **Nextflow** run; keep the Slurm strategy for now (`main.nf.example` is the starting point if moved).

## Code

- **Merged** PR #282: availability join, `expand_runs.py`, gate rulings, `curate`, site labels, panel pairs A1/A2 + Plus/Minus, per-genus span limits.
- **Open** PR #283 (branch `reads-wave3-pad`, head 953f77e at writing): `--pad` (Rhodotorula/Rhizopus/Mucor 6 kb), `BFD__<genome>` report folders (Mucoromycotina tarball), `select_runs` no-GCA genomes, `GROUPS` + `--species a|b|c`, *P. kudriavzevii* zygosity, `--parts partition:account[:qos]`, `--retry-part` (retries + summary never on preempt), `incomplete` command, defaults moved off `exfab`. Tests: `tests/test_wave_driver_gate.py` (16), `tests/test_expand_runs_script.py` (2).
- HPCC worktrees: `.claude/worktrees/run-expand300` (c07710d; first expansion batch read its worker from here) and `.claude/worktrees/run-wave3` (953f77e; everything since). Before moving `run-wave3`, check the worker is unchanged: `git -C .claude/worktrees/run-wave3 diff --quiet HEAD <new> -- datasets/reads/sra_survey/pilot/{genotype_mat_one_accession.sh,stream_type_sra.slurm,screen_reads.py}`.

## Status at handoff

Waves live in HPCC `results/sra_reads_typing/waves/<wave>/` (`out/` = per-run results; `<stage>/runs_all.tsv`, `<stage>/summary/results_all.tsv`; `gates*.tsv`; `jobs.tsv`).

| Work | Stage(s) | State |
|---|---|---|
| Expansion, wave1/wave2 | `expand300` (4,956), `expand300b` (2,782), `cleanup` (134) | **done**: 7,563 + 129 typed; 63+65 preempt-killed runs recovered by `cleanup`; 6 failed downloads remain (4 wave2, 1 basidio, 1 poryzae); 43 not at NCBI |
| Wave 3 tests + gate | `test`, `test_rhodotorula` | **done** (gates below) |
| Wave 3 full | `full` in wave3_basidio, _asco, _poryzae, _mucoro, _gattii | **done**: 542 typed (default ~60 runs per species, not yet 300) |
| Wave 4 tests | `wave4/test` (312 runs, 22 species), `wave4_groups/test` (40 runs) | **finishing**: 8 tasks on node x02 near the 2 h limit; retries 31824492/94 and 31824628, summaries **31824495** (wave4) and **31824629** (groups) |

Wave 3 gates: PASS *P. oryzae*, *R. toruloides*, *R. mucilaginosa*; PROVISIONAL_ZYGOSITY *P. kudriavzevii*; PROVISIONAL *C. deneoformans*; UNVALIDATED *C. neoformans* (0.70), *C. tetragattii* (0.57), *C. gattii* complex (0.55, 6/6 agree), *R. arrhizus* (0.60, 10/10), *R. microsporus* (0.74); HOLD *M. circinelloides* (0.41, 3 failed downloads: retest once). The Rhodotorula gate rows were merged from `gates_test_rhodotorula.tsv` into `wave3_basidio/gates.tsv`.

Wave 4 plan: 22 of 28 planned; held too_few_loci: *Blumeria graminis*, *B. hordei*, *Sporothrix brasiliensis*, *S. schenckii* (now grouped, planned in `wave4_groups`), *Hymenoscyphus fraxineus*, *Madurella mycetomatis* (no sister species to pool; stay held).

Check the queue: `ssh hpcc2 'bash -lc "squeue -u jstajich -h -o \"%i %j %P %T %r\" | grep wave"'`

## Next steps (in order)

1. **Wave 4 gate** when summaries 31824495 / 31824629 exist: `wave_driver.py gate --wave <W>/wave4` and `--wave <W>/wave4_groups`; first `incomplete --stage-name test` on each. Yeasts (*Brettanomyces, Hansenula anomala, Metschnikowia, Williopsis, Diutina*) may need a ZYGOSITY or RELEASED ruling (ask the curator; do not decide). Then `submit-full --gates gates.tsv --cap 300 --parts stajichlab:stajichlab,preempt:preempt` (retries go to stajichlab by default now).
2. **Merge expansion results into the wave tables**: for wave1 and wave2, combine `final`/existing `results_all_final.tsv` with `expand300`, `expand300b`, `cleanup` summaries (rerun `summarize_pilot.py` over the union run list, or concatenate and de-duplicate by run, keeping done over started). Then `gate --stage <merged>` or `regate`, then `curate` (keeps the C. albicans and cassette rulings).
3. **Top wave 3 and 4 up to 300**: add their `results_all_final.tsv` and `gates_final.tsv` to `DEFAULT_TYPED` / `DEFAULT_GATES` in `reads_availability.py`, `bash datasets/reads/availability/refresh.sh` (in an HPCC checkout of the branch), `expand_runs.py --cap 300 --per-study 1000000`, `submit-list ... --stage-name expand300`.
4. **Commit and publish**: copy each wave's `results_all_final.tsv` + `gates_final.tsv` (+ `species.tsv`) into `datasets/reads/sra_survey/pilot/<wave>/` (wave3_*, wave4*, updated wave1/wave2), commit the refreshed availability snapshot, `python scripts/build_reads_site.py --out web/intro/data/reads` (reads every committed `wave*` folder), update `results/2026-11-06_reads_expand300/NOTE.md`, PR. The site (mat.fungalgenomes.org/reads/) only changes when that PR merges.
5. Later: PopSketch on the expansion runs (lists: `waves/wave{1,2}/expand300/runs_all.tsv`, `results/2026-11-06_reads_expand300b/runs_wave{1,2}_extra.tsv`); consider Nextflow for the next wave.

## Pitfalls learned (do not repeat)

- **`exfab` is one node (gpu12)** with ~0.3 MB/s SRA downloads; never default to it, and never `ExcNodeList=gpu12` on an exfab job (BadConstraints, the job never starts).
- **Retry pass on preempt loses runs** (128 lost twice on 2026-11-06). Retries now go to `--retry-part` (stajichlab). After every stage run `incomplete`.
- **Submit limit 5,000 queued array tasks per account.** A refused `submit_stage` can be PARTLY submitted (the first array goes in): check `squeue` by StdOut path before resubmitting.
- **preempt AssocGrpMemLimit** starves small jobs behind big arrays: `scontrol update JobId=<big> Nice=1000`.
- `scontrol update ... ExcNodeList` on someone else's array (e.g. the deeploc job 31771400) overwrites its own exclusion list: only touch wave jobs.
- Plain `git push` in the local worktree tries to push `main` too: push the named branch.
- HPCC SSH drops under load: keep remote commands short, avoid `ls` over the big `out/` folders.
