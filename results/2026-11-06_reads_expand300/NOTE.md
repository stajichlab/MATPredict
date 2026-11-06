# Reads typing: expansion to 300 per species, released species, wave 3 (2026-11-06)

Curator decisions of 2026-11-06, after the PopSeq availability join (PR #279):

1. **Cap 300.** Every typed species whose gate is PASS, PROVISIONAL or UNVALIDATED is topped up to 300 isolates (typed + new).
2. **New wave** for the `new_candidate` species that fit a two-idiomorph model: *Cryptococcus* (6 species), *Rhodotorula* (2), Mucorales (*Rhizopus arrhizus*, *R. microsporus*, *Mucor circinelloides*), *Pichia kudriavzevii*.
3. **Count** the `needs_reference` species that fit: `needs_reference_fit.tsv` (below).
4. **Gates.** *Candida albicans* gets a zygosity gate (both = MTL a/alpha is a call); the HML/HMR cassette species are released (flagged, not held); the 50% rule (already in `wave_driver.decide`) lets the unvalidated species run in full.

## 1. Expansion (cap 300)

`expand_summary.tsv`: per species, typed before, untyped eligible, inside the 0.6-12 Gb window, selected, left out by the per-study cap, total after.

- Selector: `datasets/reads/availability/expand_runs.py` on snapshot `2026-11-06` (refreshed after the gate changes), `--cap 300 --per-study 25 --actions expand unvalidated`. Round-robin across studies, nearest 3 Gb first; panel columns copied from each species' gated `plan/<slug>/runs_full.tsv`.
- **4,956 runs, 107 species** (wave1 2,279; wave2 2,677). The per-study cap of 25 left out **2,782** further runs (a guard against clonal outbreak sets; assumption, not a curator ruling: rerun with a larger `--per-study` to add them).
- Submitted on HPCC from worktree `.claude/worktrees/run-expand300` (commit c07710d) with `wave_driver.py submit-list --stage-name expand300 --mem 8gb --env SCREEN=1` (GmG screen on): wave1 arrays 31785837-40, summary 31785841; wave2 arrays 31785842-45, summary 31785846. Results in each wave's `out/`, summaries in `<wave>/expand300/summary/`.
- Species at 300 after the run: *C. albicans*, *A. fumigatus*, *C. auris*, *F. oxysporum*, *L. parapsilosis*, *N. crassa*, *Y. lipolytica*, *Z. tritici*, *A. flavus*.

## 4. Gate changes (wave_driver.py)

- `ZYGOSITY = {Candida albicans}`: typed = one idiomorph or both; `both` agrees with an assembly's single idiomorph (an assembly of a heterozygote collapses to one haplotype); no homothallic shortcut; words `PASS_ZYGOSITY` / `PROVISIONAL_ZYGOSITY`, or `zygosity` among the reasons. Wave2 test numbers: 15 of 20 typed, 4 of 7 agree, so **UNVALIDATED:zygosity+agree_4of7** (was HOLD:typed_0.35+both_0.40+agree_2of5). The three disagreements are single calls opposite the assembly (A84, A92: `a` vs alpha; P37037: `alpha` vs a), worth a look (LOH or a mislabelled assembly).
- `RELEASED` (cassettes): *S. paradoxus*, *S. eubayanus*, *S. pastorianus*, *S. mikatae*, *Ogataea haglerorum* become `UNVALIDATED:released_cassettes+...` unless downloads failed.
- `wave_driver.py curate` applies both to an existing gates table; `regate` keeps them. Applied to wave2 (`curate_changes.tsv` in `datasets/reads/sra_survey/pilot/wave2/`; HPCC copy kept as `gates_final.tsv.pre-curate`).
- `reads_availability.py`: a species typed outside the waves but never gated is `typed_no_gate` (was `needs_reference`); *Pyricularia oryzae* (pilot4) passes the standard gate on its 19 runs (17 typed, 14 of 14 agree) but its pilot panel was cut from two rice-blast strains and missed other lineages, so it is planned again as `wave3_poryzae` with lineage panels.

## 2. Wave 3

Inputs per phylum from the v0.8.1 tables (`results/sra_reads_typing/wave_inputs_v081_{basidiomycota,mucoromycotina,ascomycota}`: `asco_genomes.tsv`/`asco_loci.tsv`/`reports_all.tar.zst` = that phylum's v0.8.1 `genomes.tsv`/`loci.tsv`/reports; survey files linked from `wave_inputs`). Waves `wave3_basidio`, `wave3_mucoro`, `wave3_asco`, `wave3_poryzae`; planned with `--borrow`. Panel builder: idiomorph pairs A1/A2 (*Rhodotorula*) and Plus/Minus (Mucorales) added; span limits per genus (*Cryptococcus* 4-250 kb, loci are 40-200 kb; *Rhodotorula* 1-60 kb, loci are about 2 kb).

Left out: *Mycosarcoma maydis* (every locus `undetermined`: the multiallelic b locus), *Phycomyces blakesleeanus* (no BFD genome), rusts, mushrooms with multiallelic loci, the *Saccharomyces* and *Ustilago* hybrids.

## 3. needs_reference species that fit (`needs_reference_fit.tsv`)

Of the 131 `needs_reference` species with at least 20 eligible isolates, counted against the v0.8.1 `detect` loci (mat_locus or partial_locus, idiomorph named):

| status | species | meaning |
|---|---|---|
| `own_both` | 15 | at least 2 loci of each idiomorph in the species' own assemblies |
| `genus_both` | 47 | both idiomorphs only with the genus (`--borrow`) |
| `one_idiomorph` | 17 | one idiomorph found (homothallic, clonal, or detect misses the other) |
| `no_loci` | 52 | nothing (mostly Agaricomycetes with multiallelic loci, chytrids, microsporidia; also *Nakaseomyces glabratus*, 99 genomes, no locus) |

**28 fit** (own or genus both, largest study under 90% of isolates, not pseudohomothallic or a cassette yeast): 4,949 eligible isolates. Largest: *Blumeria graminis* 934, *Cercospora beticola* 687, *Trichophyton indotineae* 600, *Brettanomyces bruxellensis* 561, *Cryphonectria parasitica* 410, *Exserohilum turcicum* 329, *Rhynchosporium graminicola* 196, *Sphaerulina musiva* 139. "needs_reference" (no curated db record) is therefore not a block for the wave pipeline, whose panels come from assemblies and `detect`; these can be planned as a wave 4.
