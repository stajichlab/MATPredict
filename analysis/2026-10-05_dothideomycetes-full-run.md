# Dothideomycetes: full run with the curated references (2026-10-05)

Status: results in; the polish cap options are measured (recommendation below); two decisions open (cap rule, ambiguous idiomorph calls).

## Question
The ten Dothideomycete MAT records (PR #36, from the branch `curation-mucor-dothideo`) were never in the database for the v0.6.0 Ascomycota campaign. What do they change across all
2,731 BFD Dothideomycete genomes?

## Method
- Frozen worktree `run-84361bc` (= `main` at 84361bc: code, `db/` with the ten records). 13 jobs on `short` (`results/2026-10-05_dothideomycetes_full/jobs.tsv`, job ids 29408018 to 29408030), 218
  genomes per job (115 in the last), 16 CPUs, 25 minutes to 1 hour 40 minutes each, 0 failures. Default flags (polish cap 6), lineage routing from the taxid, as in the campaign.
- Baseline: the v0.6.0 Ascomycota campaign for the same genomes (`results/2026-10-03_ascomycota_v060`). The regression panels since v0.6.0 show at most one changed Ascomycota call, so differences are attributed to the records.
- Compared with `compare_v060.py`; the table `compare_vs_v060_genome_changes.tsv` lists every genome that gained, lost or changed. Reports: `reports_all.tar.zst`.

## Results
2,728 genomes have a report in both runs (3 have none in the new run: GCA_016165855.1, GCA_986281625.1, GCA_975972555.1).

| | v0.6.0 | with records |
|---|---|---|
| Genomes with a call | 2,439 (89.4%) | 2,451 (89.8%) |
| Loci | 2,536 | 2,565 |
| Full / partial / gene only | 2,247 / 260 / 29 | 2,294 / 231 / 40 |
| Confidence high / medium / low | 1,869 / 657 / 10 | 2,120 / 442 / 3 |
| APN2 and SLA2 both in the locus | 148 | 152 |

27 genomes gained a call, 15 lost every call, 2,424 are called in both; among those the idiomorph set changed in 136, best confidence rose in 290 and fell in 50, and the locus count changed in 18.
By order (genomes, called before, called after, gained, lost, idiomorph changed): Pleosporales 1,012, 862, 875, 14, 1, 91; Mycosphaerellales 410, 348, 344, 10, 14, 18; Botryosphaeriales 807, 758, 758, 0, 0, 6;
Dothideales 204, 199, 199, 0, 0, 20; Venturiales 112, 109, 109, 0, 0, 0; Cladosporiales 90, 88, 88, 0, 0, 1.
Gains are mainly *Corynespora* 8, *Alternaria* 5 and *Cercospora* 4. The SLA2 flank share does not move, which fits the finding that SLA2 is detached from the locus in Dothideomycetes
(`2026-10-05_dothideomycetes-sla2.md`).

### The idiomorph changes in *Parastagonospora* are corrections
87 of the 91 Pleosporales idiomorph changes are *Parastagonospora*, a population set (186 genomes).
- v0.6.0: 178 MAT1-1, 4 MAT1-2, 4 none. With records: 91 MAT1-1, 91 MAT1-2, 4 none, which is the 1:1 ratio expected from a heterothallic species. A 98% MAT1-1 split was not credible.
- The new MAT1-2 calls have margins of median 343.9 and minimum 234.6 (best candidate score minus the second). The curation round-2 pilot found the same on two genomes: direct tBLASTn gives *P. nodorum* MAT1-2-1
  at 98.9% and no MAT1-1-1 hit.
- Cause: with no Pleosporales MAT1-2 reference, the HMG-box gene was named MAT1-1-3 (a known cross-match trap) and the genome called MAT1-1.
- Caveat: the new calls also list weak MAT1-1-1 and MAT1-1-4 hits in the same locus; not checked.
Table: `parastagonospora_check.txt`.

### Other idiomorph changes (49) are not resolved
- *Aureobasidium* (20): 15 change MAT1-2 to MAT1-1, 2 MAT1-1+MAT1-2 to MAT1-1, 2 to undetermined, 1 to both. *Aureobasidium* is homothallic, so a single-idiomorph label is unreliable (the pilot flagged one flip at a margin of 2.3).
- New calls of both idiomorphs, 17 genomes: *Friedmanniomyces endolithicus* 6, *Cercospora kikuchii* 3, and one each in *Aureobasidium pullulans*, *Cladosporium fusiforme*, *Neoascochyta*, *Pseudocercospora pini-densiflorae*,
  *Teratosphaeria gauchensis*, *Trematosphaeria pertusa*, a Mycosphaerellales sp., and *Nothophaeocryptopus gaeumannii* (gained). The pilot judged *C. kikuchii* a homothallism candidate and the others possibly spurious. Not checked against any truth.
- *Phyllosticta* 6 (five undetermined to MAT1-1), and single genomes elsewhere.

### The 15 lost calls are caused by the polish cap
- Genomes: 10 *Zymoseptoria* (9 *Z. tritici*, 1 *Z. pseudotritici*), 3 *Dothistroma* (2 *D. septosporum*, 1 *D. pini*), *Alternaria alternata* 1, *Hortaea werneckii* 1. All were MAT1-1, medium, partial locus in v0.6.0.
- This includes IPO323 (`GCF_000219625.1`), the genome the *Z. tritici* MAT1-1 record came from. Its locus (`NC_018206.1:616,667-623,138`, MAT1-1-1 and APN2) is still found but suppressed:
  v0.6.0 polished both genes (`polished_genes: 2`); the new run polished none, because the new references admit more clusters (suppressed 11 to 38) and the per-family cap of 6 skips the true locus.
  The round-2 pilot saw the same effect on two genomes ("the polish cap costs two loci").
- Test: the 15 genomes re-run on the same frozen tree with `--max-polished-clusters-per-family 0` (`results/2026-10-05_dothideo_cap0_lost/`). **All 15 calls return**, each identical to its v0.6.0 call. The run took 43 minutes
  on 8 CPUs for 15 genomes, so removing the cap for the whole class is costly.

## The polish cap, in plain terms
After the search, nearby hits are grouped into candidate regions. A candidate is **admitted** if it hits at least 2 different MAT genes, one of them a core MAT gene (alpha-box or HMG); there is no identity cutoff. Admitted candidates are then
**polished**: their genes are modelled with exonerate and miniprot, which is the slow step, and only polished candidates can become calls. The **cap** limits polishing to the top 6 admitted candidates per genome, because polishing everything cost
too much (measured on 561 genomes: 135 to 75 hours, six timeouts removed). "Top 6" is decided by a ranking made before any polishing: most distinct genes first, then best identity, then number of hits.
The cap does not cull noise; it only decides which candidates get a chance. A candidate that is not polished is never called.

### Why the records broke it, measured (`cap_rank_replay.py`, `cap_tiered_replay.py`)
- The median genome has 33 admitted candidates (90th percentile 47, maximum 80), so the cap binds in 99% of genomes. Across the class there are 88,371 noise candidates for 2,582 true loci.
- The ranking prefers breadth over strength. A noise region that hits three genes at 33 to 43% outranks a true locus that hits two genes at about 100%. In the 15 lost genomes the true locus ranked 7th (8 genomes), 8th (3), 10th (3) and 12th (1), just past the cap of 6.
- The replay of the current rank agrees with the candidates actually polished in 99.9% of slots, so it is a valid way to test alternatives. Truth = loci called in v0.6.0, in the full run, or in the cap-off re-run.
- Identity separates real from noise: true candidates have a median best identity of 84% (5th percentile 75%, lowest 39.5%); noise candidates have a median of 36% (95th percentile 47.5%).

### Options, with measured loss and cost (true loci lost out of 2,582; polishing work relative to today)
| option | true loci lost | polishing work |
|---|---|---|
| Today: cap 6, rank genes then identity | 17 | 1.00x (16,281 clusters) |
| Cap 10 | 2 | 1.66x |
| Cap 15 | 0 | 2.47x |
| Cap 20 / cap 30 | 0 / 0 | 3.25x / 4.56x |
| No cap | 0 | 5.59x (90,953 clusters) |
| Rank by identity first (cap 6) | 6 | 1.00x; the project's earlier out-of-sample test lost 2 Mucoromycota Minus calls with this rule |
| Rank candidates with best identity of 50% or more first, then today's rule (cap 6) | 0 | 1.00x |
| **Tiered: polish every candidate at 50% or more, then top up to 6 with today's rule** | **0** | **1.00x (16,284 clusters)** |
| Tiered at 45% / at 60% | 3 / 1 | 1.02x / 1.00x |
| Identity floor before ranking (drop candidates below the floor) at 40% / 50% | 3 / 21 | removes 72.5% / 96.1% of noise candidates |

- With the 50% tier, genomes have 0 to 8 strong candidates (801 have 1, 873 have 2, 663 have 3, 257 have 4, 50 have 5, 25 have 6 or more, 53 have none), and only 2 genomes have more than 6, so the tier costs the same as today.
- The 21 true loci below 50% identity are mostly one species, *Botryosphaeria dothidea* (44.5%), plus *Friedmanniomyces endolithicus* (39.5%), *Cladosporium* sp. BSL10, *Vermiconidia calcicola*, *Microcyclospora tardicrescens* and one *Dothideomycetes* sp. They rank in the top 3 on today's rule, so the tiered rule
  still polishes them; a hard floor would lose them.

### Recommendation
Use the tiered rule: it keeps today's rule for everything weak and only adds a guarantee that a strong hit is never skipped. On this class it recovers all 15 lost genomes at no extra cost. A hard identity floor would cut the most noise but loses divergent real loci, and raising the cap works but costs 1.7 to 5.6 times
more polishing. The 50% tier is calibrated on Dothideomycetes only; because the cap code is shared by every clade, it needs the project's cap and regression panels (Mucoromycota, Basidiomycota, Zygo 23) to show no change before it ships.

## Decisions for the curator
1. Polish cap: adopt the tiered rule (above), or one of the other options. If adopted: implement on a branch with tests, re-run the 15 lost genomes and the full Dothideomycete run, and run the regression panels.
2. *Aureobasidium* and the 17 double calls: leave as reported with the homothallism caveat, or add a rule. No ground truth is available here.
3. Treat the new run as the Dothideomycete baseline only after (1), or keep v0.6.0 for the 15 genomes.

## Limits
- The correction is checked only by the 1:1 ratio, the margins and the pilot's tBLASTn on two genomes. No mating-type assays were available.
- 277 genomes remain uncalled (reference gaps in Capnodiales, Extremaceae and others, unchanged by these records).
- The 3 genomes without a new report were not investigated.

## Files
`results/2026-10-05_dothideomycetes_full/` (list, jobs, `compare_v060.py`, `compare_vs_v060_*`, `parastagonospora_check.*`, `reports_all.tar.zst`), `results/2026-10-05_dothideo_cap0_lost/` (the 15 genomes, `cap0_recovery.txt`, reports). Ranking replays: `cap_rank_replay.py` and `cap_tiered_replay.py` with their `_summary.txt` tables, in `results/2026-10-05_dothideomycetes_full/`.
