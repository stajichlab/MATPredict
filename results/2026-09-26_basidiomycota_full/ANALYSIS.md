# Full Basidiomycota run: analysis

Run: 2026-09-26, frozen worktree `run-ad1f865` (basidio-anchors on PR #9 at
the SXI bonus slot). Predates the btbA / both-model / relaxed-pass / cap-rank /
HMM changes. Tables: `genomes.tsv`, `loci.tsv`, `by_order.tsv`,
`misses_by_order.tsv`, `size_bins.tsv`. Scripts: `analyze_full.py`,
`summarise.py`.

## Coverage

3,269 of 3,270 genomes have a report (5 suppressed). GCA_025617555.3
(Pucciniales, 1.71 Gb) timed out twice at 1 h; re-running alone as job
29117951 with a 4 h limit (`Pucciniales_long/`). The 33 gate-C pilot genomes
give identical calls. `SUBPHYLUM` is blank in samples.csv; subphylum is
derived from `CLASS`.

## Call rates

| Subphylum | Genomes | Called |
|---|---:|---:|
| Agaricomycotina | 2,417 | 83.0% |
| Ustilaginomycotina | 315 | 95.6% |
| Pucciniomycotina | 486 | 13.2% |
| Wallemiomycotina | 51 | 0% |

Routing decides the rate. Orders routed to their own curated family: 92-100%
(Agaricales 96.7%, Boletales 94.5%, Tremellales 94.4%, Polyporales 97.5%,
Russulales 93.5%, Ustilaginales 100%, Malasseziales 98.8%). Phylum-fallback
orders: mostly 0-20% (Sporidiobolales 19.8%, Cantharellales 17.2%,
Trichosporonales 4.7%, Pucciniales 1/129, Wallemiales 0/51, Sebacinales 0/12);
exceptions Hymenochaetales, Tilletiales, Microstromatales (>90%).

## What was called

3,958 loci: HD 2,738; Tremellales MAT 385; bLocus 360; aLocus 185; Aalpha 181;
Balpha/Bbeta 106. **Only 3 PR loci in the phylum**: the Agaricomycotina result
is effectively HD-only. Tremellales MAT: alpha 255, a 64 (Cryptococcus alpha
239, a 37). homothallic_candidate 0; unverified 0; assembly_gap_at_locus 56
genomes (Cantharellales 20).

## Misses (896 uncalled)

- 875 fail the bar with 0 modelled genes; 19 with 1; 1 has nothing localised.
- Assembly: N50 < 20 kb in 127 uncalled genomes; it dominates only the small
  residue in Boletales (24/28) and Agaricales (13/27).
- Anchor adjacency of withheld loci could not be tested: MIP1 is in almost
  every withheld fallback locus (Sporidiobolales 199/203), but MIP1 is in
  every genome and the withheld loci are broad. The fallback misses are a
  reference gap expressed as bar failures.

## What a re-run on current PR #9 code would change

- Relaxed-pass fix: 0-3 genomes (805/896 uncalled already have a best
  cluster at fraction 1.0; only 3 are below 0.5).
- Cap rank: not measured. 933 genomes had capped clusters; 547 of those are
  uncalled, almost all in fallback orders. Lineage orders with capped clusters
  were called anyway (Ustilaginales: 61 capped, 0 uncalled). A cap-off test on
  ~50 fallback genomes would settle it.
- btbA, both-model, HMM classifier: Mucoromycota only.
- A full re-run is not expected to change lineage-order results.

## Runtime binning

- Per-genome timeout is hardcoded: `timeout 3600` at
  `scripts/run_clade_panel.slurm:201`.
- < 500 Mb (3,212 genomes): median 91-330 s per order, max 1,737 s, no timeouts.
- > 500 Mb (58 genomes, 57 Pucciniales): median 2,500 s, 1 timeout.
- Proposal: `timeout "${GENOME_TIMEOUT:-3600}"`; bin 1 (< 500 Mb) on `short`,
  2 h jobs sized by order median; bin 2 (> 500 Mb) on `epyc`,
  `GENOME_TIMEOUT=14400`, one job of 58 at 16-way, `--time=8:00:00`.

## Open questions

1. Re-run or keep (the data do not justify a full re-run; a cap-off test on
   fallback orders answers the one open question).
2. Pucciniomycotina and Wallemiales need curated references.
3. Why PR loci are almost never called.
4. Approve the `GENOME_TIMEOUT` variable.
