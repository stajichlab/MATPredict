# Flank-carried rule: audit of the calls it changes outside Serinales

2026-09-26. Curator ruling (c): check the Ascomycota calls that the rule
(commit 3684c60, `detect/flank_carried.py`) changes before any sweep uses it.
No detect re-run. No code changed.

## Method

1. `replay.py` applies the rule exactly as implemented to existing reports.
   - Cap-6 panel (`2026-09-26_polish_cap/cap6`, 611 calls): 89 change, 52 to
     low and 37 withheld.
   - Early-diverging runs: Mucoromycota 16 of 242, Mortierellomycota 2 of 6,
     Kickxellomycota 1 of 7.
   - These counts match the implementation fork's replay exactly.
   - Output: `changed.tsv`, 108 calls.
2. `localise_e10.sh` runs tblastn on 88 genomes with each report's genetic
   code and detect's own settings (`-seg no`, e-value 10). The queries are
   the core MAT proteins from each run's `_reference.faa`. Output: `hits_e10/`.
3. `classify3.py` takes each call and each of its own core genes, and records
   the best in-locus hit's e-value and its rank among genome-wide loci.
4. `simulate.py` scores the alternative rules. Output: `simulated.tsv`.

Evidence classes use the strongest core hit in the call:
- strong: E ≤ 1e-5
- weak: 1e-5 < E ≤ 1e-2
- noise: E > 1e-2, or no hit

## Result 1: most changed calls rest on a noise-level core hit

| family | rule outcome | strong | weak | noise |
|---|---|---:|---:|---:|
| Ascomycota:MAT | low | 5 | 8 | 15 |
| Ascomycota:MAT | withheld | 1 | 8 | 26 |
| Ascomycota:MATtub | low | 2 | 2 | 5 |
| Ascomycota:MATtub | withheld | 0 | 0 | 1 |
| Ascomycota:MATyl | low | 1 | 2 | 9 |
| Ascomycota:MATyl | withheld | 0 | 0 | 1 |
| Ascomycota:MTL | low | 0 | 0 | 3 |
| Mucoromycota:MAT (all 3 runs) | low | 6 | 4 | 4 |
| Mucoromycota:MAT (all 3 runs) | withheld | 2 | 0 | 3 |

- Of the 37 Ascomycota calls the rule withholds, 28 rest on noise, 8 on weak
  hits and 1 on a strong hit.
- The withheld calls are mostly Orbiliales (20) and Dipodascales (11).
- These calls were admitted on modelled SLA2, APN2 or COX13 alone. Those are
  conserved housekeeping genes present in every genome. The only core
  evidence is a short tblastn fragment at 8–32% query coverage.
- Withholding them is mostly correct. Their mating-type call has no support.

## Result 2: the rule's geometry fits Serinales, not the SLA2/APN2 families

- In Serinales, PAP1, OBP1 and PIK1 sit inside the MTL. "Core hit inside the
  flank span ±3 kb" matches that layout.
- In Pezizomycotina and most Saccharomycotina families, SLA2, APN2 and COX13
  sit outside the idiomorph. The core genes lie next to or between them.
- When only one side is found, the span covers only that side. In Orbiliales
  the flanks found are APN2 and COX13, and 23 of 28 core hits lie on the far
  side of APN2 from COX13. That is the expected MAT position, yet the rule
  counts it as "outside".
- A single stray far hit also withholds a call, because the rule requires
  EVERY core hit to be inside. Example: *Didymobotryum rigidum* has
  MAT1-1-2, MAT1-1-3 and MAT1-2-1 between SLA2 and APN2. Two hits 40–60 kb
  away (MAT1-1-4, MAT1-2-4) push it out.

## Result 3: real loci wrongly withheld (3 calls, all with strong core hits)

- *Mucor irregularis* B50 and GCA_000697435.1: sexM at E ≈ 1.4e-13, the best
  genome-wide sexM locus. It sits 9.5–9.9 kb beyond glrA, and sexP lies inside
  the flank span. sexM was never modelled (`not_polish_candidate`, 20% query
  coverage).
- *Trigonopsis variabilis* GCA_003707065.3 (Ascomycota:MAT): MAT1-1-3 at
  E = 3.6e-10, the best genome-wide locus, 2.8 kb from the flank span. Another
  core hit falls just outside the ±3 kb padding.

## Result 4: the polish cap does not interact with the rule

- On the same genomes, 87 of the 89 changed calls also change in the cap-off
  run (`changed_off.tsv`).
- The 3 differences:
  - 2 genomes timed out in the cap-off arm and have no report.
  - 1 locus (GCA_029290875.1) was withheld by the bar under cap 6. It was
    flank-carried in the off arm too.
- No changed call overlaps a capped cluster in `evidence_diagnostics.jsonl`.
- The polish-scope cut (cut 1) skipped no core gene in any changed call.

## Simulated alternatives

Rule B keeps a call at low when the strongest core hit meets the e-value
floor and lies within W of the flank span. Otherwise it withholds.

| rule | Ascomycota low / withheld | Mucoromycota low / withheld |
|---|---|---|
| A (current) | 52 / 37 | 14 / 5 |
| B, W = 3 kb, no floor | 59 / 30 | 14 / 5 |
| B, W = 20 kb, no floor | 87 / 2 | 18 / 1 |
| B, W = 20 kb, E ≤ 1e-5 | 9 / 80 | 8 / 11 |
| B, W = 3 kb, E ≤ 1e-5 | 9 / 80 | 6 / 13 |

## Recommendation

1. Judge the strongest core hit, not all of them. This removes the
   stray-hit failure (*Didymobotryum*, the Mucor pair).
2. Add a core-evidence floor to flank-carried calls. Keep one only if its
   strongest core hit has E ≤ 1e-5; otherwise withhold it.
   - Detect already keeps the tblastn e-value (commit 895fb7f).
   - The 1e-5 cut-off is not tuned. It follows the large gap between the
     strong and noise classes in these data.
3. Make the window a per-family parameter, `flank_carried_window_bp`:
   - 3 kb for families whose flanks sit inside the idiomorph (Ascomycota:MTL).
   - 20 kb for SLA2/APN2/COX13 families and for Mucoromycota:MAT. 20 kb covers
     every strong-hit case here (the maximum distance is 9.9 kb).
4. With 1–3, Ascomycota keeps 9 of 89 changed calls at low and withholds 80.
   Mucoromycota keeps 8 of 19 and withholds 11. All 3 real loci above are
   kept.

## Limits

- The e-values come from an independent tblastn, not from the reports. They
  match detect's settings, but database size can shift them slightly.
- "Real" means a strong hit at the genome's best locus for that gene. No
  synteny re-test or phylogeny was done for these loci.
- Orbiliales was already flagged as a reference gap (5 of 40 localised). The
  noise calls there may hide real loci that need a curated reference.
- The same locus is often called by 2–3 Ascomycota families (MAT, MATtub,
  MATyl). There are 15 such duplicates among the 89. This is a separate issue.
