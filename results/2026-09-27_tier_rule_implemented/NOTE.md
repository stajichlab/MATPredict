# Allele-absent tier rule: implemented and replayed on the real code

2026-09-27. Curator ruling: variant B' with the closeness guard
(results/2026-09-27_tier_rule_replay). Code: polish-scope-cuts,
`tiering.allele_absent_genes_to_ignore` + `assign_tier(ignore_genes=...)`,
wired in `pipeline._build` (idiomorph and evidence now computed before the
tier; ignored genes also leave `_any_gene_unpolished`). Reports gain
`confidence_ignored_genes` when non-empty.

Replay: `replay_real.py` calls the real functions from the new code with each
run's own roster; post-tier caps re-applied from report fields. Reproduction
of reported confidence with the rule off: 100% in 7 runs, Mucoromycota
241/245 (the same 4 misses as the earlier replay).

## Guard scope: per call, not per gene

Applied gene by gene, the rule raised 101 calls (`pergene/`). 15 of them kept
ANOTHER other-allele gene that was strong or close, e.g. 6 Serinales alpha
calls still carrying a modelled MTLA1 45% / MTLA2 58%, and a Periconia MAT1-2
call whose kept MAT1-1-3 (39.5%) beat its own best core gene (31.8%). That is
the collapsed-locus shape the guard exists to block, so the guard is applied
to the whole call: if any other-allele gene is >= 50%, within 10 points of
the called allele's best modelled core gene, or of unknown identity, nothing
is ignored.

## Risers (medium -> high), call-level guard: 86

| Run | Family | Risers | Predicted |
|---|---|---:|---:|
| Serinales 882aa01 | MTL | 39 | 39 |
| Dothideomycetes d18ec5a + 52cd292 | MAT | 7 + 13 = 20 | 20 |
| Cap panel f7b9773 | MAT | 8 | 8 |
| Cap panel f7b9773 | MTL | 2 | 2 |
| Wallemia 3d8a755 | wallMAT | 17 | 17 |
| Basidiomycota ad1f865 | MAT (Cryptococcus) | 0 | 3 |
| Mucoromycota, early-diverging | - | 0 | 0 |

No call changes in any other direction. Cryptococcus 3 -> 0: each call also
carries the a-pheromone MFa at 52.4% (not_polish_candidate), which trips the
50% cap. MATsc: 0 risers.

## MTL / MATsc riser check (41 MTL, 0 MATsc): `mtl_matsc_risers.tsv`

- 36 supported: modelled PAP1/OBP1/PIK1 flanks and a modelled called core
  gene >= 40%.
- 2 weak core (Serinales GCA_002370695.1 MTLalpha1 38.2%; GCA_021272955.1
  MTLA1/MTLA2 37.7/36.0%), flanks modelled.
- 3 without flanks: cap-panel GCA_030583325.1 and GCA_031125185.1
  (idiomorph_gene_only; MTLalpha1 50.0% modelled; ignored MTLA1/MTLA2 38-40%),
  and Serinales GCA_030462985.1, class `homothallic_candidate`, MTLalpha1
  53.9% with the ignored MTLA2 38.75%. A homothallic-candidate call raised to
  high by ignoring one of its two alleles is inconsistent; flagged.
