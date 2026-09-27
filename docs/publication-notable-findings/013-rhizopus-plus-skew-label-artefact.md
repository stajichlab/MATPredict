# 013. The apparent Rhizopus Plus skew was a labelling artefact

- **Category:** assembly-or-annotation artefact (code); method
- **Status:** artefact-explained
- **Lineage:** Mucoromycota, Rhizopodaceae

## Summary
Thirty-three Rhizopus calls carried sexM at about 96% identity and no sexP, yet
were labelled Plus. The flank gene btbA, marked Plus-only from a single record,
outvoted sexM by bitscore. After the fix all 33 read Minus.

## Evidence
- Before: R. arrhizus 42 Plus / 1 Minus in the early-diverging scan;
  `docs/notes/2026-09-26_early-diverging-scan-and-chytrid-control.md`
  (revised, commit `f7796ae`).
- Mechanism: `idiomorph_candidates` scored each idiomorph by the best
  bitscore of any `present_in_idiomorphs` gene; btbA (long protein, ~98%)
  scored e.g. 871 vs sexM 378. Example GCA_000696915.1: sexM 95.7% vs sexP
  32.1%, per-contig resolution names sexM, label Plus.
- After: btbA `idiomorph_informative: false` (`08a2616`), then present in both
  idiomorphs (`4457c3a`); all 33 sexM-only calls read Minus
  (`results/2026-09-26_mucoro_rescore/`, committed).
- btbA sits at 66-99% in 29 Rhizopus Minus loci
  (`results/2026-09-27_tier_rule_replay/NOTE.md`).

## Method that found it
sexM/sexP tree label check (`results/2026-09-26_sexMP_phylogeny/label_check.txt`).

## Verification done / still open
Done: HMM classifier and tree agree (entry 016). Open: the true Plus/Minus
ratio in Rhizopus after the fix, and whether it reflects sampling.

## Limits
No mating-type truth for these strains.

## Related
Entry 016.
