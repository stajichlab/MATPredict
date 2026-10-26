# Mucoromycotina, v0.7.0 definitive run (`cleanrun-2026-10c`): judgement

Tables: `results/2026-10-26_v070_mucoromycotina/` (974 genomes listed, 973 reports, 1 no report).

## Result
- **Identical to the first clean run** (`cleanrun-2026-10b`): all 973 reports have the same calls and gene models. This campaign uses `detect --phylum Mucoromycota --emit-cds-fasta` and takes none of the five policy flags, and none of the code changes since 10b reaches it.
- Called 868 (BFD 258 of 289, LCG 549 of 621, Jena 61 of 64); uncalled 105; one BFD genome has no report (suppressed). Loci 914, all family `Mucoromycota:MAT`: 689 high, 212 medium, 13 low.
- The defective-gene flag marks 156 of 914 loci (149 genomes); report-only.

## Against v0.6.0 (`change_vs_baseline`)
861 identical call, 105 uncalled in both, 3 NEW calls, 3 idiomorph flips, 1 span change over 1 kb, 1 no new report.

| genome | set | change |
|---|---|---|
| GCA_900175165.2_FCH_5_7 (*Rhizomucor pusillus*) | BFD | NEW call |
| GCF_041956525.1_Rhipu1 (*R. pusillus* CBS 183.67) | LCG | NEW call |
| Mucor_pusillus_Lindt_NRRL_A-16035 (*R. pusillus*) | BFD | NEW call |
| Mucor_hiemalis_NRRL_A-16027, Mucor_subtilissimus_NRRL_6226, Mucor_varians_NRRL_1760 | LCG | idiomorph flip |
| Mucor_sp._NRRL_1454 | LCG | span changed over 1 kb |

All three new calls are *Rhizomucor pusillus*, which gained records during the clean-run preparation. The reports were not read one by one in this note.

## Not verified
- Whether the 3 idiomorph flips (all in the held-out LCG set) are corrections or errors: compare them with the curated labels in `analysis/2026-09-28_heldout-sets.md` before the held-out counts are quoted.
- Whether the one no-report genome is an intended suppression (it is the same genome as in the first run).
