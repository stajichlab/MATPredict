# 016. A profile-HMM idiomorph classifier on modelled proteins

- **Category:** method
- **Status:** verified (held-out Zygo 23)
- **Lineage:** Mucoromycota (design is general)

## Summary
Per-idiomorph HMMs score the modelled sexM/sexP proteins after gene modelling.
They call all held-out Zygo proteins correctly, and label/clade agreement across
the Mucoromycota rises from 148/195 to 204/206.

## Evidence
- `results/2026-09-26_hmm_classifier/NOTE.md` (committed `5e02dc8`): training
  sexP 76 seqs (27 genera), sexM 9 (7 genera); leave-one-genus-out 85/85, worst
  margin 22.5 bits; Zygo held out 23/23 (smallest margin 47.5); full pipeline
  23/23 on scaffolds and on contigs; 232/245 calls decided by the classifier;
  runtime 119 -> 128 s per genome.
- Concordance: `results/2026-09-27_sexMP_fasttree/NOTE.md`, `concordance.txt`
  (before 148/195, after 204/206).
- min_margin 25 bits (curator ruling 2026-09-27; changes at >= 30 bits agreed
  with the tree, those at <= 24.2 did not): `db/Mucoromycota/order.yml`, `533e266`.
- Earlier HMM typing on 6-frame ORFs failed (39% vs 51%); scoring modelled
  proteins avoids the search-mode loss (`results/2026-09-26_sexMP_hmm/NOTE.md`).

## Method that found it
pyhmmer 0.12.3; `scripts/build_idiomorph_hmms.py`; files in
`db/Mucoromycota/classifiers/MAT/` with `manifest.yaml` (commit `7277e10`, fixes
`67a122c`, `6e2f58d`, `1a00b0a`).

## Verification done / still open
Open: other families (Ascomycota MAT, MTL, Basidiomycota HD) once training sets
are diverse enough.

## Limits
sexM training is 9 records; no mating-type truth outside Zygo.

## Related
Entries 013, 014.
