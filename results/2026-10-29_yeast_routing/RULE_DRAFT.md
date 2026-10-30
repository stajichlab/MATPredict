# Draft rule: MTL outside its taxonomic scope (curator decisions recorded, 2026-10-29)

Status: decisions agreed with the curator; not implemented. Nothing in `src/` or `db/` changes until a PR does it.

## Evidence

- Routing experiment v1 (MTL not searched outside Serinales) removes the cross-call duplicates but loses
  all 22 MTL-only calls; at least 5 of those are genuine (PAP1, OBP1, PIK1 flank the MTL genes). See `README.md`.
- 181 genomes have both a MATsc and an MTL locus at the same position in v0.7.0; v1 keeps the MATsc one.

## Rule

1. **MTL stays searchable everywhere in Saccharomycotina** (v0.7.0 behaviour); no genuine call is lost.
2. **Runs after searching.** All families are searched as now; a resolution step applies rules 3 to 5.
   Searching and the database are unchanged. Same approach as the overlap dedupe (PR #209).
3. **Cross-calls**: two families calling overlapping loci (0.9 or more of the shorter locus) keep the call
   from the family whose `taxonomic_scope` contains the genome. If neither or both are scoped to it, keep the
   higher core-gene score. The other family is recorded on the locus as `cross_call` (nothing is hidden).
4. **MTL-only call outside the MTL scope** carries a flag, separate from confidence (confidence is unchanged):
   `family_unconfirmed`, reason `mtl_outside_scope`.
5. **Promotion**: the flag becomes `family_confirmed` when at least two of PAP1, OBP1 and PIK1 are found within
   15 kb of the MTL genes (same search as `flank_search.py`). On the 22 genomes: 5 confirmed, 17 unconfirmed
   (the two PIK1-only twins stay unconfirmed).
6. **Ogataea** goes to MATog (v2) only after the Ogataea reference layout is checked (#189); until then unchanged.

"Provisional" is not used here: it is already the reads-side word for a held-back gate.

## Still to do before implementing

- Check in code where the resolution step goes (`family_registry.py` routing versus the dedupe script) and how
  the flag is written to the GFF3, report and site tables.
- Re-run the 354 yeast genomes to confirm: duplicates gone, no MTL-only call lost, no-call control still empty.
