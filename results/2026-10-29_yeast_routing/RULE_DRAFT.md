# Draft rule: MTL outside its taxonomic scope (for review)

Status: draft for curator review. Nothing in `src/` or `db/` changes until it is approved.

## Evidence it rests on

- Routing experiment v1 (MTL not searched outside Serinales) removes the cross-call duplicates but loses
  all 22 MTL-only calls; at least 5 of those are genuine (PAP1, OBP1, PIK1 flank the MTL genes).
- 181 genomes have both a MATsc and an MTL locus at the same position in v0.7.0; v1 keeps the MATsc one.

## Proposed rule

1. **Keep MTL searchable everywhere in Saccharomycotina** (v0.7.0 behaviour), so no genuine call is lost.
2. **Cross-calls**: when two families call overlapping loci in one genome (overlap of 0.9 or more of the
   shorter locus, the `--match overlap` rule of PR #209), keep the call from the family whose
   `taxonomic_scope` contains the genome; keep the other as a note on that locus (`cross_call: MTL`),
   not as a second locus. If neither or both families are scoped to the genome, keep the higher-scoring call
   and note the other.
3. **MTL outside Serinales**: an MTL-only call in a genome outside the MTL scope is reported as
   **provisional** (the site word for held-back-but-shown calls), with the reason `mtl_outside_scope`.
4. **Promotion**: that call is reported normally (not provisional) when at least two of PAP1, OBP1 and PIK1
   are found within 15 kb of the MTL genes. This is measured with the same search as `flank_search.py`.
   On the 22 genomes it promotes 5 (all three flank genes) and leaves 17 provisional.
5. **Ogataea**: send to MATog (v2) only after the Ogataea reference layout is checked (#189); until then
   Ogataea routing is unchanged.

## Open points for the curator

- Does "provisional" map onto the existing `confidence` tiers, or onto a separate flag? On the reads side
  provisional is a gate word; the assembly side may need a `tier_note` value instead of a new tier.
- Threshold: two of three flank genes versus three of three. In this set the two choices give the same
  answer for the 5 genomes with all three, and the 2 PIK1-only twins stay provisional either way.
- Where it is implemented: either as a second search pass (search MTL only when the scoped family has no
  call) or as resolution after searching (rules 2 to 4 above). The second needs no change to searching and
  matches the overlap-dedupe approach in PR #209; the first saves compute. Not checked in code yet.
- The v1 setting `fallback_searchable: false` on MTL is not used by this rule.
