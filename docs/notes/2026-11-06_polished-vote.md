# Polished-evidence vote (#275)

Status: opt-in (`--polished-vote-policy on`), default off, in no preset. Decided by the curator on 2026-11-06.

## The problem
The idiomorph vote ranks each idiomorph by the best bitscore of its informative hits. A short unpolished fragment of one idiomorph's core
gene can then outscore a polished model of the other idiomorph's. In the assembly-rescue pilot, *Trichoderma viride* SRR30921444 (truth
MAT1-1) was called MAT1-2: an unpolished 47-bit MAT1-2-1 window beat a polished 39.7-bit MAT1-1-3 model. No overlap resolution is made there,
so the resolver is not involved (the bitscore-basis experiment of draft PR #285, closed, did not change this call).

## The rule
In `Ascomycota:MAT`, when the called idiomorph rests only on unpolished evidence (none of its informative core genes has a polished model)
and the other idiomorph has at least one polished informative core gene, the polished idiomorph is called. The report's
`idiomorph_vote_note` records it (`rule: polished_vote`, `decision: flipped`, the gene set aside and the polished gene with their statuses
and bitscores); as for `--partial-core-policy`, the tier is at most `medium` and the vote's margin is not reported, because it describes
the idiomorph that was set aside. The rule does not run where the partial-core rule already decided, or where a profile classifier decided.

## How often
v0.8.1 run, 12,051 `Ascomycota:MAT` loci with a called idiomorph: 8,309 have evidence for one idiomorph only; of the 3,742 with a live core
gene of both, 1,613 have both polished, 2,115 have the called idiomorph polished and the other not, and **8 (all medium) have the called
idiomorph unpolished while the other is polished**: the only loci the rule can touch in a run of that size.

## What it does not do
Other families (MATsc, MTL, MATog, Basidiomycota, Mucoromycota) are out of scope: their genes and bitscore scales were not measured. The
relaxed second pass (calls reported at tier `low`) is not changed. It does not lower tiers for weak fragments: that is
`--weak-fragment-cap-policy` (docs/notes/2026-11-06_weak-fragment-cap.md).
