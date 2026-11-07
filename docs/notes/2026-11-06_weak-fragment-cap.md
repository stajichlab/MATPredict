# Weak-fragment tier cap (#275)

Status: opt-in (`--weak-fragment-cap-policy on`), default off, in no preset. Decided by the curator on 2026-11-06.

## The problem
In some `Ascomycota:MAT` loci the call "this genome carries MAT1-1" or "MAT1-2" rests on one short fragment of a core gene with no
reliable gene model: *Microdochium nivale* (MAT1-1 low in v0.8.1; the strongest core hit is 41 bits), GCA_030077755.1 (MAT1-1 medium; no
core model by an independent `miniprot` check), Hyprub1. The tier says nothing about how thin the evidence is.

## The rule
For an `Ascomycota:MAT` call with a called idiomorph (MAT1-1 or MAT1-2) at tier `medium` or `high`: take the strongest core gene (a
`MAT1-*` gene that is not `lost_to`; best by polish status, with both-tools-agree first, then one tool or tools disagree, then unpolished,
and then by bitscore). If it is not a both-tools-agree model and scores under 50 bits, the tier becomes `low`. The idiomorph label, the locus
and the class are never changed, a `low` call is untouched, and nothing is raised. The report's `weak_fragment_cap` records the gene, its
status and bitscore, and the threshold.

## Why this trigger
Measured on the 12,051 `Ascomycota:MAT` loci with a called idiomorph of the v0.8.1 run (all 19,343 reports):

| Trigger | Loci matched | Already low | Would change tier |
|---|---|---|---|
| no core gene polished | 75 | 74 | 1 |
| strongest core single-tool only, under 50 bits | 41 | 0 | 41 |
| **strongest core not both-tools-agree, under 50 bits (this rule)** | **108** | 40 | **68** (23 medium, 45 high) |
| same, under 60 bits | 255 | 65 | 190 |
| same, identity under 40% | 1,678 | 52 | 1,626 |

99.4% of these loci have at least one polished core gene, so "no polished core gene" misses the cases of concern; an identity line
would demote ordinary divergent calls. The 50-bit line is the pipeline's own weak-fragment line (its noise floor is 39).

## What it does not do
- It does not fix the idiomorph of *T. viride* (#275): that is the separate `--polished-vote-policy` (a polished idiomorph beats an
  unpolished one in the vote), a second flag and a second pull request.
- It does not touch other families (MATsc, MTL, MATog, Basidiomycota, Mucoromycota): their core genes and bitscore scales differ and
  were not measured.
- Because the label is kept, idiomorph counts on the site do not change; only the confidence mix does (45 high calls would become low in
  a v0.8.1-sized run).

## How a locus looks
XylFL0594 (MAT1-2 medium in v0.8.1, confirmed MAT1-2 by reads): its strongest core evidence is a single-tool 30-bit model, so under the
cap it reports MAT1-2 at tier low: correct label, honest tier.
