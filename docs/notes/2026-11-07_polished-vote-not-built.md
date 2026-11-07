# Polished-evidence vote: tried, not built (#275)

Decided by the curator on 2026-11-07, after the validation below. The code is on the branch `polished-vote` (no pull request).

## The idea
The idiomorph vote ranks each idiomorph by the best bitscore of its informative hits, so a short unpolished fragment of one idiomorph's core gene
can outscore a polished model of the other's. The proposed rule (`--polished-vote-policy`, `Ascomycota:MAT` only): when the called idiomorph
rests only on unpolished evidence and the other idiomorph has a polished core gene, call the polished one. In the v0.8.1 run only 8 of
12,051 `Ascomycota:MAT` loci match (all medium).

## What validation showed (frozen v0.8.1 code and preset, flag on, against the v0.8.1 reports)
- 300 stratified Ascomycota genomes: 298 identical, 2 changed.
- 47 assembly-rescue pilot assemblies, flag off against on: 0 changed, including *Trichoderma viride* SRR30921444 (truth MAT1-1, called
  MAT1-2), the case the rule was built for. Its report lists no MAT1-1-3 among the live core genes (only APN2, COX13 and an unpolished
  MAT1-2-1), so the rule has nothing to act on there.
- **XylFL0594 (GCA_022385605.1) went from MAT1-2 to MAT1-1: wrong.** Reads confirmed MAT1-2 (2M spots: MAT1-2-1 33 spots, no MAT1-1-1, MAT1-1-2
  or MAT1-1-3). The "polished" model that won is MAT1-1-3: single tool, 35.6% identity, 30.4 bits, 50% coverage, flagged defective (internal
  stop), a broken cross-match of the HMG box.
- Leptod1 (GCA_021432625.1) went from MAT1-1 to MAT1-2: the margin was 3.1 bits (MAT1-1-1 unpolished 56.2 against MAT1-2-1 polished 53.1), so
  either call is a coin flip; no truth.

## Why a rule on these features cannot work
XylFL0594 and *T. viride* have almost the same structure, a defective polished MAT1-1-3 against an unpolished MAT1-2-1 fragment of 40-47 bits,
with opposite truth: in XylFL0594 the unpolished side is right, in *T. viride* the polished side. Excluding defective models stops the XylFL0594
flip but would not help *T. viride* either (the Fable review of the pilot found one internal stop in its polished MAT1-1-3 model; that report
no longer lists the gene). No feature I could find in the reports separates the two.

## What the merged weak-fragment cap does with them
Both have an unpolished strongest core gene under 50 bits, so `--weak-fragment-cap-policy on` caps both at tier `low` and keeps the label, which
says what the evidence allows. A possible next step, not built and needing its own design and measurement: `undetermined` when each idiomorph has
only weak evidence (one unpolished fragment and one defective or single-tool model, both under about 50 bits).
