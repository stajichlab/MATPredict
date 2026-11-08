# `nearest_reference`: a proposed definition (issue #200)

Status: **proposal for the curator to rule on; nothing is built.** Date: 2026-11-07.

## What the 2026-10-28 rulings leave open

The idiomorph labels (`<Lineage>__<System>_<Cell>`, `scripts/idiomorph_labels.py`, registry `scripts/idiomorph_label_registry.json`) name the
frozen reference organism of the system. A Pezizomycotina class with no curated reference gets the placeholder `NoRef`, and a planned
`nearest_reference` column records the closest reference "by similarity now, by tree placement later". The v0.7.0 numbers
(`results/2026-10-28_idiomorph_labels_v070/summary_*.tsv`) are 63 `Peziz__NoRef_MAT1-1` / `_MAT1-2` loci (31 / 32), and `_undet` loci
for NeuCra 473, SacCer 233, BipMay 93, CryDen 60, BotCin 34, AspFum 17, CanAlb 16, XanPar 4, GeoCan 1, PhyBla 1.

Four things were left to define. Each has a proposal here and the reason for it.

## 1. Which gene gives the identity

**The strongest live core gene of the locus**, chosen as `weak_fragment_cap` chooses it: among the core genes (`MAT1-*` for `Ascomycota:MAT`) not
`lost_to`, rank by polish status (both tools agree, then any polished model, then unpolished) and then by bitscore. One gene per locus keeps the
column a single value and means a locus is placed by the gene that carries its idiomorph, not by a flank gene that is conserved everywhere.

## 2. Which aligner

**None new.** Use the identity the pipeline already recorded for that gene: `gene_evidence.identity` against the curated record it matched
(`reference_record`), which is the result of the pipeline's own search and model (tblastn, exonerate or miniprot, whichever produced the reported
model). No second alignment is run, so the column is a pure function of the detection report and identical input gives identical output. A
tree-based placement is a later, separate column and does not replace this one.

## 3. Tie-break between equal references

Compare identity to two decimals. On a tie take the higher bitscore of the same gene; on a further tie take the record with the smaller
`record_id` in plain string order. The rule is total, so there is never an unordered choice.

## 4. What is written when there is no hit, and where the column lives

- A locus whose strongest core gene has no matched record (no `reference_record`, or no identity) is written as the explicit value `none`. It is
  never empty.
- The value is the **system code of the matched record's organism class**, taken from the registry (`pezizomycotina_class_systems`), so a `NoRef`
  locus reports which framed system it sits closest to (for example `NeuCra`) beside its placeholder label. A matched record whose class has no
  registry code is written as that record's `record_id`.
- Two further fields travel with it in the same table: the identity (two decimals) and the record id, so a reader can see how close "nearest" was.
- **Where:** the label tables and the site label data only, computed on output. Not the curated database and not the detection report YAML, in line
  with the rule that labels are output-only.

## Tests that go with it

1. Every `NoRef` locus in the label table has a `nearest_reference` that is a registry code, a record id or the explicit `none`; none is empty.
2. Identical reports give an identical column (run twice; compare).
3. A synthetic tie (equal identity, different bitscore; equal both) resolves as described above.
4. A locus with no matched record gives `none`.
5. The report YAML and `db/` are byte-identical before and after the column is added.

## Open for the curator

- Is "the strongest core gene" the right single gene, or should the column list both core genes when a locus carries two?
- Should a `_undet` locus (undetermined allele) also get a `nearest_reference`? The proposal says yes: the column is about the system, not the allele.
- Housekeeping noted in the issue: the naming note writes `Pezio` in one line; the code is `Peziz`. Fixed in the same PR as the column, not here.
