# Idiomorph names: a proposal (2026-10-28)

Status: proposal. Applied on output only for the P/M system (`scripts/idiomorph_labels.py`); no database or report change.

## Why
A bare `M` or `P` is ambiguous: Mucoromycota `sexM` / `sexP` are Minus / Plus, fission yeast P and M are cell types. The yeast names are also inconsistent: the same idiomorph is written `A`, `a`, `MATa` or `MTLa`, and the same string means different systems in different families.

## What the v0.7.0 Ascomycota loci carry (family_called: idiomorph names, loci)
| family | names | loci | system (reference) |
|---|---|---:|---|
| `Ascomycota:MAT` | MAT1-1 / MAT1-2 | 12,049 | Pezizomycotina (keep) |
| `Ascomycota:MATsc` | MATa / MATalpha | 4,813 | *S. cerevisiae* |
| `Ascomycota:MTL` | A / alpha | 2,766 | *C. albicans* (MTLa / MTLalpha) |
| `Ascomycota:MATkp` | MATa / MATalpha | 200 | *Komagataella* |
| `Ascomycota:MATog` | MATa / MATalpha | 66 | *Ogataea* |
| `Ascomycota:MATyl` | A / B | 79 | *Yarrowia lipolytica* (MATA / MATB) |
| `Ascomycota:MATdc` | A / B | 70 | reference species to be confirmed |
| `Ascomycota:MATtub` | MAT1-1 / MAT1-2 | 142 | *Tuber* (keep) |
| `Ascomycota:mat1`, `mat2`, `PM` | P / M | 195 | *S. pombe* P/M system |
(871 loci are `undetermined`.) `MATa` is used by three families and `A` by three, with different meanings.

## Grammar (curator, 2026-10-28)
`<Lineage>__<System>_<Cell>`: lineage = subphylum (`Taphr`, `Sacch`, the full word `Mucoromycotina`), system = a code for the REFERENCE ORGANISM (three letters of the genus, three of the species), cell last: `Taphr__SchPom_M`, `Taphr__SchPom_P`. The label names the reference system, not the species the gene was found in, and not a clade or a genus (genus names change; the Candida genera were split and do not predict the system). Genes carry the system code and the gene only (`SchPom_Mc`); the lineage is its own column. A gene that is only similar to a reference gene, not in a locus the pipeline reported, gets the suffix `_like` (`SchPom_Mc_like`; column `label_basis`: `pipeline_detected` or `similarity_only`).

## Applied now (output only)
`scripts/idiomorph_labels.py`: candidate tables and FASTA headers of `results/2026-10-28_taphrinomycotina_mat_candidates/` (columns `ref_gene_label`, `idiomorph_label`, `lineage`, `label_basis`; raw `ref_gene` kept) and the reads pages (branch `wave2`).

## Review by a second model (Opus, 2026-10-28) and what I checked
- Use a reference-organism code, not a clade code. Checked: MTL calls are 2,541 of 2,766 in Serinales but also 98 in Pichiales, 66 in Phaffomycetales, 38 in Ascoideales; MATsc covers 4,145 in Saccharomycetales plus 382 in Pichiales and 188 in Phaffomycetales. An order code (`Ser`) would be wrong for the off-clade calls.
- The Saccharomycotina "systems" are not separate homology groups (MTLalpha1 and MATalpha1 are orthologs); they differ in gene content and in which panel scored best. `Sacch__SacCer_a` on a *Pichia kudriavzevii* locus means "called by the S. cerevisiae-type panel", not a clade or homology claim; say so in the About text or add a `label_basis = panel` column.
- Mixed genera look like cross-calls: of 109 *Ogataea* genomes with loci, 29 have both MATsc and MTL; of 136 *Pichia*, 27 (the reviewer said 35 for *Ogataea*; I count 29). Their MTL loci are nearly all `idiomorph_gene_only`. Label them with the genome's primary system plus a flag (`secondary_crosscall`), or re-route; do not give them a CanAlb label.
- Formerly-Candida genera: label by the system called, never the genus: *Candidozyma*, *Clavispora*, *Meyerozyma*, *Lodderomyces*, *Metschnikowia*, *Debaryomyces* -> CanAlb (MTL); *Nakaseomyces* (*C. glabrata*) -> SacCer; *Pichia* mostly SacCer.
- Cell token: lower case `a` / `alpha` for a/alpha systems (the stored capital `A` in MTL is the inconsistency to fix); `P` / `M` upper case; Yarrowia `YarLip_A` / `YarLip_B` with the a/alpha equivalence recorded in the registry; an explicit `undet` token for `undetermined`.
- Parseability: `__` then `_` splits cleanly; under strict Newick an unquoted `_` is a space, so quote tip labels or test the viewer; avoid `+ - ? ( )` in tokens; the viewer should group by the registry's lineage field, not by parsing the string. The reviewer suggested fixed-width 5-letter lineage codes (`Pezio`, `Sacch`, `Taphr`, `Mucor`); the curator chose the full word `Mucoromycotina`: open.
- Mucoromycotina: reference *Phycomyces blakesleeanus* (`PhyBla`; sexP/sexM first defined as a locus there), cell token `plus` / `minus` (`sexP` / `sexM` are gene names, the same P/M collision). Basidiomycota HD and P/R loci are multiallelic, not two idiomorphs; the grammar needs a `<locus>-<allele>` slot there and should not be promised to generalise yet.

## Proposed labels (not applied)
| now | proposed | genes |
|---|---|---|
| MATsc MATa / MATalpha | `Sacch__SacCer_a` / `Sacch__SacCer_alpha` | `SacCer_a1`, `SacCer_alpha1`, `SacCer_alpha2` |
| MTL A / alpha | `Sacch__CanAlb_a` / `Sacch__CanAlb_alpha` | `CanAlb_a1`, `CanAlb_a2`, `CanAlb_alpha1`, `CanAlb_alpha2` |
| MATkp MATa / MATalpha | `Sacch__KomPha_a` / `Sacch__KomPha_alpha` | |
| MATog MATa / MATalpha | `Sacch__OgaPol_a` / `Sacch__OgaPol_alpha` | |
| MATyl A / B | `Sacch__YarLip_A` / `Sacch__YarLip_B` | |
| MATdc A / B | confirm the reference species first | |
| Mucoromycotina sexP / sexM | `Mucoromycotina__PhyBla_plus` / `_minus` (reference to confirm) | `PhyBla_sexP`, `PhyBla_sexM` |
| mat1 / mat2 / PM P / M | `Taphr__SchPom_P` / `Taphr__SchPom_M` (applied on output) | `SchPom_Pc`, `SchPom_Pi`, `SchPom_Mc`, `SchPom_Mi` |

## Database (needs curator rulings)
Do not store labels in records. Compute them at output time from (family_called, raw name) with a registry mapping (`system_code`, `lineage`, raw name -> cell token, raw gene -> gene label, the a/alpha equivalence for A/B systems); no `record_version` bump and no schema break. If records must carry labels later, add an optional `label` next to `name` in each `idiomorphs` entry. Rulings needed: each system code; the MATdc reference; Yarrowia A/B -> a/alpha; the Mucoromycotina reference and lineage token; whether MAT1-1 / MAT1-2 stay unprefixed or become `Pezizo__MAT1-1`; the policy for secondary cross-calls; five-letter against full-word lineage codes.
