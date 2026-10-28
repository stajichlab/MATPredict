# Idiomorph names: a proposal (2026-10-28)

Status: proposal. Applied on output only for the P/M system (see below); no database or report change.

## Why
A bare `M` or `P` is ambiguous: Mucoromycota `sexM` / `sexP` are Minus / Plus, fission yeast P and M are cell types. The yeast names are also inconsistent: the same idiomorph is written `A`, `a`, `MATa` or `MTLa`, and the same string means different systems in different families.

## What the v0.7.0 Ascomycota loci carry (family_called: idiomorph names, loci)
| family | names | loci | system (reference) |
|---|---|---:|---|
| `Ascomycota:MAT` | MAT1-1 / MAT1-2 | 12,049 | Pezizomycotina (Turgeon-Yoder names, keep) |
| `Ascomycota:MATsc` | MATa / MATalpha | 4,580 | *S. cerevisiae* |
| `Ascomycota:MTL` | A / alpha | 2,750 | *C. albicans* (MTLa / MTLalpha) |
| `Ascomycota:MATkp` | MATa / MATalpha | 200 | *Komagataella* |
| `Ascomycota:MATog` | MATa / MATalpha | 66 | *Ogataea* |
| `Ascomycota:MATyl` | A / B | 79 | *Yarrowia lipolytica* (MATA / MATB) |
| `Ascomycota:MATdc` | A / B | 70 | reference species to be confirmed |
| `Ascomycota:MATtub` | MAT1-1 / MAT1-2 | 142 | *Tuber* (standard names, keep) |
| `Ascomycota:mat1`, `mat2`, `PM` | P / M | 195 | *S. pombe* P/M system (Taphrinomycotina) |

`MATa` is used by three families and `A` by three, with different meanings; a table that mixes them (the landing page, the explorer, the reads pages) cannot tell them apart without the family.

## Rule
The label names the reference system the idiomorph is defined in, not the species it was found in, like the family name `MATsc`:
`<three letters of the genus><three of the species><cell type>`; the genes follow (`SchPomMc`, `SchPomPi`, ...). A homologue in another species takes the same label. Pezizomycotina `MAT1-1` / `MAT1-2` stay as they are.

## Applied now (output only): `SchPomP`, `SchPomM`
`scripts/idiomorph_labels.py`: candidate tables and FASTA headers (`results/2026-10-28_taphrinomycotina_mat_candidates/`, columns `ref_gene_label`, `idiomorph_label`; `best=SchPomMc(PM)` in the headers), the reads pages (species with the P/M panel). Raw names stay in the tables, the report YAML and the database.

## Suggested next (not applied)
| now | proposed label | genes |
|---|---|---|
| MATsc MATa / MATalpha | `SacCerA` / `SacCerAlpha` | `SacCerA1`, `SacCerA2`, `SacCerAlpha1`, `SacCerAlpha2` |
| MTL A / alpha | `CanAlbA` / `CanAlbAlpha` | `CanAlbA1`, `CanAlbA2`, `CanAlbAlpha1`, `CanAlbAlpha2` |
| MATkp MATa / MATalpha | `KomPasA` / `KomPasAlpha` | |
| MATog MATa / MATalpha | `OgaPolA` / `OgaPolAlpha` | |
| MATyl A / B | `YarLipA` / `YarLipB` | |
| MATdc A / B | (confirm the reference species first) | |
| mat1 / mat2 / PM P / M | `SchPomP` / `SchPomM` (done on output) | |

## Database change (needs a curator ruling)
Changing the stored names touches the `idiomorphs` vocabulary of the metadata schema, the idiomorph values and gene names of the records of these families (a `record_version` bump each), the family registry, the report text and tests, the site data and explorer, and the About page. Output-level labels avoid all of that; a table that keeps both (`idiomorph`, `idiomorph_label`) is the cheapest step for the other families.
