# Genomes listed as "(unassigned)" on the results page: name rule and search

2026-11-12. Code: `scripts/build_site_data.py` (`curated_name`), `web/intro/results.html`, `web/intro/charts.js`.
Data: BFD `samples.csv` (23,679 genomes) and the v0.9.0 tables.

## Problem
The site took the first word of the recorded species name as the genus. A name that starts with a phylum
name or a qualifier gave `(unassigned)`. `Ascomycota sp. Antarctolithica obscura` (8 genomes, MAGs) was grouped
under `(unassigned)`, and a genus search for "Antarctolithica" found nothing. The search box matched only the
name at the selected rank.

## Rule (`curated_name(recorded) -> (genus, species, name_status)`)
1. `<higher taxon> sp. <Genus-like> <epithet>` (for example `Ascomycota sp. Antarctolithica obscura`): the
   submitter's name is used. `name_status = provisional`. It is not a validated genus.
2. A leading qualifier (`uncultured`, `cf.`, `aff.`, `Candidatus`, `unidentified`, `fungal`) is dropped. If a
   genus-like word is left, it is the genus. `name_status = stripped`. A rank name left after a qualifier
   (`uncultured Hypocreales`) is not a genus.
3. Otherwise the old reading stands: `ok` for a normal name, `unassigned` for a placeholder or an identifier.
   `Ascomycota sp.` and `Ascomycota sp. Tripp 6053` stay unassigned.

`name_status` is a new last column of the per-genome tables. The recorded name stays in `species_recorded`.

## Measured effect (all of `samples.csv`)
22 of 23,679 genomes change. All 22 were `(unassigned)` before. No other genome changes.

| Genomes | Before | After |
|---:|---|---|
| 8 | (unassigned) | Antarctolithica obscura (provisional) |
| 3 | (unassigned) | Caloplaca sp. |
| 2 | (unassigned) | Porphyrellus sp. |
| 1 each | (unassigned) | Muellerella sp., Tylopilus sp. |
| 1 each | (unassigned) | genus only: Alectoria, Ascochyta, Malassezia, Penicillium, Pseudogymnoascus, Saccharomyces, Sordaria |

In the v0.9.0 site data, genomes under `(unassigned)` drop from 52 to 30 (26 Ascomycota, 4 Basidiomycota).

## Search
A row matches when the text is in the name at the selected rank, or (3 or more characters) in another rank or in
the recorded species name. The table counts only matching genomes and shows how many matched through another
field. Drilling into a group keeps the search on. The per-genome download, the locus maps and a new list
"Genomes without a usable genus" use the same rule. `taxa.json` has a new field `recorded` (names that were
changed or not used for the genus: 27 entries; the file grew from 642,513 to 661,701 bytes).

## Open (curator)
- 235 genomes in 76 "genera" are rank names, for example Boletaceae 19, Xylariaceae 16, Saccharomycetaceae 14,
  Debaryomycetaceae 13, Sordariomycetes 10. They come from names such as `Boletaceae sp. KM206666`. The old
  reading (genus = the family name) is kept. They are findable today under that word. Moving them to
  `(unassigned)` is a separate decision.
- The provisional genus rests on the submitter's name. It could be a collector name or a placeholder.
- The page was not tested in a browser; the page script passes `node --check` and the Python tests pass.
