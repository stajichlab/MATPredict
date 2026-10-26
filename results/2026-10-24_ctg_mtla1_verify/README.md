# CTG-clade MTLA1 candidate verification (2026-10-24)

Working evidence for `db/candidates/Ascomycota/43968_nrrl-yb-4298_MTL_A`.

- `verify.tsv`: 14 candidate MTLA1 proteins (BFD annotated genomes of CTG-clade genera that lost a polished MTLA1
  model) with their gene-order neighbours (offset:label:percent identity to the curated C. lusitaniae / C. auris
  PAP1, OBP1, PIK1, MTLA2 proteins).
- `pfam_candidates.txt` / `pfam_candidates.py`: Pfam-A hmmscan (pyhmmer, E <= 1e-3) of the candidates and their
  neighbours. The five *Scheffersomyces*/*Priceomyces* candidates hit `MATalpha_HMGbox`, so they are the alpha-type gene,
  not MTLA1, and were not turned into records. The *Aciculoconidium* (43968) and *Hyphopichia* (717740) candidates hit
  `Homeodomain`, like all five curated MTLA1 reference proteins; the curated MTLA2 proteins have no Pfam hit.

## Serinales regression (records added to a copy of the v0.6.0-era database, clean-run settings)

Compared with the clean-run baseline on the 737 genomes that changed in a seven-record trial or are CTG/source genomes:

| records added | confidence up | confidence down | alpha->A flips | CTG-clade genomes fixed (of 16) |
|---|---|---|---|---|
| all seven candidates | 76 | 157 | 274 | 8 (4 flipped idiomorph) |
| Hyphopichia + Aciculoconidium | 85 | 14 | 0 | 14 |
| Hyphopichia only | 36 | 12 | 0 | 6 |
| Aciculoconidium only | 69 | 3 | 0 | 14 |

Only *Aciculoconidium* is proposed. The *Hyphopichia* record adds MTLA1/MTLA2 calls inside alpha loci in six genomes and
drops five A loci to medium. The three remaining *Aciculoconidium* drops are A loci whose MTLA1 hit was already a
weak fragment (17-25% coverage, 33-36 bits) and scores 42-44 bits with the new record; the gene set is unchanged.
The other ~1,630 Serinales genomes were not re-run with the single record.
