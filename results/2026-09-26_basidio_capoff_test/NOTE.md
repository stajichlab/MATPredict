# Cap-off test on phylum-fallback Basidiomycota genomes

2026-09-27. Same frozen code as the full run (`run-ad1f865`), cap off
(`DETECT_ARGS="--max-polished-clusters-per-family 0"`). Jobs 29118039-41,
`short`, 16-way, 22-25 min each.

## Sample

50 genomes (`selected.tsv`, chosen by `select.py`): uncalled in the full run,
phylum-fallback routing, >= 1 polish-capped cluster, < 750 Mb, stratified by
order and taking the most-capped genomes first. Sporidiobolales 6,
Trichosporonales 6, Wallemiales 5, Microbotryales 5, Filobasidiales 5,
Cystofilobasidiales 5, Cystobasidiales 4, Auriculariales 4, Cantharellales 4,
Holtermanniales 2, Kriegeriales 2, Dacrymycetales 2. Eligible pool: 544
genomes. No Pucciniales genome qualified (none had a capped cluster).

## Result

- 50/50 reports. **2 of 50 genomes gain a call** (4%). 0 timeouts.
- Runtime: median 440 s capped vs 596 s cap-off (x1.25 median, x2.39 max);
  total 7.3 h vs 9.1 h.

| Genome | Family | Call | Evidence (identity / status) |
|---|---|---|---|
| *Rhodotorula sphaerocarpa* (Sporidiobolales) | bLocus | undetermined, high, idiomorph_gene_only | bE 60.6% polished; bW 32.3% polished |
| *Microbotryum violaceum* (Microbotryales) | bLocus | undetermined, high, idiomorph_gene_only | bE 26.6% polished; bW 23.6% polished |

Both gained calls sit on a locus the capped run had withheld
(`overlaps_capped_withheld`). Rhodotorula looks like a real bE/bW (HD)
locus: bE at 60.6% with a modelled bW partner. Microbotryum is weak (both
genes under 27%); Microbotryum has known HD genes, but at these identities
the call cannot be told from a paralog without synteny or a tree.

## Reading

- Extrapolated to the 544-genome pool: about 22 genomes might gain a call,
  roughly half of them weak like Microbotryum. Not measured beyond 50.
- The cap is not what keeps the fallback orders uncalled. The references are.
- Out of scope, noted: a call with both core genes at 24-27% identity gets
  confidence `high`. That tier may be too generous for `idiomorph_gene_only`
  calls in fallback routing.

## Recommendation

Do not re-run the fallback orders for the cap. The new cap rank in PR #9 may
already recover the Rhodotorula-type case; curation of Sporidiobolales /
Microbotryales references would help more.

Files: `select.py`, `selected.tsv`, `lists/`, `jobs.tsv`, `compare.py`,
`per_genome.tsv`, `gained_calls.tsv`, `capoff_0*/` (runs, rollout).
