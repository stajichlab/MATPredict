# Offline taxonomy check (2026-10-04)

`detect --taxid 669874` on P. tannophilus GCA_001661245.1 with
`MATPREDICT_TAXONOMY` (table from NCBI taxdmp_2026-10-01.zip, 3,016,750 taxa,
23.7 MB), `MATPREDICT_OFFLINE=1` and an empty NCBI cache (`run.slurm`, job
29396722, 8 min 47 s, 306 MB peak).

- NCBI cache after the run: 0 files (no E-utilities call).
- Report: `genetic_code: 26` (from the table), `routing_mode: phylum_fallback`,
  `routing_error: null`, `taxonomy_source: local NCBI taxonomy snapshot 2026-10-01`.
- Calls identical to the online run of the same commit's code
  (`../2026-10-04_alaninales_v061_exostring/`): MAT MAT1-1 high, MATsc MATa
  medium, MATyl A medium, MTL alpha medium, same coordinates.
- Table vs 8,185 cached efetch answers: see `tests/db/test_local_taxonomy.py`.
