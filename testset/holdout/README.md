# Ascomycota held-out manifest

`ascomycota_holdout.tsv` lists every curated record under `db/Ascomycota` (72) with its subphylum, family
(locus name), mating system (`mating_type.system`), idiomorphs and assembly accession, and 26 records have an
assembly, so only those can be held-out genomes. For the genus, family and order radius it gives how many
curated records the holdout withholds and how many references of the same family and of the same idiomorph
remain. Rebuild with `python scripts/build_asco_holdout_manifest.py --db db --out testset/holdout/ascomycota_holdout.tsv`.

Run it with `scripts/run_holdout_benchmark.py` (it withholds by radius, runs detect on the record's genome and
scores the idiomorph against the record's own). A radius that leaves no reference of the record's idiomorph is
`no_reference_idiomorph`, a gap in the database, not a miss.

Testable records: Pezizomycotina `MAT` 6 heterothallic and 3 homothallic; Saccharomycotina `MATsc` 5, `MATyl` 1,
`MTL` 6; Taphrinomycotina `PM` 3 and `mat1` 2. In Saccharomycotina `MATsc`, `MATyl` and Taphrinomycotina `mat1`
the family and order radius leave no reference of the same idiomorph, so only the Pezizomycotina `MAT` set (and
`MTL` at genus radius) can be tested at the wider radii.
