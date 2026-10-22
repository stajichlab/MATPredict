# Pilot list for the clean full-BFD run (phase 3 of `docs/notes/2026-10-22_clean-run-plan.md`)

Built 2026-10-22 by `build_pilot_list.py` (deterministic, seed 20261022; run on the HPCC, it reads `samples.csv`, `tables/asm_stats.parquet`, the v0.6.0 `genomes.tsv` of the Ascomycota and Basidiomycota campaigns, and the *Mucoromycotina* `inputs.tsv`/`bfd_kept.tsv`).
`build.log` is the builder's output. A genome counts only if it has a taxid, a genome file in `input_clean_genomes/` and a v0.6.0 row; the suppress list is applied by the worker at run time (`detect suppress-filter`), not here.

## Files
- `pilot_list.tsv`: **569 genomes** (462 Ascomycota, 107 Basidiomycota), `ASMID<TAB>TAXID`, for the per-genome `detect` path (`scripts/submit_binned_panel.py --sample-list`).
- `pilot_meta.tsv`: the same genomes with phylum, class, order, family, genus, species, assembly size and bin, the v0.6.0 status and routing, and the strata each was drawn for.
- `pilot_mucoromycotina.tsv`: **20 genomes** for the separate *Mucoromycotina* path (`detect --phylum Mucoromycota --emit-cds-fasta`, no taxid): 10 BFD plus 5 LCG and 5 Jena (the last two are the held-out sets, `source` column). **The BFD share is 10, not 15**: `bfd_kept.tsv` contains fewer eligible rows than the target after the suppress list.

## Design
1. **Base, stratified by class** like the 246-genome v0.6.0 pilot: about 180 x class share for Ascomycota and 60 for Basidiomycota, at least 3 per class, round-robin over orders inside a class. The 3-per-class floor rounds the totals up (228 and 87 drawn; all 30 Ascomycota classes and 20 Basidiomycota classes are covered).
2. **Large genomes (>= 500 Mb)** to exercise the 14,400 s path: 10 Basidiomycota (Pucciniales). **No Ascomycota genome at or above 500 Mb was available**, so the planned 5 were not drawn.
3. **Known problem clades**, each with a reason to expect changes: *Bipolaris*/*Curvularia*/*Exserohilum* 40 (30 uncalled in v0.6.0; the scoped-roster question), *Pyrenophora* 25 (homothallic candidates), *Aureobasidium* 25 (fragment models), Sclerotiniaceae (*Botrytis*, *Botryotinia*, *Sclerotinia*, *Monilinia*, *Stromatinia*) 30 (remnant and fused loci), Cordycipitaceae/Clavicipitaceae/Ophiocordycipitaceae 25 (three-exon `MAT1-1-1`), Dipodascales/Lipomycetales 25 (the guard that removes calls), *Ceratocystis* 10, *Komagataella* 8, *Ogataea* 8, Xylariales 15 (absent MAT), *Alternaria* 10, *Parastagonospora* 5, *Fusarium* 10, *Hemileia* 5 and *Puccinia* 5 (`rustHD`).
4. Included on purpose: 10 genomes that failed and 1 with no report in v0.6.0 (genetic code 26), to confirm the fix.
Totals: 569 + 20 = 589 genomes; 403 called, 155 uncalled, 10 failed, 1 no_report in v0.6.0; 2 genomes belong to two strata. At the v0.6.0 pilot mean of 158 s per genome this is about 25 CPU-hours (about 25 minutes at 64-way); the new rules may be slower, which is what the pilot measures.

## Use
`DETECT_ARGS="--partial-core-policy coverage_profile --confidence-cap-policy off --strong-single-gene-policy off --homothallic-policy both_genes" python3 scripts/submit_binned_panel.py --sample-list results/2026-10-22_clean_run_pilot/pilot_list.tsv --out <results dir> --worktree <frozen worktree> --name pilot --dry-run` (check the printed `sbatch` commands and the job log's `detect_args=`, then submit without `--dry-run`). The *Mucoromycotina* list goes through `results/2026-10-03_mucoromycotina_mat/run_chunk.slurm`.
