# STE3-like receptor inventory, Agaricomycetes (2026-10-06 study)

Protein sequences, gene-model status, taxonomy and call status for every STE3-like locus found by the array study
(`analysis/2026-10-06_agaricomycetes-pr-arrays.md`). Locus coordinates alone were saved before; this makes the genes revisitable.

## Set definitions (which file holds which set)

| Set | Definition | Genomes | Loci | Where |
|---|---|---|---|---|
| scanned | Agaricomycetes of the v0.6.0 run, assembly < 500 Mb (`scan==True` in `../agari_genomes.tsv`) | 1,853 (1,847 with >= 1 locus) | 17,129 | `../loci_all.tsv.gz` (one row per locus), `genomes_scanned.tsv.gz` |
| qpass (BUSCO, N50 only) | BUSCO complete >= 70 and N50 >= 20 kb (the `qpass` column of `agari_genomes.tsv` / `genome_table.tsv`) | 1,292 | 12,576 | derivable; NOT the headline set |
| **qpass (headline)** | the above **and contigs <= 5,000** (removes 5 genomes, 58 loci) | **1,287 (705 species)** | **12,518** (9,436 arrays) | `qpass==True` in `ste3_loci_table.tsv.gz`; `ste3_loci_proteins_qpass.faa.gz` |
| not qpass | scanned minus headline qpass | 566 | 4,611 | `qpass==False` |

`loci_all.tsv.gz` is the all-scanned set (17,129 rows), not the 12,518 quoted in the study; it has no `qpass` column and was left unchanged.
The `qpass` column in `ste3_loci_table.tsv.gz` (headline definition) is the join. `set_reconciliation.tsv` has the counts.
(`genome_table.tsv` `qpass` is BUSCO+N50 only, 1,293 genomes; the contig limit was applied inside `analyze.py` downstream.)

## Files

| File | Content |
|---|---|
| `ste3_loci_proteins.faa.gz` | 17,129 translated best gene models, one per locus, all scanned genomes; 80-column FASTA |
| `ste3_loci_proteins_qpass.faa.gz` | the 12,518 loci of the headline qpass set |
| `ste3_loci_table.tsv.gz` | one row per locus (below); `locus_id` is the FASTA identifier |
| `genomes_scanned.tsv.gz` | the 1,853 scanned genomes: qpass flags, taxonomy, assembly stats, pipeline status, `n_loci` |
| `tax_loci_arrays_genomes_species_by_{class,order,family}.tsv` | per taxon: genomes/species scanned and qpass; loci, arrays, genomes and species with loci, complete models, CAAX-flagged loci, loci in a pipeline call, loci in a withheld cluster; prefixes `all_` and `qpass_` |
| `per_genome_receptor_count_by_{order,family}.tsv` | per-genome STE3-like locus count distribution (n genomes, species, mean, median, quartiles, min, max, genomes with 0,1,2,3,4,5-6,7-10,11+ loci); sets `qpass` and `all_scanned` (denominator = all genomes of the set, zero-locus genomes included) |
| `call_status_by_family.tsv` | per family and set: genomes by pipeline `status`, loci in a pipeline call / only in a withheld cluster / not called, CAAX-flagged loci and those not in a call |
| `set_reconciliation.tsv` | row counts per set |

FASTA header: `>locus_id q=<model query> complete=<True|False>`. `locus_id` = `<genome>|<contig>:<start>-<end><strand>`. Genome FASTAs are
`/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes/<genome>.fa.gz`. Internal stop codons are written as `X`; the terminal stop is removed.

## Columns of `ste3_loci_table.tsv.gz`

- Identity and taxonomy: `locus_id`, `genome` (BFD ASMID), `qpass` (headline), `phylum`, `class`, `order`, `family` (BFD `samples.parquet`; blank family -> `unclassified`, 128 loci), `genus`, `species` (binomial, as used in the study), `species_bfd` (BFD SPECIES), `taxid` (NCBI).
- Locus (from the scan, as in `loci_all.tsv.gz`): `contig`, `contig_len`, `start`, `end` (0-based half-open, PAF), `strand`, `best_query` (query of the highest-scoring alignment at the locus), `best_qcov` (that alignment's query coverage), `locus_ident` (its identity, matches/block length), `locus_ref_ident` (best curated REF query, 0 if none), `n_queries`.
- Gene model (miniprot `--gff`, best query): `model_query` (equals `best_query` for all loci with a model), `model_start`, `model_end` (1-based GFF), `model_identity`, `model_positive`, `model_qcov` (aligned query residues / query length), `n_cds`, `prot_len` (aa, stop excluded).
- Completeness: `complete` = starts with Met AND ends in a stop codon AND no internal stop AND no frameshift. `partial_reason` = comma list of `no_start`, `no_stop`, `internal_stop`, `frameshift` (empty when complete; `no_model` if miniprot made none, which did not occur). `start_met`, `has_stop`, `n_internal_stop`, `frameshift` are the components. Partial models are kept and translated (translation of the best model, genetic code 1). Frameshift is miniprot's `Frameshift` tag or a CDS phase break; models with frameshift typically carry internal stops.
- Array and call context: `array_id` (`<genome>|<array number>`, 50 kb single-linkage arrays of `analyze.py`), `array_size`, `caax_flag` (strict CAAX precursor ORF within 10 kb), `precursor_homology_near` (`hx`, tblastn precursor homology E <= 1 within 10 kb), `d_T`, `nT_10kb`, `d_Hx`, `nHx_10kb`, `pipeline_call` (locus inside a pipeline PR call), `withheld_cluster` (inside a withheld cluster), `call_status` (`in_pipeline_call` | `in_withheld_cluster` | `not_called`; the two flags never co-occur).

## Results

- Models for 17,129 / 17,129 loci. Complete models: 2,878 (16.8%) overall, 1,991 (15.9%) of the 12,518 qpass loci. Most partial models are `no_stop` (alignment of the best query ends before a stop; the loci are the aligned region, not extended gene calls) or `no_start,no_stop` (low-identity fragments); 4,264 (3,811 with a frameshift) carry a frameshift or internal stop (pseudogene-like or assembly errors, or a mis-spliced model). Treat `complete` as conservative: it means "miniprot reproduced a full ORF from one query", not "the gene is partial".
- Complete-model length median 476 aa (range 84-1,472).

## Provenance and exact commands

Study scripts (this directory): `scan_genome.py`, `array_scan.py`, `run_scan.sh`, `make_genome_list.py`; run on HPCC in `/bigdata/stajichlab/jstajich/agari_work`.
miniprot 0.18-r281 (pixi default env of `/bigdata/stajichlab/jstajich/projects/MATPredict`).

```
# 1. models (Slurm array, 62 tasks x 30 genomes, 8 cpu, 16 GB each; 1,853 genomes, 0 failures); script re-derives the loci with scan_genome.miniprot_loci
sbatch run_extract_proteins.sh        # -> /bigdata/stajichlab/jstajich/agari_work/prot_out/<genome>.prot.{tsv,faa}
# extract_loci_proteins.py: miniprot -t8 -I --outn=100 --gff <genome.fa.gz> <best queries of the genome>; per locus choose the model of best_query
# overlapping the locus (most overlap, then identity); CDS from the genome (miniprot CDS features include the stop); translate.
# 2. taxonomy (BFD tables), on HPCC with /usr/bin/python3.12:
#    pyarrow.parquet.read_table(".../Fungi_BFD/tables/samples.parquet", columns=[ASMID,NCBI_TAXONID,PHYLUM,SUBPHYLUM,CLASS,ORDER,FAMILY,GENUS,SPECIES]), rows with CLASS == Agaricomycetes -> tax.tsv (tab, no header)
# 3. assemble (pandas, numpy), from this directory:
rsync -a --include='*.prot.tsv' --include='*.prot.faa' --exclude='*' hpcc2:/bigdata/stajichlab/jstajich/agari_work/prot_out/ prot_out/
python build_inventory.py prot_out tax.tsv
```

All 17,129 locus identifiers of the protein run match `loci_all.tsv.gz` (contig, start, end, strand) one to one; the loci were re-derived deterministically by the same miniprot call as the scan.
The per-genome raw outputs stay on HPCC (`prot_out/`, not committed). Checksums (sha256): see `SHA256SUMS` in this directory.
