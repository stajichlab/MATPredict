# 018. Non-genome records deposited as genome assemblies

- **Category:** assembly-or-annotation artefact
- **Status:** verified (NCBI Datasets)
- **Lineage:** various

## Summary
Nine NCBI "Complete Genome" assemblies from one BioProject are single Sanger
amplicons of 356-655 bp. Two other BFD entries are not genomes either. All are
suppressed in the BFD list.

## Evidence
- PRJEB104476, assembly method "Chromas", one sequence each:
  GCA_986280975.1 (DAH1005FM, 655 bp; Lichtheimia ramosa), GCA_986280995.1,
  GCA_986281615.1, GCA_986280865.1, GCA_986280875.1, GCA_986281355.1,
  GCA_986281505.1, GCA_986281625.1, GCA_986281645.1 (356-655 bp).
  NCBI Datasets dataset_report (queried 2026-09-26 and 2026-09-27).
- GCA_046252445.1 (C. auris 27-CA): 12 kb of repeats (entry 012).
- GCA_975972555.1 (Tubeufia hainanensis metagenome bin): 86,093 bp in 4 contigs.
- All in `/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/data/curation/suppress.txt`.

## Method that found it
An empty genome file (55-byte gzip) in the early-diverging scan, traced to the
NCBI source and dataset report.

## Verification done / still open
Done: NCBI metadata. Open: report to NCBI.

## Limits
Only records encountered in these runs were checked.

## Related
MATPredict honours the BFD suppress list (`detect/suppress.py`).
