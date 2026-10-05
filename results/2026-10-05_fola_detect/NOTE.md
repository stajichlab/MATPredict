# Fola (F. oxysporum f. sp. lactucae): MATPredict assembly calls vs read-based calls (2026-10-05)

Code: run-94c1a3b (adds the F. oxysporum MAT1-1/MAT1-2 records, AB011379.2 / AB011378.1), `detect --taxid 5507` (`run.slurm`).
Assemblies: /bigdata/stajichlab/nicolel/Fola/genomes/ (19 files; 15 not readable by this account, so no report; see `runs_*.log`).
Read-based calls: N. L.'s samtools coverage of the same two references (06_Align/Mating_Types/Coverage/).
The read mapper she used is not recorded in the CRAM headers. Rule used here: an idiomorph is present when
breadth >= 90%. This rule reproduces her table (MAT2 129, MAT1 17, Both 1, Neither 1; n = 148). See `reads_vs_assembly.txt`.

| Strain | MATPredict (assembly) | Reads: MAT-1 breadth / depth | Reads: MAT-2 breadth / depth | Read call | Agree |
|---|---|---|---|---|---|
| VSP-0916 (flye subset assembly) | MAT1-2 high, mat_locus | 11.6% / 23.2x | 99.6% / 201.4x | MAT2 | yes |
| VSP-0980 (AVITI assembly) | MAT1-1 high, mat_locus | 100.0% / 44.3x | 24.5% / 10.5x | MAT1 | yes |
| JCP043 | MAT1-2 high, mat_locus | 9.6% / 2.5x | 99.6% / 69.5x | MAT2 | yes |
| AT141 | MAT1-2 high, mat_locus | - | - | no reads locally | - |

JCP043 reads (SRR28734937) are not in N. L.'s coverage set. They were typed here with minimap2 -ax sr and
samtools coverage (`jcp043_reads.slurm`, `JCP043.MAT.coverage.txt`).

Notes:
- 3 of 3 strains with both an assembly call and reads agree.
- The idiomorph that is not carried still has 10-25% breadth. These are the parts of each reference that are
  shared with the flanks. Tool A must mask these parts (see the A. fumigatus ~270-bp shared region).
- Rerun the other 15 assemblies after N. L. gives group read permission (`run.slurm`).
- Tracked here: scripts, logs, `detection_report.yaml` and `detected_loci.gff3` for the 4 typed assemblies.
  Not tracked: `runs/*/_reference.faa` (a copy of the database proteins, about 0.6 MB each).
