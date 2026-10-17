# Raw-read datasets for MATPredict reads typing

One folder per genus. Each holds a `manifest.tsv` (every read set found), notes on it, and, for the SRA genera, `selected.tsv` and `download.sh`.

| Folder | Source | Status |
|---|---|---|
| `coccidioides/` | Reads already on hpcc (`2025_All_Cocci/Assembly/input`), joined to the assembly freeze | manifest built; no reads copied |
| `cryptococcus/` | ENA/SRA, Illumina paired WGS | manifest under review |
| `fusarium/` | ENA/SRA, Illumina paired WGS (*F. graminearum* complex and other species) | manifest under review |

## Coccidioides

`build_manifest.py` joins three hpcc files (read-only copies; the commands are in the script header):

- `.../2025_All_Cocci/Assembly/samples.csv` - one row per read set: strain, centre, locus tag, and the species label of that file
- `.../NovInvenio_Investigations/studies/fungi/coccidioides_pangenome/species.csv` - corrected species and the freeze assembly path for 529 strains
- `ls -L` of `.../Assembly/input` - file names and sizes

Species corrections (`coccidioides_pangenome/CORRECTIONS.md`, 2026-09-21), shown with `relabelled = yes`:
`485B-1_L_OLD_CPA0023` immitis to posadasii, `B3476` posadasii to immitis. `samples.csv` still carries the old labels and the freeze file names (`Coccidioides_<species>_<strain>.*`) also still carry the wrong species word; `species.csv` has the corrected labels.

Counts: 558 read sets (383 *C. posadasii*, 175 *C. immitis*), about 1.23 TB of gzipped FASTQ. 529 have an assembly in the freeze; 29 do not.
`TX13` is labelled *Fusarium solani* in `samples.csv`; it is dropped (decided 2026-10-17).
