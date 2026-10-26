# Pilot plan: raw-read MAT typing of public fungal read sets on the HPCC

Draft 2026-10-26 for issue #173. Test 1 (20 *N. crassa* runs) has been run: see `pilot/ncrassa20/README.md`. The transfer section below was rewritten after it: the source is NCBI's copy of the SRA data on AWS open data, not ENA. Numbers marked (est.) are estimates, not measurements.

## What the pilot answers
1. Does the pipeline work end to end on public runs: stream, subsample, type, delete, keep the observations?
2. How long and how much transfer does one run take, so the full survey can be budgeted?
3. How well do the read calls agree with an independent call (assembly `detect`, or a mating type stated in the strain name)?
4. How many strains of a species are really one genotype, and can the rest be skipped?

## Species (about 50 strains each, one run per strain)
| Species | Tier | Why |
|---|---|---|
| *Neurospora crassa* | A | Both idiomorphs recorded. About 1,300 strains, many single-end (about 320) and many morphological mutants from a few ancestral genotypes. A test of the identical-strain rule. Louisiana wild isolates (61, PRJNA1108904) are the clean natural set; mating type is often in the strain name (e.g. FGSC 2489 mat A). |
| *Zymoseptoria tritici* | A | Large natural population panels (hundreds per study). Tetrads and mutation-accumulation sets are left out by the filter. |
| *Botrytis cinerea* | A | MAT1-1 / MAT1-2 records, field isolates. |
| *Aspergillus flavus* | B | Genus records only: tests building a panel without a species reference. |
| *Verticillium dahliae* | A | Records exist (curated in 2026-10). Mostly clonal lineages: tests the identity rule again. |
| *Sclerotinia sclerotiorum* | A | Homothallic: every strain should read `both` (a control). 270 strains. |
| ***Pyricularia oryzae*** | C | No records for *Pyricularia* or *Magnaporthe*. 1,474 strains; BFD has 495 assemblies with 471 `detect` MAT calls (MAT1-1 227, MAT1-2 244), so the panel can be cut from assemblies as for *N. crassa* and the read calls checked against `detect`. Clonal lineages and the wheat-blast pathotype to watch. |
| *Cryphonectria parasitica* | C | Second no-reference species (curator choice): 409 strains, 87 BFD assemblies called (MAT1-1 46, MAT1-2 41). |

## Panels (before any download)
`reads-type` needs one FASTA per idiomorph and, for the single-copy depth control, a shared flank. The *C. posadasii* panel failed because its two records had no shared flank (every strain `low_depth`). For each pilot species: build `panels/<species>/<prefix>_<IDIOMORPH>.fasta` from the accepted records, trim to the locus plus shared flank, and check on 5 strains with an assembly before the array job (as the Coccidioides pilot did). The idiomorph names follow the species (MAT1-1 / MAT1-2, MATa / MATalpha). Panel building is the manual step; everything below is mechanical.

## What runs on a compute node (one array task = one run)
`pilot/stream_type_sra.slurm` (tested on 20 runs). Per task:
1. **Download the run's `.sra` from NCBI's copy on AWS open data** (`https://sra-pub-run-odp.s3.amazonaws.com/sra/<run>/<run>`; resumable `curl`, 3 tries) into `$SCRATCH`. No login; ENA FASTQ was 100 times slower from the HPCC (see below).
2. **Cut evenly spaced windows of spots** with `fastq-dump -N start -X end -Z --split-spot` on the local `.sra` (8 windows, 2 million spots in all; a run with fewer spots is read whole). Random access reads along the whole run, so a position-sorted old run is no longer sampled from its start only. Paired and single-end runs use the same code. Then the `.sra` is deleted.
3. **Type** against the species panel (`matpredict reads-type`), one result row.
4. **Keep, then delete.** Kept per run: `<run>.<panel>.tsv` (call, breadth, depth, shared depth, flags) and `<run>.run.json` (accession, layout, spots, windows, bytes, seconds per stage, state). The `.sra` and the FASTQ are removed by the exit trap; nothing but these small files goes to `/bigdata`.
5. **Not done yet:** the identity sketch (below), an md5 check (the `.sra` is not checked against a published sum; a truncated download fails in `fastq-dump`), and a retry ledger beyond the `run.json` state.
6. **Identity sketch (optional, for the identical-strain question).** From the same subsample write a small k-mer sketch so strains can be clustered afterwards; the pilot's result decides the skip rule below.

## Speed and transfer (measured 2026-10-26, from the HPCC `exfab` partition)
| Source | Result |
|---|---|
| ENA FASTQ over HTTPS/FTP | 0.13 MB/s (5 MB in 40 s) and 0.9 MB/s on two files: not usable at scale |
| NCBI `sra-instant` FTP path | 404 for the run tested |
| NCBI copy on AWS open data (`.sra`) | 20 to 26 MB/s per stream; 4 streams 76 MB/s; the two ERR (European) runs are there too |
- **Test 1:** 20 runs, 22.7 GB, 17 MB/s per stream on average, median 108 s per run, the 19-run array 10.7 minutes at 4 tasks at a time (16 CPUs), peak memory 4.9 GB.
- **Whole survey (est.):** about 1.1 GB of `.sra` per run; 16,500 runs at 50 strains per species is about 18 TB, about 2 minutes a run: at 4 tasks at a time about 5.7 days; at 8 tasks (32 CPUs, inside the under-60 limit) about 3 days if the link scales. The pilot species (about 300 runs) are about 2.5 hours.
- **Cost:** nothing is downloaded from a paid cloud. The AWS bucket is NCBI's open-data copy (no account; the objects we fetched are not requester-pays); GCP was not used. Moving the work into a cloud region only makes sense with a budget.
- **ENA** remains the way to see runs not yet in the SRA copy; the slowness may be the route from this campus at this hour, so it can be retried.
- **Failures.** A broken download is retried 3 times; a task that still fails writes `state: failed_*` in its `run.json`; re-running the same array does only the missing runs (the result file is the done-marker).
- **Be polite.** One run per strain, no re-downloads, 4 to 8 streams at a time.

## Choosing the runs
- One run per strain; the biggest paired run, else the biggest single-end run.
- At most 50 per species, spread over studies (no more than 5 per BioProject), so one project does not stand for the species. For the *N. crassa* mutant set, at most 2 per project and no more than 5 per strain-number prefix until the sketches show how many ancestors there are.
- The curator's rule (2026-10-26): strains from several studies are pooled when they are the same species; identical strain sets need only one member. Metadata can only catch an identical strain string. Sketches catch the rest, after the fact; once the pilot shows the clusters, later runs skip all but one member before downloading.

## Calls, truth and what is reported
- Read call per strain: one idiomorph, `both`, `none` or `low_depth` (existing triage: `typed`, `both`, `review`). For a diploid or dikaryotic species `both` is expected; the report states the rule per species.
- **Truth.** (a) `detect` on the strain's assembly where one exists (as in Coccidioides 487/495 and *A. fumigatus* 293/296); (b) a mating type in the strain name or BioSample (FGSC strains); (c) simulated mixtures for sensitivity. Where there is none the call is reported as unverified.
- **Species-label check.** The shared-flank depth tells a wrong species or a failed run from a real `none`: a run with no hits on either idiomorph and no shared-flank depth is reported as `review: wrong species or low depth`. This also catches mislabelled records.
- **Per-species report.** Counts of each call, the number of independent genotypes (clusters), the studies, and what was left out and why. Observations are not frequencies for clonal sets.

## After the pilot
Decide from the numbers: typed rate (aim for 90% or more where a panel exists), agreement with truth (95% or more), seconds and MB per run, failure rate, the depth target and early-stop rule, and the identical-strain skip rule. Then write the species list and panel queue for the full run and ask for the go-ahead.

## Open
- Genome size per species for the depth target.
- Whether early-stop is safe for the older runs (test).
- A synonym table so a renamed species finds its records (*Pyricularia oryzae* / *Magnaporthe oryzae*, *Nakaseomyces glabratus* / *Candida glabrata*).
- Legacy ENA records with no study accession (381 runs; 33 are *Neurospora*, e.g. "*Neurospora crassa* 7035 Project", old per-strain JGI/BI submissions): they are real records without a project, not mislabels, but the survey groups them as one blank study; handle by sample.
