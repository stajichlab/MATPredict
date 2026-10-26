# Pilot plan: raw-read MAT typing of public fungal read sets on the HPCC

Draft 2026-10-26 for issue #173. Nothing here has been run. Numbers marked (est.) are estimates from the ENA metadata (read and base counts), not measurements; the pilot's first job is to measure them.

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
| one tier C species | C | Tests the case with no reference (expected: needs curation first; the pilot may stop at the panel step). |

## Panels (before any download)
`reads-type` needs one FASTA per idiomorph and, for the single-copy depth control, a shared flank. The *C. posadasii* panel failed because its two records had no shared flank (every strain `low_depth`). For each pilot species: build `panels/<species>/<prefix>_<IDIOMORPH>.fasta` from the accepted records, trim to the locus plus shared flank, and check on 5 strains with an assembly before the array job (as the Coccidioides pilot did). The idiomorph names follow the species (MAT1-1 / MAT1-2, MATa / MATalpha). Panel building is the manual step; everything below is mechanical.

## What runs on a compute node (one array task = one run)
Existing code: `datasets/reads/fetch_type.slurm` (download, subsample, type, delete) and `accumulate.py` (merge, triage). Changes proposed:
1. **Stream instead of download.** Today the whole FASTQ is downloaded to `$SCRATCH` and then subsampled. Instead `curl` the ENA `fastq_ftp`/HTTPS URL, decompress with `pigz -dc`, keep every Nth read (N from the run's `read_count`, so a position-sorted old file is sampled along its whole length), and write only the subsample to `$SCRATCH`. Paired mates are streamed separately with the same N, which keeps them in step. `$SCRATCH` then holds about 0.5 to 1 GB per task at the peak, not the whole run.
2. **Single-end runs** use the same path with one file (`--samples` already takes a comma-separated list; the *A. fumigatus* run had 20 single-end sets).
3. **Integrity.** In full-stream mode tee the compressed stream through `md5sum` and compare with the ENA md5; a mismatch fails the task. In early-stop mode (below) only gzip integrity of the part read can be checked.
4. **Type** against the species panel (`matpredict reads-type`), writing one result row.
5. **Keep, then delete.** Kept per run, in `<out>/<run>.<panel>.tsv` and a small `<run>.run.json`: the reads-type row (call, depths, breadth, flags), run accession, species, strain and study, layout, `read_count`, bases streamed, subsample step and reads used, wall seconds per stage (transfer, decompress, type), md5 status, and the exit state. Deleted by the exit trap: all reads. Nothing but these small files goes to `/bigdata`.
6. **Identity sketch (optional, for the identical-strain question).** From the same subsample write a small k-mer sketch (for example a 1,000-hash MinHash, a few KB) so strains can be clustered afterwards. The pilot types everything and stores sketches; the pilot's result decides the skip rule below.

## Speed and transfer
- **Transfer is the bottleneck, not compute.** One typing needs minutes; the stream of a run is hundreds of MB to several GB. Download from ENA (fastq.gz over FTP or HTTPS) or NCBI to the cluster is free; running on the HPCC avoids any cloud egress charge (the GCP option is held back until a budget exists).
- **Volume (est.).** About 0.3 gzip bytes per base. With one run per strain and a cap of 50 strains per species, the whole survey is about 16,500 runs and about 20 TB if every run were streamed in full. The five pilot species at 50 runs each are about 24 to 66 GB per species (about 200 GB in all).
- **Read only what is needed.** To type a 40 Mb genome 15 to 20x genome-wide is plenty (about 2 to 4 million pairs), far below the 8 million pairs used for Coccidioides. For runs whose file is in random order (modern platforms), stop reading after the target (`curl` is closed when the pipe ends): this cuts the transfer for a typical 1.4 to 4 Gb run to about a fifth to a half (est.). Older runs may be position-sorted, so they use the full-stream sampling; the rule is chosen from platform and date and tested in the pilot by typing the same strains both ways.
- **Per-species depth target.** Take the genome size of each species from its NCBI assembly (open item) and set the target pairs to about 20x; do not use one number for all species.
- **Concurrency.** `%8` array tasks with 4 CPUs each is 32 CPUs, inside the under-60 limit; ENA and NCBI are shared servers, so no more than about 8 simultaneous streams from the cluster. Throughput per stream (the first 20 runs measure it) sets the pace: at 20 to 50 MB/s a 1 GB stream is under a minute, so the pilot (about 300 runs) should take hours, not days (est.).
- **Alternatives to test, not assume.** NCBI "SRA Lite" files via `prefetch` (smaller, but the conversion needs the SRA toolkit and CPU); Aspera if installed (faster on long-haul links); a cloud-region copy of the SRA data if a budget is ever set (transfer inside the same cloud region is free).
- **Failures.** A broken stream cannot be resumed: retry the whole run up to 3 times with back-off, then log `failed` in the ledger with the reason. Re-running the same array only does the missing runs (the result file is the done-marker).
- **Be polite.** One run per strain, no re-downloads, a ledger so a run is never fetched twice.

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
