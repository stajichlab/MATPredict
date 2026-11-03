# Contamination spike-in for the raw-reads screen

Question: how little of a second organism can GmG see in a read set, and is 2M read pairs enough?

The raw-reads typing streams up to 2M pairs of each run. A run that is not the expected species, or a mix of two
organisms, would still be typed (the MAT k-mers of the expected species are found in a mixed set). A screen before typing
should flag it. [GmG](https://github.com/stajichlab/GmG) classifies reads against a 23,677-genome fungal index and reports a
mixture block (`qc.mixture`: the second genus, its hashes, how many sit in one genome, `suspected`). Its own rule fires only
for a different genus that explains at least 1% of the hashes and at least 20 hashes, and its documentation says there is no
reliable mixture signal beyond that. This measures it on real reads.

## Design

Real reads are mixed at known fractions, so any prefix of a mixed file is a mix at about that fraction (`spikein.py mix`):

| | |
|---|---|
| backgrounds | *C. posadasii* GT-134 (22.3M pairs, 93 bp), *C. immitis* 3M3 (7.7M pairs, 126 bp) |
| contaminants | *Histoplasma* NM_7898 (3.2M pairs, 150 bp; a read set submitted as *C. posadasii* that GmG calls *Histoplasma*), *Lecanicillium* CA10 (labelled *C. immitis*, GmG calls it *Lecanicillium*), and the other *Coccidioides* species (a same-genus negative: a genus-level screen should not flag it) |
| fractions of contaminant pairs | 0.5, 1, 2, 5, 10, 25% |
| controls | each background alone (fraction 0) |
| depths | the first 100k, 500k, 2M and 6M pairs (`gmg run --max-reads N` counts reads per file, so N pairs) |

38 conditions. Each is one 6M-pair file pair; the realised contaminant pairs at each depth are in its `.mix.json` sidecar.

## Run

    H=datasets/reads/contamination_spikein
    python3 $H/spikein.py sheet --sources $H/sources.cocci.tsv --outdir OUT        # conditions.tsv and the GmG sample sheet
    sbatch --array=1-38%12 --export=ALL,HERE=$PWD/$H,SRC=$PWD/$H/sources.cocci.tsv,OUT=OUT $H/run_mix.slurm
    sbatch --dependency=afterok:<array job> --export=ALL,OUT=OUT $H/run_gmg.slurm
    python3 $H/spikein.py summarize --conditions OUT/conditions.tsv --sources $H/sources.cocci.tsv \
        --runs 100000=OUT/gmg/depth_100000/results.jsonl 500000=... --clean-jsonl CLEAN.jsonl --out OUT/spikein

`--clean-jsonl` is a GmG `results.jsonl` of read sets known to be clean (here the 552 *Coccidioides* sets that passed the
TE project's species audit, screened at 2M reads); it gives the false-positive reference at the depth the pipeline uses.

## Output

`spikein.rows.tsv`: one row per condition and depth with GmG's call and mixture block. `spikein.detection.tsv`: per
contaminant, fraction and depth, how many conditions GmG's own flag sees and how many a (minimum hashes, minimum
concentration) rule sees, for a grid of both. `spikein.false_positives.json`: what the same rules flag in the controls, the
clean sets and the same-genus mixtures.

## Limits

One run per condition and one background pair per species. The mix takes each source's reads in file order, so it does not
model a contaminant that is spread unevenly across lanes. GmG 0.0.1, index 2026.10.5. Read length varies between read sets
(93 to 150 bp), so the same number of pairs is a different number of bases.

## Read-level confirmation of flagged sets

`confirm_flagged.slurm` checks a set GmG flagged by looking at the reads themselves, with no guess about which genome the contaminant is: it
maps the first 8M pairs to two *Coccidioides* references (`confirm_setup.sh` builds the bowtie2 index), keeps the pairs where neither mate maps,
and classifies those with GmG. `confirm_summary.py` calls a set `confirmed` when GmG names the expected second genus in the unmapped reads
(as the primary genus or as a mixture member) and the unmapped fraction is above that of the negative controls, `not_confirmed` when it does not,
and `no_unmapped` when too few reads are left. `confirm_ids.tsv` lists the sets: the 13 low-level flagged *Coccidioides* sets, 4 strongly flagged
ones as positive controls, and 6 clean or borderline sets as negatives.
