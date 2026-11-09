# Clone groups by solid containment (2026-11-08)

Which runs of a species are the same strain sequenced again, so the MAT counts can be given once per genotype as well as once per run.
Written by `../clone_solid.py`, read by `scripts/build_reads_site.py`. It replaces the older pair tables (`../clone_pilot/`, `../wave2/clone_pairs/`)
for every species that has a table here.

## Files

| File | What |
|---|---|
| `<slug>.clones.tsv` | One row per sketched run: hashes, solid hashes, median k-mer abundance (a coverage proxy for the subsample), nearest run and its ANI, clone id at 0.9999, 0.9995, 0.999 and 0.998 (`c1` the largest group), and the least similar pair inside the run's group at each threshold |
| `summary.tsv` | One row per species group: runs typed and sketched, median abundance, low-coverage runs, between-run ANI median and p99, clones at each threshold, the same-strain check |
| `same_strain_pairs.tsv.gz` | Every pair of runs of one species whose strain names match (after normalising case and punctuation; placeholders such as `missing` skipped), with its ANI |

## Method

- **Sketches:** the worker's sketch of each run's subsample (`genotype_mat_one_accession.sh` with `SKETCH=1`: up to 2 million read pairs, k=21, scaled 500, with abundances). The sketches were written during typing and in the sketch pass of 2026-11-08 (`waves/<wave>/{out,sketch_out}` on the HPCC).
- **Before the distance:** k-mers seen once are dropped, as the sketch store does when it sketches (`../../sketch_store/sketch_one.sh`). Then the sketch is downsampled to scaled 1000.
- **Distance:** solid containment, as in `../../sketch_store/sketch_query.py`, which was validated on 136 *A. fumigatus* strains against SNP calls.
  - A hash is solid in run a if its abundance is at least max(3, 0.25 × a's median). Only solid hashes count in a's denominator. A hash counts as shared if run b has it at all.
  - ANI(a,b) is the mean over both directions of containment^(1/21).
  - So sequencing errors and thin coverage do not count as differences.
- **Groups:** average linkage (scipy) on 1−ANI, cut at each threshold.

**A first run left k-mers seen once in the sketch, and it was wrong.** In these subsamples, error k-mers outnumber the genome's k-mers. That pulled the median abundance to 1 to 5, so 167 of 203 species looked low coverage, and the solid threshold fell to its floor of 3. That run is not used; its output is kept on the HPCC in `results/2026-11-08_clone_solid/` for comparison.

## Calibration: runs filed under the same strain name

There are no split-half sketches left from the October pilot, so the check uses the data itself. It takes 30,045 pairs of runs whose strain names match. Some of these are true resequencing of one strain; others are a name reused for a different line or a lab-derived strain.

| ANI band | same-name pairs |
|---|---|
| ≥ 0.9999 | 1,235 |
| 0.9995–0.9999 | 5,057 |
| **0.999–0.9995** | **881** |
| 0.998–0.999 | 5,533 |
| 0.995–0.998 | 5,977 |
| < 0.995 | 11,362 |

The distribution has two modes, with a valley between 0.999 and 0.9995: copies sit at or above 0.9995, and different lines with one name sit at or below 0.999. This agrees with the sketch store's validation at 30x (same strain median 0.99996; between strains p99 0.999). Hence:

- **Clone threshold 0.9995** (average linkage).
- **Plateau check against 0.999.** Counts per clone are given only if the number of clones moves by at most 10% of the runs between the two thresholds.

## Rules on the site (`build_reads_site.py`)

Counts per clone are given only when all of these hold:

1. **At least 10 sketched runs.**
2. **The plateau holds**, as above.
3. **Coverage is enough:** the species' median k-mer abundance is at least 8. Below that, a thin subsample cannot separate a copy from a close relative. 23 of 202 groups fall below.
4. **The species is not uniform:** the median ANI between its runs is below 0.999.

Rule 4 covers species whose strains differ by fewer variants than the sketch resolves, at roughly one in 10 kb or less. There, a group means "the sketch cannot tell these runs apart", not "these are copies". This applies to 13 groups:

- *Exophiala dermatitidis*: between-run median 0.99990.
- *Trichophyton indotineae*: 0.99984, a clonal epidemic.
- *Lodderomyces parapsilosis*: 0.99957.
- Several *Candidozyma* species.
- Small groups.

Groups are still shown for every species. Only the counts are withheld.

## Limits

- **Biology is not known here.** Sibling progeny of one cross fall below the threshold, and a clonal lineage with a few hundred SNPs can fall above it.
- **The same-strain names are a guide, not truth.**
- **A few species have very many same-name pairs:** *P. kudriavzevii* (4,191), *S. uvarum* (9,311), *C. albicans* (5,135) and the *C. gattii* complex (4,529). These are mostly one study naming many runs alike. They weigh heavily in the table above, but the two modes are also seen in the other species.
- **Diploids and heterozygous species (*C. albicans*):** the k-mers of both haplotypes are in the sketch. Two runs of one strain that lost heterozygosity differently can fall below the threshold.

## Reproduce (HPCC)

```
W=results/sra_reads_typing/waves
A=; for d in $W/wave*/sketch_out $W/wave*/out; do A="$A --sigs $d"; done
sbatch -p short --qos short -c 16 --mem 96gb -t 2:00:00 --wrap \
  "/opt/linux/rocky/8.x/x86_64/pkgs/sourmash_branchwater/0.9.13/bin/python datasets/reads/sra_survey/pilot/clone_solid.py \
   --pilot datasets/reads/sra_survey/pilot --out results/2026-11-08_clone_solid_v2 --threads 16 $A"
```

This run was job 31874107, at commit cd5648d, and took about 25 minutes on 16 cores. It sketched 21,618 of the 22,138 typed runs across 202 groups. A run without a sketch is shown without a clone id.
