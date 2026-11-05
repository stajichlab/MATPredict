# #148: the Hyphopichia burtonii MTLA1 record on the v0.8.0 database (2026-11-04)

**Question (#148):** why does adding `717740_makgeolli_MTL_A` (*H. burtonii*, GCA_003856795.1; candidate on the HPCC, never in `db/`) lower confidence in some Serinales genomes? **Setup:** the v0.8.0 database with and without the record, the same v0.8.1 code, 389 Serinales genomes (the 12 named in #148 plus a random 377 with an MTL call; `setup_and_run.sh`, `compare.py`).

## Result
354 of 389 genomes unchanged; 9 confidence up, 11 down, 15 other (gene scores or models change, tier the same).

**The five A-locus drops are all the 39-bit floor cliff of #149** (`a_locus_genes.py`):
| genome | gene that crossed the floor | bits before -> after |
|---|---|---|
| GCA_030556945.1 | MTLA2 (36% id, 39% cov), while MTLA1 improved to a 40/100% polished_agree model | 34 -> 45 |
| GCA_030563265.1 | MTLA2 | 38 -> 45 |
| GCA_030567015.1 | MTLA1 | 39 -> 43 |
| GCA_030674215.1 | MTLA2 | 31 -> 42 |
| GCA_030444895.1 | **OBP1, a flanking gene** (16% coverage); the MTLA genes are unchanged | 38 -> 42 |
So the cliff also bites on a flanking gene, not only on MTLA1/MTLA2.

**The six alpha-locus drops are new weak MTLA1 or MTLA2 hits** (`alpha_locus_hits.py`): 16-28% coverage, 36-39% identity, 26-28 bits, lying 28-30 kb from the locus's MTLalpha genes at the ends of 31-40 kb spans (GCA_030557965.1, GCA_030566005.1, GCA_030571075.1, GCA_030578995.1, GCA_030580355.1, GCF_019202705.1). They are other-allele cross-hits, not a second locus. The cross-hit guard (`allele_absent_genes_to_ignore`: identity under 50% and at least 10 points under the called allele's best core gene) does not ignore them because the called alpha genes are themselves only 35-45% identical, so the gap is under 10.

## Reading
The record is not wrong: it gives a better MTLA1 model in some genomes (MTLA1 51/25% -> 40/100%). The drops come from two threshold rules acting on weak fragments: the 39-bit noise floor (#149) and the 10-point cross-hit gap. #149's variant 1b (a gene is noise when every hit is under 39 bits or under 30% coverage with a score under 50 bits) removed the #147 drops without any down in 2,368 genomes; it reverses part of ruling F6 and needs a curator ruling. Not yet tested with this record.

## Variant 1b of #149 does not clear the drops (2026-11-05; `run_1b.sh`, `compare_1b.py`)
Same 389 genomes, same code, with the experimental noise rule of #149 variant 1b (branch `floor-1b-experiment`, an environment switch, not for merge): a hit under 30% coverage and under 50 bits counts as noise.
| comparison | same | up | down |
|---|---|---|---|
| record added, floor as is | 354 | 9 | 11 |
| 1b alone (no record) | 386 | 3 | 0 |
| record added under 1b | 356 | 8 | 9 |
1b clears two of the eleven drops (GCA_030444895.1, the OBP1 flank gene at 16% coverage, and GCA_030567015.1, an MTLA1 fragment at 27%). Nine remain: the six alpha-locus cross-hits (rest on the 10-point cross-hit gap, not the floor) and three A loci (GCA_030556945.1, GCA_030563265.1, GCA_030674215.1) whose MTLA2 hit covers 39-46% of its reference, above the 30% cut. So variant 1b alone is safe (3 up, 0 down on 389 genomes) but is not the fix for this record. **Recommendation: hold `717740_makgeolli_MTL_A` and move #148 to v0.9.0 with #149**; the source genome (GCA_003856795.1) is already called `A`, medium, without it.

## Outcome (2026-11-05)
The curator decided to add `717740_makgeolli_MTL_A` now, although the recommendation above was to hold it for #149 (PR #273, merged). The 11 downward changes on the 389-genome sample stand and stay open under #149.
