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
