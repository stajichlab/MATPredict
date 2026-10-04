# R. toruloides HD: generic-HD call (2026-09-26) versus redHD (v0.6.0)

Status: open (curator ruling pending). Recommendation below.

## Question
In 26 Sporidiobolales genomes, the 2026-09-26 run called a generic
`Basidiomycota:HD` locus, and v0.6.0 calls `Basidiomycota:redHD` at another
place. Which one is the mating-type HD locus?

## Data
- Old: `results/2026-09-26_basidiomycota_full/` (run-ad1f865). Sporidiobolales
  had no curated record then: every genome routed `phylum_fallback` and was
  searched with the Agaricomycete-centred HD family (HD1, HD2, MIP1, beta_fg).
- New: `results/2026-10-03_basidiomycota_v060/` (run-7c7ed99, = v0.6.0).
  Every genome routes `lineage` to the Rhodotorula records (5 redHD + 6 redPR
  from the group's preprint, doi 10.1101/2025.09.11.675505; tier 2, signed off
  2026-09-27).
- 253 Sporidiobolales genomes in both runs. Table:
  `results/2026-10-04_rtoruloides_hd/per_genome.tsv` (`compare_hd.py`).

## Results

### 1. Known positions (curated records; not independent for redHD)
| Assembly | Known HD locus | Old generic HD | New redHD |
|---|---|---|---|
| GCA_000988875.2 R. toruloides NBRC 0880 | LCTV02000005.1:1,270,620-1,273,396 | LCTV02000013.1:669,326-734,979 (PR contig) | LCTV02000005.1:1,256,423-1,276,797 |
| GCA_921037615.3 R. toruloides CBS 14 | CAKLCE030000016.1:73,460-76,222 | CAKLCE030000013.1:273,532-301,593 | CAKLCE030000016.1:73,463-76,219 |
| GCA_920103745.3 R. glutinis CBS 20 | CAKKSX030000026.1:389,350-392,484 | none | CAKKSX030000026.1:389,353-412,378 |
| GCA_024748845.1 R. mucilaginosa JY1105 | JANBVD010000009.1:486,366-489,433 | none | JANBVD010000009.1:486,369-489,430 |
| GCA_002917965.1 R. kratochvilovae LS11 | PQDI01000056.1:305,841-308,574 | PQDI01000056.1:306,413-308,614 | PQDI01000056.1:305,710-316,513 |

The old call misses the known locus in both R. toruloides assemblies. The new
call matches all five, but these are its own training records.

### 2. Sequence check at the old locus (independent)
tblastn of the curated R. toruloides HD1 and HD2 proteins (CBS 14 and
NBRC 0880) against the NBRC 0880 genome:
- every hit is at LCTV02000005.1:1,270,623-1,273,393 (the redHD locus);
- the CBS 14 alleles hit there at 69% (HD1) and 91% (HD2) identity, as expected
  for different MAT alleles;
- **no hit at all in the old locus** LCTV02000013.1:669-735 kb, even at
  e-value 10.

The old call's own evidence in NBRC 0880: HD1 27.6% identity with a gene
model spanning 33.5 kb; HD2 32% identity, 216 bp, unpolished; MIP1 50.6%.

### 3. All 253 Sporidiobolales genomes
| | Old generic HD | New redHD |
|---|---|---|
| Genomes called | 41 | 247 |
| Routing | phylum_fallback (253) | lineage (253) |
| HD1 identity to its reference, median (range) | 34.9% (25.9-50.0), n=28 | 65.5% (32.3-100), n=247 |
| HD2 identity, median (range) | 28.8% (25.6-46.5) | 82.9% (34.8-100) |
| HD1 model | polished_single 19, unpolished 9, absent 13; span median 192 bp, max 33,533 bp | polished_agree 70, polished_disagree 172, polished_single 5 |
| On the same contig as the P/R call | 8/41 | 0/241 |

Where both exist (41 genomes): same locus 15 (all 11 R. kratochvilovae, plus
4 others), other contig 24, same contig but not overlapping 2.

R. toruloides (35 genomes):
- Old: 24 called; genes `HD2|beta_fg` (13) or `HD1|HD2|MIP1` (11); HD1
  identity exactly 27.6% or 40.4% in every genome. 0 of 24 at the redHD locus.
- New: 34 called. HD1 identity to the closer reference: 100% in 25 genomes
  (they carry the NBRC 0880 B39 or CBS 14 B41 allele; many BFD R. toruloides
  are related lab strains, so these are not independent), 64.5-88.2% in 9
  genomes (other alleles).
- A MAT HD gene is multiallelic, so identity should vary between strains, as
  it does at the redHD locus. The fixed 27.6% / 40.4% at the old locus is what
  one conserved gene, scored against distant Agaricomycete references, gives.

### 4. Linkage
In all four curated assemblies the HD and P/R loci are on different contigs,
and in the 2026-09-27 run 0 of 221 Sporidiobolales genomes had redHD and redPR
on one contig. The old note "generic HD on the PR contig, about 200 kb from
PR" held for 8 of 41 old calls only; it is not evidence of linkage.

## Recommendation
Keep `redHD` (v0.6.0 behaviour). The old generic-HD call in R. toruloides is a
weak homology hit to Agaricomycete HD references at a locus with no
R. toruloides HD sequence; it is not the MAT HD locus. The old call agrees with
redHD only in R. kratochvilovae.

## Limits
- The known positions come from the 5 records redHD was built from, so check 1
  is not independent; checks 2 and 3 are.
- Allele diversity is shown for HD1 only, and the R. toruloides set includes
  related lab strains.
- What the old locus encodes was not identified (the HD-like hits are 26-50%
  identity to Agaricomycete HD1/HD2).

## Decision
Pending curator.
