# C. auris no-calls, flank-carried calls, and C. albicans a/alpha collapse

2026-09-26. Three open questions from `docs/HANDOFF-2026-09-26.md`.

## 1. The six uncalled C. auris genomes: explained, and now called

Genomes: GCA_035771965.1 (CA220), GCA_048414545.1 (CBS15605),
GCA_053223135.1 / 053223155.1 / 057661105.1 (LACPHL series),
GCA_059629915.1 (MRSN134489). On the pre-flank roster (`8dec9e1`) their best
locus was withheld, labelled alpha, 0 modelled genes.

* **Assembly quality is not the cause.** All ~12.2 Mb with 0% N; three are
  near chromosome level (7-21 contigs, N50 0.84-2.6 Mb), three fragmented but
  ordinary (229-579 contigs, N50 33-98 kb).
* **They are alpha-idiomorph genomes.** On the current code (`882aa01`, PAP1/
  OBP1/PIK1 flanks, miniprot-polished) all six are called **alpha, high,
  `mat_locus`, 4 modelled genes (MTLalpha1 + the three flanks)**. C. auris
  alpha carries ONE core gene (alpha2 is absent in the clade), so before the
  flanks an alpha genome could not reach the two-modelled-gene bar. The
  "MTLA1" in their old withheld cluster was a cross-matching fragment.
* Not explained: why alpha1 itself got 0 models on the old run despite a
  same-species reference. The flanks make it moot; noted in case it recurs.

Evidence: scratch runs of the six genomes on `run-882aa01` (reports reproduced
by the Serinales-wide scan `results/2026-09-26_serinales_all_882aa01/`).

## 2. Flank-carried calls: a recommendation

A flank-carried call has no modelled MTL core gene; it cleared the bar on the
flank models alone. In the 910 Serinales genomes scanned on the flank roster,
**7 of 1,096 calls (0.6%)** are flank-carried, all in genera with no curated
reference (Australozyma x4, Kurtzmaniella, Limtongozyma, Candida argentea).

| genome | species | call | core hit(s), all unmodelled | core inside flank block |
|---|---|---|---|---|
| GCA_030581555.1 | Candida argentea | A | MTLA2 33.8%, MTLA1 26.7% | yes |
| GCA_019775655.1 | Candida anglica | A | MTLA2 43.6%, MTLA1 27.9% | yes |
| GCA_030566735.1 | Australozyma touchengensis | A | MTLA2 31.7% | yes |
| GCA_030563905.1 | Australozyma succicola | alpha | MTLalpha1 26.6% | yes |
| GCA_030563825.1 | Australozyma nongkhaiensis | A | MTLA1 33.3% | yes |
| GCA_030569455.1 | Australozyma saccharicola | alpha | MTLalpha1 28.3% | yes |
| GCA_030558265.1 | Candida pseudocylindracea | A | MTLA2 27.7% | **no** |

All were reported `medium`, as `mat_locus` or `partial_locus`.

**Recommendation (not implemented; curator decision):** when a call has no
modelled core gene,
1. core hit inside the PAP-OBP-PIK span (+-3 kb): keep, cap at `low`, class
   `partial_locus`, and flag `idiomorph_unmodelled: true` -- a real MTL position
   in an uncurated lineage whose idiomorph rests on one weak hit;
2. core hit outside the flank span: withhold, like any bar failure.
On this data that downgrades 6 calls and withholds 1.

## 3. C. albicans assemblies hide a/alpha heterozygosity

Question: 134 C. albicans genomes on the flank roster are called a+alpha 77,
alpha only 35, a only 17, none 5 -- far fewer heterozygotes than wild isolates
show. Test: WGS reads (ENA, library_strategy WGS), R1 randomly subsampled to
~3 M reads across the WHOLE file, mapped with minimap2 (-ax sr, MAPQ >= 20) to
both SC5314 idiomorphs (AF167162.1 a, AF167163.1 alpha) and a 20-kb single-copy
control (NC_032089.1:1,000,001-1,020,000). Depth ratios near 0.5/0.5 of the
control mean heterozygous; ~1/0 homozygous. `results/2026-09-26_calbicans_mtl_reads/`
(`depth2.tsv`, `map2.slurm`).

| assembly genotype | reads say | isolates |
|---|---|---:|
| alpha only | **a/alpha** (a 0.36-0.53, alpha 0.42-0.61 of control) | **7 / 7** (A48, A67, A92, A123, A203, CHN1, ATCC 36802) |
| a only | a only (a 0.68-0.97, alpha 0.00) | 5 / 6 |
| a only | **a/alpha** (P37037) | 1 / 6 |
| a + alpha | a/alpha | 3 / 3 |

* Every alpha-only assembly tested is a collapsed heterozygote; the collapse
  kept the alpha haplotype. Six of the seven share one assembly series
  (Broad), so this may be a pipeline effect, not general.
* a-only assemblies are mostly genuine MTLa homozygotes.
* **Detection is correct about the assembly; the assembly is wrong about the
  genome.** For diploid heterozygous species, genotype needs reads.

Method trap, recorded: the first round took the FIRST 2 M reads; some SRA
submissions are ordered by genome position, so that sampled only the start of
chr1 (control 5-8x too deep, MTL at zero). Always subsample randomly across the
whole file.

**Possible feature (not built):** an optional read-based zygosity check that
maps reads to both idiomorphs of the detected family and reports depth ratios.
