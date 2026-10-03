# Split-locus rule measured (code 4174440 vs 8d80bed)

Rule (curator's ruling 2026-09-27; src/MATPredict/detect/split_locus.py): for a
family with no reported call, one MODELLED core gene >= 95% within 500 bp of a
contig end, >= 2 roster flanks at >= core identity - 10 (best genome-wide hit
of each, outside any reported locus, >= 1 on another contig), and the core hit
the genome's best for its gene -> `partial_locus`, `low`, `split_locus` block.
Idiomorph from the classifier on the modelled core protein.

| panel | genomes | called before -> after | gained | lost | changed | split calls |
|---|---|---|---|---|---|---|
| Mucoromycota | 293 | 235 -> 247 | 12 | 0 | 0 | 12 |
| Cryptococcus | 243 | 237 -> 237 | 0 | 0 | 0 | 0 |
| Serinales (random 200) | 200 | 191 -> 191 | 0 | 0 | 0 | 0 |
| Dothideomycetes pilot | 99 | 89 -> 89 | 0 | 0 | 0 | 0 |

Zygo 23: 23/23 locus and idiomorph on scaffolds and on contigs.

The 12 new calls: 11 of the 12 GL genomes (all Plus, sexP 98.4%, edge 1-199 bp,
classifier on the modelled sexP, margin ~306 bits) plus R. microsporus 56028
GCA_039881115.1 (sexP 98.82%, edge 78, margin 248.5). The 12th GL genome,
GCA_011764265.1, has no sexP hit in the diagnosis either (best_hits.tsv) and
stays uncalled. Flanks: tptA 100%, btbA 96-98.6%, rnhA 92.5-93.0% on 2-3 other
contigs. The "flank>5below" flag in compare_output.txt is rnhA at ~93% vs sexP
98.4%, inside the 10-point rule; no call looks like a paralog (no genome had
another call; every core is the genome's best sexP hit).

Note: R. microsporus 56028 has tptA 100% to the R. arrhizus reference; its
species label may be wrong (not checked).

Not recovered, and not split cases: R. microsporus MLY36 (best cluster tptA/sexM/
algA on one contig, 1 modelled gene) and R. stolonifer PG92-21 (sexM+rnhA, 0
modelled). Neither has a >= 95% modelled core gene at a contig end.

Runtime: unchanged within noise (Cryptococcus 15 vs 21 min, Serinales 22 vs 22,
Dothideomycetes 35 vs 35 min per 16-way job).
