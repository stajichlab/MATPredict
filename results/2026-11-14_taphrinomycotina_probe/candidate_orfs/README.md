# Candidate short ORFs (for later search against nr, other genomes and proteomes)

Created 2026-11-16. `candidate_orfs.faa` holds protein sequences from the Taphrinomycotina MAT probe. None is a confirmed gene.
Header format: `set|region or transcript|frame|length|note`. Coordinates are 1-based, inclusive, on the forward strand of the contig.

| Set | Source | What it is |
|---|---|---|
| `Tdef_A2_minusgene` | StringTie on PYCC 5710 24H1 reads mapped to *T. deformans* A2, RHGF01000006.1:52,192-52,790, minus strand, 5 exons | Spliced, highly expressed gene that overlaps the matPi region by 45 bp. Longest ORF 89 aa. No HMG-box or homeodomain hit at E <= 1 in six frames. Function unknown; the first search to run is blastp/DIAMOND against nr. |
| `Tdef_A2_plus_HDregion` | RHGF01000006.1:51,694-52,236, plus strand | Translation behind the curated matPi record (5011), 181 aa. A tblastn region, not a gene model. Low, flat RNA-seq depth, no intron inside. |
| `Tdef_A2_plus_HDregion_ORF` | Stop-free stretches >= 40 aa in the same region (no intron model) | For searches of the homeodomain-like region. |
| `Plact_12-1054_HDregion_ORF` | *P. lactucae-debilis* 12-1054 NW_024468021.1:38,146-38,599 | Stop-free stretches >= 40 aa (no intron model) in the homeodomain-like region; the annotation has no gene model there. |

Also here: `tx_6frame.faa` (six-frame translation of the two StringTie transcripts) for the same region.
Transcript nucleotides: `../rnaseq/stringtie_pycc/pycc.tx.fa.zst` (all PYCC transcripts; the two at this locus are `ST_PYCC.1275.1` and `ST_PYCC.1275.2`).
Limits: ORFs come from stop-codon spans, not from gene models, except the `minusgene` set (spliced transcripts). Hits to nr are not evidence of MAT function.
Files not committed (size): BAMs, `.sra` files and downloaded reads under `rnaseq/`. Regenerate with `scripts/rnaseq_locus_batch.slurm`.

## nr search (2026-11-16, job 32114380)
Local DIAMOND 2.1.21 blastp `--more-sensitive --evalue 1e-3 --max-target-seqs 10` against the NCBI nr cluster db (build 2026-01-28). Table: `results/2026-11-16_decisions_followup/nr_diamond/nr_hits.tsv` (41 rows, 5 of the 17 queries with hits). A hit says nothing about MAT function.
- `Tdef_A2_minusgene` (89 aa and 52 aa ORFs): best hits are 33-42% identical to hypothetical proteins and to "protein Asterix" / PAT-complex subunit Pat10 (*S. osmophilus*, 34.5%). The gene looks like an Asterix (Get/PAT complex) homolog, not a transcription factor. This agrees with the missing HMG/homeodomain hit.
- `Tdef_A2_plus_HDregion` (181 aa) and its ORF: hits to *Pneumocystis* hypothetical proteins (29-44%), *Neolecta irregularis* "Homeobox protein knotted-1-like 7" (33.2%) and *P. lactucae-debilis* "homeobox KN domain-containing protein, partial" (54.0%). The region is homeobox-like (KN/TALE type).
- `Plact_12-1054_HDregion_ORF` (frame -2, 151 aa): hits to "homeobox KN domain-containing protein" of *Aspergillus recurvatus* (70.2%), *Neolecta* KNOX-like 7 (51.5%) and *Pneumocystis* proteins. The other 12 ORFs had no hit at E <= 1e-3.
Limit: a KN-type homeobox gene is common in fungi; the hits do not show a mating-type role.
