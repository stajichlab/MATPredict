# Zygomycete (Mucoromycota) detection validation — 23 genomes with prior ground truth

Run 2026-09-20 on branch `detect-targeting` (commits `efdfac8..11f650a`, 377 tests passing).
All numbers below are measured from real runs, not estimated.

## What was run

`matpredict detect --genome <g>.fna --proteins <g>.faa --phylum Mucoromycota
--evidence-diagnostics <...> --emit-cds-fasta`, once per genome, 23 genomes.

Inputs: `testset/Zygo/query/` (813 symlinks to annotated ZyGoLife GenBank genomes).
The 23 are those with prior absence/presence ground truth in
`testset/Zygo/plus.abspres.csv` (16) and `minus.abspres.csv` (7), produced by an
earlier cblaster analysis. GenBank was converted to nucleotide FASTA plus a
proteome whose deflines carry `contig:start-end:strand` (required by the
diamond fast path). Genomes are heavily fragmented drafts: median 1851 contigs,
range 590–3349, 20/23 over 1000 contigs.

## Performance

| | Value |
|---|---|
| Genomes completed | 23/23, 0 errors |
| Mean wall time | 41.8 s |
| Max wall time | 55.1 s |
| Total | 835 s (~14 min) |
| Reference proteins queried | **19** (Mucoromycota only) vs 181 unrestricted |

For contrast, the same pipeline before this branch's routing fix was killed after
**>24 minutes on one small yeast genome without completing**, because routing fell
through to all 19 families across 3 phyla and polished spurious cross-phylum
candidates serially.

## Locus-level recall: 23/23

Every genome: exactly one locus detected, on the ground-truth scaffold.

Per-gene recall, denominator = organisms whose ground truth marks the gene present:

| Gene | Recall |
|---|---|
| `tptA` (flanking_conserved) | 23/23 (100%) |
| `rnhA` (flanking_conserved) | 23/23 (100%) |
| `sexP` (core_MAT) | 16/16 (100%) |
| `sexM` (core_MAT) | 7/7 (100%) |

**Caveat on this number.** The ground truth was produced by cblaster, which finds
gene *clusters* by construction, so a genome whose locus was split across contigs
could never have entered the CSV. These 23 are therefore the cases where the locus
is intact on one scaffold — confirmed: all 23 ground-truth loci have all three genes
co-located, spanning 6795–13089 bp. Genuinely fragmented loci would hide among the
790 genomes without ground truth. 100% here is real but is not evidence about the
fragmented case.

## Finding 1 — idiomorph calling fails in every genome, and the fix is measured

**`idiomorph=undetermined` in 23/23.** In every genome BOTH `sexM` and `sexP` were
reported, at ~100% overlapping coordinates (92–100% overlap of the shorter feature).
`sexM` and `sexP` are mutually exclusive alternative idiomorph genes sharing an HMG
domain, so both curated references hit the same locus gene.

**Taking the higher-identity of the two as the idiomorph is correct 23/23 (100%).**

| Truth | n | `sexP` identity | `sexM` identity |
|---|---|---|---|
| plus | 16 | 45.3–100.0 | 27.3–35.8 |
| minus | 7 | 25.9–32.9 | 30.6–43.5 |

Margins are wide for plus strains and **thin for minus (2.3–13.0 points)**. The
narrowest is `Cunninghamella_bertholletiae_NRRL_1376`: `sexP` 28.3 vs `sexM` 30.6.

This also explains an observed confidence pattern: **all 7 minus genomes came out
`confidence: medium`, 14 of 16 plus genomes `high`** — a direct consequence of the
lower minus identities.

Root cause is reference coverage, not algorithm: the curated DB holds 3 Minus and 3
Plus Mucoromycota records (*Mucor circinelloides*, *Phycomyces blakesleeanus*,
*Rhizopus arrhizus*), and none is close to *Cunninghamella*, *Absidia*,
*Chaetocladium* or *Actinomucor* — 5 of the 7 genera in this test set.

## Finding 2 — the evidence floor works, and identity would sharpen it

Measured on `Absidia_cuneospora_RSA_623_Plus`: 60 candidate (cluster, family) pairs
evaluated, **38 rejected** by the new floor (`min_hits=2`, `require_core_role=True`),
22 admitted.

| | Best identity |
|---|---|
| True locus | **75.9%** |
| 20 spurious admissions (10 scaffolds) | 30.6–50.0% |

All spurious admissions were core_MAT-only 2-gene clusters — HMG-domain proteins are
common in fungal genomes, so `sexM`/`sexP` references hit unrelated HMG genes. A
`min_identity` of ~55–60% would remove all 20 and keep the true locus, a 25.9-point
margin. Downstream tiering already filtered them out of the report (1 locus reported),
but each still paid for two polish subprocesses.

`min_identity` remains `None`. It should not be set from one genome, and the minus
identities above (25.9–43.5%) show a global cutoff near 55% would destroy minus-strain
detection entirely. Any identity floor must be role- or gene-aware, not global.

## Finding 3 — evidence diagnostics rows are duplicated

Every `(cluster, family)` pair appears **twice** in `evidence_diagnostics.jsonl`
(60 rows = 30 unique pairs), inflating the calibration dataset 2×. Not investigated.

## Validated loci available for curation

All 23 detections carry real `CDS` features with `translation=` attributes and a
companion FASTA (60-column wrapped, contig names matching the GFF3 seqids). They are
candidates for promotion into the curated reference DB, which would directly address
Finding 1's root cause — especially the 5 genera with no curated representative.

Outputs are in `$SCRATCH/zygo/runs/<organism>/` and are NOT committed (per-genome
GFF3 + FASTA + report + diagnostics).
