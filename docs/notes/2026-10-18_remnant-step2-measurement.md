# Step 2 of the remnant status: measured effect of removing remnant proteins from the references (2026-10-18)

Change (branch `remnant-status-step2`, on top of `remnant-status-step1`): `gff_export.write_proteins_fasta` omits genes with a `status`, and the three
remnant/pseudogene proteins that were `present: true` leave `proteins.faa`: `746128_a1163` MAT1-2-1 (250 aa), `230073_uamh-11059` MAT1-1-1 (224 aa),
`107463_unknown-1` MAT1-2-1 (104 aa). `reference_fasta` and `_curated_protein_lengths` both read `proteins.faa`, so these three stop being search references
and stop being the denominator of the coverage rule.

## Method
`detect --taxid <the genome's own taxid>`, `MATPREDICT_HTML=0`, own `MATPREDICT_DB_ROOT` and `PYTHONPATH` per variant: A = origin/main (v0.6.3), B = step 2.
Default settings (`--partial-core-policy off`), and, for a subset, `--partial-core-policy coverage_profile`. No records withheld. A report differs when its
`Ascomycota:MAT` loci differ (contig, start, end, idiomorph, confidence, class, genes found); db path, db hash and timings are ignored. HPCC `rem-*` arrays,
at most 15 tasks. Genome ids and taxids: `pcv` job lists and `assembly_summary_GCA/GCF.txt` column 6.

## Result, default policy (627 genomes, A vs B)
| set | genomes | loci changed | what changed |
|---|---|---|---|
| *A. fumigatus*: all 180 genomes with a `homothallic_candidate` MAT1-1 locus in v0.6.0, plus 30 other *A. fumigatus* | 210 | 180 | `homothallic_candidate` -> `mat_locus`; idiomorph MAT1-1 and confidence high unchanged; genes found unchanged |
| Ophiostomatales (all in the campaign) | 126 | 1 | *Leptographium procerum* GCA_000806385.1: the MAT1-2 locus start moves 836 -> 1436 (same call, genes, end) |
| Teloschistales (all in the campaign, incl. Xanthoria) | 88 | 2 | *Flavoplaca citrina* and *F.* cf. *flavocitrina*: a 102-bp, 14%-coverage, 27.7-bit MAT1-2-1 "hit" 4 kb past SLA2, built on the 107463 protein, disappears; MAT1-1 / medium / `mat_locus` unchanged |
| Botrytis | 34 | 0 | all identical |
| BioSample truth set (203 with a stated MAT1-1/MAT1-2, plus Botrytis T4 and B05.10 counted in the Botrytis row) | 203 | 0 | no model in any of them was built on a remnant protein |

No call (idiomorph, confidence) changed anywhere. 12 more loci changed only the reference record of a model (see below).

## Was any model built from a remnant protein as its best reference?
Yes, often, in the clades the remnants come from.
- a1163 MAT1-2-1 (250 aa) was the reference of the MAT1-2-1 model in 186 of the 420 clade loci (all MAT1-1 *A. fumigatus* loci but one). Its model covered 100% (median; min 52.8%) of the
  reference: "full length". After step 2 the same model is built on Af293's 322-aa protein (`746128_af293`, 181 of 181) and covers a median 37.6% (14.3-45.3%).
- uamh-11059 MAT1-1-1 (224 aa) was the reference of the MAT1-1-1 model in 10 Ophiostomatales loci (MAT1-2 loci that carry the truncated MAT1-1-1; *O. montium*
  97.8% identity at 100% coverage). After step 2 they are built on `230073_uamh-1363` MAT1-1-1 (711 aa) and cover 10.8-31.8%. Same loci, same calls.
- 107463 MAT1-2-1 (104 aa): the reference of the weak 14% MAT1-2-1 hits in the two Flavoplaca loci only.
- Truth set: none. Botrytis: none.

## The `homothallic_candidate` label of *A. fumigatus*
It clears in all 180 loci, and the evidence says that is the right result. The label came from the relaxed rule (`HOMOTHALLIC_MIN_MODEL_COVERAGE` 0.5): the
3' MAT1-2-1 remnant of the MAT1-1 idiomorph, modelled on a1163's own 250-aa remnant protein, read 100% "full length", so the locus looked like it carried a full
alpha-box gene and a full HMG gene. Against a real MAT1-2-1 (Af293) the same model covers 14-45%. This agrees with `analysis/2026-10-07_afum-reads-and-assembly.md`
(a degenerate remnant 117 bp from MAT1-1-1 in all MAT1-1 assemblies) and with *A. fumigatus* being heterothallic. The loci are not reclassified as anything new:
they become plain `mat_locus` MAT1-1. This is a consequence the curator should confirm (question 1).

## The edge-kept contig (curator request)
Genome GCA_020498755.1 (*A. fumigatus* E-1-75s-2), contig JAIBPP010000651.1 (14.5 kb), locus 8,366-14,446. Main, default policy: MAT1-1, high, `homothallic_candidate`.

| gene | span | identity | coverage | reference record (main) | reference record (step 2) |
|---|---|---|---|---|---|
| COX13 | 8366-9103 | 100 | 100 | a1163 | a1163 |
| APN2 | 9464-11649 | 99.8 | 100 | a1163 | a1163 |
| MAT1-2-1 | 12567-13601 (main) / 12567-12929 (step 2) | 100 / 93.4 | 100 / **37.6** | **a1163 (the remnant protein)** | af293 |
| MAT1-1-1 | 13718-14437 (ends 30 bp from the contig end) | 100 | 60.6 | a1163 | a1163 |

Profile scores (`scripts/score_core_profiles.py`, `MAT_core`): alpha-on-alpha 301.5, alpha-on-HMG 12.9; HMG-on-HMG 142.6 (main model) / 113.2 (step 2 model),
HMG-on-alpha 2.3 / 3.2; both call MAT1-1 (margin 158.9 / 188.3). So the "complete MAT1-2-1, 100% identical to Af293" in the flip review is the a1163 remnant read against its own
protein: it is 37.6% of a real MAT1-2-1, not a full-length gene, and the contig is an ordinary MAT1-1 locus with its usual remnant (not a chimeric contig); the truncated MAT1-1-1 is an
edge truncation of the 720-bp-modelled gene. Step 2: idiomorph MAT1-1 unchanged, confidence high unchanged, class `homothallic_candidate` -> `mat_locus`. With
`--partial-core-policy coverage_profile` main fires the rule here (`kept_edge_artefact`: MAT1-1-1 61% but 30 bp from the edge, "other" MAT1-2-1 at 100%); step 2 gives no vote note
because the other gene is no longer full length. Same label.

## `coverage_profile` policy (298 genomes: 50 *A. fumigatus* incl. the edge contig, all 214 Ophiostomatales/Teloschistales, 34 Botrytis)
Identical findings: 43 *A. fumigatus* class changes, the same Leptographium and Flavoplaca changes, Botrytis 0 changed; the 25 flips and 1 `undetermined` vote notes are the same in A and B.
The step has no effect on the partial-coverage rule's decisions in these genomes.

## Held-out at species and genus radius
`run_asco_holdout.py` with a manifest extended by species-radius counts (`holdout_manifest_s2.tsv` style; only 3 of the 10 relevant records have an assembly, so the
records of uamh-11059 and 107463 cannot be held out on their own genome). Run: a1163 and Af293 on their own genomes, plus four record-less genomes as proxies
(*Ophiostoma montium*, *O. novo-ulmi*, *Xanthoria mediterranea*, *X. calcicola*) with the species or genus records withheld. 12 runs, A vs B: all 12 identical
(same locus, idiomorph, confidence, class, genes, profile call; a1163 and Af293 found at the truth span). By construction a record's own protein is withheld in its own
radius, so this check cannot show an effect of step 2 for that record; the full runs above are the informative measurement.

## Clades that lose their only reference of an idiomorph
None. Ophiostomatales keep `uamh-1363` MAT1-1-1 (711 aa, an unusually long protein, which is why remnant models now read 10-32%) and `uamh-11059` MAT1-2-1;
Teloschistales keep MAT1-2-1 in 115814, 1301501_voucher-267 and the two Xanthoria sp. records (351-384 aa); Eurotiales keep a1163's MAT1-1-1 (368 aa) and Af293's MAT1-2-1 (322 aa).

## Cases it makes worse
None found. Two caveats: (1) an Ophiostomatales MAT1-1-1 now has a 711-aa denominator, so a short-but-real MAT1-1-1 would read as partial under `coverage_profile`;
none of the 126 genomes changed call, but the 711-aa record is worth a look; (2) the *Leptographium* locus start shifts by 600 bp (a weak remnant-protein hit no longer extends the cluster).

## Not verified
Species/genus held-out for uamh-11059 and 107463 on their own genomes (no assemblies); the truth set was not scored against stated types (no call changed); the 417 non-sampled *A. fumigatus* MAT1-2 genomes were not run;
`coverage_profile` was not run on the truth set.
