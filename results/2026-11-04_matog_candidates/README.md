# Candidate MATog records for divergent Ogataea (2026-11-04)

Curator decision (2026-11-04): draft candidate `MATog` records, tier 3, `needs_review`, for *O. pignaliae*, *O. wickerhamii*, *O. dorogensis*, *O. siamensis* and *Candida boidinii*, to test whether more references close the recall gap of #228 (divergent *Ogataea* score 67-142 against 269-409 for the polymorpha group). **Nothing is written to `db/`; each record needs curator sign-off.**

## What was built (six records for five species)
`records/*.yaml` (+ `.proteins.faa`), made by `make_records.py` from `models.py`. Genes are the miniprot projection of the DL-1 and CBS 4732 MATog proteins (`run_miniprot_relaxed.sh` and the default run in `data/`), kept only if the projection has no frameshift and no internal stop. No RNA-seq, no annotation: tier 3 throughout; every model is partial.

| record | assembly | genes recorded |
|---|---|---|
| `51923_nrrl-y-17664_MATog_MATa-MATalpha` | *O. pignaliae* GCA_030570775.1 | sla2 (253 aa), MATa1 (58 aa), MATalpha2 (53 aa): fragments on 2 contigs |
| `1249076_nrrl-yb-4943_MATog_MATa-MATalpha` | *O. wickerhamii* GCA_030559035.1 | sla2, MATa1 (83), MATalpha1 (150), MATalpha2 (54), dic1 |
| `169506_nrrl-y-27599_MATog_MATa-MATalpha` | *O. dorogensis* GCA_030568015.1 | sla2, MATa1 (185), dic1, MATalpha2 (161), MATalpha1 (161) |
| `246385_nrrl-y-27963_MATog_MATa-MATalpha` | *O. siamensis* GCA_030574315.1 | sla2, MATa1 (171), MATalpha1 (169), MATalpha2 (148), dic1 |
| `5477_ndk27a1_MATog_MATalpha` | *C. boidinii* NDK27A1 GCA_002008015.1 | sla2, MATalpha1 (151), MATalpha2 (41), dic1 |
| `5477_uniss-cb60_MATog_MATa` | *C. boidinii* UNISS-Cb60 GCA_002008025.1 | sla2, MATa1 (59), dic1 |

Deviations from the proposal, and why:
- **Six records, not five.** Each *C. boidinii* genome carries one idiomorph (some genomes call MATa, some MATalpha; my earlier "all call only MATalpha" was wrong), so there is one record per strain, not one hand-merged record.
- **The four *Ogataea* records have two segments** (the MATa and MATalpha cassette ends sit on different contigs, as in the assemblies).
- **Several are small fragments**: all four *Ogataea* assemblies hold the cassette ends on contigs of 2-20 kb. The low score of *O. pignaliae* (67-77) is partly fragment length, not only divergence; its MATa1 and MATalpha2 models are 58 and 53 aa. *O. siamensis* MATalpha2 uses the CBS 4732 projection because the DL-1 one has a frameshift.

All six pass `validate_metadata`; `curate-db validate` (2026-11-04): `taxonomy_current` and `accession_resolved` true for all six; sequence match not applicable (no protein accessions).

## Leave-one-out check (`loo_setup.sh`, `loo_run.sh`, `compare_loo.py`, `loo/compare.tsv`)
v0.8.1 on the 134 *Ogataea* genomes with the base DB (`base`), and with the six records added (`cand`: 128 other genomes; `self0-5`: each source genome with its own record excluded by `--exclude-records`; the *boidinii* genomes still see the other strain's record).

| | same calls | calls differ | idiomorph gained | idiomorph lost |
|---|---|---|---|---|
| two cassettes (71) | 52 | 19 | 0 | 0 |
| one cassette (56) | 40 | 16 | **0** | 0 |
| none (7) | 7 | 0 | 0 | 0 |

- **No missing cassette was recovered.** The 56 one-cassette genomes (including the 22 *C. boidinii*) stay at one cassette.
- **Scores and confidence rise in the divergent species**: e.g. *O. glucozyma* MATalpha 108 -> 134, *O. xyloterini* 131 -> 186, *O. nitratoaversa* and the four *O. sonorensis* genomes medium -> high with a fourth gene, *O. methanolica* medium -> high. The 49 polymorpha-group genomes (269-409) are unchanged.
- **One new call to look at**: *C. boidinii* NDK27A1 gains a second MATalpha call, MSSE01000034.1:49,406-85,483 (2 genes, score 38.5, 36 kb span), which is probably noise.
- The self-excluded source genomes also gain score from the other species' records (e.g. *O. wickerhamii* MATa 78.6 -> 80.5), so the records generalise a little across species.

Conclusion: the candidates raise confidence for the divergent species and change nothing for the polymorpha group, but they do **not** solve the one-cassette gap of #228; that gap needs the missing cassette ends to be present in the assemblies (most sit on contig fragments) or a different search (#254-style review of the withheld loci).
