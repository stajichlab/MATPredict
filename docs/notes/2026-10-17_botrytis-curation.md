# Botrytis cinerea MAT1-2 record curated (2026-10-17)

Branch `botrytis-curation`, from origin/main (624048d). Adds `db/Ascomycota/Helotiales/40559_sas405_MAT_MAT1-2/` (metadata.yaml, locus.gbk, locus.gff3, proteins.faa) from GenBank KF944386.1 (strain SAS405, `/note="type: MAT1-2"`, 4,501 bp). No change to `src/`.

## Outcome

The record alone does NOT flip the 17 stated-MAT1-2 genomes. They still call MAT1-1 (high). The MAT1-2 score rose (222 -> 326-331 bits) but the MAT1-1-1 remnant still scores 410-421. The stated-MAT1-1 genomes are unchanged (MAT1-1, high). The record is correct and useful, and is a prerequisite for any vote fix, but a `src/` change (below) is needed to correct the calls.

## Record decisions

- Naming follows `40559_sas56_MAT_MAT1-1`: `<taxid>_<strain lower>_MAT_<idiomorph>`.
- Genes: MAT1-2-1 = AHX22634.1 (complement(join(423..898,954..1327,1378..1658)), HMG); the deposit's MAT1-2-4 = AHX22635.1 (join(2404..2594,2644..2882,2931..3115)) recorded as **MAT1-2-10**, the name `db/Ascomycota/order.yml` declares (MAT1-2-4 in order.yml is a distinct Eurotiomycete gene). Exons are recorded, so validate re-derives both proteins from the nucleotide record at 100% identity / 100% coverage.
- MAT1-1-1 remnant: KF944386.1 carries an unannotated 5'-truncated MAT1-1-1 at about nt 3267..3960 (+ strand; miniprot of AHX22632.1: reference residues 140-353, 214 of 353 aa, 96.3% identity; the N-terminal 139 aa incl. the alpha-box region are absent). The schema has no remnant/pseudogene flag; the only mechanism is `present: false`. It is recorded as gene_index 2, MAT1-1-1, `present: false`, no coordinates, no protein. `present: false` genes are skipped by the vocabulary check, GFF/GenBank/FASTA export and the reference FASTA, so the remnant contributes no sequence and is not a searchable alpha-box gene. The remnant span and nature are documented in `definition_note` only. Caveat: `present: false` elsewhere in the DB means "documented absence" (e.g. Fusarium/Teloschistales records); here it means "no functional gene", explained in the note.
- Evidence: all three axes tier 2 (GenBank REFERENCE is "Unpublished", as for KF944385.1). Cited, with abstracts opened via PubMed: Abate et al. 2018, PMID 29767553 (Sclerotiniaceae MAT1-2 = MAT1-2-1 + MAT1-2-10; dMAT1-1-1 fragments near the flanks), and Wilson et al. 2023, PMID 37824487 (abstract only). Not cited: Wilken et al. 2017 (could not find it in PubMed, so the MAT1-2-4 -> MAT1-2-10 rename rests on order.yml's comment and on Abate 2018); Mercier et al. 2021 (PMID 33829853) only appears here as context, abstract read, MAT calls not in it.
- No same-strain assembly anchor exists for SAS405 (deposit only), so no assembly locus tags. No genome in BioProject PRJNA624742 is SAS405, so there is no self-hit risk.
- `validate` printed `taxonomy_current: false` because `taxonomy.resolve_lineage` runs `taxonkit reformat -I 1` on a bare taxid, which returns empty with taxonkit 0.20.0 on the HPCC (`taxonkit lineage` resolves 40559 fine to Botrytis cinerea). I set `taxonomy_current: true` by hand after that check. Environment/tool issue to be fixed separately (not edited).
- `curate-db accept` and `build-gff` ran on the HPCC; records pass the db tests.

## Before / after detect

`detect --taxid 40559`, MATPREDICT_HTML=0, v0.6.2; before = main db (166 records), after = main db plus the new record. Stated = BioSample (T4 MAT1-2 and B05.10 MAT1-1 from the literature/known layout). Scores are `idiomorph_candidates` (best bitscore per idiomorph).

| genome | stated | before (conf; scores) | after (conf; scores) |
|---|---|---|---|
| Hm1 | MAT1-1 | MAT1-1 (high; MAT1-1=532,MAT1-2=70) | MAT1-1 (high; MAT1-1=532,MAT1-2=85) |
| Rf1 | MAT1-2 | MAT1-1 (high; MAT1-1=416,MAT1-2=222) | MAT1-1 (high; MAT1-1=416,MAT1-2=326) |
| Sl12 | MAT1-1 | MAT1-1 (high; MAT1-1=540,MAT1-2=74) | MAT1-1 (high; MAT1-1=540,MAT1-2=91) |
| Sl11 | MAT1-2 | MAT1-1 (high; MAT1-1=414,MAT1-2=222) | MAT1-1 (high; MAT1-1=414,MAT1-2=326) |
| Sl10 | MAT1-2 | MAT1-1 (high; MAT1-1=414,MAT1-2=222) | MAT1-1 (high; MAT1-1=414,MAT1-2=326) |
| Sl9 | MAT1-1 | MAT1-1 (high; MAT1-1=540,MAT1-2=70) | MAT1-1 (high; MAT1-1=540,MAT1-2=87) |
| Rf2 | MAT1-2 | MAT1-1 (medium; MAT1-1=152) | MAT1-1 (medium; MAT1-1=152) |
| Sl8 | MAT1-1 | MAT1-1 | MAT1-1 (medium | medium; MAT1-1=540 | MAT1-1=389,MAT1-2=70) | MAT1-1 | MAT1-1 (medium | medium; MAT1-1=540 | MAT1-1=389,MAT1-2=87) |
| Sl7 | MAT1-2 | MAT1-1 (high; MAT1-1=414,MAT1-2=222) | MAT1-1 (high; MAT1-1=414,MAT1-2=326) |
| Sl6 | MAT1-1 | MAT1-1 (high; MAT1-1=540,MAT1-2=70) | MAT1-1 (high; MAT1-1=540,MAT1-2=87) |
| Sl4 | MAT1-1 | MAT1-1 (high; MAT1-1=540,MAT1-2=74) | MAT1-1 (high; MAT1-1=540,MAT1-2=91) |
| Sl5 | MAT1-2 | MAT1-1 (high; MAT1-1=414,MAT1-2=222) | MAT1-1 (high; MAT1-1=414,MAT1-2=326) |
| Vv16 | MAT1-2 | MAT1-1 (high; MAT1-1=421,MAT1-2=234) | MAT1-1 (high; MAT1-1=421,MAT1-2=331) |
| Sl2 | MAT1-1 | MAT1-1 (high; MAT1-1=540,MAT1-2=74) | MAT1-1 (high; MAT1-1=540,MAT1-2=91) |
| Sl1 | MAT1-1 | MAT1-1 (high; MAT1-1=540,MAT1-2=74) | MAT1-1 (high; MAT1-1=540,MAT1-2=91) |
| Vv15 | MAT1-2 | MAT1-1 (high; MAT1-1=412,MAT1-2=222) | MAT1-1 (high; MAT1-1=412,MAT1-2=326) |
| Vv14 | MAT1-2 | MAT1-1 (high; MAT1-1=421,MAT1-2=234) | MAT1-1 (high; MAT1-1=421,MAT1-2=331) |
| Vv13 | MAT1-2 | MAT1-1 (high; MAT1-1=421,MAT1-2=234) | MAT1-1 (high; MAT1-1=421,MAT1-2=331) |
| Vv12 | MAT1-1 | MAT1-1 (high; MAT1-1=540,MAT1-2=70) | MAT1-1 (high; MAT1-1=540,MAT1-2=87) |
| Vv11 | MAT1-2 | MAT1-1 (high; MAT1-1=410,MAT1-2=222) | MAT1-1 (high; MAT1-1=410,MAT1-2=326) |
| Vv10 | MAT1-2 | MAT1-1 (high; MAT1-1=421,MAT1-2=234) | MAT1-1 (high; MAT1-1=421,MAT1-2=331) |
| Vv9 | MAT1-2 | MAT1-1 (high; MAT1-1=414,MAT1-2=222) | MAT1-1 (high; MAT1-1=414,MAT1-2=326) |
| Vv8 | MAT1-1 | MAT1-1 (high; MAT1-1=540,MAT1-2=70) | MAT1-1 (high; MAT1-1=540,MAT1-2=87) |
| Vv7 | MAT1-1 | MAT1-1 (high; MAT1-1=540,MAT1-2=70) | MAT1-1 (high; MAT1-1=540,MAT1-2=87) |
| Vv6 | MAT1-1 | MAT1-1 (high; MAT1-1=540,MAT1-2=70) | MAT1-1 (high; MAT1-1=540,MAT1-2=87) |
| Vv5 | MAT1-2 | MAT1-1 (high; MAT1-1=421,MAT1-2=234) | MAT1-1 (high; MAT1-1=421,MAT1-2=331) |
| Sl3 | MAT1-2 | MAT1-1 (high; MAT1-1=414,MAT1-2=222) | MAT1-1 (high; MAT1-1=414,MAT1-2=326) |
| Vv4 | MAT1-1 | MAT1-1 (high; MAT1-1=540,MAT1-2=70) | MAT1-1 (high; MAT1-1=540,MAT1-2=87) |
| Vv3 | MAT1-1 | MAT1-1 (high; MAT1-1=540,MAT1-2=70) | MAT1-1 (high; MAT1-1=540,MAT1-2=87) |
| Vv2 | MAT1-1 | MAT1-1 (high; MAT1-1=540,MAT1-2=70) | MAT1-1 (high; MAT1-1=540,MAT1-2=87) |
| Vv1 | MAT1-2 | MAT1-1 (high; MAT1-1=421,MAT1-2=234) | MAT1-1 (high; MAT1-1=421,MAT1-2=331) |
| Pv1 | MAT1-2 | MAT1-1 (high; MAT1-1=396,MAT1-2=224) | MAT1-1 (high; MAT1-1=396,MAT1-2=322) |
| T4 | MAT1-2 | MAT1-1 (high; MAT1-1=414,MAT1-2=222) | MAT1-1 (high; MAT1-1=414,MAT1-2=326) |
| B05.10 | MAT1-1 | MAT1-1 (high; MAT1-1=540,MAT1-2=70) | MAT1-1 (high; MAT1-1=540,MAT1-2=87) |

Rf2 (GCA_015148645.1): only a 222 bp MAT1-1-1 remnant on the contig with SLA2 (APN2 and the rest of the locus on other contigs), so MAT1-1 (medium) with no MAT1-2 candidate; unchanged by the new record. Sl8's two reported loci are both MAT1-1 (medium).

## Diagnosis (Rf1 after, same pattern in all 17)

Gene evidence, locus JACVFJ010000135.1: MAT1-2-1 99.73% to the new record (exonerate_refine, polished_agree, bitscore 326), MAT1-2-10 99.01% (251), MAT1-1-1 remnant 92.69% to sas56 (416, polished_disagree). Resolutions list only MAT1-2-1 beating MAT1-1-3 (a 105 bp tblastn bit); there is no resolution between the remnant and MAT1-2-1. `idiomorph_candidates` takes the best single `bitscore` per idiomorph (idiomorph.py), and that bitscore is the first-pass hit's, not the polished full-length gene's: a 3-exon, 376 aa MAT1-2-1 at 99.7% scores only 326, while an intron-free 215 aa remnant at 92.7% scores 416. That is consistent with a per-HSP score on an intron-split gene (not verified at the hit level), and the earlier agent's exonerate raw scores (1914 vs 1085) invert the ranking. Margin is now 90 bits MAT1-1 (was 194), still `high`.

Needed fix (src, not done): weight or filter votes by reference coverage (the remnant covers 61% of its reference) and/or score by the polished or summed alignment rather than best single HSP, or tie-break with the HMG vs alpha profiles (branch ascomycota-core-hmms). See docs/notes/2026-10-17_botrytis-mislabels.md on branch botrytis-mislabels.

## Tests

HPCC, `pixi run pytest tests/db tests/detect -q` (PYTHONPATH and MATPREDICT_DB_ROOT set to this worktree): tests/db 138 passed (test_draw.py and test_synteny.py excluded: no `pygenomeviz` in the environment). tests/detect 1029 passed, 31 failed, all `ModuleNotFoundError: pyhmmer` (test_idiomorph_classifier, test_paralog_class, test_classifier_build_deterministic): environment, unrelated to this record; the same module is absent without the record. Not run: the baseline for those 31 on main.
