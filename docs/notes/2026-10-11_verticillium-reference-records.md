# Verticillium reference records: withheld-locus rise and regression panels

Date: 2026-10-11. Written from the Verticillium MAT run in
`/bigdata/stajichlab/shared/projects/KassonLab/Verticillium/nf/pipeline/mating_type/`.
Records: 3 accepted *V. dahliae* records in `db/Ascomycota/Glomerellales/` and 8 proposed records in `db/candidates/Ascomycota/`.

## Why the withheld count rises (four Verticillium genomes)

The final run has 16 to 20 withheld loci per genome, against 12 to 15 in the baseline. All 17 new loci in the four genomes have the same pattern. This was checked on 2026-10-11 by searching each genome with `tblastn` as `matpredict detect` does (default e-value 10, `-seg no`).

- **What the new loci are.** Each is a cluster of 2 or 3 weak hits on a contig away from the real MAT locus. Examples: MAT1-2-1 and COX13 5 kb apart (MD4 `MTK05_134`); APN2, SLA2 and MAT1-1-3 over 30 kb (MD4 `MTK05_213`). All are in the `Ascomycota:MAT` family. All are withheld as `modelled_gene_bar`, and no gene was polished.
- **Cause.** The six new reference proteins (`AB505214`, `AB505215` and the four VdLs.17 proteins) add weak hits. In all 17 loci, at least one weak hit from a new protein lies inside the locus span, next to weak hits from older references. I did not extract the exact hit set that `detect` clustered. The likely mechanism is that the new hit completes a cluster of two distinct genes, which is enough to be reported as a withheld candidate. The baseline did not have these loci.
- **The hits are not significant.** Across the 17 loci, the best new-protein hit in each locus scores 26 to 35 bits, with e-values of 0.2 to 9.4. In the four MD4 loci, those hits are 28 to 43% identity and 34 to 72 aa long. This is the kind of hit that `tblastn` returns at e-value 10.
- **The reference proteins do not hit the new loci strongly.** At e-value 1e-5, the four VdLs.17 proteins hit the real MAT locus in MD4 (85 to 100% identity). Their one other hit is a weak MAT1-2-1 match (34%, 62 aa) on `MTK05_104`, which is not one of the new regions. None of them hits a new region at that cutoff.
- **Effect on calls.** None in these four genomes. The real locus is unchanged, and the rejected loci stay withheld.
- **Not tested.** The same cause is likely for the 273 withheld-only changes in the Ascomycota regression panel, but I checked only the four *Verticillium* genomes. I also did not test whether a stricter e-value cutoff in MATPredict would remove these loci.

## Regression check (MATPredict)

Before sign-off, the panel in `docs/regression-check.md` was run on a baseline worktree and on a candidate worktree that held the TV103 and NBRC 6119 records only. The VdLs.17 record was accepted later and **is not covered by this panel**.

- Output: `~/projects/MATPredict/results/2026-10-11_regression_vdahliae/diff/summary.md` (jobs 29858226 to 29858236).
- Panels: 30 Ascomycota, 33 Basidiomycota, 80 Mucoromycota, 23 record-source genomes, Zygo 23 (46 inputs).
- Basidiomycota, Mucoromycota and Zygo 23: no change.
- Record-source genomes: 5 loci changed, all withheld on both sides. No call changed.
- Ascomycota: 42 loci changed. 40 are withheld on both sides (listed in `regression_withheld_changes.md`). Two touch a call and need curator review:

| Genome | Species (order) | Change |
|---|---|---|
| `GCA_001557505.1` | *Reticulascus tulasneorum* (Glomerellales) | MAT1-1 call unchanged (medium). The MAT1-1-1 model changed from 432 bp, 39.6% identity to 495 bp, 35.6%. |
| `GCA_025802115.2` | *Bertia moriformis* (Coronophorales) | A locus (`JAOYNC020017679.1:1278-5158`) went from withheld (`below_fraction_floor`) to called MAT1-2 (medium). This genome already had called MAT1-1 and MAT1-2 loci. |

I have not judged whether these two changes are improvements or errors.

### Second panel: all 11 new records

The second panel tests the VdLs.17 record and the eight proposed records together with the two Usami records. The candidate worktree is `run-e9060daf32-vd11`: main plus 11 new records, accepted in the worktree only. The baseline is the same as in the first panel (its outputs were reused). The eight proposed records are not accepted in the live database.

- Output: `~/projects/MATPredict/results/2026-10-11_regression_vd11/diff/summary.md` (jobs 29860946 to 29860951).
- Basidiomycota (33), Mucoromycota (80) and Zygo 23 (46 inputs): no change.
- Record-source genomes (23): 34 loci changed, all withheld on both sides. No call changed.
- Ascomycota (30): 278 loci changed. 273 are withheld on both sides (up from 40 in the first panel; total loci rose from 1,115 to 1,239). Five touch a call:

| Genome | Species (order) | Change |
|---|---|---|
| `GCA_001557505.1` | *Reticulascus tulasneorum* (Glomerellales) | MAT1-1 unchanged; MAT1-1-1 model 432 bp, 39.6% → 495 bp, 35.6%. Also in the first panel. |
| `GCA_025802115.2` | *Bertia moriformis* (Coronophorales) | Withheld → called MAT1-2 (medium); gene set APN2, MAT1-1-3 → APN2, MAT1-2-1. Also in the first panel. |
| `GCA_004918105.1` | *Aureobasidium pullulans* (Dothideales) | MAT1-2 confidence medium → high; span changed; MAT1-2-1 model 454 bp, 47.0% → 459 bp, 49.3%. |
| `GCA_026023515.1` | *Claussenomyces* sp. (Leotiales) | MAT1-1 unchanged; a MAT1-1-2 model is added. |
| `GCA_026023555.1` | *Chrysothrix* sp. (Arthoniales) | MAT1-2 high unchanged; MAT1-2-1 model 264 bp, 43.2% → 252 bp, 50.0%. |

- The last three changes are new in this panel. I did not trace them to a specific record. Because all 11 records were added together, a change cannot be assigned to one record without checking the genes at that locus.
- The VdLs.17 genome is not in the panel, so the self-call check reports `no_report` there. A direct run called MAT1-2 (high) at the record's own genome. The other ten records have no source assembly (`not_possible`).

## Open items

1. Review the five call changes in the second panel (two of them also appear in the first).
2. Decide whether the rise in withheld loci is acceptable. Check whether the same cause explains the 273 withheld-only changes in the Ascomycota panel. Test whether a stricter tblastn e-value removes the new loci.
3. Decide whether to accept the eight proposed records. The *F. sacchari* MAT1-2 record leaves out MAT1-2-3, which is not in the MAT vocabulary of `order.yml`.
4. Nothing from this reference work is committed.
