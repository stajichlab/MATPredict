# 014. sexP-type HMG genes detected in Lichtheimiaceae and Syncephalastraceae

- **Category:** biology; curation-first (no deposit exists)
- **Status:** candidate
- **Lineage:** Mucoromycota, Mucorales: Lichtheimiaceae, Syncephalastraceae

## Summary
No MAT locus deposit exists for these families, and their flank order differs
from other Mucorales. HMG genes at two loci fall in the supported sexP clade
and are classified sexP by HMMs and blastp.

## Evidence
- Syncephalastrum racemosum GCA_000696955.1 (called, Plus): sexP by full HMM
  (margin 83.1), HMG-box HMM (23.5), blastp (52.4).
- Rhizomucor pusillus GCA_023512895.1 (withheld): 100.8, 32.4, 78.2.
- Lichtheimia ornata GCF_029851405.1 (withheld): unresolved (-8.9, -5.7, 0.8).
- Source: `results/2026-09-26_sexMP_hmm/disputed_scores.tsv`; tree
  `results/2026-09-26_sexMP_phylogeny/` (sexP clade UFBoot 99), both committed.
- Flank order: tptA and rnhA never within 100 kb in Lichtheimiaceae (0/30) or
  Syncephalastraceae (0/9) (`results/2026-09-26_flank_synteny_ED/NOTE.md`).
- No deposit: NCBI and PubMed checked (`docs/notes/2026-09-26_mucorales-dothideomycetes-curation.md`, curation-mucor-dothideo).

## Method that found it
hmmalign (PF00505) + IQ-TREE; leave-one-genus-out sexM/sexP HMMs; blastp.

## Verification done / still open
Open (curator ruling): propose a tier-2 record only when a candidate sits in
the sexP clade at UFBoot >= 95 in the final ML tree.

## Limits
One 69-column domain; per-leaf support not reported.

## Related
Entries 015, 016.
