# Read-level confirmation of the flagged Coccidioides sets (2026-11-03)

Question: are the 13 low-level sets that GmG flagged (second genus at roughly 1-3% of the reads) real contamination, or a chance match?

Method (`datasets/reads/contamination_spikein/confirm_flagged.slurm`): map each set's reads to Coccidioides references (bowtie2), keep the read pairs that do not map,
and run GmG on those. If the contamination is real, the unmapped reads should be classified as the contaminant genus.

Result (`confirm_table.tsv`, 23 sets):
- All 13 low-level sets are confirmed: GmG names the expected contaminant genus as the second organism of the unmapped reads, with its own mixture flag raised, on 31-439 hashes
  (Candidozyma, Aspergillus, Malassezia, Nakaseomyces).
- The 4 positive controls with a known contaminant (M192 Histoplasma, TX7 Lecanicillium, M223 Calonectria, B12526 Candidozyma) are confirmed.
- The 6 negative controls are never flagged: second genus at 2-26 hashes, `mixed` False.

The unmapped fraction does not decide. A first version compared it with the largest control value, but one clean control (GT-134) leaves 19% unmapped because its strain maps
poorly to the references, so that baseline was dropped and the rule uses GmG's call on the unmapped reads instead.

Caveats: this shows that reads of the named genus are present, not their source (a lab contaminant, a mixed culture, or a labelling error). B10919 has 37% unmapped, so its
Nakaseomyces evidence sits among many other unmapped reads. Six negative controls are a small set for the false-positive rate; the 552 clean records in
`results/2026-11-02_contam_spikein/` give the larger estimate.
