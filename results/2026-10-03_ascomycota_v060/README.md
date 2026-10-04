# Full Ascomycota run on v0.6.0 (2026-10-03)

All BFD Ascomycota genomes with a taxid and a genome file (samples.csv), lineage routing
from the taxid. Code: frozen worktree run-7c7ed99 (= v0.6.0). Gate: regression 0 changed loci;
pilot ../2026-10-03_ascomycota_pilot (246 genomes, mean 158 s/genome, max 624 s).
9 waves (classes/orders interleaved), 64-way each, partition exfab, GENOME_TIMEOUT 3600.

Known defect in v0.6.0: genetic code 26 (Alaninales, CUG=Ala) makes exonerate exit
("No built in genetic code corresponding to id [26]") and the genome gets no report.
Exonerate has no table for codes >= 24. These genomes are re-run after a fix, from a
new frozen worktree, into a separate folder; the re-run is recorded here.
