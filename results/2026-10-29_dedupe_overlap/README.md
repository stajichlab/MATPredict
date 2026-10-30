# Overlap mode of the family-call dedupe (2026-10-29)

`scripts/dedupe_family_calls.py --match overlap` (opt-in; the default is still `identical`, the rule you approved on 2026-10-28) groups loci of different families on one contig when the overlap covers at least 0.9 of the shorter locus, then keeps the best-scoring family as before. Output only; no call, report or database changes.

On the v0.7.0 Ascomycota loci (`results/2026-10-26_v070_ascomycota/loci.tsv`):
- identical coordinates (default): 185 loci marked in 164 genomes;
- overlap mode: 264 loci marked in 233 genomes; dropped by family: MTL 149, mat2 69, MATsc 41, MAT 3, MATtub 2. The 79 additional loci are in `newly_marked_by_overlap.tsv` (63 MTL calls kept as MATsc, 14 MATsc kept as MTL, one MAT/MATtub pair each way). Orders: Phaffomycetales 35, Pichiales 21, Ascoideales 20.

Why: of the 80 yeast genomes that `results/2026-10-28_idiomorph_crosscall_flags/` sends to curation, in 73 one family's locus lies entirely inside the other's (the shorter locus is at least 90% covered; the median overlap is 5.8 kb): one locus called with different boundaries. Seven are second loci (6 on another contig, 1 more than 50 kb away). After the overlap rule, 12 of the 80 still show two undeduplicated yeast-family loci (7 genuine plus others with several loci per family; not examined).

Not checked: whether the kept family is the right one for each of the 79 (the score is the same locus-class, confidence, polished genes and genes-found order as before); a threshold other than 0.9 (73 of the 80 have a fraction of 0.9 or more and the other 7 have none, no genome falls between, so any threshold from just above 0 to 0.9 gives the same set among these 80 genomes); other taxa outside Saccharomycotina (the 2 MAT/MATtub rows are the only non-yeast additions).
