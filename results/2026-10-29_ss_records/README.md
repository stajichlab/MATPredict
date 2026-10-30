# S. schoenii chromosome IV MAT loci: gene models rebuilt from public RNA-seq (2026-10-30)

**Why:** Kriti et al., G3 2026 ([10.1093/g3journal/jkag067](https://doi.org/10.1093/g3journal/jkag067), PMID 41849659) describe SLA2-MAT cassettes on every chromosome of *S. schoenii* (GDD12, GCA_054953955.1, WGS JBTFFQ01) and say the BRAKER2 GTF is on GitHub and Zenodo. Neither holds it (GitHub `dkriti/Saccharomycopsis-schoenii-Genome-Assembly` and Zenodo 10.5281/zenodo.18950794 carry scripts only; checked 2026-10-29), and the NCBI assembly has no annotation. The models below are rebuilt here and are **ours, not the authors'**.

**Method** (`ss_rnaseq.slurm`, `ss_models.slurm`, `chr4_orfs.py`, `chr4_models.py`, `make_refs.sh`): the three paired non-stranded RNA-seq runs of PRJNA1364152 (SRR37453601-03; 86% overall alignment) mapped with hisat2; StringTie transcripts; TransDecoder long ORFs; ORFs identified by phmmer against the 63 curated Saccharomycotina MAT proteins (`refs.faa`) and Pfam-A (HMG_box, Homeodomain, MATalpha_HMGbox, I_LWEQ). The library is not stranded, so transcript and ORF strands are inferred from the ORF alone.

## What agrees with the paper
- Chromosome IV is CM143129.1. Splice junction **536,635-536,692 (58 bp, 218 reads)** = the paper's conserved 58 bp intron at 536,635. The paper's mirrored junction at 529,128 is not seen (the two arms are an inverted repeat, so reads map to one copy).
- RNA depth in the left cluster ends near 474,500 (paper: SLA2 to 475,000); the right cassette has two blocks, 530-532 kb and 534-536 kb, with no depth at 532.5-533.5 kb (the paper's spacer).
- SLA2 (I_LWEQ domain, 546-712 bits to curated SLA2) is found as one ORF in each arm of the palindrome (CDS 529,992-532,490 and 533,981-535,828, mirrored) and, diverged, at 472,708-474,606: the paper's "SLA2 anchors".

## Models (ours)
| locus | ORF | CDS (CM143129.1) | evidence |
|---|---|---|---|
| left | MATalpha1 | 471,020-471,595, minus | 191 aa, MATalpha_HMGbox 54 bits, 73 bits to *Ogataea* MATalpha1 |
| left | HMG-like (a2-type) | 470,349-470,951 plus, or 470,176-470,243 + 470,336-470,951 with a 92 bp intron (97 reads) | 200/227 aa, 44 bits to *S. cerevisiae* MATA2, no Pfam hit |
| left | homeodomain | 471,961-472,227 + 472,445-472,489 plus | 103 aa (partial), Homeodomain 57 bits; the 217 bp intron (198 reads) joins it toward SLA2, the paper's "MATa1-SLA2 splice junction" |
| left | SLA2 | 472,708-473,478 and 473,680-474,606 | two diverged pieces, I_LWEQ 186 bits |
| right (IR copy 2) | homeodomain | 536,046-536,090 + 536,296-536,634 + 536,693-536,944, minus | 211 aa, Homeodomain 54 bits, introns of 205 bp (221 reads) and 58 bp (218 reads) |
| right (IR copy 1) | homeodomain | 529,258-529,749 | 163 aa, 34 bits; the 58 bp intron is not seen here |
Raw table: `data/chr4_orfs.tsv`; junctions `data/junctions_chrIV.tsv`; transcripts `data/stringtie.gtf`.

## What is not settled
- **Gene names.** A homeodomain protein is MATa1 or MATalpha2; Pfam cannot say which. The paper calls the left cluster alpha (MATalpha2, MATalpha1) and the right cassette a (MATa1, MATa2). I use those names only by position, as a hypothesis.
- **No HMG-box (MATa2) ORF in the right cassette** passed the filters; the left "HMG-like" ORF has no Pfam hit, as for *Komagataella* MATa2 (13.6 bits, below GA).
- Palindrome copies are collapsed by short-read mapping; strand of the unstranded transcripts is inferred.
- These are our own ORF calls, not an annotation of the deposit, so a record built on them is **tier 3**, not tier 2, until the authors' GTF is obtained or the models are checked another way.
