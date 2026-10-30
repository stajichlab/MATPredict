# Petasospora (Cyberlindnera) fabianii: our annotation of the MAT locus (2026-10-30)

**Question.** The only public lead was ONH67106 "Mating-type-like protein ALPHA1" (strain 65, MPUK01000005.1, a depositor annotation with no paper on the locus), and the v0.7.0 replay called a 55 kb MATsc locus on GCA_022641835.1 (CP048715.1:630,345-685,843) with DIC1 beside it. The curator decided a manual-annotation record at **tier 3**.

**Answer in one line.** The MAT locus is a small block between DIC1 and SLA2 at the same place in all four assemblies; JOY008 and JCM 3601 carry **MATa** (MATa1 homeodomain + MATa2 HMG-box), TBRC 4498 and strain 65 carry **MATalpha** (MATalpha1 alpha-box). It is not 55 kb: the genes span 16 kb (DIC1 to SLA2) and the MAT genes themselves 2.3 kb (MATa) or 0.9 kb (MATalpha).

## Evidence available
- **No RNA-seq exists.** ENA has 24 runs for tax 36022 (checked 2026-10-30, `data/ena_runs.tsv`): all WGS or amplicon, none RNA-Seq. The models are therefore homology and ab initio gene prediction, not transcript-supported.
- Four assemblies, all aligned with miniprot (`pf_miniprot.slurm`; the 63 curated Saccharomycotina MAT proteins plus ONH67106 and the *Williopsis jadinii* HMG-box protein XP_020071302): JOY008 GCA_022641835.1 (complete genome, ecotype Missouri, SAMN14051719), TBRC 4498 GCA_050575385.1 (chromosome; xylitol strain, PubMed 40558965, [DOI](https://doi.org/10.3390/jof11060453)), strain 65 GCA_001983305.1 (draft from fermented masau fruit, PubMed 28385833, [DOI](https://doi.org/10.1128/genomeA.00064-17)), JCM 3601 GCA_001599195.1 (type strain, scaffold). According to PubMed, no paper describes the MAT locus of this species; the YJS4271 genome (PubMed 25103752, [DOI](https://doi.org/10.1128/genomeA.00638-14)) is not annotated for MAT.
- AUGUSTUS 3.5.0 (`pf_augustus.slurm`, `classify_aug.py`) with four yeast species models on the four 19-22 kb windows (`make_windows.py`), classified by phmmer against the references, Pfam-A and the Saccharomycotina MATa2 profile. `scan_window.py` independently scans all six frames. Models checked in `build_models.py` (start codon, stop, introns).

## Layout and models (CDS coordinates include the stop codon; minus-strand exons listed in transcript order)
| allele (strain, contig) | gene | strand | CDS | aa | evidence |
|---|---|---|---|---|---|
| MATa (JOY008, CP048715.1) | dic1 | - | 628,449-629,375 | 308 | 374 bits to curated DIC1; 4 of 4 species models agree |
| | **MATa2** | - | 630,863-630,492 + 630,419-629,850 | 313 | HMG-box, Saccharomycotina MATa2 profile **76.8 bits** (P. kudriavzevii 52.4); intron 630,420-630,491 (72 bp) GTAGGT..TACTAAC..TAG; P. stipitis model |
| | **MATa1** | + | 631,437-631,737 + 631,835-632,136 | 200 | homeodomain: Pfam Homeodomain 64 bits, 66 bits to Komagataella MATa1; intron 631,738-631,834 (97 bp) GTACCA..AAAG (weak, see limits) |
| | sla2 | + | 641,746-644,871 | 1,041 | I_LWEQ + ANTH, 1,233 bits to curated SLA2; single exon, 4 of 4 models |
| MATalpha (TBRC 4498, CM115473.1) | dic1 | - | 631,646-632,572 | 308 | as above |
| | **MATalpha1** | + | 637,540-638,457 | 305 | alpha box PF04769 36 bits; 407 bits to ONH67106; single exon, 4 of 4 models, identical protein in strain 65 (MPUK01000005.1:617,300-618,217) |
| | sla2 | + | 648,035-651,160 | 1,041 | as above |

Between the MAT genes and SLA2 lie six conserved non-MAT ORFs (409, 173, 589, 219, 465, 622 aa) in the same order in both alleles; SLA2 is therefore about 9.5 kb from the MAT genes. DIC1 is 0.5 kb from MATa2 in the MATa allele and 5.0 kb from MATalpha1 in the MATalpha allele (the 800 aa ORF of finding 2 lies between). `data/models.tsv`, `data/models_proteins.faa`, `data/aug_classified.tsv` hold everything; the JCM 3601 (MATa) and strain 65 (MATalpha) copies are in the same table. The two MATa copies differ by one amino acid each in MATa1 and MATa2.

## Findings that matter for the database
1. **Heterothallic-looking, no second cassette.** No MAT gene other than the locus above is found by miniprot in any of the four assemblies (only unrelated TGL1 hits elsewhere), and each genome carries one idiomorph.
2. **The two idiomorphs are at one position** (DIC1 side) with different gene content: MATa = MATa2 + MATa1 (two genes, 2.3 kb); MATalpha = MATalpha1 (one gene). The alpha allele also has an 800 aa MutS_III-domain ORF (TBRC 4498 CM115473.1:633,233-635,626) between DIC1 and MATalpha1 that has no counterpart in the MATa allele; not recorded.
3. **The depositor's ONH67106 model is probably wrong at the N-terminus.** It joins a 69 aa exon (616,190-616,396) to the 185 aa alpha-box exon across a 1,263 bp intron with a poor donor (GTCCTT); the open frame upstream of the alpha box is continuous, giving a 305 aa single-exon protein that all four species models produce. The exon-1 sequence matches no reference. No RNA can settle it.
4. **MATalpha2 not found** in either alpha allele: only a 57 aa Homeobox_KN fragment at 12-18 bits about 1.2 kb upstream of MATalpha1 (TBRC 4498 636,200-636,370 minus strand). Treated as absent; the draft family makes MATalpha2 optional.
5. **MATsc vs MTL.** This locus has DIC1 0.5-5 kb away and SLA2 9.5 kb away, no PAP1/OBP1/PIK1 in the window: the SLA2-type layout, not MTL, and MATa1/MATa2 are not MATsc cassettes either. A scoped family is the natural home (`order_MATpf_proposal.yml`, scope species 36022, not applied).

## Limits (tier 3 for these reasons)
- No transcript evidence: every intron, start and stop is prediction. MATa2's intron has the full yeast signal set and the profile is strong; **MATa1's intron has a weak donor (GTACCA)** and an upstream exon is only required because the homeodomain open frame (631,726-632,133) has no start codon; its N-terminal 100 aa match no reference. The K. lactis model puts the acceptor at 631,864 (127 bp), the S. cerevisiae model at 631,834 (97 bp, same donor); I used the latter.
- **Gene names are by homology and position**: homeodomain = MATa1 (it could in principle be an alpha2 in the other allele; none was found there), HMG-box = MATa2, alpha box = MATalpha1.
- MATa2 is a 313 aa protein with a divergent HMG box (no Pfam HMG_box hit above GA, as for Komagataella and Pichia).
- Taxon: NCBI Taxonomy lists the species as *Petasospora fabianii* (taxid 36022; synonyms *Cyberlindnera fabianii*, *Hansenula fabianii*), class Saccharomycetes, order Phaffomycetales, family Wickerhamomycetaceae.
- Strain provenance: JOY008 has no publication that I found; TBRC 4498 and strain 65 are cited; the type strain JCM 3601 is the independent MATa replicate.

## Files
- Drafts (not in db/): `36022_joy008_MATpf_MATa.yaml`, `36022_tbrc-4498_MATpf_MATalpha.yaml`, written by `make_records.py` (metadata schema, gene-status and exon-bound checks pass; the gene vocabulary and idiomorph checks pass against `order_MATpf_proposal.yml` appended to `db/Ascomycota/order.yml`). `curate-db validate` (taxonkit, HPCC) has not been run.
- Scripts: `pf_miniprot.slurm`, `pf_augustus.slurm`, `summarize_miniprot.py`, `make_windows.py`, `scan_window.py`, `classify_aug.py`, `build_models.py`. Genome FASTAs and `work/` are not committed (`.gitignore`); `data/` is 780 KB.
