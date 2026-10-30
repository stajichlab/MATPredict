# Barnettozyma californica: our own annotation of the MAT locus (2026-10-30)

Draft only. Nothing is written to `db/`. Everything here was produced from sequence; **there is no RNA-seq for the genus** (ENA, 2026-10-30: *Barnettozyma* runs are WGS and Hi-C only), so the models are predictions, **tier 3**.

## Literature
According to PubMed, Hession, Byrne, Wolfe & Butler, PLoS Genet 2025 ([10.1371/journal.pgen.1011814](https://doi.org/10.1371/journal.pgen.1011814), PMID 41359681) report that a centromere formed at the MAT locus of *B. californica*; its 460 bp centromere sequence lies in a **2.8 kb intergenic region between MATa1 and MATa2 on one side and MATalpha1 and MATalpha2 on the other**; all three *Barnettozyma* species appear homothallic (both MATa and MATalpha genes present). The paper gives no gene models or coordinates for the locus; its UCD09 ONT assembly is not in the NCBI assembly database; the S1 Table was not available (the supplement bundle holds figures only).

## Assemblies
Replay MTL calls of v0.7.0 (`results/2026-10-29_near_tie/flank_synteny.tsv`): UCD09 GCA_003123585.1, UCD89 GCA_003123885.1, NRRL Y-17395 GCA_003706655.3, PB4207 GCA_023627795.1 (*B. californica*, taxid 36038) and *B. discipulorum* GCA_049863505.1 (comparator). All are scaffold-level drafts except *B. discipulorum* (complete). **Primary: PB4207** (scaffold N50 421,778; the locus lies on a 694 kb contig, 63 kb from one end); the other three californica assemblies carry the same proteins at 99.6-100% identity (`data/*.lift.gff`), i.e. they are the same strain family and not independent evidence.

## Method (all scripts here; window = call +/- 25 kb, `fetch_windows.sh`)
1. `run_miniprot.sh`: curated Saccharomycotina MAT/flank proteins (`../2026-10-29_ss_records/refs.faa`) against each window. `run_miniprot_a2.sh`: the seven curated MATa2 proteins (+ *P. kudriavzevii*) at relaxed settings: **no hit** in any window, so MATa2 is not detectable by protein homology.
2. `orf_scan.py` / `run_orf_scans.sh`: six-frame stop-to-stop segments, Pfam (Homeodomain, MATalpha_HMGbox, ANTH, I_LWEQ) and phmmer against the references.
3. `scan_one_intron.py`, `scan_two_introns.py`: exhaustive GT..AG intron searches (one, or up to two) in the region of each gene, scored with the curated Saccharomycotina MATa2 profile (`db/Ascomycota/classifiers/MATa2_saccharomycotina`, supporting evidence only), a Pfam profile or phmmer against curated proteins.
4. `models.py`: the chosen models, checked for start codon, single stop, canonical introns -> `models.gff3`, `models_proteins.faa`, `models.tsv`. `a2_benchmark.py`: strength of the MATa2 call. `make_records.py`: the draft record. `MATbz_family_proposal.yml`: the proposed family block.

## Models (PB4207, `JAMAIE010000005.1`; 1-based, strand on the contig)
| gene | strand | CDS | exons | aa | evidence |
|---|---|---|---|---|---|
| DIC1 | - | 605,238-606,188 | 1 | 316 | mitochondrial carrier; phmmer 362 bits to curated DIC1 |
| **MATa2** | - | 613,590-614,455 | 2 (intron 614,112-614,167, 56 bp GTATGT..CTTAG) | 269 | HMG-like; 82 bits to the Saccharomycotina MATa2 profile vs a maximum of 32 over all 124 six-frame segments >= 80 aa of the window (curated MATa2 proteins score 144-256; *P. kudriavzevii* 52); no Pfam HMG_box hit; the 56 bp intron beat the next donors (64 and 50 bits) |
| MATa1 | + | 615,080-615,568 | 1 | 162 | Pfam Homeodomain 53 bits, helix 3 `WFTNKR`; 56% identity to *Ogataea* MATa1 over the homeodomain; divergent partner of MATa2 (5' ends 625 bp apart) |
| MATalpha1 | + | 618,402-619,103 | 1 | 233 | MATalpha_HMGbox 33 bits, 62-64 bits to curated MATalpha1; zero-intron ORF; **C-terminus uncertain** (one-intron alternatives with introns of 45-202 bp score within 1-2 bits) |
| SLA2 | - | 627,731-630,871 | 1 | 1,046 | ANTH + I_LWEQ, 1217 bits to curated SLA2 |
| MATalpha2 | | **not modelled** | | | a homeodomain region follows MATalpha1 (genome 620,370-620,700; `WISNKR` helix-3-like, introns of 45-61 bp GTATGT..TAG in the best two-intron scans), but no complete gene could be built (best partial: 104 aa, 22-26 Pfam bits) |

Agreement with the paper: MATa1 ends at 615,568 and MATalpha1 begins at 618,402, an **intergenic region of 2,833 bp (paper: 2.8 kb)**, with the a-side genes (MATa2/MATa1, divergent) on one side and the alpha-side on the other. DIC1 lies 7.4 kb and SLA2 8.6 kb outside the cassette block (the replay saw DIC1 about 8.2 kb from the call).

## Limits
- **No RNA, no annotation to check against**: every intron, start and stop is predicted. MATalpha1's end and the whole of MATalpha2 are unresolved.
- **Gene names are by position**: Pfam cannot tell a MATa1 from a MATalpha2 homeodomain. The homeodomain next to the HMG-like gene is named MATa1; the one after MATalpha1 would be MATalpha2.
- **MATa2 identity rests on a profile built from seven proteins** (thin; the profile is not used by any call) and on position; the same was true for *Komagataella*.
- *B. discipulorum* (`GCA_049863505.1`) has SLA2, MATa1, MATalpha1 and DIC1 at 45-86% identity but arranged differently (MATa1 and MATalpha1 adjacent, both minus strand) and no MATa2 hit: the arrangement is not conserved in the genus.
- The centromere-associated 460 bp sequence (S3 Fig) was not checked.
- Genetic code assumed standard.

## Blocks to a record
The draft `36038_pb4207_MATbz_MATa-MATalpha.yaml` passes the JSON-schema check (`validate_metadata`) and exon bounds, but has no family (`MATbz` is only proposed), is tier 3 throughout, omits MATalpha2, and its assembly anchor is the PB4207 assembly itself (the paper's UCD09 assembly is not public). `curate-db validate` (taxonkit) was not run. A tier-2 record would need the authors' annotation (S1 Table accessions, the YGAP annotation used in the paper) or RNA-seq.
