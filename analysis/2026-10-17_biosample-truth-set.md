# A mating-type truth set from BioSample (2026-10-17)

Code: `scripts/harvest_biosample_mating_types.py` (assembly accession to BioSample through the NCBI fungal assembly
summaries, one E-utilities fetch per BioSample, standard library only) and `scripts/biosample_truth_set.py`
(normalise the stated value, join with a campaign). Data: `results/2026-10-17_biosample_truth/`
(`truth_set.tsv`, `agreement.tsv`, `biosample_mating_attributes.tsv`).

The stated mating type is the submitter's free text and unverified. It is independent of this pipeline, but it is
not ground truth, so a disagreement is something to read, not an error by assumption.

## Coverage
Of the 19,415 Ascomycota genomes in the v0.6.0 campaign, 19,363 map to a BioSample in the assembly summaries, and
E-utilities returned a record for 17,697 of them (the other 1,666 came back empty and were not investigated, so the
truth set may be missing some); 339 of the 17,697 (1.9%) carry a mating-type attribute (`mating_type`, `mating-type`, `mating type`). `samp_mat_type` is the MIxS sample
material type (values such as `solid` and `petroleum`) and is excluded. After normalisation 296 genomes remain:
MAT1-2 103, MAT1-1 100, MATa 58, MATalpha 8, homothallic or diploid 27. 43 values are not a mating type.
None of the 296 is a curated record's assembly or strain, so none is in the database.

## Agreement with the v0.6.0 campaign (not a withheld-record test)
Pezizomycotina with a stated MAT1-1 or MAT1-2: 203 genomes. The campaign label equals the stated one in 178 (88%):
Eurotiomycetes 98 of 101, Sordariomycetes 62 of 66, Leotiomycetes 16 of 33, Dothideomycetes 2 of 2. Excluding the
17 *Botrytis* genomes below it is 178 of 186 (96%).
The Neurospora genomes stated `A` or `a` agree with MAT1-1 and MAT1-2 in 13 of 13 when A is read as MAT1-1 and
a as MAT1-2 (case-sensitive; the normaliser leaves them out of the MAT1-1/MAT1-2 set).

The 25 that do not agree:
- 17 *Botrytis* (16 *B. cinerea*, 1 *B. fabae*): stated MAT1-2, called MAT1-1 with high confidence. In 16 the locus
  models both core genes: MAT1-1-1 at 91 to 94% identity (about 415 bits) and MAT1-2-1 at 77 to 79% (222 to 234
  bits), plus MAT1-2-10. The curated *Botrytis* record is MAT1-1, so a 93% MAT1-1-1 hit is a real alpha-box gene.
  These genomes carry both idiomorphs' genes at one locus (mixed or heterokaryotic isolates, a mislabel, or a
  stated value that means something else); not resolved.
- *Penicillium roqueforti* 3 (stated MAT1-2: called MAT1-1 in 2, not called in 1), *Cryphonectria parasitica* 2
  (stated MAT1-1, not called), and one each *Morchella conica*, *Podospora tetraspora* and a Magnaporthales sp.
  (called the opposite idiomorph).

## Open
- Whether the 16 *Botrytis* genomes are mixed-idiomorph isolates; the assemblies and their read data would show it.
- Running the truth set with the withheld-record held-out runner (genus and order radius) rather than the
  campaign's all-records run. The 296 genomes are not curated records, so the genomes are held out already; the
  radius removes their relatives from the references.
- 58 MATa, 8 MATalpha and 27 homothallic or diploid values are Saccharomycotina, Sordariomycetes `a/A` (Neurospora)
  and diploid yeasts; they are in the truth set but not compared here beyond Neurospora.
