# Pucciniomycotina follow-up (2026-09-27)

Branch `curation-puccinio` (local). Curator rulings of 2026-09-27.

## Step A: rustHD on Sporidiobolales (option c). Result: reverted

Trial commit `8a49b00` extended `rustHD` (Pucciniales bE/bW, tier 2) to
Sporidiobolales (231213). Routing unions lineage matches, so both `redPR` and
`rustHD` ran. Run: frozen `run-8a49b00`, job 29118549, the same 22-genome list as
`results/2026-09-26_puccinio_curation/pilot_small.tsv` (15 Sporidiobolales + 7
Pucciniales as control). Tables: `stepA_summary.txt`, `stepA_sporid/`.

- rustHD calls in Sporidiobolales: **0 of 15**.
- 10 of 15 genomes have a co-located bE+bW cluster; best rustHD cluster identity
  per genome 41.9-56.7% (median 50.0). No gene in any rustHD cluster could be
  modelled, so every cluster was withheld by the 2-modelled-gene bar.
- redPR and Pucciniales calls: identical to `run-293640d` in 22 of 22 genomes.
- Median wall per Sporidiobolales genome: 195 s with rustHD vs 77 s without
  (2.5x).
- Per the ruling ("both genes < 30% or unmodelled" -> revert), the scope change
  was reverted in `a6770ca`. Sporidiobolales HD stays queued for curation.
- Not tested: whether the withheld bE+bW clusters are the Sporidiobolales HD
  locus. They could be, given 42-57% identity, but the rust models do not polish
  across that distance.

## Step B: putative tier-2 Wallemia record

Record `671144_cbs-633-66_wallMAT_v1` (commit `6590929`), family `wallMAT`,
scope Wallemiales (431958, order, verified). Build: `wallemia/build_wallemia.py`.

Locus on *W. mellicola* CBS 633.66 WALSEscaffold_2, `JH668224.1:17,113-42,263`:

| gene | protein | coordinates | role |
|---|---|---|---|
| BAP31 | EIM23475.1 (WALSEDRAFT_59214) | 17,113-18,142 + | flanking_conserved |
| STE3 | EIM23479.1 (WALSEDRAFT_59221) | 21,336-22,572 - | core (pheromone_receptor) |
| CAF1 | EIM23480.1 (WALSEDRAFT_14687) | 23,330-24,258 - | flanking_conserved |
| HMG | EIM23483.1 (WALSEDRAFT_62153) | 25,382-26,300 - | core (HMG_box) |
| SXI1 | none (unannotated) | 41,384-41,863, 41,907-42,263 + | core (no gene_class) |

How the genes were identified:
- The papers give no gene IDs. Gostincar et al. 2019 names the locus genes
  (BAP31, SXI1, HMG, STE3, CAF1). In CBS 633.66, BAP31, STE3 and CAF1 are
  adjacent on scaffold_2; the HMG gene is the Pfam HMG_box hit next to CAF1
  (E=1.9e-6). blastp to the *W. ichthyophaga* EXF-994 locus on NW_008806275.1:
  BAP31 87.5%, CAF1 96.4%, HMG 51.2%.
- SXI1: the homeodomain gene XP_009266116.1 on the same *W. ichthyophaga* contig
  hits CBS 633.66 scaffold_2 at 41.4-42.2 kb, where no CDS is annotated.
  miniprot: 2 exons, 48.7% identity over all 263 aa, GT-AG intron; extended in
  frame to the first stop (TAA ending at 42,263): 278 aa. Curated with
  `protein_accession: null`; `build-gff` translates it from the exons.
- The other RefSeq homeodomain protein (XP_006956117.1, scaffold_1) sits in a
  conserved block (topo II, laccase, RRS1) that is also away from the locus in
  *W. ichthyophaga*; it is not SXI1.
- Validation: genes 0-3 match their accessions at 100%; SXI1 has no accession to
  check. Tests: 738 pass. Logged as B5.4.

### Detect on all 51 Wallemiales genomes (run-a6770ca, job 29118717)

- Called: **51 / 51** (0 before). Routing `lineage` for all.
- v1-type genomes (33): all `high`, `mat_locus`, all five genes modelled.
- Other-version genomes (18): all `medium`, `mat_locus`, on the divergent STE3
  (34.7-37.9%, modelled) plus BAP31/CAF1 (87-100%); HMG weak (33-60%) or
  unmodelled; no SXI1.
- Idiomorph is `undetermined` for all (pattern vocabulary, by design).
- Wall: the whole job took 1.5 min.

### Check 1 (synteny) -- recorded, not a gate

Independent of detect: tblastn and miniprot of the five record proteins per
genome (`wallemia/checks.py`, `checks_miniprot.py`, `versions.tsv`).
- v1-type genomes: SXI1, HMG and STE3 on one contig within 50 kb in **33/33**
  (span about 20-21 kb); all five record genes on one contig within 60 kb in
  **32/33**.
- Other version, *W. ichthyophaga* (4) and *W. canadensis* (1): STE3 and HMG on
  one contig 3.4-3.8 kb apart; no SXI1.
- Other version, *W. mellicola* (13): only STE3 is found by tblastn (29%); HMG
  and SXI1 are not found (detect finds a weak 33-37% HMG-like hit, unmodelled).

### Check 2 (a mix of versions) -- recorded, not a gate

Version = SXI1 modelled by miniprot (v1, the CBS 633.66 type) or not (other).

| species | n | v1 | other | STE3 identity to CBS 633.66 allele |
|---|---:|---:|---:|---|
| *W. mellicola* | 27 | 14 | 13 | v1 100.0%; other 28.7-29.3% (tblastn HSP) |
| *W. ichthyophaga* | 22 | 18 | 4 | v1 64.4-65.0%; other 33.1% |
| *W. hederae* | 1 | 1 | 0 | 70.6% |
| *W. canadensis* | 1 | 0 | 1 | 30.3% |

- The 4 *W. ichthyophaga* other-version genomes are EXF-759, EXF-3555, EXF-8622
  and EXF-8623: exactly the four "inverted" strains of Gostincar et al. 2019.
- *W. mellicola* splits 14:13, as Sun et al. 2019 reported (about half each).
- Within each species the STE3 of the two versions is highly divergent (about
  29-33% vs 65-100%), as expected for idiomorph-specific receptors.

### Limits

- The locus is still putative: no mating, meiosis or cross has been observed in
  *Wallemia*. Two versions at near 1:1 in *W. mellicola* fit a heterothallic
  locus but do not prove it.
- The record carries only the v1 version. The other version is called through
  a 35-38% STE3 and conserved flanks; a second record from an other-version
  genome (e.g. *W. mellicola* EXF-1262) would let detect label versions. Not
  built (not asked).
- The v1 SXI1 is this project's model, not a deposited gene.
