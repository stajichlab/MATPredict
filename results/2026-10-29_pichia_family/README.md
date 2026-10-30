# MATpk: the Pichia kudriavzevii family (#189)

Curator decision of 2026-10-29: add a Pichia family, with the CBS 573 MATa record and the MATalpha record MH260578, using the 202 aa MATa2.

## What was added

| Item | Where |
|---|---|
| Family `MATpk` (scope: *P. kudriavzevii*, taxid 4909; heterothallic; genes MATa1, MATa2, MATalpha1, MATalpha2, sla2, tgl1) | `db/Ascomycota/order.yml` |
| Record `4909_cbs-573_MATpk_MATa`: CP028774.1 (= RefSeq NC_042507.1) 307,183-313,684, assembly GCA_003054445.1 (BioSample match, 100% identity) | `db/Ascomycota/Pichiales/` |
| Record `4909_cbs-573_MATpk_MATalpha`: MH260578.1 780-7,515, no assembly anchor (the 11.8 kb contig is not in the assembly) | `db/Ascomycota/Pichiales/` |
| Label registry entry for `Ascomycota:MATpk` with the **proposed** system code `PicKud` | `scripts/idiomorph_label_registry.json` |
| Tests | `tests/detect/test_matpk_family.py`, `tests/detect/test_fallback_searchable.py` |

Records were made with `curate-db propose`, validated under Slurm (`validate.slurm`, taxonkit on the node; job 31408439) and accepted: all 8 genes re-translate from their recorded coordinates to the recorded proteins at 100% identity and coverage, `taxonomy_current` true. `make_records.py` builds the proposals from the committed GenBank files (`cbs573_NC_042507.1_305000-320000.gb`, `MH260578.gb`).

Gene models: sla2, MATa1, MATa2 (two exons, 311,575-311,841 and 311,904-312,245, 202 aa) and tgl1 are the RefSeq annotation of GCF_003054445.1, whose notes name MATa1 and MATa2; MATalpha1, MATalpha2 (three exons) and the flanks are the deposit's. SUI1 and SEC3, one gene further out, are not recorded (SUI1 and SEC3 are partial in MH260578). Evidence for the MATa2 model and the open N-terminus question: `../2026-10-29_yeast_routing/ogataea_pichia_reference_check.md`.

Why two alleles and not cassettes: the species is heterothallic with no silent loci (Douglass et al. 2018, PMID 30024981). A heterozygous diploid carries MATa on one chromosome and MATalpha on the other (or on a separate contig). So `cassette_report` is off for MATpk.

## Replay: the 158 Pichia genomes of the v0.7.0 Ascomycota run

`--preset v0.8` throughout; `before` = tree_v5 (no MATpk), `after` = the current tree. `replay_*_scope.txt` list every genome before and after; `cmp_replay.py` and `identity_by_species.py` make them.

**Scope = the genus Pichia (4919), first attempt** (job 31409451):

| | before | after |
|---|---|---|
| genomes called | 136 | 150 (none lost, 14 gained) |
| *P. kudriavzevii* (56 genomes) | 42 called (MATsc) | 54 called (MATpk); core identity median 100% |
| other Pichia (102 genomes) | 94 called, MATsc `medium` or `high` | 96 called, MATpk, mostly `high` at 31-60% core identity (median about 50%) |

Calling the other species `high` from a one-species reference overstates the evidence, so the scope was narrowed.

**Scope = *P. kudriavzevii* (4909), as shipped** (job 31410832):

- The 56 *P. kudriavzevii* genomes give exactly the same calls as with the genus in scope: 54 genomes with a MATpk call (42 called before; 84 high and 4 medium loci), 2 uncalled (GCA_001879685.1, GCA_029290715.1; not examined).
- The 102 other Pichia genomes are unchanged (the same 94 called with the same MATsc calls as before; the only difference is the flag the v0.8 policy adds).
- The four complete assemblies: GCF_003054445.1 (CBS 573) MATa 307,183-313,684, high (exactly the record span); GCA_003054405.1 (CBS 5147) MATa high; GCA_014873065.1 (CY902) MATa high; **GCA_003033855.1 (SJP, an independent isolate) MATalpha 305,672-312,398, high**, called `MATsc medium` before. The SJP alpha call is a test the records were not built from.
- Of the 54 called *P. kudriavzevii* genomes, 32 carry both idiomorphs (MATa and MATalpha loci, heterozygous diploids), 13 only MATa and 9 only MATalpha.

## For the curator

1. **Family code `MATpk`** and the registry label code **`PicKud`** are proposals (the grammar says genus3+species3 of the reference organism); `PicKud` is flagged PROPOSED in the registry.
2. **Scope.** Species (shipped) or genus. The genus needs a second *Pichia* record or an identity floor; the 101 other species matched at 31-60% core identity.
3. **MATa2 N-terminus.** The record carries ATG 311,575 (homology and RefSeq); the RNA transcript starts about 311,660. Recorded in the record's note.
4. **The MATalpha record has no assembly anchor.** `assembly_accession` is null, as for the first MATog record.
5. **Evidence tiers are 2** for all three claims (no functional assay in this species).
6. The two uncalled *P. kudriavzevii* genomes were not looked at.
