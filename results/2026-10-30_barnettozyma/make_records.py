"""Write the draft proposal for the B. californica PB4207 MAT locus (DRAFT, in this results directory only; nothing goes to db/).

One combined record (both idiomorphs at one locus; Hession et al. 2025 describe MATa1+MATa2 on one side and MATalpha1+MATalpha2 on the other
of a 2.8 kb intergenic region, and consider the species homothallic). Gene models are OURS (models.gff3, tier 3).
Run from this directory: python3 make_records.py
"""
import re, yaml, collections
CONTIG, ACC, ASM = "JAMAIE010000005.1", "JAMAIE010000005.1", "GCA_023627795.1"
exons, meta = collections.defaultdict(list), {}
for l in open("models.gff3"):
    if l.startswith("#"): continue
    c = l.rstrip("\n").split("\t")
    gid = re.search(r"(?:ID|Parent)=(BC_[^;.]+)", c[8]).group(1)
    if c[2] == "gene": meta[gid] = (c[6], int(c[3]), int(c[4]))
    if c[2] == "CDS": exons[gid].append((int(c[3]), int(c[4])))
# MATalpha2 (added 2026-10-30 from results/2026-10-30_mat_hmms): a partial alpha2 homeodomain ORF, 620,080-620,421 plus strand, 113 aa to the stop;
# the ORF begins right after an in-frame stop (620,077-620,079), so the N-terminal exon is not found. alpha2 profile 37 bits; same ORF in all four assemblies.
meta["BC_MATalpha2"] = ("+", 620080, 620421); exons["BC_MATalpha2"] = [(620080, 620421)]
ROLE = {"dic1": "flanking_conserved", "sla2": "flanking_conserved"}
IDIO = {"MATa1": ["MATa"], "MATa2": ["MATa"], "MATalpha1": ["MATalpha"], "MATalpha2": ["MATalpha"]}
GCLASS = {"MATalpha1": "alpha_box"}
order = sorted(meta, key=lambda g: meta[g][1])
genes = []
for i, gid in enumerate(order):
    name = gid[3:]; strand, lo, hi = meta[gid]
    ex = sorted(exons[gid]); ex = ex if strand == "+" else ex[::-1]
    g = {"gene_index": i, "name": name, "role": ROLE.get(name, "core_MAT"), "present": True, "segment_index": 0,
         "start": lo, "end": hi, "strand": strand, "order_in_locus": i, "exons": [{"start": a, "end": b} for a, b in ex]}
    if name in GCLASS: g["gene_class"] = GCLASS[name]
    genes.append(g)
lo = min(g["start"] for g in genes); hi = max(g["end"] for g in genes)

NOTE = ("DRAFT, GENE MODELS ARE OURS (tier 3). B. californica strain PB4207 (Dalian Medical University, GCA_023627795.1, ctg_5 = JAMAIE010000005.1, "
        "694,392 bp; the locus is 63 kb from one contig end). The one published description is Hession, Byrne, Wolfe & Butler, PLoS Genet 2025 "
        "(PMID 41359681, doi 10.1371/journal.pgen.1011814): a centromere formed at the MAT locus of B. californica, whose 460 bp centromere sequence lies in a "
        "2.8 kb intergenic region between MATa1 and MATa2 on one side and MATalpha1 and MATalpha2 on the other; all three Barnettozyma species appear homothallic. "
        "The paper gives no gene models, coordinates or accessions for the locus (its UCD09 ONT assembly is not in the NCBI assembly database as of 2026-10-30; "
        "the supplement S1 Table was not available). NO RNA-seq exists for the genus (ENA, 2026-10-30: only WGS and Hi-C), so every model is a prediction from "
        "sequence: stop-to-stop six-frame ORFs, Pfam (Homeodomain, MATalpha_HMGbox, ANTH/I_LWEQ), phmmer against the curated Saccharomycotina MAT proteins, an "
        "exhaustive one-intron GT..AG search scored with the curated Saccharomycotina MATa2 profile. WHAT AGREES WITH THE PAPER: MATa1 ends at 615,568 and MATalpha1 "
        "starts at 618,402, a 2,833 bp intergenic region (paper: 2.8 kb); MATa2 and MATa1 are the divergent pair on one side. SLA2 and DIC1 flank the locus "
        "(DIC1 7.4 kb from MATa2; SLA2 8.6 kb from MATalpha1). MODELS: SLA2 (1,046 aa, ANTH + I_LWEQ, single exon), DIC1 (316 aa, single exon), MATa1 (162 aa, "
        "homeodomain, single exon), MATalpha1 (233 aa, alpha-box, zero-intron ORF whose C-terminus is uncertain: one-intron alternatives score the same), MATa2 (269 aa, "
        "one 56 bp GTATGT..CTTAG intron at 614,112-614,167, 82 bits to the MATa2 profile against a background maximum of 32 over all 124 six-frame segments of the "
        "window; no Pfam HMG_box hit, as for Komagataella). MATalpha2 IS RECORDED AS PARTIAL (added 2026-10-30): a 113 aa alpha2-type homeodomain ORF at 620,080-620,421 (plus strand; the alpha2 profile of results/2026-10-30_mat_hmms scores 37 bits), the same ORF in all four californica assemblies; it begins right after an in-frame stop, so the N-terminal exon is not found and the model is a fragment. NAMING: the a1 and alpha2 profiles of results/2026-10-30_mat_hmms (supporting evidence, a score difference) call the MATa1 homeodomain a1 (104 bits against 23 for alpha2). The other three californica assemblies (UCD09 "
        "GCA_003123585.1, UCD89 GCA_003123885.1, NRRL Y-17395 GCA_003706655.3) carry the same proteins at 99.6-100% identity, which is one strain family, not independent "
        "evidence. B. discipulorum (GCA_049863505.1) has SLA2, MATa1, MATalpha1 and DIC1 at 45-86% identity but arranged differently and no MATa2 hit. "
        "GENETIC CODE: standard (table 1), not tested. Details: results/2026-10-30_barnettozyma/README.md.")
tax = {"taxid": 36038, "lineage": "k__Fungi;p__Ascomycota;subphylum__Saccharomycotina;c__Saccharomycetes;o__Phaffomycetales;f__Phaffomycetaceae;g__Barnettozyma;s__Barnettozyma_californica",
       "lineage_resolved_date": "2026-10-30"}
CIT = [{"pmid": "41359681", "doi": "10.1371/journal.pgen.1011814"}]
rec = {"record_id": "36038_pb4207_MATbz_MATa-MATalpha", "record_version": 1, "taxonomy": tax,
       "organism": {"species": "Barnettozyma californica", "strain": {"name": "PB4207", "known": True, "culture_collection_ids": ["PB4207"], "differs_from_sequenced": False}},
       "mating_type": {"locus_name": "MATbz", "idiomorphs": ["MATa", "MATalpha"], "system": "homothallic"},
       "locus": {"assembly_accession": ASM, "assembly_accession_basis": "strain_match_verified", "coordinate_provenance": "curator_derived",
                 "excluded_from_coordinate_benchmark": False,
                 "assembly_location": {"contig": CONTIG, "start": lo, "end": hi, "identity": 100.0},
                 "core": {"completeness": "partial", "reference_orientation": "DIC1->MATa2->MATa1->MATalpha1->MATalpha2 (partial)->SLA2", "definition_note": NOTE,
                          "segments": [{"segment_index": 0, "sequence_source": {"type": "insdc_nucleotide", "accession": ACC, "seq_region": ACC},
                                        "start": lo, "end": hi, "contig_edge_distance": None, "sequence_checksum": None}]},
                 "extended_flank": []},
       "genes": genes,
       "evidence": {"locus_existence": {"tier": 3, "experimental_method": "published description of the locus (PMID 41359681) plus our own homology/ORF analysis of the PB4207 assembly; no RNA-seq", "citations": CIT},
                    "boundaries": {"tier": 3, "experimental_method": "gene models predicted here from sequence (ORF, Pfam, phmmer, one-intron search); not an annotation of the deposit", "citations": []},
                    "idiomorph_assignment": {"tier": 3, "experimental_method": "gene content and position: HMG-like MATa2 + homeodomain MATa1 on one side, alpha-box MATalpha1 on the other, separated by 2.8 kb (PMID 41359681); the a1 profile of results/2026-10-30_mat_hmms calls the MATa1 homeodomain a1 (104 vs 23 bits); MATa2 scores 81-82 bits on the a2 profiles", "citations": CIT}},
       "validation": {"status": "needs_review", "rejection_reason": None, "accession_resolved": False, "accession_resolved_date": None,
                      "accession_resolved_version": None, "sequence_match": {"status": "pass", "per_gene": [], "notes": "models derived from the assembly itself"}, "taxonomy_current": None},
       "curation": {"proposed_by": "claude-code-curation-agent", "proposal_dedupe_key": "41359681|36038|MATbz|MATa-MATalpha|PB4207", "reviewed_by": None, "reviewed_date": None},
       "model_provenance": None}
open(f"{rec['record_id']}.yaml", "w").write(yaml.safe_dump(rec, sort_keys=False, width=100))
print(rec["record_id"], [(g["name"], g["start"], g["end"], g["strand"], len(g["exons"])) for g in genes])
