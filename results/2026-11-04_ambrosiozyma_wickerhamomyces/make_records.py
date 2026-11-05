"""Draft candidate records: MATam (Ambrosiozyma philentoma NRRL Y-7523) and MATwe (Wickerhamomyces edaphicus TBRC-BCC18952), tier 3.
Models: models.py (miniprot projections). Run from this folder with the repo's test python: python make_records.py"""
import yaml
from models import build
CIT = [{"pmid": "41849659", "doi": "10.1093/g3journal/jkag067"}]
LIN_AMB = "k__Fungi;p__Ascomycota;subphylum__Saccharomycotina;c__Pichiomycetes;o__Pichiales;f__Pichiaceae;g__Ambrosiozyma;s__Ambrosiozyma_philentoma"
LIN_WIC = "k__Fungi;p__Ascomycota;subphylum__Saccharomycotina;c__Saccharomycetes;o__Phaffomycetales;f__Wickerhamomycetaceae;g__Hansenula;s__Wickerhamomyces_edaphicus"
RECS = [
 dict(rid="43984_nrrl-y-7523_MATam_MATa-MATalpha", fam="MATam", taxid=43984, lin=LIN_AMB, species="Ambrosiozyma philentoma", strain="NRRL Y-7523", cc=["NRRL Y-7523"], asm="GCA_003707675.3",
      wins=["philentoma__PPKZ03000004.1", "philentoma__PPKZ03000012.1"],
      note=("CANDIDATE MATam record of Ambrosiozyma philentoma NRRL Y-7523 (GCA_003707675.3, contig-level assembly). TWO CASSETTE ENDS AT TWO UNRELATED POSITIONS: a MATa side on PPKZ03000004.1 "
            "(DIC1 - MATa2 - MATa1 - SLA2, with SUI1 beside SLA2; the Komagataella-type MATa arrangement) and a MATalpha side on PPKZ03000012.1, 211 kb from a contig edge (MATalpha2 and MATalpha1 on one strand 5 kb apart, "
            "and a 35 aa SLA2 3-prime fragment 2 kb beyond MATalpha1 that is not recorded). The two regions share no nucleotide sequence (blastn, 2026-11-04), so this is not an allele pair: a complementary pair of "
            "cassettes in one genome, as in the Ogataea two-cassette genomes; which one is expressed is not determined. In the same genus A. platypodis and A. monospora show the same two-locus pattern in the v0.8.0 run "
            "(fragments, not recorded). OUR OWN PROJECTION (tier 3): miniprot alignments of curated Saccharomycotina and Ogataea MAT proteins, kept only without frameshift or internal stop; no RNA-seq or annotation. "
            "GENE NAMES: MTLA1 of the Candida references is MATa1 and MTLA2 the HMG MATa2; the gene identities rest on homology at 28-70% identity. Every MAT model is partial (MATa1 58 aa homeodomain fragment, MATalpha2 107 aa). "
            "Independent review (Opus, 2026-11-04): the strongest second-cassette case of the 80 genomes of #189. GENETIC CODE: standard (NCBI code 1).")),
 dict(rid="640698_tbrc-bcc18952_MATwe_MATa-MATalpha", fam="MATwe", taxid=640698, lin=LIN_WIC, species="Wickerhamomyces edaphicus", strain="TBRC-BCC18952", cc=["TBRC-BCC18952"], asm="GCA_054949465.1",
      wins=["edaphicus__JBSXOT010000020.1", "edaphicus__JBSXOT010000274.1"],
      note=("CANDIDATE MATwe record of Wickerhamomyces edaphicus TBRC-BCC18952 (GCA_054949465.1; NCBI Taxonomy lists the species as '[Wickerhamomyces] edaphicus (nom. inval.)' under Hansenula). ONE HETEROZYGOUS "
            "MAT LOCUS WITH BOTH ALLELES ASSEMBLED AS SEPARATE CONTIGS: the MATalpha allele on JBSXOT010000020.1 (SUI1 - SLA2 - MATalpha1 - MATalpha2 - DIC1) and the MATa allele on the 3.5 kb contig JBSXOT010000274.1 "
            "(SLA2 3-prime end - MATa1 - MATa2 - DIC1 5-prime end). Both ends of contig 274 match contig 20 at 97.2% (533 bp, contig 20 41,292-41,824) and 97.3% (147 bp, 45,203-45,347) on either side of the "
            "idiomorph, so the two are the a and alpha alleles at the same position, NOT two cassettes (blastn, 2026-11-04; independent review, Opus). Only 680 bp of flank is shared; read mapping would confirm half coverage "
            "for each allele. OUR OWN PROJECTION (tier 3): miniprot alignments of curated Saccharomycotina and Ogataea MAT proteins, kept only without frameshift or internal stop; no RNA-seq or annotation. "
            "GENE NAMES: the MATa2 HMG gene (170 aa at 28% identity to the Pichia/Komagataella MATa2) and MATa1 (55 aa homeodomain fragment) rest on homology. GENETIC CODE: standard (NCBI code 1).")),
]
for R in RECS:
    genes, segs = [], []
    for si, w in enumerate(R["wins"]):
        ms = [(g, m) for g, m in build(w) if m]
        assert all(m["stops"] == 0 and not m["fs"] and m["mod3"] == 0 for g, m in ms)
        for g, m in ms:
            ex = sorted(m["exons"], reverse=(m["strand"] == "-"))
            genes.append({"name": g, "protein_accession": None, "role": m["role"], "present": True, "locus_tag": None, "segment_index": si,
                          "start": min(a for a, b in m["exons"]), "end": max(b for a, b in m["exons"]), "strand": m["strand"],
                          "exons": [{"start": a, "end": b} for a, b in ex], "_prot": m["prot"], "_contig": m["contig"]})
        gs = [x for x in genes if x["segment_index"] == si]
        segs.append({"segment_index": si, "sequence_source": {"type": "insdc_nucleotide", "accession": gs[0]["_contig"], "seq_region": gs[0]["_contig"]},
                     "start": min(x["start"] for x in gs), "end": max(x["end"] for x in gs), "contig_edge_distance": None, "sequence_checksum": None})
    genes.sort(key=lambda g: (g["segment_index"], g["start"]))
    for i, g in enumerate(genes): g["gene_index"] = i; g["order_in_locus"] = i
    out_genes = [{"gene_index": g["gene_index"], **{k: v for k, v in g.items() if k not in ("gene_index", "_prot", "_contig", "order_in_locus")}, "order_in_locus": g["order_in_locus"]} for g in genes]
    rec = {"record_id": R["rid"], "record_version": 1, "taxonomy": {"taxid": R["taxid"], "lineage": R["lin"], "lineage_resolved_date": "2026-11-04"},
           "organism": {"species": R["species"], "strain": {"name": R["strain"], "known": True, "culture_collection_ids": R["cc"], "differs_from_sequenced": False}},
           "mating_type": {"locus_name": R["fam"], "idiomorphs": ["MATa", "MATalpha"], "system": "unknown"},
           "locus": {"assembly_accession": R["asm"], "assembly_accession_basis": "strain_match_verified", "coordinate_provenance": "curator_derived", "excluded_from_coordinate_benchmark": False,
                     "core": {"completeness": "partial", "reference_orientation": "->".join(g["name"] for g in genes), "definition_note": R["note"], "segments": segs}, "extended_flank": []},
           "genes": out_genes,
           "evidence": {"locus_existence": {"tier": 3, "experimental_method": "homology to curated SLA2-type MAT cassettes (Ogataea, Komagataella; SLA2-MAT architecture, Kriti et al. 2026); the v0.8.0 detection called both loci; no study of this species", "citations": CIT},
                        "boundaries": {"tier": 3, "experimental_method": "OUR OWN miniprot projection of curated MAT and flank proteins onto the assembly contigs; no transcript support", "citations": CIT},
                        "idiomorph_assignment": {"tier": 3, "experimental_method": "gene identity by best-scoring curated reference (MATa1 and MATa2 versus MATalpha1 and MATalpha2); the idiomorph follows the gene content", "citations": CIT}},
           "validation": {"status": "needs_review", "rejection_reason": None, "accession_resolved": False, "accession_resolved_date": None, "accession_resolved_version": None,
                          "sequence_match": {"status": "pass", "per_gene": [], "notes": ""}, "taxonomy_current": None},
           "curation": {"proposed_by": "claude-code-curation-agent", "proposal_dedupe_key": f"41849659|{R['taxid']}|{R['fam']}|MATa-MATalpha", "reviewed_by": None, "reviewed_date": None},
           "model_provenance": None}
    import os; os.makedirs("records", exist_ok=True)
    yaml.safe_dump(rec, open(f"records/{R['rid']}.yaml", "w"), sort_keys=False, width=110)
    with open(f"records/{R['rid']}.proteins.faa", "w") as fh:
        for g in genes: fh.write(f">{R['rid']}|gene_index={g['gene_index']}|name={g['name']}|role={g['role']}\n{g['_prot']}\n")
    print(R["rid"], [(g["name"], g["segment_index"], g["start"], g["end"], g["strand"]) for g in genes])
