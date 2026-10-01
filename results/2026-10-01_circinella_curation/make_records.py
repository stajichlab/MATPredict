"""Write the two Circinella-group tier-2 record metadata.yaml files from
verify_sources.tsv (exon coordinates verified there at 100% identity).

Usage: python make_records.py verify_sources.tsv DB_ROOT
"""
import csv, sys
from pathlib import Path
import yaml

rows = {r["protein"]: r for r in csv.DictReader(open(sys.argv[1]), delimiter="\t")}
db = Path(sys.argv[2])
NOSRCH = "TIER 2 (genome-derived). No Circinella, Thamnostylum, Fennellomyces or Zychaea MAT locus has a locus-specific deposit or a publication (NCBI nuccore/protein and PubMed searched 2026-10-01)."


def gene(i, name, role, prot, order):
    r = rows[prot]
    ex = [tuple(map(int, e.split("-"))) for e in r["exons"].split(",")]
    if r["strand"] == "-":
        ex = ex[::-1]  # transcript order, as the existing records
    return {"gene_index": i, "name": name, "protein_accession": f"ncbi_protein:{prot}", "role": role,
            "present": True, "locus_tag": None, "segment_index": 0, "start": int(r["start"]),
            "end": int(r["end"]), "strand": r["strand"], "order_in_locus": order,
            "exons": [{"start": s, "end": e} for s, e in ex],
            "codon_start": int(r["phase0"]) + 1, "transl_table": int(r["transl_table"])}


def absent(i, name, role):
    return {"gene_index": i, "name": name, "protein_accession": None, "role": role, "present": False,
            "locus_tag": None, "segment_index": None, "start": None, "end": None, "strand": None,
            "order_in_locus": None}


SPECS = [
    dict(record_id="64656_rsa-1403_MAT_Minus", taxid=64656, species="Zychaea mexicana",
         genus="Zychaea", strain="RSA 1403", ccids=["RSA 1403"], idiomorph="Minus", sex="sexM",
         asm="GCF_025766255.1", contig="NW_026516701.1",
         genes=[("VMA1", "flanking_conserved", "XP_052979470.1"),
                ("pntA", "flanking_conserved", "XP_052979471.1"),
                ("sexM", "core_MAT", "XP_052979473.1"),
                ("gpmI", "flanking_conserved", "XP_052979474.1"),
                ("rnhA", "flanking_conserved", "XP_052979475.1")],
         note=("Zychaea mexicana RSA 1403 (RefSeq GCF_025766255.1 = GenBank GCA_025766255.1, JGI "
               "Zycmex1; BioProject PRJNA331829, BioSample SAMN05444742; contig N50 307.7 kb) MAT Minus "
               "locus on NW_026516701.1 (identical to GenBank JAIXMT010000038.1; the RefSeq copy is used "
               "because it is the sequence the BFD genome library carries). Gene order: VMA1 "
               "(XP_052979470.1, -) -- pntA (XP_052979471.1, -) -- [XP_052979472.1, uncharacterized, -, "
               "not recorded] -- sexM (XP_052979473.1, +, 191 aa, single exon) -- gpmI (XP_052979474.1, -) "
               "-- rnhA (XP_052979475.1, +). GenBank copies of the proteins: KAI9493205.1-KAI9493210.1 "
               "(byte-identical). sexM carries a PF00505 HMG box (aa 46-113, E 5.0e-26); best curated "
               "match 4837_nrrl1555 sexM (blastp 43.4% over 203 aa); HMM classifier (a8863f3 build) sexM "
               "132.9 vs sexP 70.7 bits (margin +62.2). sexM ends 112 bp from gpmI and 2.1 kb from rnhA. "
               "rnhA is 54.9% identical to 29922_cbs109-16 rnhA. tptA (XP_052979568.1, NW_026516700.1), "
               "algA (XP_052976134.1, NW_026516743.1) and glrA (XP_052972719.1, NW_026516849.1) are on "
               "other contigs, so they are recorded present: false. VMA1, pntA and gpmI are recorded in "
               "locus.extended_flank, not as roster genes, so detection never searches them. Every protein matches the "
               "translation of its annotated exons from NW_026516701.1 at 100% (table 1, codon_start 1). "
               "TIER 2 (genome-derived). Not a published idiomorph.")),
    dict(record_id="101103_nrrl1351_MAT_Plus", taxid=101103, species="Circinella umbellata",
         genus="Circinella", strain="NRRL1351", ccids=["NRRL 1351"], idiomorph="Plus", sex="sexP",
         asm="GCA_025093555.1", contig="JAIWNF010000086.1",
         genes=[("VMA1", "flanking_conserved", "KAI7847719.1"),
                ("pntA", "flanking_conserved", "KAI7847720.1"),
                ("sexP", "core_MAT", "KAI7847721.1"),
                ("gpmI", "flanking_conserved", "KAI7847722.1"),
                ("rnhA", "flanking_conserved", "KAI7847723.1")],
         note=("Circinella umbellata NRRL 1351 (GCA_025093555.1, JGI Cirumb1; BioProject PRJNA331974, "
               "BioSample SAMN05444945; contig N50 688.3 kb) MAT Plus locus on JAIWNF010000086.1. Gene "
               "order: VMA1 (KAI7847719.1, -) -- pntA (KAI7847720.1, -) -- sexP (KAI7847721.1, +, 310 aa, "
               "single exon) -- gpmI (KAI7847722.1, -) -- rnhA (KAI7847723.1, +). KAI7847721.1 is "
               "annotated only as 'hypothetical protein' (ANNOTATION_ERRORS_FIXED_REPORT.md). sexP "
               "carries a PF00505 HMG box (aa 117-177, E 1.5e-15); best curated match 4837_ubc21 sexP "
               "(blastp 30.1% over 166 aa); HMM classifier (a8863f3 build) sexP 115.8 vs sexM 32.0 bits "
               "(margin +83.8). sexP ends 98 bp from gpmI and 2.1 kb from rnhA. rnhA is 55.8% identical "
               "to 29922_cbs109-16 rnhA. tptA (KAI7852361.1) and algA (KAI7852358.1) sit together on "
               "JAIWNF010000034.1 and glrA (KAI7860806.1) on JAIWNF010000001.1, so they are recorded "
               "present: false. VMA1, pntA and gpmI are recorded in locus.extended_flank, not as roster "
               "genes, so detection never searches them. The same neighbourhood (VMA1, pntA, HMG gene, gpmI, rnhA) is in the Minus "
               "record 64656_rsa-1403 (blastp VMA1 92.2%, pntA 77.3%, gpmI 87.5%, rnhA 77.6%). The Plus "
               "assignment of this sexP clade rests on the classifier and on the two Fennellomyces "
               "strains labelled Plus (results/2026-10-01_circinella_label_tree/); the HMG-box-only tree "
               "does not support the clade (IQ-TREE 0/19; full-length 83/93). Every protein matches the "
               "translation of its annotated exons from JAIWNF010000086.1 at 100% (table 1, codon_start "
               "1). TIER 2 (genome-derived). Not a published idiomorph.")),
]

#: Conserved neighbours recorded in `locus.extended_flank`, NOT as genes: a
#: roster gene (even unsearched) enters other Mucorales' scoring denominators.
NEIGHBOURS = {"VMA1": "V-type proton ATPase catalytic subunit A",
              "pntA": "P-loop containing nucleoside triphosphate hydrolase protein",
              "gpmI": "2,3-bisphosphoglycerate-independent phosphoglycerate mutase"}


def flank_entry(name, prot):
    r = rows[prot]
    ex = [tuple(map(int, e.split("-"))) for e in r["exons"].split(",")]
    if r["strand"] == "-":
        ex = ex[::-1]
    return {"name": name, "product": NEIGHBOURS[name], "protein_accession": f"ncbi_protein:{prot}",
            "searched": False, "contig": r["contig"], "start": int(r["start"]), "end": int(r["end"]),
            "strand": r["strand"], "exons": [{"start": a, "end": b} for a, b in ex],
            "codon_start": int(r["phase0"]) + 1, "transl_table": int(r["transl_table"]),
            "protein_match": "100% to the translation of its annotated exons"}


for sp in SPECS:
    genes = []
    present = sorted((g for g in sp["genes"] if g[0] not in NEIGHBOURS), key=lambda g: int(rows[g[2]]["start"]))
    flank = sorted((flank_entry(n, p) for n, _, p in sp["genes"] if n in NEIGHBOURS), key=lambda e: e["start"])
    for order, (name, role, prot) in enumerate(present):
        genes.append(gene(order, name, role, prot, order))
    n = len(genes)
    for k, (name, role) in enumerate([("tptA", "flanking_conserved"), ("algA", "flanking_variable"),
                                      ("glrA", "flanking_variable")]):
        genes.append(absent(n + k, name, role))
    starts = [g["start"] for g in genes if g["present"]]
    ends = [g["end"] for g in genes if g["present"]]
    orient = "|".join(f"{g['name']}{'->' if g['strand'] == '+' else '<-'}" for g in genes if g["present"])
    rec = {
        "record_id": sp["record_id"], "record_version": 1,
        "taxonomy": {"taxid": sp["taxid"],
                     "lineage": f"k__Fungi;p__Mucoromycota;c__Mucoromycetes;o__Mucorales;f__Lichtheimiaceae;g__{sp['genus']};s__{sp['species'].replace(' ', '_')}",
                     "lineage_resolved_date": "2026-10-01"},
        "organism": {"species": sp["species"], "strain": {"name": sp["strain"], "known": True,
                     "culture_collection_ids": sp["ccids"], "differs_from_sequenced": False}},
        "mating_type": {"locus_name": "MAT", "idiomorphs": [sp["idiomorph"]], "system": "heterothallic"},
        "locus": {
            "assembly_accession": sp["asm"],
            "coordinate_provenance": "curator_derived", "excluded_from_coordinate_benchmark": True,
            "core": {"completeness": "complete", "reference_orientation": orient,
                     "definition_note": sp["note"],
                     "segments": [{"segment_index": 0,
                                   "sequence_source": {"type": "insdc_nucleotide", "accession": sp["contig"],
                                                       "seq_region": sp["contig"]},
                                   "start": min(starts), "end": max(ends),
                                   "contig_edge_distance": None, "sequence_checksum": None}]},
            "extended_flank": flank},
        "genes": genes,
        "evidence": {
            "locus_existence": {"tier": 2, "citations": [], "experimental_method": NOSRCH + (
                " Locus identified by the HMG-box gene typed by the MATPredict HMM idiomorph classifier "
                "(db/Mucoromycota/classifiers/MAT) and by the Circinella-group neighbourhood "
                "VMA1-pntA-[sex]-gpmI-rnhA, conserved between Zychaea and Circinella.")},
            "boundaries": {"tier": 2, "citations": [], "experimental_method": NOSRCH + (
                f" Record segment spans the sex gene to rnhA on {sp['contig']} ({sp['asm']}); the "
                "idiomorph boundaries themselves are not mapped. VMA1, pntA and gpmI are in "
                "extended_flank (recorded, not searched).")},
            "idiomorph_assignment": {"tier": 2, "citations": [], "experimental_method": NOSRCH + (
                f" Idiomorph from the HMG-box gene ({sp['sex']}): classifier and PF00505 scores in "
                "definition_note; Plus = sexP, Minus = sexM, as the curated Mucorales records.")},
        },
        "validation": {"status": "accepted", "rejection_reason": None, "accession_resolved": True,
                       "accession_resolved_date": "2026-10-01", "accession_resolved_version": sp["contig"],
                       "sequence_match": {"status": "pass", "per_gene": [
                           {"gene_index": g["gene_index"], "percent_identity": 100.0, "coverage": 100.0,
                            "status": "pass"} for g in genes if g["present"]], "notes": (
                           "Deposited proteins compared to the translation of their annotated exons "
                           "(results/2026-10-01_circinella_curation/verify_sources.tsv).")},
                       "taxonomy_current": True},
        "curation": {"proposed_by": "claude-genome-curation",
                     "proposal_dedupe_key": f"{sp['asm']}|{sp['taxid']}|MAT|{sp['idiomorph']}",
                     "reviewed_by": None, "reviewed_date": None,
                     "notes": ("Proposed 2026-10-01 on curator ruling (Circinella-group curation, tier-2 "
                               "records for both idiomorphs). PENDING curator sign-off.")},
        "model_provenance": None,
    }
    d = db / "Mucoromycota" / "Mucorales" / sp["record_id"]
    d.mkdir(parents=True, exist_ok=True)
    (d / "metadata.yaml").write_text(yaml.safe_dump(rec, sort_keys=False, width=100))
    print(d)
