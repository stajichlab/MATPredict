"""Build the tier-2 Syncephalastrum racemosum NRRL 2496 MAT Plus record (pending sign-off).

Curator ruling 2026-09-27 (revised rule for Lichtheimiaceae/Syncephalastraceae):
build a tier-2 record where the HMM idiomorph classifier assigns the core gene
with margin > 25 bits AND the locus shows the Mucorales gene order.

sexP is NOT annotated in the JGI assembly (GCA_002105135.1, Synrac1); it is the
single-exon ORF MCGN01000004.1:1,755,515-1,756,447 (+), 310 aa, 69 bp upstream
of rnhA (ORY97819.1, +). The same HMG protein is 100% identical over 228 aa in
S. racemosum B6101 (GCA_000696955.1), where detection called it Plus beside rnhA.
Usage: build_synrac.py DB_ORDER_DIR   (.../db/Mucoromycota/Mucorales)
"""
import os, sys, yaml

S = "/scratch/jstajich/28984668/claude-1181/-bigdata-stajichlab-jstajich-projects-MATPredict/48b08d1a-4f03-410b-bf7c-431cd980f588/scratchpad/lichsync"
out_root = sys.argv[1]
A, C = "GCA_002105135.1", "MCGN01000004.1"
rnhA_exons = [(1756516, 1756535), (1756587, 1757253), (1757301, 1757612), (1757660, 1757687),
              (1757737, 1758015), (1758071, 1758181), (1758231, 1759411), (1759462, 1759589),
              (1759638, 1760241), (1760296, 1760526), (1760578, 1760805)]
TIER2 = ("TIER 2 (genome-derived). No Syncephalastraceae MAT locus has a locus-specific deposit or a "
         "publication (NCBI and PubMed searched 2026-09-26/27). ")
genes = [
    {"gene_index": 0, "name": "sexP", "protein_accession": None, "role": "core_MAT", "present": True,
     "locus_tag": None, "segment_index": 0, "start": 1755515, "end": 1756447, "strand": "+",
     "order_in_locus": 0, "exons": [{"start": 1755515, "end": 1756447}], "codon_start": 1, "transl_table": 1},
    {"gene_index": 1, "name": "rnhA", "protein_accession": "ncbi_protein:ORY97819.1", "role": "flanking_conserved",
     "present": True, "locus_tag": None, "segment_index": 0, "start": 1756516, "end": 1760805, "strand": "+",
     "order_in_locus": 1, "exons": [{"start": a, "end": b} for a, b in rnhA_exons], "codon_start": 1,
     "transl_table": 1},
]
for i, (nm, role) in enumerate([("tptA", "flanking_conserved"), ("algA", "flanking_variable"),
                                 ("glrA", "flanking_variable"), ("btbA", "flanking_variable")], start=2):
    genes.append({"gene_index": i, "name": nm, "protein_accession": None, "role": role, "present": False,
                  "locus_tag": None, "segment_index": None, "start": None, "end": None, "strand": None,
                  "order_in_locus": None})
definition = (
    "Syncephalastrum racemosum NRRL 2496 (GCA_002105135.1, JGI Synrac1; BioProject PRJNA330704) MAT Plus "
    "locus on MCGN01000004.1. Gene order: sexP (+) -- 69 bp -- rnhA (ORY97819.1, +, 1262 aa, 11 exons), the "
    "canonical Mucorales sexP->rnhA adjacency and orientation (cf. tptA->btbA->sexP->rnhA in R. arrhizus "
    "CBS 346-36). tptA, algA, glrA and btbA are NOT at this locus: the genome's best tptA (ORZ02716.1) and "
    "glrA (ORZ02423.1) lie on MCGN01000001.1, algA hits are weak (ORY97962.1, 32.9%); recorded present: "
    "false. So the Mucorales order is PARTIAL here (sexP-rnhA only). sexP is not annotated in the JGI "
    "assembly (annotation gap); recorded as the single-exon ORF 1,755,515-1,756,447 (+), 310 aa, from the "
    "first Met after an upstream in-frame stop to the stop codon (stop-to-stop 334 aa). PF00505 HMG box aa "
    "112-175 (E 1.2e-16); best tblastn to curated sexP: 4837_ubc21 sexP (35.2% over 165 aa, 92 bits). HMM "
    "classifier (db/Mucoromycota/classifiers/MAT) on the ORF: sexP 179.1 vs sexM 38.9 bits (margin +140). "
    "The same HMG protein is 100% identical (228 aa) in S. racemosum B6101 (GCA_000696955.1, "
    "JNDN01000979.1:849,656), where MATPredict called Plus beside rnhA. Final ML tree "
    "(results/2026-09-27_sexMP_final_ml): the called S. racemosum loci sit outside the sexP clade (other "
    "HMG; the 69-column HMG box gives UFBoot 62 even for the sexP references). " + TIER2)
rec = {
    "record_id": "13706_nrrl-2496_MAT_Plus", "record_version": 1,
    "taxonomy": {"taxid": 13706, "lineage": "k__Fungi;p__Mucoromycota;subphylum__Mucoromycotina;c__Mucoromycetes;"
                 "o__Mucorales;f__Syncephalastraceae;g__Syncephalastrum;s__Syncephalastrum_racemosum",
                 "lineage_resolved_date": "2026-09-27"},
    "organism": {"species": "Syncephalastrum racemosum",
                 "strain": {"name": "NRRL 2496", "known": True, "culture_collection_ids": ["NRRL 2496"],
                            "differs_from_sequenced": False}},
    "mating_type": {"locus_name": "MAT", "idiomorphs": ["Plus"], "system": "heterothallic"},
    "locus": {"coordinate_provenance": "curator_derived", "excluded_from_coordinate_benchmark": True,
              "core": {"completeness": "partial", "reference_orientation": "sexP->|rnhA->",
                       "definition_note": definition,
                       "segments": [{"segment_index": 0,
                                     "sequence_source": {"type": "insdc_nucleotide", "accession": C, "seq_region": C},
                                     "start": 1755515, "end": 1760805, "contig_edge_distance": None,
                                     "sequence_checksum": None}]},
              "extended_flank": []},
    "genes": genes,
    "evidence": {
        "locus_existence": {"tier": 2, "citations": [], "experimental_method": TIER2 + (
            "Locus identified by the MATPredict HMM idiomorph classifier on the HMG-box protein (margin +140 "
            "bits) and by the Mucorales sexP-rnhA adjacency (partial gene order; see definition_note).")},
        "boundaries": {"tier": 2, "citations": [], "experimental_method": TIER2 + (
            f"Record segment spans sexP and rnhA on {C} ({A}).")},
        "idiomorph_assignment": {"tier": 2, "citations": [], "experimental_method": TIER2 + (
            "Idiomorph from the HMG-box gene: HMM classifier sexP 179.1 vs sexM 38.9 bits.")},
    },
    "validation": {"status": "needs_review", "rejection_reason": None, "accession_resolved": None,
                   "accession_resolved_date": None, "accession_resolved_version": None,
                   "sequence_match": None, "taxonomy_current": None},
    "curation": {"proposed_by": "claude-genome-curation", "proposal_dedupe_key": f"{A}|13706|MAT|Plus",
                 "reviewed_by": None, "reviewed_date": None,
                 "notes": ("Proposed 2026-09-27 on curator ruling (revised Lichtheimiaceae/Syncephalastraceae "
                           "rule: classifier margin > 25 bits and Mucorales gene order). PENDING CURATOR "
                           "SIGN-OFF: accepted on the unmerged branch curation-umbelopsis only so detect can "
                           "be measured with it.")},
    "model_provenance": None,
}
d = os.path.join(out_root, rec["record_id"])
os.makedirs(d, exist_ok=True)
with open(os.path.join(d, "metadata.yaml"), "w") as fo:
    yaml.safe_dump(rec, fo, sort_keys=False, width=100)
print("wrote", d)
