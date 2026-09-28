"""Build the tier-2 Umbelopsis MAT records (Plus and Minus), pending curator sign-off.

Curator ruling 2026-09-27: tier-2 Umbelopsis records for both idiomorphs so
Umbelopsis sexM/sexP can be modelled. No Umbelopsidales MAT deposit exists.

Plus  = U. ramanniana gzUmbRama1, GCA_977110945.1 (canu/PromethION, 23 contigs,
        N50 1.65 Mb), CDSBDH010000020.1. sexP = CAO3688886.1 (annotated).
Minus = U. vinacea WA0000051536, GCA_016758895.1 (SPAdes, 150 scaffolds,
        N50 1.33 Mb), JAEPRA010000016.1. sexM = 3-exon model built here; the
        assembly's annotated KAG2174494.1 is a 3'-part-only partial of it.
Usage: build_umbelopsis.py DB_ORDER_DIR   (…/db/Mucoromycota/Umbelopsidales)
"""
import os
import sys
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
DL = os.path.join(HERE, "dl")
out_root = sys.argv[1]


def cds_from_gff(asm, pid):
    """Exons (1-based, inclusive) and strand of an annotated protein's CDS."""
    gff = f"{DL}/{asm}/ncbi_dataset/data/{asm}/genomic.gff"
    parts, strand, contig, phase0 = [], None, None, None
    for line in open(gff):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS" or f"protein_id={pid}" not in f[8]:
            continue
        parts.append((int(f[3]), int(f[4]), int(f[7]) if f[7] in "012" else 0))
        strand, contig = f[6], f[0]
    assert parts, pid
    parts.sort(key=lambda p: -p[0] if strand == "-" else p[0])
    first_phase = parts[0][2]
    exons = [{"start": a, "end": b} for a, b, _ in parts]
    return contig, strand, exons, first_phase + 1


def gene(asm, pid, name, role):
    contig, strand, exons, codon_start = cds_from_gff(asm, pid)
    return contig, {"name": name, "protein_accession": "ncbi_protein:" + pid, "role": role, "present": True,
                    "locus_tag": None, "segment_index": 0,
                    "start": min(e["start"] for e in exons), "end": max(e["end"] for e in exons),
                    "strand": strand, "exons": exons, "codon_start": codon_start, "transl_table": 1}


def absent(name, role):
    return {"name": name, "protein_accession": None, "role": role, "present": False, "locus_tag": None,
            "segment_index": None, "start": None, "end": None, "strand": None}


KEYS = ("gene_index", "name", "protein_accession", "role", "present", "locus_tag", "segment_index",
        "start", "end", "strand", "order_in_locus", "exons", "codon_start", "transl_table")


def finish(genes_present, genes_absent):
    genes_present.sort(key=lambda g: g["start"])
    genes = genes_present + genes_absent
    out = []
    for i, g in enumerate(genes):
        g["gene_index"] = i
        g["order_in_locus"] = i if g["present"] else None
        out.append({k: g[k] for k in KEYS if k in g})
    return out, genes_present[0]["start"], genes_present[-1]["end"]


TIER2 = ("TIER 2 (genome-derived). No Umbelopsidales MAT locus has a locus-specific deposit "
         "or a publication (NCBI and PubMed searched 2026-09-27). ")
RNHA_NOTE = ("rnhA is recorded present: false: it does not sit at this locus in Umbelopsis "
             "(its ortholog is elsewhere, see definition_note), unlike the Mucorales "
             "algA-tptA-[sex]-rnhA-glrA arrangement. ")


def record(rid, taxid, lineage, species, strain, known, idiom, asm, contig, genes, s, e, definition,
           dedupe):
    return {
        "record_id": rid, "record_version": 1,
        "taxonomy": {"taxid": taxid, "lineage": lineage, "lineage_resolved_date": "2026-09-27"},
        "organism": {"species": species,
                     "strain": {"name": strain, "known": known,
                                "culture_collection_ids": [strain] if known else [],
                                "differs_from_sequenced": False}},
        "mating_type": {"locus_name": "MAT", "idiomorphs": [idiom], "system": "heterothallic"},
        "locus": {"coordinate_provenance": "curator_derived", "excluded_from_coordinate_benchmark": True,
                  "core": {"completeness": "complete",
                           "reference_orientation": "glrA<-|sexP<-|tptA->|algA->" if idiom == "Plus" else "glrA<-|sexM<-|tptA->|algA->",
                           "definition_note": definition,
                           "segments": [{"segment_index": 0,
                                         "sequence_source": {"type": "insdc_nucleotide", "accession": contig,
                                                             "seq_region": contig},
                                         "start": s, "end": e, "contig_edge_distance": None,
                                         "sequence_checksum": None}]},
                  "extended_flank": []},
        "genes": genes,
        "evidence": {
            "locus_existence": {"tier": 2, "citations": [], "experimental_method": TIER2 + (
                "Locus identified by synteny with the Mucorales MAT locus (tptA, algA, glrA at the "
                "locus; an HMG-box gene 1.5 kb from tptA) and by the MATPredict HMM idiomorph "
                "classifier (db/Mucoromycota/classifiers/MAT) on the HMG-box protein.")},
            "boundaries": {"tier": 2, "citations": [], "experimental_method": TIER2 + (
                f"Record segment spans the recorded genes on {contig} ({asm}).")},
            "idiomorph_assignment": {"tier": 2, "citations": [], "experimental_method": TIER2 + (
                "Idiomorph from the HMG-box gene: classifier and PF00505 scores in definition_note.")},
        },
        "validation": {"status": "needs_review", "rejection_reason": None, "accession_resolved": None,
                       "accession_resolved_date": None, "accession_resolved_version": None,
                       "sequence_match": None, "taxonomy_current": None},
        "curation": {"proposed_by": "claude-genome-curation", "proposal_dedupe_key": dedupe,
                     "reviewed_by": None, "reviewed_date": None,
                     "notes": ("Proposed 2026-09-27 on curator ruling (tier-2 Umbelopsis records for both "
                               "idiomorphs). PENDING CURATOR SIGN-OFF: accepted on the unmerged branch "
                               "curation-umbelopsis only so detect can be measured with it.")},
        "model_provenance": None,
    }


LIN = "k__Fungi;p__Mucoromycota;c__Umbelopsidomycetes;o__Umbelopsidales;f__Umbelopsidaceae;g__Umbelopsis;s__"

# ---- Plus: U. ramanniana gzUmbRama1 ----
A = "GCA_977110945.1"
gp = []
for pid, nm, role in [("CAO3688830.1", "glrA", "flanking_variable"), ("CAO3688886.1", "sexP", "core_MAT"),
                      ("CAO3688890.1", "tptA", "flanking_conserved"), ("CAO3688894.1", "algA", "flanking_variable")]:
    contig, g = gene(A, pid, nm, role)
    gp.append(g)
genes, s, e = finish(gp, [absent("rnhA", "flanking_conserved")])
plus = record("41833_gzumbrama1_MAT_Plus", 41833, LIN + "Umbelopsis_ramanniana", "Umbelopsis ramanniana",
              "gzUmbRama1", False, "Plus", A, contig, genes, s, e, (
    "Umbelopsis ramanniana gzUmbRama1 (GCA_977110945.1, canu v2.2 on PromethION, 23 contigs, contig N50 "
    "1.65 Mb; BioProject PRJEB96040, BioSample SAMEA118965535; the assembly name is used as the strain "
    "label because the BioSample gives no strain) MAT Plus locus on CDSBDH010000020.1. Gene order: glrA "
    "(CAO3688830.1, -) ... sexP (CAO3688886.1, -, 320 aa, single exon) -- tptA (CAO3688890.1, +) -- algA "
    "(CAO3688894.1, +). sexP carries a PF00505 HMG box (aa 108-172, E 6.8e-14); best curated match "
    "36080_nrrl-3631 sexP (blastp 25.5% over 153 aa); HMM classifier sexP 145.2 vs sexM 34.2 bits "
    "(margin +111). The detection scan modelled only an 88-aa sexP fragment here (29.6%); this record "
    "supplies the full-length protein. " + RNHA_NOTE +
    "rnhA ortholog (CAO3689817.1, 48.2% to 3108396 rnhA) lies at CDSBDH010000020.1:683,887-690,792, "
    "~513 kb away. glrA has several annotated isoforms (CAO3688830/34/38/42/46.1); the longest "
    "(CAO3688830.1, 465 aa) is recorded. " + TIER2),
    f"{A}|41833|MAT|Plus")

# ---- Minus: U. vinacea WA0000051536 ----
B = "GCA_016758895.1"
gm = []
for pid, nm, role in [("KAG2174489.1", "glrA", "flanking_variable"), ("KAG2174464.1", "tptA", "flanking_conserved"),
                      ("KAG2174465.1", "algA", "flanking_variable")]:
    contig, g = gene(B, pid, nm, role)
    gm.append(g)
sexM_exons = [{"start": 154336, "end": 154788}]
gm.append({"name": "sexM", "protein_accession": None, "role": "core_MAT", "present": True, "locus_tag": None,
           "segment_index": 0, "start": 154336, "end": 154788, "strand": "-", "exons": sexM_exons,
           "codon_start": 1, "transl_table": 1})
genes, s, e = finish(gm, [absent("rnhA", "flanking_conserved")])
minus = record("44442_wa0000051536_MAT_Minus", 44442, LIN + "Umbelopsis_vinacea", "Umbelopsis vinacea",
               "WA0000051536", True, "Minus", B, contig, genes, s, e, (
    "Umbelopsis vinacea WA0000051536 (GCA_016758895.1, SPAdes 3.10.5 on Illumina and PacBio RSII, 150 "
    "scaffolds, N50 1.33 Mb; BioProject PRJNA668042) MAT Minus locus on JAEPRA010000016.1. Gene order: "
    "glrA (KAG2174489.1, -) ... sexM (-) -- tptA (KAG2174464.1, +) -- algA (KAG2174465.1, +). sexM is not "
    "annotated; recorded here as the single-exon ORF 154,336-154,788 (minus strand, 150 aa, from the first "
    "Met after an upstream in-frame stop to the TAA stop). Its 3' end is uncertain: curated Mucorales "
    "sexM proteins are 185-249 aa, and no intron could be supported, so completeness is taken as the ORF. "
    "The HMG-box exon was found by tblastn of the curated sexM set (HSP 44.3 bits); PF00505 HMG box at "
    "aa 33-94. The ORF has no match in the Plus genome gzUmbRama1 proteome (idiomorph-specific). HMM "
    "classifier on the ORF: sexM 39.8 vs sexP 26.6 bits. The adjacent annotated KAG2174494.1 (131 aa, "
    "'partial') is NOT part of sexM: it is a conserved neighbouring gene (89% identical to CAO3688882.1 "
    "beside sexP in the Plus genome); a first draft of this record wrongly fused it to the HMG exon. "
    "The better-assembled U. vinacea gzUmbVina2 (GCA_977110975.1) was not used: its sexM region maps only "
    "with frameshifts and is not annotated. " + RNHA_NOTE +
    "rnhA ortholog (KAG2186781.1) is on JAEPRA010000004.1. " + TIER2),
    f"{B}|44442|MAT|Minus")

for rec in (plus, minus):
    d = os.path.join(out_root, rec["record_id"])
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "metadata.yaml"), "w") as fo:
        yaml.safe_dump(rec, fo, sort_keys=False, width=100)
    print("wrote", d, rec["locus"]["core"]["segments"][0]["start"], rec["locus"]["core"]["segments"][0]["end"])
