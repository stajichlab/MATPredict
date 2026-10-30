"""Write the two MATpk proposals (db/candidates) from the committed GenBank files.

MATa  allele: CBS 573 chromosome 2, CP028774.1 (= RefSeq NC_042507.1, identical sequence; RefSeq gene models and
              notes, window 305,000-320,000 in cbs573_NC_042507.1_305000-320000.gb).
MATalpha allele: MH260578.1 (the same strain, BioSample SAMN08558011).
Run from the repo root: python3 results/2026-10-29_pichia_family/make_records.py
"""
import re, yaml
from Bio import SeqIO

HERE = "results/2026-10-29_pichia_family/"
OFFSET = 304999   # window index 0 -> genome position 305,000 (1-based = offset + 1-based index in window)

def exons_of(loc, off, plus_first=True):
    parts = sorted(((int(p.start) + 1 + off, int(p.end) + off) for p in loc.parts), key=lambda t: t[0])
    strand = "+" if loc.strand == 1 else "-"
    if strand == "-":
        parts = parts[::-1]          # transcript order, descending
    return strand, [{"start": s, "end": e} for s, e in parts]

def gene(idx, name, prot, role, tag, strand, exons, order):
    lo = min(x["start"] for x in exons); hi = max(x["end"] for x in exons)
    g = {"gene_index": idx, "name": name, "protein_accession": f"ncbi_protein:{prot}", "role": role, "present": True,
         "locus_tag": tag, "segment_index": 0, "start": lo, "end": hi, "strand": strand, "order_in_locus": order}
    g["exons"] = exons
    return g

def cds_by_tag(rec, off):
    out = {}
    for f in rec.features:
        if f.type == "CDS":
            q = f.qualifiers
            out[q.get("locus_tag", [q.get("note", [""])[0]])[0]] = (f, off)
    return out

# ---- MATa allele (RefSeq window)
rs = next(SeqIO.parse(HERE + "cbs573_NC_042507.1_305000-320000.gb", "genbank"))
a = cds_by_tag(rs, OFFSET)
def pick(c, tag, name, role, order, idx):
    f, off = c[tag]
    strand, ex = exons_of(f.location, off)
    return gene(idx, name, f.qualifiers["protein_id"][0], role, tag, strand, ex, order)
genes_a = [pick(a, "C5L36_0B01570", "sla2", "flanking_conserved", 0, 0),
           pick(a, "C5L36_0B01580", "MATa1", "core_MAT", 1, 1),
           pick(a, "C5L36_0B01585", "MATa2", "core_MAT", 2, 2),
           pick(a, "C5L36_0B01590", "tgl1", "flanking_conserved", 3, 3)]
# ---- MATalpha allele (MH260578.1; coordinates are the deposit's own)
mh = next(SeqIO.parse(HERE + "MH260578.gb", "genbank"))
m = {}
for f in mh.features:
    if f.type == "CDS":
        m[f.qualifiers["gene"][0]] = f
def pickm(g, name, role, order, idx):
    f = m[g]; strand, ex = exons_of(f.location, 0)
    return gene(idx, name, f.qualifiers["protein_id"][0], role, f.qualifiers.get("note", [None])[0], strand, ex, order)
genes_m = [pickm("SLA2", "sla2", "flanking_conserved", 0, 0), pickm("MATalpha1", "MATalpha1", "core_MAT", 1, 1),
           pickm("MATalpha2", "MATalpha2", "core_MAT", 2, 2), pickm("TGL1", "tgl1", "flanking_conserved", 3, 3)]

COMMON_TAX = {"taxid": 4909, "lineage": "k__Fungi;p__Ascomycota;subphylum__Saccharomycotina;c__Pichiomycetes;o__Pichiales;f__Pichiaceae;g__Pichia;s__Pichia_kudriavzevii",
              "lineage_resolved_date": "2026-10-29"}
ORG = {"species": "Pichia kudriavzevii", "strain": {"name": "CBS 573", "known": True, "culture_collection_ids": ["CBS 573"], "differs_from_sequenced": False}}
CIT = [{"pmid": "30024981", "doi": "10.1371/journal.ppat.1007138"}]
def record(rid, idiomorph, genes, accession, orient, note, anchor, ev_note):
    lo = min(g["start"] for g in genes); hi = max(g["end"] for g in genes)
    loc = {"coordinate_provenance": "curator_derived" if accession.startswith("CP") else "published_explicit",
           "excluded_from_coordinate_benchmark": False,
           "core": {"completeness": "complete", "reference_orientation": orient, "definition_note": note,
                    "segments": [{"segment_index": 0, "sequence_source": {"type": "insdc_nucleotide", "accession": accession, "seq_region": accession},
                                  "start": lo, "end": hi, "contig_edge_distance": None, "sequence_checksum": None}]},
           "extended_flank": []}
    loc.update(anchor)
    return {"record_id": rid, "record_version": 1, "taxonomy": COMMON_TAX, "organism": ORG,
            "mating_type": {"locus_name": "MATpk", "idiomorphs": [idiomorph], "system": "heterothallic"},
            "locus": loc, "genes": genes,
            "evidence": {"locus_existence": {"tier": 2, "experimental_method": ev_note, "citations": CIT},
                         "boundaries": {"tier": 2, "experimental_method": "genome annotation (RefSeq or the depositor's); gene extents are annotation-derived", "citations": CIT},
                         "idiomorph_assignment": {"tier": 2, "experimental_method": "gene content and homology: MATa1 and MATa2 with the MATa allele, MATalpha1 (alpha box, MAT1-1-1 profile) and MATalpha2 with the MATalpha allele, in a heterozygous strain", "citations": CIT}},
            "validation": {"status": "needs_review", "rejection_reason": None, "accession_resolved": False, "accession_resolved_date": None,
                           "accession_resolved_version": None, "sequence_match": {"status": "pass", "per_gene": [], "notes": ""}, "taxonomy_current": None},
            "curation": {"proposed_by": "claude-code-curation-agent", "proposal_dedupe_key": f"30024981|4909|MATpk|{idiomorph}", "reviewed_by": None, "reviewed_date": None},
            "model_provenance": None}

NOTE_A = ("MATa allele of P. kudriavzevii CBS 573 on chromosome 2 (CP028774.1; RefSeq NC_042507.1, the same sequence), about 310 kb from the "
          "left end (Douglass et al. 2018, PMID 30024981). The species is heterothallic with no silent loci; CBS 573 is a heterozygote, and the other allele "
          "is the record 4909_cbs-573_MATpk_MATalpha (MH260578.1). GENE MODELS are the RefSeq annotation of GCF_003054445.1: SLA2 C5L36_0B01570 "
          "(XP_029320369.1), MATa1 C5L36_0B01580 (XP_029320370.1, 810 bp, single exon, minus strand), MATa2 C5L36_0B01585 (XP_029320371.1, plus strand, "
          "two exons 311,575-311,841 and 311,904-312,245, 202 aa) and TGL1 C5L36_0B01590 (XP_029320372.1). The RefSeq /note fields name MATa1 and MATa2 and "
          "a 'mating type locus Ya region' feature at 310,398-312,401. RNA-SEQ CHECK (SRR6966973, the same strain, 72.1 M reverse-stranded reads, 2026-10-29): "
          "the 62 bp intron at 311,842-311,903 is GT..AG with 27 spliced reads; the plus-strand transcript starts about 311,660, so the annotated start codon "
          "(ATG 311,575) is supported by homology and annotation, not by this RNA; ATG 311,788 would remove the HMG-box core. The MATa2 protein is divergent "
          "(Saccharomycotina MATa2 profile 52.4 bits, supporting evidence only). Flank synteny with the alpha allele: SLA2 and TGL1 (100% and 98% identical "
          "proteins by miniprot). SUI1 (306,404-306,947) and SEC3 (313,903-317,976) lie one gene further out and are not recorded. GENETIC CODE: standard (table 1). "
          "Details: results/2026-10-29_yeast_routing/ogataea_pichia_reference_check.md and results/2026-10-29_pichia_family/README.md.")
NOTE_M = ("MATalpha allele of P. kudriavzevii CBS 573, the deposit MH260578.1 (11,807 bp; SUI1, SLA2, MATalpha1, MATalpha2, TGL1, SEC3; Douglass et al. 2018, "
          "PMID 30024981). The deposit is the PacBio contig for the alternative allele of the heterozygous CBS 573; the chromosome-2 assembly carries the MATa allele "
          "(record 4909_cbs-573_MATpk_MATa). Same BioSample as the assembly (SAMN08558011, strain CBS 573). The record covers SLA2 to TGL1; SUI1 (5-prime end of "
          "the deposit, partial) and SEC3 (partial) are not recorded. MATalpha1 (180 aa) is single-exon on the minus strand; MATalpha2 has three exons on the minus "
          "strand (the last exon is 14 bp, as deposited). The MAT1-1-1 profile hits MATalpha1 (E 3.2e-7). GENETIC CODE: standard (table 1). No assembly anchor "
          "(assembly_accession null): the 11.8 kb contig is not part of GCA_003054445.1; the alpha allele of the independent isolate SJP "
          "(GCA_003033855.1, CP021089.1) is the test case in results/2026-10-29_pichia_family/README.md.")
ancho_a = {"assembly_accession": "GCA_003054445.1", "assembly_accession_basis": "strain_match_verified",
           "assembly_location": {"contig": "CP028774.1", "start": min(g["start"] for g in genes_a), "end": max(g["end"] for g in genes_a), "identity": 100.0}}
ra = record("4909_cbs-573_MATpk_MATa", "MATa", genes_a, "CP028774.1", "SLA2->MATa1->MATa2->TGL1", NOTE_A, ancho_a,
            "genome sequencing and annotation of the heterozygous CBS 573 (PMID 30024981); RNA-seq splice check of the MATa2 intron")
rm_ = record("4909_cbs-573_MATpk_MATalpha", "MATalpha", genes_m, "MH260578.1", "SLA2->MATalpha1->MATalpha2->TGL1", NOTE_M, {"assembly_accession": None},
             "targeted sequencing and annotation of the alpha allele of the heterozygous CBS 573 (MH260578.1; PMID 30024981)")
for r in (ra, rm_):
    open(f"/tmp/{r['record_id']}.yaml", "w").write(yaml.safe_dump(r, sort_keys=False))
    print(r["record_id"], [(g["name"], g["start"], g["end"], g["strand"], len(g["exons"])) for g in r["genes"]])
