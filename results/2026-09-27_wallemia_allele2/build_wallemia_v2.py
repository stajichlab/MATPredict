"""Build the putative tier-2 Wallemia canadensis EXF-10342 MAT record (version v2).

Second Wallemia record (curator ruling 2026-09-27): a genome carrying the OTHER
locus version, so detect can name both versions. Coordinates, exons,
codon_start and protein accessions of BAP31, STE3v2 and CAF1 are read from the
contig's own GenBank CDS features (JBHFOR010000004.1). SXI1 and HMG are recorded
present: false (see definition_note). Usage: build_wallemia_v2.py OUTDIR
"""
import sys
import yaml
from Bio import SeqIO

HERE = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-27_wallemia_allele2/work"
out_dir = sys.argv[1]

gb = SeqIO.read(f"{HERE}/JBHFOR010000004.1.gbk", "genbank")
cds = {f.qualifiers["protein_id"][0]: f for f in gb.features
       if f.type == "CDS" and "protein_id" in f.qualifiers}

CIT = [
    {"doi": "10.1016/j.fgb.2012.01.007", "pmid": "22326418"},
    {"doi": "10.3390/genes10060427", "pmid": "31167502"},
    {"doi": "10.3389/fmicb.2019.02019", "pmid": "31551960"},
]


def annotated(pid, name, role):
    f = cds[pid]
    loc = f.location
    strand = "+" if loc.strand == 1 else "-"
    exons = [{"start": int(p.start) + 1, "end": int(p.end)} for p in loc.parts]
    exons.sort(key=lambda e: -e["start"] if strand == "-" else e["start"])
    q = f.qualifiers
    return {"name": name, "protein_accession": "ncbi_protein:" + pid, "role": role, "present": True,
            "locus_tag": q.get("locus_tag", [None])[0], "segment_index": 0,
            "start": int(loc.start) + 1, "end": int(loc.end), "strand": strand,
            "exons": exons, "codon_start": int(q.get("codon_start", ["1"])[0]),
            "transl_table": int(q.get("transl_table", ["1"])[0])}


def absent(name):
    return {"name": name, "protein_accession": None, "role": "core_MAT", "present": False,
            "locus_tag": None, "segment_index": None, "start": None, "end": None, "strand": None,
            "exons": None, "codon_start": None, "transl_table": None}


present = [
    annotated("KAN3001147.1", "STE3v2", "core_MAT"),
    annotated("KAN3001148.1", "CAF1", "flanking_conserved"),
    annotated("KAN3001158.1", "BAP31", "flanking_conserved"),
]
present.sort(key=lambda g: g["start"])
genes = present + [absent("HMG"), absent("SXI1")]
for i, g in enumerate(genes):
    g["gene_index"] = i
    g["order_in_locus"] = i if g["present"] else None
KEYS = ("gene_index", "name", "protein_accession", "role", "present", "locus_tag", "segment_index",
        "start", "end", "strand", "order_in_locus", "exons", "codon_start", "transl_table")
genes = [{k: g[k] for k in KEYS} for g in genes]
for g in genes:
    if not g["present"]:
        for k in ("exons", "codon_start", "transl_table"):
            g.pop(k)

start, end = present[0]["start"], present[-1]["end"]
PUT = ("PUTATIVE MAT LOCUS, TIER 2 (genome-derived). No MAT gene of Wallemia has a "
       "locus-specific GenBank deposit, and no mating or meiosis has been observed in the genus. ")
rec = {
    "record_id": "1708542_exf-10342_wallMAT_v2",
    "record_version": 1,
    "taxonomy": {"taxid": 1708542,
                 "lineage": "k__Fungi;p__Basidiomycota;sc__Wallemiomycotina;c__Wallemiomycetes;"
                            "o__Wallemiales;f__Wallemiaceae;g__Wallemia;s__Wallemia_canadensis",
                 "lineage_resolved_date": "2026-09-27"},
    "organism": {"species": "Wallemia canadensis",
                 "strain": {"name": "EXF-10342", "known": True, "culture_collection_ids": ["EXF-10342"],
                            "differs_from_sequenced": False}},
    "mating_type": {"locus_name": "wallMAT", "idiomorphs": ["v2"], "system": "heterothallic"},
    "locus": {
        "coordinate_provenance": "curator_derived",
        "excluded_from_coordinate_benchmark": True,
        "core": {
            "completeness": "complete",
            "reference_orientation": "STE3v2<-|CAF1<-|BAP31<-",
            "definition_note": (
                "PUTATIVE. Wallemia canadensis EXF-10342 (GCA_056320075.1, SPAdes, 159 contigs, "
                "contig N50 186 kb) putative mating-type locus, the second of the two versions "
                "described by Sun et al. 2019 (W. mellicola) and Gostincar et al. 2019 (W. "
                "ichthyophaga), on JBHFOR010000004.1. Chosen as the best-assembled genome of this "
                "version among the 18 BFD Wallemiales genomes that carry it, and the only one whose "
                "assembly annotates the locus receptor. STE3v2 = KAN3001147.1 (ACAZ38_004777, "
                "'Pheromone B beta 1 receptor'), CAF1 = KAN3001148.1 (ACAZ38_004780, 'CCR4-NOT "
                "transcription complex subunit 7'; 99.6% to the v1 CAF1), BAP31 = KAN3001158.1 "
                "(ACAZ38_004824, 'hypothetical protein'; 98.4% to the v1 BAP31). STE3v2 is 34.8% "
                "to the v1 receptor EIM23479.1 over 155 aa, and 80% (tblastn, 136 aa) to the "
                "receptor of this version in W. mellicola, 65% in W. ichthyophaga. The annotated "
                "STE3v2 (192 aa) is probably 3'-truncated: tblastn and exonerate of the v1 STE3 "
                "(366 aa) continue about 260 bp past its annotated stop (logged B5.5); the record "
                "keeps the annotated protein. HMG: recorded present: false. Only a 43-aa HMG-box "
                "fragment is found (exonerate of the v1 HMG, 60.5% identity, 331,565-331,693, "
                "unannotated), consistent with Gostincar et al. 2019 who report a truncated HMG "
                "gene in the inverted version; no intact ORF is asserted. SXI1: recorded present: "
                "false; Gostincar et al. 2019 report it absent from the inverted version and no "
                "homeodomain hit lies in the locus. The pheromone-processing gene STE14 "
                "(KAN3001145.1) and an NDUFA6-like gene (KAN3001146.1) lie next to STE3v2, as in "
                "v1, and are not curated. BAP31 sits on the far side of CAF1 here (about 20 kb), "
                "the reverse of v1, consistent with the published inversion. Idiomorph 'v2' is a "
                "placeholder for this version; which version is which mating type is not known. "
                "TIER 2: genome-derived."),
            "segments": [{"segment_index": 0,
                          "sequence_source": {"type": "insdc_nucleotide", "accession": "JBHFOR010000004.1",
                                              "seq_region": "JBHFOR010000004.1"},
                          "start": start, "end": end, "contig_edge_distance": None,
                          "sequence_checksum": None}],
        },
        "extended_flank": [],
    },
    "genes": genes,
    "evidence": {
        "locus_existence": {"tier": 2, "citations": CIT, "experimental_method": PUT + (
            "The two locus versions are described as putative by Sun et al. 2019 and Gostincar "
            "et al. 2019 from population genomes; this record is the version without SXI1.")},
        "boundaries": {"tier": 2, "citations": CIT, "experimental_method": PUT + (
            "Record segment spans STE3v2 to BAP31 on JBHFOR010000004.1; three genes from the "
            "assembly's own CDS features.")},
        "idiomorph_assignment": {"tier": 2, "citations": CIT, "experimental_method": PUT + (
            "'v2' is a placeholder for the second locus version; which version is which mating "
            "type is not known. mating_type.system 'heterothallic' is itself putative: it rests on "
            "two locus versions each carried by part of the sequenced strains.")},
    },
    "validation": {"status": "needs_review", "rejection_reason": None, "accession_resolved": None,
                   "accession_resolved_date": None, "accession_resolved_version": None,
                   "sequence_match": None, "taxonomy_current": None},
    "curation": {"proposed_by": "claude-literature-search",
                 "proposal_dedupe_key": "10.3389/fmicb.2019.02019|1708542|wallMAT|v2",
                 "reviewed_by": None, "reviewed_date": None,
                 "notes": ("Proposed 2026-09-27 on curator ruling (a second record from a genome "
                           "carrying the other version). PENDING CURATOR SIGN-OFF: accepted on the "
                           "unmerged branch curation-puccinio only so detect can be measured with it.")},
    "model_provenance": None,
}
path = f"{out_dir}/{rec['record_id']}.yaml"
open(path, "w").write(yaml.safe_dump(rec, sort_keys=False, width=100))
print("wrote", path, start, end)
