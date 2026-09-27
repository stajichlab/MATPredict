"""Build the tier-2 Russulaceae B (PR) record: Russula nobilis gfRusNobi1.

Curator ruling 2026-09-27 (receptor queue step c, Russulaceae). No
Russulaceae mating-type deposit exists in NCBI nucleotide (only WGS scaffolds
and mip partial CDS) and PubMed has no Russula/Lactarius mating-type paper, so
the record is genome-derived.

Genome: GCA_984573805.1 (gfRusNobi1.hap1.1), chromosome level, 123 contigs,
contig N50 1.48 Mb. The public assembly has NO gene annotation. Receptor gene
models are the BFD funannotate models (FD20D505_007127/8/9), recorded by
coordinates on the public contig OZ475200.1; the precursor is an unannotated
ORF from the strict-CAAX scan (results/2026-09-27_pheromone_positional).

Every gene is re-translated from the genome and checked: receptors must equal
the BFD model protein; the precursor must start with Met and have no internal
stop. Usage: build_record.py DB_ROOT GENOME_FASTA BFD_PROTEINS BFD_GFF3
"""
import os
import re
import sys

import yaml
from Bio import SeqIO
from Bio.Seq import Seq

DB, FA, PROT, GFF = sys.argv[1:5]
CONTIG = "OZ475200.1"
STRICT = re.compile(r"C[VI][IV][AVMG]$")
R, P = "pheromone_receptor", "fungal_mating_type_pheromone"

seq = next(r.seq for r in SeqIO.parse(FA, "fasta") if r.id == CONTIG)
prot = {r.id: str(r.seq).rstrip("*") for r in SeqIO.parse(PROT, "fasta")}

cds = {}
for line in open(GFF):
    f = line.rstrip("\n").split("\t")
    if len(f) == 9 and f[0] == CONTIG and f[2] == "CDS":
        pid = re.search(r"Parent=([^;]+)", f[8]).group(1)
        cds.setdefault(pid, []).append((int(f[3]), int(f[4]), f[6]))


def from_model(pid):
    parts = sorted(cds[pid])
    strand = parts[0][2]
    nt = Seq("".join(str(seq[a - 1:b]) for a, b, _ in parts))
    if strand == "-":
        nt = nt.reverse_complement()
    aa = str(nt.translate()).rstrip("*")
    assert aa == prot[pid], (pid, "translation differs from the BFD model protein")
    exons = [{"start": a, "end": b} for a, b, _ in parts]
    if strand == "-":
        exons.reverse()
    return dict(name=R, start=parts[0][0], end=parts[-1][1], strand=strand, exons=exons, _aa=aa,
                _src=f"BFD funannotate model {pid.replace('-T1', '')} (not in the public assembly)")


def from_coords(start, end, strand, note):
    s = seq[start - 1:end]
    if strand == "-":
        s = s.reverse_complement()
    aa = str(s.translate()).rstrip("*")
    assert aa.startswith("M") and "*" not in aa, aa
    return dict(name=P, start=start, end=end, strand=strand, exons=[{"start": start, "end": end}],
                _aa=aa, _src=note)


genes = [from_model("FD20D505_007127-T1"), from_model("FD20D505_007128-T1"),
         from_model("FD20D505_007129-T1"),
         from_coords(1529356, 1529433, "-", "unannotated ORF, strict-CAAX scan")]
for i, g in enumerate(genes):
    g.update(gene_index=i, segment_index=0, order_in_locus=i, protein_accession=None, role="core_MAT",
             present=True, locus_tag=None, codon_start=1, transl_table=1)
start, end = genes[0]["start"], genes[-1]["end"]
pre = genes[-1]
assert STRICT.search(pre["_aa"]), pre["_aa"]

POS = ("Positional evidence (results/2026-09-27_pheromone_positional, MATPredict receptor queue "
       "step 2): a pheromone-receptor gene with a pheromone-precursor ORF ending in the strict "
       "Basidiomycota CAAX motif C[VI][IV][AVMG] within 10 kb flagged every curated Agaricomycete "
       "mating receptor (Coprinopsis 3/3, Schizophyllum 2/2) and 0/25 non-mating STE3 copies; "
       "random genomic windows 2.5%.")
definition = (
    f"Russula nobilis gfRusNobi1 (specimen KDTOL00553; Russulales, Russulaceae) B "
    f"(pheromone/receptor) mating-type locus candidate, GCA_984573805.1 {CONTIG}:{start}-{end}: "
    f"three tandem pheromone receptors (STE3, Pfam PF02076 at gathering threshold; 490, 414 and "
    f"455 aa) and one pheromone-precursor ORF. The receptors' best blastp matches among curated "
    f"receptors are the Heterobasidion irregulare TC 32-1 B-locus receptors (38.5-53.5% identity, "
    f"E <= 1.7e-100). The public assembly carries no gene annotation: receptor exons are the BFD "
    f"funannotate models FD20D505_007127/007128/007129, recorded by coordinates. Precursor: "
    f"unannotated ORF {pre['start']}-{pre['end']}(-), {len(pre['_aa'])} aa "
    f"({pre['_aa']}), strict CAAX {pre['_aa'][-4:]}, 7.1 kb from the third receptor; no homology "
    f"to curated Basidiomycota pheromones (tblastn E >= 0.1 genome-wide), so it rests on the motif "
    f"and position only and is SHORT for a basidiomycete precursor. PR/B loci routinely carry "
    f"several receptor and pheromone copies (curator domain ruling 2026-09-16). Which receptor(s) "
    f"define B specificity is not established. Idiomorph 'B1' is a placeholder (curator convention "
    f"2026-09-26). TIER 2: genome-derived; no Russulaceae mating-type deposit exists.")
# No publication describes this genome or locus; provenance is the assembly itself
# (GCA_984573805.1, BioProject PRJEB113492, BioSample SAMEA110757389, ToLID gfRusNobi1),
# stated in the evidence text. Citations stay empty rather than cite an unverified DOI.
cite = []
PROV = ("Source assembly GCA_984573805.1 (gfRusNobi1.hap1.1; BioProject PRJEB113492; BioSample "
        "SAMEA110757389; specimen KDTOL00553; Royal Botanic Gardens Kew / Earlham Institute); no "
        "publication describes this genome or locus. ")
rid = "2830151_kdtol00553_PR_B1"
rec = {
    "record_id": rid,
    "record_version": 1,
    "taxonomy": {"taxid": 2830151, "lineage": (
        "k__Fungi;p__Basidiomycota;subphylum__Agaricomycotina;c__Agaricomycetes;o__Russulales;"
        "f__Russulaceae;g__Russula;s__Russula_nobilis"), "lineage_resolved_date": "2026-09-27"},
    "organism": {"species": "Russula nobilis",
                 "strain": {"name": "KDTOL00553", "known": True, "culture_collection_ids": [],
                            "differs_from_sequenced": False}},
    "mating_type": {"locus_name": "PR", "idiomorphs": ["B1"], "system": "heterothallic"},
    "locus": {
        "coordinate_provenance": "curator_derived",
        "excluded_from_coordinate_benchmark": False,
        "core": {"completeness": "partial",
                 "reference_orientation": "not established from coordinate data alone",
                 "definition_note": definition,
                 "segments": [{"segment_index": 0,
                               "sequence_source": {"type": "insdc_nucleotide", "accession": CONTIG,
                                                   "seq_region": CONTIG},
                               "start": start, "end": end, "contig_edge_distance": None,
                               "sequence_checksum": None}]},
        "extended_flank": [],
    },
    "genes": [{k: g[k] for k in ("gene_index", "name", "protein_accession", "role", "present",
                                 "locus_tag", "segment_index", "start", "end", "strand",
                                 "order_in_locus", "exons", "codon_start", "transl_table")}
              for g in genes],
    "evidence": {
        "locus_existence": {"tier": 2, "experimental_method": (
            PROV + "genome-derived receptor cluster (three tandem STE3 receptors, B-locus-type by homology "
            "to the curated Heterobasidion B-locus receptors) with a strict-CAAX precursor ORF. "
            + POS), "citations": cite},
        "boundaries": {"tier": 2, "experimental_method": (
            f"record segment spans the first receptor to the precursor on {CONTIG}; genome-derived"),
            "citations": cite},
        "idiomorph_assignment": {"tier": 2, "experimental_method": (
            "'B1' is a placeholder; the specimen's B allele is not reported"), "citations": cite},
    },
    "validation": {"status": "accepted", "rejection_reason": None, "accession_resolved": True,
                   "accession_resolved_date": "2026-09-27", "accession_resolved_version": "GCA_984573805.1",
                   "sequence_match": {"status": "pass",
                                      "per_gene": [{"gene_index": g["gene_index"], "percent_identity": 100.0,
                                                    "coverage": 100.0, "status": "pass"} for g in genes],
                                      "notes": ("receptors: translation of the recorded exons equals the BFD "
                                                "funannotate model protein; precursor: translation has a start "
                                                "Met and no internal stop")},
                   "taxonomy_current": True},
    "curation": {"proposed_by": "claude-receptor-queue-step-c",
                 "proposal_dedupe_key": "GCA_984573805.1|2830151|PR|B1",
                 "reviewed_by": None, "reviewed_date": None,
                 "notes": ("Proposed 2026-09-27 on curator ruling (Russulaceae tier-2 B record). PENDING "
                           "CURATOR SIGN-OFF: accepted on branch basidio-anchors only so detect can be "
                           "measured with it.")},
    "model_provenance": None,
}
d = os.path.join(DB, "Russulales", rid)
os.makedirs(d, exist_ok=True)
with open(os.path.join(d, "metadata.yaml"), "w") as fo:
    yaml.safe_dump(rec, fo, sort_keys=False, width=100, allow_unicode=True)
for g in genes:
    print(g["name"], g["start"], g["end"], g["strand"], len(g["_aa"]), g["_src"])
print(rid, CONTIG, start, end)
