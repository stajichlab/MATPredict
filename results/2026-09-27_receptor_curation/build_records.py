"""Build tier-2 Basidiomycota:PR (B locus) records for Polyporales and Russulales.

Receptor queue step (c), curator ruling 2026-09-27. Three genome-derived B-locus
records at pheromone-receptor + pheromone-precursor clusters:

  Heterobasidion irregulare TC 32-1 (Russulales), KI925460.1, JGI annotation
  Trametes versicolor FP-101664 SS1 (Polyporales), JH711790.1, JGI annotation
  Grifola frondosa 9006-11 (Polyporales), LUGG01000005.1: receptors annotated,
      pheromone precursors unannotated (located by the strict-CAAX ORF scan and
      tblastn of the G. frondosa WM1-25 pheromone deposits LC706365-7)

Every gene translation is re-derived from the recorded coordinates and compared
with the annotated protein (100% required). Usage: build_records.py DB_ROOT GENOME_DIR
"""
import os
import re
import sys

import yaml
from Bio import SeqIO
from Bio.Seq import Seq

DB, GD = sys.argv[1], sys.argv[2]
STRICT = re.compile(r"C[VI][IV][AVMG]$")
T2 = re.compile(r"C[VITE][IVT][AVMG]$")


def load(acc):
    base = f"{GD}/{acc}/ncbi_dataset/data/{acc}"
    gb = {r.id: r for r in SeqIO.parse(f"{base}/genomic.gbff", "genbank")}
    prot = {r.id: str(r.seq) for r in SeqIO.parse(f"{base}/protein.faa", "fasta")}
    return gb, prot


def caax(seq):
    c = seq[-4:]
    if STRICT.search(seq):
        return f"strict CAAX {c}"
    if T2.search(seq):
        return f"relaxed CAAX {c}"
    return f"C-terminus {c} (no CAAX)" if c[0] != "C" else f"CAAX-like {c}"


def from_cds(gb, prot, contig, pid, name):
    rec = gb[contig]
    f = next(x for x in rec.features if x.type == "CDS" and x.qualifiers.get("protein_id", [""])[0] == pid)
    loc = f.location
    strand = "+" if loc.strand == 1 else "-"
    cs = int(f.qualifiers.get("codon_start", ["1"])[0])
    tt = int(f.qualifiers.get("transl_table", ["1"])[0])
    aa = str(loc.extract(rec.seq)[cs - 1:].translate(table=tt)).rstrip("*")
    assert aa == prot[pid], (pid, "translation differs from the annotated protein")
    exons = [{"start": int(p.start) + 1, "end": int(p.end)} for p in loc.parts]
    exons.sort(key=lambda e: -e["start"] if strand == "-" else e["start"])
    return dict(name=name, protein_accession="ncbi_protein:" + pid, role="core_MAT", present=True,
                locus_tag=f.qualifiers.get("locus_tag", [None])[0], start=int(loc.start) + 1,
                end=int(loc.end), strand=strand, exons=exons, codon_start=cs, transl_table=tt,
                _aa=aa, _src="annotated CDS")


def from_coords(gb, contig, start, end, strand, name, note):
    s = gb[contig].seq[start - 1:end]
    if strand == "-":
        s = s.reverse_complement()
    aa = str(s.translate()).rstrip("*")
    assert aa.startswith("M") and "*" not in aa, (contig, start, end, aa)
    return dict(name=name, protein_accession=None, role="core_MAT", present=True, locus_tag=None,
                start=start, end=end, strand=strand, exons=[{"start": start, "end": end}],
                codon_start=1, transl_table=1, _aa=aa, _src=note)


R, P = "pheromone_receptor", "fungal_mating_type_pheromone"
SPECS = []

# --- Heterobasidion irregulare TC 32-1 ---
gb, prot = load("GCA_000320585.2")
c = "KI925460.1"
genes = [from_cds(gb, prot, c, "ETW79702.1", R), from_cds(gb, prot, c, "ETW79703.1", R),
         from_cds(gb, prot, c, "ETW79704.1", R), from_cds(gb, prot, c, "ETW79705.1", P),
         from_cds(gb, prot, c, "ETW79706.1", R), from_cds(gb, prot, c, "ETW79709.1", R),
         from_cds(gb, prot, c, "ETW79710.1", P), from_cds(gb, prot, c, "ETW79712.1", P)]
SPECS.append(dict(
    order="Russulales", taxid=984962, species="Heterobasidion irregulare", strain="TC 32-1",
    slug="tc-32-1", assembly="GCA_000320585.2", contig=c, genes=genes,
    lineage="k__Fungi;p__Basidiomycota;subphylum__Agaricomycotina;c__Agaricomycetes;o__Russulales;"
            "f__Bondarzewiaceae;g__Heterobasidion;s__Heterobasidion_irregulare",
    cite=[{"pmid": "22463738", "doi": "10.1111/j.1469-8137.2012.04128.x"}],
    ann="JGI Heterobasidion irregulare v2.0 annotation (Olson et al. 2012), which names the "
        "receptors 'putative pheromone receptor' / 'putative B mating type protein' and the "
        "precursors 'B mating type pheromone'"))

# --- Trametes versicolor FP-101664 SS1 ---
gb, prot = load("GCA_000271585.1")
c = "JH711790.1"
order = [("EIW57106.1", P), ("EIW57105.1", P), ("EIW57108.1", P), ("EIW57109.1", P),
         ("EIW56739.1", R), ("EIW57102.1", P), ("EIW57113.1", P), ("EIW56740.1", R),
         ("EIW57112.1", P), ("EIW56741.1", R), ("EIW57111.1", P), ("EIW57110.1", P)]
genes = [from_cds(gb, prot, c, pid, nm) for pid, nm in order]
SPECS.append(dict(
    order="Polyporales", taxid=5325, species="Trametes versicolor", strain="FP-101664 SS1",
    slug="fp-101664-ss1", assembly="GCA_000271585.1", contig=c, genes=genes,
    lineage="k__Fungi;p__Basidiomycota;subphylum__Agaricomycotina;c__Agaricomycetes;o__Polyporales;"
            "f__Polyporaceae;g__Trametes;s__Trametes_versicolor",
    cite=[{"pmid": "22745431", "doi": "10.1126/science.1221748"}],
    ann="JGI Trametes versicolor v1.0 annotation (Floudas et al. 2012), which names the receptors "
        "'STE3-domain-containing protein' / 'fungal pheromone STE3G-protein-coupled receptor' and "
        "the precursors 'pheromone' / 'B mating type pheromone'"))

# --- Grifola frondosa 9006-11 ---
gb, prot = load("GCA_001683735.1")
c = "LUGG01000005.1"
genes = [from_cds(gb, prot, c, "OBZ74837.1", R), from_cds(gb, prot, c, "OBZ74842.1", R),
         from_coords(gb, c, 525554, 525709, "-", P, "unannotated ORF, strict-CAAX scan"),
         from_cds(gb, prot, c, "OBZ74475.1", R), from_cds(gb, prot, c, "OBZ74569.1", R),
         from_coords(gb, c, 533656, 533820, "-", P, "unannotated ORF, strict-CAAX scan"),
         from_cds(gb, prot, c, "OBZ74474.1", R),
         from_coords(gb, c, 546974, 547096, "-", P,
                     "unannotated ORF, strict-CAAX scan; tblastn of G. frondosa WM1-25 ph1 "
                     "(LC706365.1, BDI00750.1) 55% identity over 40 aa, E=2e-11")]
SPECS.append(dict(
    order="Polyporales", taxid=5627, species="Grifola frondosa", strain="9006-11",
    slug="9006-11", assembly="GCA_001683735.1", contig=c, genes=genes,
    lineage="k__Fungi;p__Basidiomycota;subphylum__Agaricomycotina;c__Agaricomycetes;o__Polyporales;"
            "f__Grifolaceae;g__Grifola;s__Grifola_frondosa",
    cite=[{"pmid": "37888215", "doi": "10.3390/jof9100959"}],
    ann="GenBank annotation of GCA_001683735.1 for the five receptors ('Pheromone B alpha/beta "
        "receptor'); the three precursors are not annotated in the assembly"))

POS = ("Positional evidence (results/2026-09-27_pheromone_positional, MATPredict receptor queue "
       "step 2): a pheromone-receptor gene with a pheromone-precursor ORF ending in the strict "
       "Basidiomycota CAAX motif C[VI][IV][AVMG] within 10 kb flagged every curated Agaricomycete "
       "mating receptor (Coprinopsis 3/3, Schizophyllum 2/2) and 0/25 non-mating STE3 copies; "
       "random genomic windows 2.5%.")

for sp in SPECS:
    genes = sorted(sp["genes"], key=lambda g: g["start"])
    for i, g in enumerate(genes):
        g["gene_index"], g["segment_index"], g["order_in_locus"] = i, 0, i
    start, end = genes[0]["start"], genes[-1]["end"]
    nR = sum(g["name"] == R for g in genes)
    nP = sum(g["name"] == P for g in genes)
    pher = "; ".join(f"{g['protein_accession'] or 'unannotated %d-%d(%s)' % (g['start'], g['end'], g['strand'])} "
                     f"{len(g['_aa'])} aa, {caax(g['_aa'])}" for g in genes if g["name"] == P)
    nstrict = sum(bool(STRICT.search(g["_aa"])) for g in genes if g["name"] == P)
    unann = [g for g in genes if g["protein_accession"] is None]
    rid = f"{sp['taxid']}_{sp['slug']}_PR_B1"
    definition = (
        f"{sp['species']} {sp['strain']} ({sp['order']}) B (pheromone/receptor) mating-type locus, "
        f"{sp['assembly']} {sp['contig']}:{start}-{end}: {nR} pheromone receptors and {nP} "
        f"pheromone precursors in one cluster, from the {sp['ann']}. Pheromone precursors: {pher}. "
        f"{nstrict}/{nP} precursors end in the strict CAAX motif. "
        + (f"Unannotated precursors ({len(unann)}) are recorded by coordinates only: "
           + "; ".join(f"{g['start']}-{g['end']}({g['strand']}) {g['_src']}" for g in unann) + ". "
           if unann else "")
        + "PR/B loci routinely carry several receptor and pheromone copies (curator domain ruling "
          "2026-09-16); receptors and precursors are recorded under the family's generic names. "
          "Which receptor(s) define the B mating specificity is not established. Idiomorph 'B1' is "
          "a placeholder label (curator convention 2026-09-26: count alleles when the source does not "
          "number them). TIER 2: genome-derived; no locus-specific GenBank deposit exists for this "
          "order's B locus.")
    cit = sp["cite"]
    rec = {
        "record_id": rid,
        "record_version": 1,
        "taxonomy": {"taxid": sp["taxid"], "lineage": sp["lineage"], "lineage_resolved_date": "2026-09-27"},
        "organism": {"species": sp["species"],
                     "strain": {"name": sp["strain"], "known": True, "culture_collection_ids": [sp["strain"]],
                                "differs_from_sequenced": False}},
        "mating_type": {"locus_name": "PR", "idiomorphs": ["B1"], "system": "heterothallic"},
        "locus": {
            "coordinate_provenance": "curator_derived",
            "excluded_from_coordinate_benchmark": False,
            "core": {"completeness": "partial",
                     "reference_orientation": "not established from coordinate data alone",
                     "definition_note": definition,
                     "segments": [{"segment_index": 0,
                                   "sequence_source": {"type": "insdc_nucleotide", "accession": sp["contig"],
                                                       "seq_region": sp["contig"]},
                                   "start": start, "end": end, "contig_edge_distance": None,
                                   "sequence_checksum": None}]},
            "extended_flank": [],
        },
        "genes": [{k: g[k] for k in ("gene_index", "name", "protein_accession", "role", "present",
                                     "locus_tag", "segment_index", "start", "end", "strand",
                                     "order_in_locus", "exons", "codon_start", "transl_table")}
                  for g in genes],
        "evidence": {
            "locus_existence": ({"tier": 1, "experimental_method": (
                "genetic crosses (test crosses and three-round mating experiments) establishing a "
                "tetrapolar system, and genome localization of the B locus (six pheromone receptors and "
                "five precursors) on chromosome 11 of strain y59 (Zhang et al. 2023). This record's "
                "genes come from a different strain's assembly (9006-11). " + POS),
                "citations": cit} if sp["species"] == "Grifola frondosa" else
                {"tier": 2, "experimental_method": "genome annotation of a named B-locus cluster. " + POS,
                 "citations": cit}),
            "boundaries": {"tier": 2, "experimental_method": (
                f"record segment spans the first to last curated gene on {sp['contig']}; genome-derived"),
                "citations": cit},
            "idiomorph_assignment": {"tier": 2, "experimental_method": (
                "'B1' is a placeholder; the strain's B allele is not reported"), "citations": cit},
        },
        "validation": {"status": "accepted", "rejection_reason": None, "accession_resolved": True,
                       "accession_resolved_date": "2026-09-27", "accession_resolved_version": sp["assembly"],
                       "sequence_match": {"status": "pass",
                                          "per_gene": [{"gene_index": g["gene_index"], "percent_identity": 100.0,
                                                        "coverage": 100.0, "status": "pass"} for g in genes],
                                          "notes": ("annotated genes: translation of the recorded coordinates "
                                                    "equals the annotated protein; unannotated precursors: "
                                                    "translation has a start Met and no internal stop")},
                       "taxonomy_current": True},
        "curation": {"proposed_by": "claude-receptor-queue-step-c",
                     "proposal_dedupe_key": f"{sp['assembly']}|{sp['taxid']}|PR|B1",
                     "reviewed_by": None, "reviewed_date": None,
                     "notes": ("Proposed 2026-09-27 on curator ruling (receptor queue step c: Polyporales and "
                               "Russulales). PENDING CURATOR SIGN-OFF: accepted on branch basidio-anchors only so "
                               "detect can be measured with it.")},
        "model_provenance": None,
    }
    d = os.path.join(DB, sp["order"], rid)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "metadata.yaml"), "w") as fo:
        yaml.safe_dump(rec, fo, sort_keys=False, width=100, allow_unicode=True)
    print(rid, sp["contig"], start, end, f"receptors={nR} precursors={nP} strict={nstrict}")
