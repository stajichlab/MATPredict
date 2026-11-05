"""Draft MATog candidate records for six divergent Ogataea genomes from our own miniprot projections (tier 3).

Models: models.py (miniprot projection of the DL-1 and CBS 4732 MATog proteins onto the called contigs; a gene is kept only when the
projection has no frameshift and no internal stop). Windows are the contig, or the called locus +/- 20 kb; coordinates are written in contig
coordinates. Run from the repo root: python3 results/2026-11-04_matog_candidates/make_records.py [outdir]
"""
import glob, os, sys, yaml
HERE = "results/2026-11-04_matog_candidates/"
os.chdir(HERE); sys.path.insert(0, ".")
from models import models
OUT = sys.argv[1] if len(sys.argv) > 1 else "."
LIN = "k__Fungi;p__Ascomycota;subphylum__Saccharomycotina;c__Pichiomycetes;o__Pichiales;f__Pichiaceae;g__Ogataea;s__"
CIT = [{"pmid": "25412462", "doi": "10.1371/journal.pgen.1004796"}, {"pmid": "25349420", "doi": "10.1073/pnas.1416014111"}]
CORE = {"MATa1", "MATalpha1", "MATalpha2"}

# slug, record id, taxid, species, species lineage tail, strain, culture ids, assembly, [(window glob, contig)], idiomorphs, note
RECORDS = [
 ("pignaliae", "51923_nrrl-y-17664_MATog_MATa-MATalpha", 51923, "Ogataea pignaliae", "Ogataea_pignaliae", "NRRL Y-17664", ["NRRL Y-17664"], "GCA_030570775.1",
  ["pignaliae__JAJLWO010000194.1", "pignaliae__JAJLWO010000212.1"], ["MATa", "MATalpha"]),
 ("wickerhamii", "1249076_nrrl-yb-4943_MATog_MATa-MATalpha", 1249076, "Ogataea wickerhamii", "Ogataea_wickerhamii", "NRRL YB-4943", ["NRRL YB-4943"], "GCA_030559035.1",
  ["wickerhamii__JAKTQX010000069.1", "wickerhamii__JAKTQX010000058.1"], ["MATa", "MATalpha"]),
 ("dorogensis", "169506_nrrl-y-27599_MATog_MATa-MATalpha", 169506, "Ogataea dorogensis", "Ogataea_dorogensis", "NRRL Y-27599", ["NRRL Y-27599"], "GCA_030568015.1",
  ["dorogensis__JAJLTE010000037.1", "dorogensis__JAJLTE010000018.1"], ["MATa", "MATalpha"]),
 ("siamensis", "246385_nrrl-y-27963_MATog_MATa-MATalpha", 246385, "Ogataea siamensis", "Ogataea_siamensis", "NRRL Y-27963", ["NRRL Y-27963"], "GCA_030574315.1",
  ["siamensis__JAJMET010000031.1", "siamensis__JAJMET010000024.1"], ["MATa", "MATalpha"]),
 ("boidinii_alpha", "5477_ndk27a1_MATog_MATalpha", 5477, "Candida boidinii", "Candida_boidinii", "NDK27A1", ["NDK27A1"], "GCA_002008015.1",
  ["boidinii_alpha__MSSE01000022.1"], ["MATalpha"]),
 ("boidinii_a", "5477_uniss-cb60_MATog_MATa", 5477, "Candida boidinii", "Candida_boidinii", "UNISS-Cb60", ["UNISS-Cb60"], "GCA_002008025.1",
  ["boidinii_a__MSRY01000009.1"], ["MATa"]),
]
ORDER = {"sla2": 0, "MATa1": 1, "MATalpha1": 2, "MATalpha2": 3, "dic1": 4}


def window_of(prefix):
    fa = glob.glob(f"data/{prefix}_*.fa")[0]
    contig, rng = fa.split("__")[1].rsplit(".fa", 1)[0].rsplit("_", 1)
    return contig, int(rng.split("-")[0]) - 1, fa


def build(rec):
    slug, rid, taxid, species, sptail, strain, cc, asm, wins, idios = rec
    genes, segs, notes = [], [], []
    for si, w in enumerate(wins):
        contig, off, fa = window_of(w)
        ms = models(w)
        rej = ms.pop("_rejected")
        sel = sorted(ms.items(), key=lambda kv: kv[1]["start"])
        for name, m in sel:
            exons = sorted(m["exons"], reverse=(m["strand"] == "-"))
            genes.append({"name": name, "role": "core_MAT" if name in CORE else "flanking_conserved", "present": True, "locus_tag": None,
                          "segment_index": si, "start": m["start"] + off, "end": m["end"] + off, "strand": m["strand"],
                          "exons": [{"start": a + off, "end": b + off} for a, b in exons],
                          "_ref": m["ref"], "_aa": len(m["prot"].rstrip("*")), "_reflen": m["reflen"], "_t0": m["t0"], "_t1": m["t1"], "_id": m["ident"],
                          "_contig_len": m["len"] + off if off == 0 else None, "_prot": m["prot"].rstrip("*")})
        lo = min(g["start"] for g in genes if g["segment_index"] == si); hi = max(g["end"] for g in genes if g["segment_index"] == si)
        segs.append({"segment_index": si, "sequence_source": {"type": "insdc_nucleotide", "accession": contig, "seq_region": contig},
                     "start": lo, "end": hi, "contig_edge_distance": None, "sequence_checksum": None})
        for m in rej:
            notes.append(f"{m['name']} on {contig}: projection rejected (frameshift or internal stop)")
    for si in sorted({g["segment_index"] for g in genes}):
        pass
    genes.sort(key=lambda g: (g["segment_index"], g["start"]))
    for i, g in enumerate(genes):
        g["gene_index"] = i; g["order_in_locus"] = i
    return genes, segs, notes


def note_for(rec, genes, segs, notes):
    slug, rid, taxid, species, sptail, strain, cc, asm, wins, idios = rec
    parts = [f"{g['name']} {g['segment_index'] and 'seg1' or 'seg0'} {g['start']}-{g['end']}({g['strand']}) {g['_aa']} aa, {g['_t0']}-{g['_t1']} of {g['_reflen']} aa of the {'DL-1' if g['_ref'].startswith('1005962') else 'CBS 4732'} protein at {100*g['_id']:.0f}% identity" for g in genes]
    shape = ("two cassette ends on two contigs of the assembly (SLA2-MATa1 end on one, MATalpha1-MATalpha2-DIC1 end on the other)" if len(wins) == 2
             else f"one cassette end (MAT{'alpha' if 'MATalpha' in idios else 'a'} idiomorph; the strain carries one idiomorph)")
    return (f"CANDIDATE MATog record of {species} {strain} ({asm}), a divergent member of the genus (MATog score 67-142 in the v0.8.0 run against the two polymorpha-clade references). "
            f"OUR OWN PROJECTION (tier 3): each gene is the miniprot alignment of the O. parapolymorpha DL-1 or O. polymorpha CBS 4732 MATog proteins to the contig, kept only when it has no frameshift and no internal stop; "
            f"no RNA-seq or annotation was used (results/2026-11-04_matog_candidates/README.md). The record holds {shape}. Genes: " + "; ".join(parts) + ". "
            f"PARTIAL BY CONSTRUCTION: every model starts or stops at the alignment end or at a contig end, so none has a verified start or stop codon, and no MAT gene is recorded as absent. "
            + ("; ".join(notes) + ". " if notes else "") + "GENETIC CODE: standard (NCBI Taxonomy code 1 for the species).")


def record(rec):
    slug, rid, taxid, species, sptail, strain, cc, asm, wins, idios = rec
    genes, segs, notes = build(rec)
    out_genes = [{k: v for k, v in g.items() if not k.startswith("_")} for g in genes]
    out_genes = [{"gene_index": g["gene_index"], "name": g["name"], "protein_accession": None, **{k: v for k, v in g.items() if k not in ("gene_index", "name")}} for g in out_genes]
    core_names = [g["name"] for g in genes if g["role"] == "core_MAT"]
    orient = "->".join(g["name"] for g in genes)
    return {"record_id": rid, "record_version": 1,
            "taxonomy": {"taxid": taxid, "lineage": LIN + sptail, "lineage_resolved_date": "2026-11-04"},
            "organism": {"species": species, "strain": {"name": strain, "known": True, "culture_collection_ids": cc, "differs_from_sequenced": False}},
            "mating_type": {"locus_name": "MATog", "idiomorphs": idios, "system": "unknown"},
            "locus": {"assembly_accession": asm, "assembly_accession_basis": "strain_match_verified", "coordinate_provenance": "curator_derived",
                      "excluded_from_coordinate_benchmark": False,
                      "core": {"completeness": "partial", "reference_orientation": orient, "definition_note": note_for(rec, genes, segs, notes), "segments": segs},
                      "extended_flank": []},
            "genes": out_genes,
            "evidence": {"locus_existence": {"tier": 3, "experimental_method": "homology to the two MATog cassette records of the genus (DL-1, CBS 4732; two-locus architecture in H. polymorpha, PMIDs 25412462 and 25349420); the v0.8.0 detection called this locus; no study of this species",
                                             "citations": CIT},
                         "boundaries": {"tier": 3, "experimental_method": "OUR OWN miniprot projection of the DL-1 and CBS 4732 MATog proteins onto the assembly contigs; no transcript support", "citations": CIT},
                         "idiomorph_assignment": {"tier": 3, "experimental_method": "gene identity by best-scoring MATog reference protein (MATa1 versus MATalpha1/MATalpha2); the idiomorph follows the gene content", "citations": CIT}},
            "validation": {"status": "needs_review", "rejection_reason": None, "accession_resolved": False, "accession_resolved_date": None,
                           "accession_resolved_version": None, "sequence_match": {"status": "pass", "per_gene": [], "notes": ""}, "taxonomy_current": None},
            "curation": {"proposed_by": "claude-code-curation-agent", "proposal_dedupe_key": f"25412462|{taxid}|MATog|{'-'.join(idios)}|{strain}", "reviewed_by": None, "reviewed_date": None},
            "model_provenance": None}, genes


if __name__ == "__main__":
    for rec in RECORDS:
        r, genes = record(rec)
        with open(f"{OUT}/{r['record_id']}.yaml", "w") as fh:
            yaml.safe_dump(r, fh, sort_keys=False, width=100)
        with open(f"{OUT}/{r['record_id']}.proteins.faa", "w") as fh:
            for g in genes:
                fh.write(f">{r['record_id']}|gene_index={g['gene_index']}|name={g['name']}|role={g['role']}\n{g['_prot']}\n")
        print(r["record_id"], [(g["name"], g["segment_index"], g["start"], g["end"], g["strand"], g["_aa"]) for g in genes])
