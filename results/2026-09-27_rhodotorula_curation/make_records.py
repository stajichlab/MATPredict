"""Write Rhodotorula redPR / redHD metadata.yaml records (tier 2, pending sign-off)."""
import json, os, sys
import yaml

DB = sys.argv[1]  # .../db/Basidiomycota/Sporidiobolales
genes = json.load(open("final_genes.json"))
CONTIG_FIX = {"CAKKSX020000013.1": "CAKKSX030000020.1", "CAKKSX020000026.1": "CAKKSX030000026.1"}
for g in genes:
    g["contig"] = CONTIG_FIX.get(g["contig"], g["contig"])

LIN = "k__Fungi;p__Basidiomycota;sc__Pucciniomycotina;c__Microbotryomycetes;o__Sporidiobolales;f__Sporidiobolaceae;g__Rhodotorula;s__"
STRAINS = {  # key: taxid, NCBI species, strain name, slug, assembly, manuscript species, clade, P/R, HD allele
    "CBS14": (5286, "Rhodotorula toruloides", "CBS 14", "cbs-14", "GCA_921037615.3", "R. toruloides", "B", "A1", "B41"),
    "NBRC_0880": (5286, "Rhodotorula toruloides", "NBRC 0880", "nbrc-0880", "GCA_000988875.2", "R. toruloides", "B", "A2", "B39"),
    "CBS_20": (5535, "Rhodotorula glutinis", "CBS 20", "cbs-20", "GCA_920103745.3", "R. glutinis", "C", "A1", "B63"),
    "JY1105": (5537, "Rhodotorula mucilaginosa", "JY1105", "jy1105", "GCA_024748845.1", "R. aff. mucilaginosa", "A", "A1", "B8"),
    "Y-7192": (86839, "Rhodotorula sphaerocarpa", "NRRL Y-7192", "nrrl-y-7192", "GCA_060419215.1", "R. sphaerocarpa", "A", "A2", "B36"),
    "LS11": (86836, "Rhodotorula kratochvilovae", "LS11", "ls11", "GCA_002917965.1", "R. kratochvilovae", "C", "A2", "B74"),
    "JJ10.1": (29898, "Rhodotorula graminis", "JJ10.1", "jj10-1", "GCA_026119225.1", "Rhodotorula sp. (clade I; R. aff. babjevae)", "C", "A2", "B54"),
}
PREPRINT = {"doi": "10.1101/2025.09.11.675505", "pmid": None}
ZENODO = {"doi": "10.5281/zenodo.18230104", "pmid": None}
COELHO = {"doi": "10.1186/1471-2148-11-249", "pmid": "21880139"}
NAME = {"STE3_A1": "STE3a1", "STE3_A2": "STE3a2", "STE20": "STE20", "HD1": "HD1", "HD2": "HD2"}
ROLE = {"STE3a1": "core_MAT", "STE3a2": "core_MAT", "STE20": "flanking_variable", "HD1": "core_MAT", "HD2": "core_MAT"}

by = {}
for g in genes:
    by.setdefault((g["strain"], g["locus"]), []).append(g)

written = []
for (strain, locus), gs in sorted(by.items()):
    taxid, sp, sname, slug, asm, ms_sp, clade, pr, hd = STRAINS[strain]
    if locus == "PR" and not any(g["gene"].startswith("STE3") for g in gs):
        continue
    if locus == "HD" and len({g["gene"] for g in gs}) < 2:
        continue
    fam = "redPR" if locus == "PR" else "redHD"
    idio = pr if locus == "PR" else hd
    rid = f"{taxid}_{slug}_{fam}_{idio}"
    contigs = {g["contig"] for g in gs}
    assert len(contigs) == 1, (strain, locus, contigs)
    contig = contigs.pop()
    gs = sorted(gs, key=lambda g: g["start"])
    seg_start, seg_end = min(g["start"] for g in gs), max(g["end"] for g in gs)
    gene_rows, per_gene, model_notes = [], [], []
    for i, g in enumerate(gs):
        nm = NAME[g["gene"]]
        gene_rows.append(dict(gene_index=i, name=nm, protein_accession=None, role=ROLE[nm], present=True,
                              locus_tag=None, segment_index=0, start=g["start"], end=g["end"], strand=g["strand"],
                              order_in_locus=i, exons=[dict(start=a, end=b) for a, b in g["exons"]],
                              codon_start=g["codon_start"], transl_table=1))
        per_gene.append(dict(gene_index=i, percent_identity=100.0, coverage=100.0, status="pass"))
        src = ("the group's funannotate/miniprot model (" + g["id"] + ")" if g["model_source"] == "group_gff"
               else "a miniprot 0.18 model built here (" + g["model_source"].split(" from ", 1)[1] + ") because the group's model did not translate cleanly on the public contig")
        model_notes.append(f"{nm} ({g['aa_len']} aa{'' if g['starts_M'] else ', no start Met: 5-prime partial'}) is {src}")
    if locus == "PR":
        ste3 = next(g for g in gs if g["gene"].startswith("STE3"))
        definition = (f"{ms_sp} {sname} (manuscript clade {clade}), P/R MAT locus allele {pr}: "
                      f"{'; '.join(model_notes)}. Assignment {pr} = the manuscript's P/R allele for this strain, "
                      f"and its STE3 groups with the curated Sporidiobolales STE3{pr.lower()} receptors (Coelho et al. 2011). "
                      f"Pheromone precursor (RHA) genes in the locus are not curated: they are too short for the current search "
                      f"(queued with the receptor work).")
        ref_orient = "->".join(NAME[g["gene"]] for g in gs)
        idio_method = (f"{pr} by the manuscript's allele assignment (STE3.{pr} receptor clade) and by blastp of the "
                       f"receptor to the curated Coelho 2011 STE3{pr.lower()} records; genome-derived, not by crosses")
    else:
        definition = (f"{ms_sp} {sname} (manuscript clade {clade}), HD MAT locus, divergently transcribed HD1/HD2 pair; "
                      f"allele {hd} is the manuscript's HD allele name for this strain. {'; '.join(model_notes)}. "
                      f"The HD locus sits on a different contig from the P/R locus, as the manuscript reports for the genus.")
        ref_orient = "HD2<->HD1" if gs[0]["gene"] == "HD2" else "HD1<->HD2"
        idio_method = (f"allele label {hd} from the manuscript's HD allele catalogue; alleles are not callable by homology "
                       f"(pattern vocabulary)")
    evid_text = ("TIER 2 (genome-derived, preprint). Locus located and genes re-annotated by the curator's group "
                 "(Liu, Tsai, Coelho, ... Stajich) with Leucosporidium scottii CBS 5931 MAT genes as reference, manual "
                 "re-annotation and miniprot for short or hypervariable genes; reported in a bioRxiv preprint that is NOT "
                 "yet peer-reviewed. No locus-specific GenBank deposit exists. The public assembly " + asm + " carries the "
                 "sequence; the gene coordinates are the group's models, verified here to translate on the public contig.")
    rec = {
        "record_id": rid, "record_version": 1,
        "taxonomy": {"taxid": taxid, "lineage": LIN + sp.replace(" ", "_"), "lineage_resolved_date": "2026-09-27"},
        "organism": {"species": sp, "strain": {"name": sname, "known": True,
                                                "culture_collection_ids": [sname] if strain not in ("JY1105", "LS11", "JJ10.1") else [],
                                                "differs_from_sequenced": False}},
        "mating_type": {"locus_name": fam, "idiomorphs": [idio], "system": "heterothallic"},
        "locus": {"coordinate_provenance": "published_explicit", "excluded_from_coordinate_benchmark": True,
                  "core": {"completeness": "partial", "reference_orientation": ref_orient, "definition_note": definition,
                           "segments": [{"segment_index": 0,
                                         "sequence_source": {"type": "insdc_nucleotide", "accession": contig, "seq_region": contig},
                                         "start": seg_start, "end": seg_end, "contig_edge_distance": None,
                                         "sequence_checksum": None}]},
                  "extended_flank": []},
        "genes": gene_rows,
        "evidence": {
            "locus_existence": {"tier": 2, "experimental_method": evid_text, "citations": [PREPRINT, ZENODO]},
            "boundaries": {"tier": 2, "experimental_method": "record segment spans the curated gene models on " + contig,
                           "citations": [PREPRINT]},
            "idiomorph_assignment": {"tier": 2, "experimental_method": idio_method,
                                     "citations": [PREPRINT] + ([COELHO] if locus == "PR" else [])},
        },
        "validation": {"status": "accepted", "rejection_reason": None, "accession_resolved": True,
                       "accession_resolved_date": "2026-09-27", "accession_resolved_version": contig,
                       "sequence_match": {"status": "pass", "per_gene": per_gene,
                                          "notes": "No protein accession: sequences are translations of the recorded "
                                                   "coordinates on the public contig (checked: no internal stop)."},
                       "taxonomy_current": True},
        "curation": {"proposed_by": "claude-curation-from-group-preprint",
                     "proposal_dedupe_key": f"{PREPRINT['doi']}|{taxid}|{fam}|{idio}|{slug}",
                     "reviewed_by": None, "reviewed_date": None,
                     "notes": ("Proposed 2026-09-27 on curator request (Rhodotorula MAT from the group's preprint, incl. HD). "
                               "PENDING CURATOR SIGN-OFF: accepted on the unmerged branch curation-puccinio only so detect "
                               "can be measured with it." + (f" NCBI names this assembly {sp}; the manuscript places it in "
                                                             f"{ms_sp}." if sp.split()[-1] not in ms_sp else ""))},
        "model_provenance": None,
    }
    d = os.path.join(DB, rid)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "metadata.yaml"), "w") as fo:
        yaml.safe_dump(rec, fo, sort_keys=False, width=100, allow_unicode=True)
    written.append(rid)
print("\n".join(written))
