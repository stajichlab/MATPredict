"""Write the synthetic detect reports in fixtures/ (shapes copied from detect/report.py).

Run: python tests/report/make_fixtures.py tests/report/fixtures. fola_50a_old_format.yaml is a real run
(results/2026-10-06_fola_50a), copied, not generated.
"""
import copy, yaml, sys
from pathlib import Path
out = Path(sys.argv[1])

def gene(name, role, contig, start, end, strand, ident, status="polished_agree", rec="rec", exons=None, method="exonerate_refine", cov=None, alt=None, ev=1e-50, bs=300.0):
    return {"gene": name, "role": role, "contig": contig, "start": start, "end": end, "strand": strand,
            "identity": ident, "coverage": cov, "reference_record": rec, "method": method, "status": status,
            "alternate_model": alt, "exons": exons or [{"start": start, "end": end}], "evalue": ev, "bitscore": bs}

def base(run=None, **kw):
    d = {"routing_mode": "lineage", "routing_error": None, "not_searched_reason": None, "genetic_code": 1,
         "genetic_code_error": None, "taxonomy_source": "local NCBI taxonomy snapshot 2026-10-01",
         "suppressed_unpolished": 0, "suppressed_flank_carried": 0, "suppressed_mat_gene_gate": 0,
         "suppressed_paralog_class": 0, "suppressed_loci": [], "assembly_gap_at_locus": [], "zygosity": None,
         "two_idiomorphs": [], "receptor_arrays_note": "Report only: receptor arrays occur for mating and non-mating receptors alike.",
         "receptor_arrays": [], "families_attempted": [], "detected": [], "not_detected": []}
    if run: d = {"run": run, **d}
    d.update(kw); return d

def run(sample, organism, taxid, phylum, contigs=812, length=48_213_377, n50=1_203_321):
    return {"sample": sample, "organism": organism, "taxid": taxid, "phylum": phylum, "matpredict_version": "0.6.1",
            "database": {"root": "/app/db", "content_sha256": "3f2a9c41d07be5aa0c1e9d4f6b2a7781e3d05c6b9a1f4e2d8c7b6a5f4e3d2c1b", "records": 150},
            "taxonomy_source": "local NCBI taxonomy snapshot 2026-10-01",
            "genome": {"file": f"{sample}.fna", "sha256": "9d1e4c0b7a2f3e5d6c8b9a0f1e2d3c4b5a6f7e8d9c0b1a2f3e4d5c6b7a8f9e0d", "contigs": contigs, "length_bp": length, "n50": n50},
            "proteins": None, "parameters": {"phylum": None, "exhaustive": False, "genetic_code": None, "min_hits": 2, "require_core_role": True, "max_polished_clusters_per_family": 6, "polish_strong_identity": 50.0},
            "started": "2026-10-07T14:03:11Z", "wall_seconds": 96}

def call(fam, contig, start, end, idio, conf, genes, cls="mat_locus", **kw):
    core_bits = max([g["bitscore"] for g in genes if g["role"] == "core_MAT" and g["bitscore"]] or [0.0])
    d = {"family": fam, "contig": contig, "start": start, "end": end, "confidence": conf, "idiomorph": idio,
         "idiomorph_class": None, "idiomorph_candidates": [{"idiomorph": idio, "score": core_bits}], "detection_pass": "strict",
         "locus_class": cls, "idiomorph_unmodelled": False, "verification": None, "merged_from": [], "subloci": [],
         "polished_genes": sum(1 for g in genes if g["status"].startswith("polished")), "distinct_intervals": len(genes),
         "span_exceeds_plausible_bound": False, "idiomorph_margin": None, "idiomorph_classifier": None,
         "idiomorph_resolutions": [], "ambiguous_with": [], "genes_found": [g["gene"] for g in genes],
         "genes_missing": [], "genes_not_searchable": [], "fragmented": False,
         "reference_records": sorted({g["reference_record"] for g in genes if g["reference_record"]}),
         "segments": [{"contig": contig, "start": start, "end": end, "contig_edge_distance": 210_455}], "gene_evidence": genes}
    d.update(kw); return d

# 1. Phycomyces: Mucoromycota with the HMM classifier.
c = "NW_017265134.1"
genes = [gene("tptA", "flanking_conserved", c, 3978283, 3980101, "+", 98.2, rec="4837_nrrl1555_MAT_Minus", exons=[{"start": 3978283, "end": 3978390}, {"start": 3978460, "end": 3980101}]),
         gene("sexM", "core_MAT", c, 3984120, 3984918, "-", 100.0, rec="4837_nrrl1555_MAT_Minus", exons=[{"start": 3984120, "end": 3984602}, {"start": 3984680, "end": 3984918}]),
         gene("rnhA", "flanking_conserved", c, 3987022, 3988410, "+", 99.1, rec="4837_nrrl1555_MAT_Minus"),
         gene("algA", "flanking_variable", c, 3989877, 3991589, "-", 97.5, status="polished_single", rec="4837_nrrl1555_MAT_Minus")]
phy = base(run("Phybl2", "Phycomyces blakesleeanus NRRL 1555", 4837, "Mucoromycota", contigs=80, length=53_939_167, n50=4_521_222),
           families_attempted=["Mucoromycota:MAT"],
           detected=[call("Mucoromycota:MAT", c, 3978283, 3991589, "Minus", "high", genes,
                          idiomorph_candidates=[{"idiomorph": "Minus", "score": 402.3, "basis": "hmm_classifier"},
                                                {"idiomorph": "Plus", "score": 61.9, "basis": "hmm_classifier"}],
                          idiomorph_margin=340.4,
                          idiomorph_classifier={"method": "profile_hmm", "idiomorph": "Minus", "scores": {"Minus": 402.3, "Plus": 61.9},
                                                "margin": 340.4, "classifier_input": "polished_model", "min_margin": 25.0,
                                                "proteins_scored": 1, "genes_scored": ["sexM"], "manifest_sha256": "ab12"})])
(out / "phycomyces_classifier.yaml").write_text(yaml.safe_dump(phy, sort_keys=False))

# 2. Two idiomorphs, unlinked, Rhizopus-like.
c1, c2 = "scaffold_12", "scaffold_88"
g1 = [gene("tptA", "flanking_conserved", c1, 50100, 51800, "+", 91.0), gene("sexP", "core_MAT", c1, 53210, 53900, "+", 88.4), gene("rnhA", "flanking_conserved", c1, 55020, 56300, "-", 93.2)]
g2 = [gene("sexM", "core_MAT", c2, 1200, 1980, "-", 71.3, status="polished_disagree", alt={"contig": c2, "start": 1150, "end": 1980, "strand": "-", "exons": [{"start": 1150, "end": 1500}, {"start": 1560, "end": 1980}], "identity": 69.0, "method": "miniprot_refine"})]
two = base(run("RS-114", "Rhizopus stolonifer", 4846, "Mucoromycota"), families_attempted=["Mucoromycota:MAT"],
           detected=[call("Mucoromycota:MAT", c1, 50100, 56300, "Plus", "high", g1),
                     call("Mucoromycota:MAT", c2, 1200, 1980, "Minus", "low", g2, cls="idiomorph_gene_only",
                          segments=[{"contig": c2, "start": 1200, "end": 1980, "contig_edge_distance": 1199}])],
           two_idiomorphs=[{"family": "Mucoromycota:MAT", "arrangement": "unlinked",
                            "plus_call": {"idiomorph": "Plus", "contig": c1, "start": 50100, "end": 56300, "confidence": "high", "locus_class": "mat_locus", "classifier_input": "polished_model", "margin": 210.0, "flanks": []},
                            "minus_call": {"idiomorph": "Minus", "contig": c2, "start": 1200, "end": 1980, "confidence": "low", "locus_class": "idiomorph_gene_only", "classifier_input": "hsp_fragment", "margin": 12.0, "flanks": []},
                            "evidence": {"distance_bp": None, "gc_plus_contig": 38.1, "gc_minus_contig": 46.9, "gc_difference_pct": 8.8, "weak_calls": ["Minus"], "shared_flanks": []},
                            "possible_causes": [], "supported_causes": ["duplication", "mixed_culture_or_heterokaryon"],
                            "not_assessed": ["read depth"]}])
(out / "two_idiomorphs.yaml").write_text(yaml.safe_dump(two, sort_keys=False))

# 3. Not searched.
ns = base(run("Bd-JEL423", "Batrachochytrium dendrobatidis JEL423", 109871, "Chytridiomycota"), routing_mode="not_searched",
          not_searched_reason="no curated MAT family covers phylum Chytridiomycota; run with --exhaustive to search every family anyway")
(out / "not_searched.yaml").write_text(yaml.safe_dump(ns, sort_keys=False))

# 4. Nothing called, assembly gap, withheld loci.
nc = base(run("Xylaria-G536", "Xylaria cubensis G536", 2512241, "Ascomycota", contigs=3120, length=56_100_200, n50=38_400),
          families_attempted=["Ascomycota:MAT"],
          suppressed_unpolished=2, suppressed_loci=[
              {"family": "Ascomycota:MAT", "contig": "VFLP01000025.1", "start": 244382, "end": 252058, "idiomorph": None, "polished_genes": 1, "genes_found": ["SLA2", "APN2"], "withheld_reason": "mat_gene_gate", "idiomorph_classifier": None},
              {"family": "Ascomycota:MAT", "contig": "VFLP01000310.1", "start": 1880, "end": 2410, "idiomorph": "MAT1-2", "polished_genes": 0, "genes_found": ["MAT1-2-1"], "withheld_reason": "modelled_gene_bar", "idiomorph_classifier": None}],
          assembly_gap_at_locus=[{"family": "Ascomycota:MAT", "contig": "VFLP01000025.1", "start": 248000, "end": 249600, "n_bases": 1600, "anchors": ["SLA2", "APN2"]}],
          not_detected=[{"family": "Ascomycota:MAT", "reason": "best cluster carried 1 modelled gene(s), below the 2 required to report a locus; modelled: SLA2", "best_fraction_found": 0.5, "genes_found": ["SLA2", "APN2"], "genes_missing": ["MAT1-1-1", "MAT1-2-1"], "genes_not_searchable": []}])
(out / "none_called_gap.yaml").write_text(yaml.safe_dump(nc, sort_keys=False))

# 5. Basidiomycota: merged HD call with subloci, and an unverified PR call in an array.
c = "scaffold_3"
hd = [gene("HD1", "core_MAT", c, 410200, 412300, "-", 64.0, rec="5334_h48_HD_A43", exons=[{"start": 410200, "end": 410900}, {"start": 410960, "end": 412300}]),
      gene("HD2", "core_MAT", c, 412900, 415100, "+", 58.2, rec="5334_h48_HD_A43"),
      gene("MIP", "flanking_conserved", c, 405000, 407800, "+", 81.0, rec="5334_h48_HD_A43"),
      gene("beta-fg", "flanking_variable", c, 417000, 418400, "-", 70.2, rec="5334_h48_HD_A43", status="unpolished", method="tblastn_genome", cov=61.0)]
pr = [gene("STE3.1", "core_MAT", "scaffold_9", 88000, 89400, "+", 47.0, rec="5334_h48_PR_B43"),
      {**gene("phb1", "core_MAT", "scaffold_9", 90100, 90220, "-", None, status="not_polish_candidate", method="caax_scan", rec=None), "orf_length_aa": 40, "caax_motif": "CVIA", "orf_count": 1}]
bas = base(run("Schco-H4-8", "Schizophyllum commune H4-8", 5334, "Basidiomycota", contigs=36, length=38_500_000, n50=3_100_000),
           families_attempted=["Basidiomycota:HD", "Basidiomycota:PR"],
           detected=[call("Basidiomycota:HD", c, 405000, 418400, "A43-like", "medium", hd,
                          subloci=[{"sublocus": "Aalpha", "generic": False, "idiomorph": "Aa1", "contig": c, "start": 410200, "end": 415100, "genes": ["HD1", "HD2"], "genes_missing": [], "completeness": "complete"},
                                   {"sublocus": "Abeta", "generic": False, "idiomorph": "Ab3", "contig": c, "start": 416000, "end": 418400, "genes": ["beta-fg"], "genes_missing": ["Abeta HD1"], "completeness": "partial"}],
                          genes_missing=["Abeta HD1"], span_exceeds_plausible_bound=True),
                     call("Basidiomycota:PR", "scaffold_9", 88000, 90220, "B43-like", "low", pr, cls="mat_locus",
                          verification={"status": "unverified", "reason": "admitted only through a strict-CAAX scan precursor; the scan's false-positive rate is not yet bounded"},
                          receptor_array_id="RA1", receptor_array_size=3, receptor_array_members=["scaffold_9:88000-89400(+)"], receptor_array_support="supported", receptor_array_support_reasons=["array_size>=2"],
                          genes_not_searchable=["pheromone_receptor"])],
           receptor_arrays=[{"receptor_array_id": "RA1", "family": "Basidiomycota:PR", "contig": "scaffold_9", "start": 61000, "end": 89400, "receptor_array_size": 3, "receptor_array_members": [], "strict_caax_orfs": 2, "precursor_homology_hits": 0, "receptor_array_support": "supported", "receptor_array_support_reasons": ["array_size>=2", "strict_caax_orfs=2"], "calls": 1}])
(out / "basidio_hd_pr.yaml").write_text(yaml.safe_dump(bas, sort_keys=False))

# 6. Hostile strings (escaping).
ev = copy.deepcopy(phy)
ev["run"]["sample"] = "<script>alert(1)</script>"
ev["detected"][0]["contig"] = 'ctg"><img src=x onerror=alert(1)>'
for g in ev["detected"][0]["gene_evidence"]:
    g["contig"] = ev["detected"][0]["contig"]
ev["detected"][0]["segments"][0]["contig"] = ev["detected"][0]["contig"]
(out / "hostile_strings.yaml").write_text(yaml.safe_dump(ev, sort_keys=False))
