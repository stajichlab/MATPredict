"""Write the two draft MATpf proposals for Petasospora (Cyberlindnera) fabianii into this directory (NOT db/).

MATa    allele: JOY008 (St. Louis College of Pharmacy; ecotype Missouri), complete genome GCA_022641835.1, CP048715.1.
MATalpha allele: TBRC 4498, chromosome-level GCA_050575385.1, CM115473.1 (the strain-65 draft, GCA_001983305.1, carries the same locus and
                 the depositor's ONH67106 model; JCM 3601 GCA_001599195.1 carries the MATa locus again).
Gene models are our own: AUGUSTUS ab initio (no RNA-seq exists for the species) selected by agreement of several yeast species models and by
homology/profile evidence; coordinates include the stop codon. Run from the repo root:
    python3 results/2026-10-30_petasospora/make_records.py
"""
import yaml

HERE = "results/2026-10-30_petasospora/"
# (gene, role, strand, exons in genomic order, stop-codon-excluded as AUGUSTUS gives them)
A = ("CP048715.1", [("dic1", "flanking_conserved", "-", [(628452, 629375)]),
                    ("MATa2", "core_MAT", "-", [(629853, 630419), (630492, 630863)]),
                    ("MATa1", "core_MAT", "+", [(631437, 631737), (631835, 632133)]),
                    ("sla2", "flanking_conserved", "+", [(641746, 644868)])])
AL = ("CM115473.1", [("dic1", "flanking_conserved", "-", [(631649, 632572)]),
                     ("MATalpha1", "core_MAT", "+", [(637540, 638454)]),
                     ("sla2", "flanking_conserved", "+", [(648035, 651157)])])

def gene(idx, name, role, strand, ex):
    ex = sorted(ex)
    if strand == "+": ex[-1] = (ex[-1][0], ex[-1][1] + 3)       # add the stop codon
    else: ex[0] = (ex[0][0] - 3, ex[0][1])
    tx = ex if strand == "+" else ex[::-1]
    return {"gene_index": idx, "name": name, "protein_accession": None, "role": role, "present": True, "locus_tag": None,
            "segment_index": 0, "start": min(a for a, b in ex), "end": max(b for a, b in ex), "strand": strand, "order_in_locus": idx,
            "exons": [{"start": a, "end": b} for a, b in tx]}

TAX = {"taxid": 36022, "lineage": "k__Fungi;p__Ascomycota;subphylum__Saccharomycotina;c__Saccharomycetes;o__Phaffomycetales;f__Wickerhamomycetaceae;g__Petasospora;s__Petasospora_fabianii",
       "lineage_resolved_date": "2026-10-30"}
CIT_65 = [{"pmid": "28385833", "doi": "10.1128/genomeA.00064-17"}]
CIT_TBRC = [{"pmid": "40558965", "doi": "10.3390/jof11060453"}]
COMMON = ("GENE MODELS ARE OURS (tier 3, manual annotation; curator decision 2026-10-30). No RNA-seq exists for the species in ENA/SRA (checked 2026-10-30), "
          "so the models are AUGUSTUS ab initio predictions (species models S. cerevisiae S288C, C. guilliermondii, P. stipitis, K. lactis) chosen where several "
          "agree and where homology supports them; MATa2 is also supported by the Saccharomycotina MATa2 profile (76.8 bits, above the 52.4 bits of P. kudriavzevii). "
          "Coordinates include the stop codon. Intron signals: MATa2 72 bp GTAGGT..TACTAAC..TAG (strong); MATa1 97 bp GTACCA..AAAG (weak donor; an upstream exon is required because the homeodomain open frame has no start codon, but the "
          "splice sites are unconfirmed). The MATalpha1 model is a single exon of 305 aa; the depositor's ONH67106 (254 aa, strain 65) is the same C-terminal 185 aa plus a "
          "69 aa exon joined across a 1,263 bp intron (GTCCTT..CCCCAG) that has no support here: the N-terminus is unresolved. GENE NAMES by homology and position: "
          "MATa1 is the homeodomain protein (Pfam Homeodomain 64 bits, 66 bits to Komagataella MATa1), MATa2 the HMG-box protein (divergent; no Pfam HMG_box hit above GA), "
          "MATalpha1 the alpha-box protein (PF04769 36 bits, 407 bits to ONH67106). No MATalpha2 homeodomain gene was found in either alpha allele (only a 57 aa Homeobox_KN fragment at "
          "12-18 bits about 1.2 kb upstream of MATalpha1). LAYOUT: DIC1 (minus) - 0.5 kb (MATa) or 5.0 kb (MATalpha) - MAT genes - six conserved non-MAT ORFs (409, 173, 589, 219, 465, 622 aa; identical in both alleles, not recorded) - SLA2 (plus); "
          "SLA2 is therefore about 9.5 kb from the MAT genes. The alpha allele also carries an 800 aa MutS_III-domain ORF between DIC1 and MATalpha1 that has no counterpart in the MATa allele (not recorded). "
          "No second MAT cassette: miniprot of the MAT proteins finds no other MAT gene in any of the four assemblies, so the species looks heterothallic with one idiomorph per genome. "
          "GENETIC CODE: standard (table 1). Details: results/2026-10-30_petasospora/README.md.")

def record(rid, idiomorph, contig, genes, strain, acc, cit, note, defn):
    lo = min(g["start"] for g in genes); hi = max(g["end"] for g in genes)
    return {"record_id": rid, "record_version": 1, "taxonomy": TAX,
            "organism": {"species": "Petasospora fabianii", "strain": {"name": strain[0], "known": True, "culture_collection_ids": strain[1], "differs_from_sequenced": False}},
            "mating_type": {"locus_name": "MATpf", "idiomorphs": [idiomorph], "system": "heterothallic"},
            "locus": {"coordinate_provenance": "curator_derived", "excluded_from_coordinate_benchmark": False,
                      "assembly_accession": acc, "assembly_accession_basis": "sequence_link",
                      "assembly_location": {"contig": contig, "start": lo, "end": hi, "identity": 100.0},
                      "core": {"completeness": "complete", "reference_orientation": defn,
                               "definition_note": note + " " + COMMON,
                               "segments": [{"segment_index": 0, "sequence_source": {"type": "insdc_nucleotide", "accession": contig, "seq_region": contig},
                                             "start": lo, "end": hi, "contig_edge_distance": None, "sequence_checksum": None}]},
                      "extended_flank": []},
            "genes": genes,
            "evidence": {"locus_existence": {"tier": 3, "experimental_method": "manual annotation: homology (curated Saccharomycotina MAT proteins, Pfam, MATa2 profile) and AUGUSTUS ab initio models on the deposited assembly; no RNA-seq", "citations": cit},
                         "boundaries": {"tier": 3, "experimental_method": "AUGUSTUS ab initio gene models checked against homology; DIC1 and SLA2 are the flanking anchors, MAT genes lie between them", "citations": cit},
                         "idiomorph_assignment": {"tier": 3, "experimental_method": "gene content: MATa1 (homeodomain) with MATa2 (HMG-box) versus MATalpha1 (alpha box), at the same DIC1-SLA2 position in different strains", "citations": cit}},
            "validation": {"status": "needs_review", "rejection_reason": None, "accession_resolved": False, "accession_resolved_date": None,
                           "accession_resolved_version": None, "sequence_match": {"status": "pass", "per_gene": [], "notes": ""}, "taxonomy_current": None},
            "curation": {"proposed_by": "claude-code-curation-agent", "proposal_dedupe_key": f"manual|36022|MATpf|{idiomorph}", "reviewed_by": None, "reviewed_date": None},
            "model_provenance": None}

ga = [gene(i, n, r, s, e) for i, (n, r, s, e) in enumerate(A[1])]
gl = [gene(i, n, r, s, e) for i, (n, r, s, e) in enumerate(AL[1])]
NOTE_A = ("MATa idiomorph of P. fabianii isolate JOY008 (ecotype Missouri, St. Louis College of Pharmacy), complete genome GCA_022641835.1, CP048715.1 (BioSample SAMN14051719). "
          "The type strain JCM 3601 (GCA_001599195.1, BCGI01000001.1:641,150-657,563) carries the same MATa locus with 99-100% identical proteins (MATa1 and MATa2 one substitution each). "
          "The v0.7.0 MATsc call (CP048715.1:630,345-685,843, 55 kb) is wider than the locus: the genes are DIC1 628,452-629,375 to SLA2 641,746-644,871.")
NOTE_AL = ("MATalpha idiomorph of P. fabianii TBRC 4498 (xylitol-producing strain; PMID 40558965), chromosome-level GCA_050575385.1, CM115473.1. The draft genome of strain 65 from fermented masau fruit "
           "(GCA_001983305.1, MPUK01000005.1:611,410-630,920; PMID 28385833) carries the same locus with an identical MATalpha1 protein and holds the depositor's ONH67106 "
           "('Mating-type-like protein ALPHA1', 254 aa, MPUK01000005.1:616,190-616,396 + 617,660-618,217), which our model replaces.")
recs = [record("36022_joy008_MATpf_MATa", "MATa", A[0], ga, ("JOY008", ["JOY008"]), "GCA_022641835.1", [], NOTE_A, "DIC1->MATa2->MATa1->(six conserved genes)->SLA2"),
        record("36022_tbrc-4498_MATpf_MATalpha", "MATalpha", AL[0], gl, ("TBRC 4498", ["TBRC 4498"]), "GCA_050575385.1", CIT_TBRC + CIT_65, NOTE_AL, "DIC1->MATalpha1->(six conserved genes)->SLA2")]
for r in recs:
    open(f"{HERE}{r['record_id']}.yaml", "w").write(yaml.safe_dump(r, sort_keys=False, width=120))
    print(r["record_id"], [(g["name"], g["start"], g["end"], g["strand"], len(g["exons"])) for g in r["genes"]])
