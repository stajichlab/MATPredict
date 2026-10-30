"""Draft MATss proposals for S. schoenii (GDD12, GCA_054953955.1) from our own chromosome IV gene models.

Models: results/2026-10-29_ss_records (StringTie + TransDecoder on PRJNA1364152 RNA-seq; the authors' GTF is not public).
Each CDS is re-extracted from the window fasta and checked (ATG start, no internal stop, stop at the end, length) before writing.
Run from the repo root: python3 results/2026-10-30_saccharomycopsis/make_records.py [outdir]
"""
import sys, yaml
from Bio.Seq import Seq

HERE = "results/2026-10-30_saccharomycopsis/"
OUT = sys.argv[1] if len(sys.argv) > 1 else HERE
WIN_START = 463000
seq = "".join(l.strip() for l in open(HERE + "chrIV_463k-542k.fa") if not l.startswith(">"))

def cds(exons, strand):
    s = "".join(seq[a - WIN_START:b - WIN_START + 1] for a, b in sorted(exons))
    s = s if strand == "+" else str(Seq(s).reverse_complement())
    return s

def protein(exons, strand):
    return str(Seq(cds(exons, strand)).translate(table=12))   # CUG = Ser (NCBI table 12)

# name: (role, strand, exons ascending, expected aa, note)
LEFT = {
 "MATalpha2": ("core_MAT", "+", [(470349, 470951)], 200, "alpha2 profile 73 bits (a1 12); a 92 bp intron (97 reads) gives a 227 aa isoform"),
 "MATalpha1": ("core_MAT", "-", [(471020, 471595)], 191, "alpha box (MATalpha_HMGbox 54 bits)"),
 "MATa1": ("core_MAT", "+", [(471961, 472227), (472445, 472489)], 103, "a1 profile 86 bits (alpha2 29); partial: 217 bp intron (198 reads) joins it toward SLA2"),
 "sla2": ("flanking_conserved", "-", [(472708, 473478)], 256, "diverged SLA2 (I_LWEQ 186 bits); a second piece 473,680-474,606 (309 aa) is not recorded"),
}
RIGHT = {
 "sla2": ("flanking_conserved", "+", [(533981, 535828)], 615, "SLA2 anchor (I_LWEQ 183 bits), SS.2002.2 ORF p1; mirrored copy 529,992-532,490"),
 "MATa1": ("core_MAT", "-", [(536046, 536090), (536296, 536634), (536693, 536944)], 211, "a1 profile 85 bits (alpha2 27); introns 205 bp (221 reads) and 58 bp (218 reads)"),
 "MATa2": ("core_MAT", "+", [(537379, 538170)], 263, "a2 profile hit, head-to-head with MATa1; ATG to stop, single exon; mirror copy 527,654-528,454 (minus)"),
}
for nm, genes in (("left", LEFT), ("right", RIGHT)):
    for g, (role, strand, ex, aa, _) in genes.items():
        p = protein(ex, strand)
        print(nm, g, len(p) - (1 if p.endswith("*") else 0), "aa expected", aa, "start", p[0], "stop", p.endswith("*"), "internal stops", p[:-1].count("*"))
        assert p.endswith("*") and p[:-1].count("*") == 0 and len(p) - 1 == aa, (g, len(p))

CIT = [{"pmid": "41849659", "doi": "10.1093/g3journal/jkag067"}, {"pmid": "38492364", "doi": "10.1016/j.micres.2024.127691"}]
TAX = {"taxid": 44080, "lineage": "k__Fungi;p__Ascomycota;subphylum__Saccharomycotina;c__Saccharomycetes;o__Ascoideales;f__Saccharomycopsidaceae;g__Saccharomycopsis;s__Saccharomycopsis_schoenii",
       "lineage_resolved_date": "2026-10-30"}
ORG = {"species": "Saccharomycopsis schoenii", "strain": {"name": "GDD12", "known": True, "culture_collection_ids": [], "differs_from_sequenced": False}}

def genes_of(spec):
    out = []
    for i, (name, (role, strand, ex, aa, note)) in enumerate(sorted(spec.items(), key=lambda kv: min(a for a, b in kv[1][2]))):
        exons = [{"start": a, "end": b} for a, b in (sorted(ex) if strand == "+" else sorted(ex, reverse=True))]
        out.append({"gene_index": i, "name": name, "role": role, "present": True, "locus_tag": None, "segment_index": 0,
                    "start": min(a for a, b in ex), "end": max(b for a, b in ex), "strand": strand, "order_in_locus": i,
                    "transl_table": 12, "exons": exons})
    return out

def record(rid, idiomorph, spec, orient, note, span):
    genes = genes_of(spec)
    lo = min(g["start"] for g in genes); hi = max(g["end"] for g in genes)
    loc = {"coordinate_provenance": "curator_derived", "excluded_from_coordinate_benchmark": False,
           "core": {"completeness": "partial", "reference_orientation": orient, "definition_note": note,
                    "segments": [{"segment_index": 0, "sequence_source": {"type": "insdc_nucleotide", "accession": "CM143129.1", "seq_region": "CM143129.1"},
                                  "start": lo, "end": hi, "contig_edge_distance": None, "sequence_checksum": None}]},
           "extended_flank": [], "assembly_accession": "GCA_054953955.1", "assembly_accession_basis": "strain_match_verified",
           "assembly_location": {"contig": "CM143129.1", "start": lo, "end": hi, "identity": 100.0}}
    return {"record_id": rid, "record_version": 1, "taxonomy": TAX, "organism": ORG,
            "mating_type": {"locus_name": "MATss", "idiomorphs": idiomorph if isinstance(idiomorph, list) else [idiomorph], "system": "homothallic"},
            "locus": loc, "genes": genes,
            "evidence": {"locus_existence": {"tier": 2, "experimental_method": "published description of the SLA2-MAT cassettes and RNA-seq of the active chromosome IV locus (Kriti et al. 2026)", "citations": CIT},
                         "boundaries": {"tier": 3, "experimental_method": "OUR OWN gene models from the public RNA-seq (PRJNA1364152; StringTie, TransDecoder, phmmer/Pfam) because the authors' GTF is not public; coordinates are ours, the locus span is the paper's", "citations": CIT},
                         "idiomorph_assignment": {"tier": 3, "experimental_method": "new MATa1, MATalpha2 and MATa2 profiles (a score difference, supporting evidence: results/2026-10-30_mat_hmms); the sub-clusters are the paper's", "citations": CIT}},
            "validation": {"status": "needs_review", "rejection_reason": None, "accession_resolved": False, "accession_resolved_date": None,
                           "accession_resolved_version": None, "sequence_match": {"status": "pass", "per_gene": [], "notes": ""}, "taxonomy_current": None},
            "curation": {"proposed_by": "claude-code-curation-agent", "proposal_dedupe_key": f"41849659|44080|MATss|{'-'.join(idiomorph) if isinstance(idiomorph, list) else idiomorph}", "reviewed_by": None, "reviewed_date": None},
            "model_provenance": None}

NOTE_L = ("Left (alpha-dominant) sub-cluster of the active MAT locus on chromosome IV of S. schoenii GDD12 (CM143129.1; the paper's span 470,610-475,000, "
          "Kriti et al. G3 2026, PMID 41849659). HOMOTHALLIC (Dost et al. 2024, PMID 38492364); one of several SLA2-MAT cassettes (one per chromosome plus a palindromic one on chromosome IV). "
          "GENE MODELS ARE OURS (tier 3): the BRAKER2 GTF the paper cites is not on its GitHub or Zenodo, so the models were rebuilt from the three public paired RNA-seq runs "
          "(PRJNA1364152; hisat2, StringTie, TransDecoder, phmmer against the curated MAT proteins, Pfam-A); see results/2026-10-29_ss_records/README.md. The library is not stranded; "
          "strand is the ORF's. GENE NAMES come from new MATa1, MATalpha2 and MATa2 profiles (results/2026-10-30_mat_hmms; supporting evidence, a score difference, not a gate): the ORF at 470,349-470,951 "
          "(200 aa, single exon, 227 aa with a 92 bp intron of 97 reads) is MATalpha2 (alpha2 profile 73 bits, a1 12) and not an HMG-box MATa2; the homeodomain at 471,961-472,489 is MATa1 (a1 86, alpha2 29), "
          "so the cluster holds both an alpha pair and an a1 beside SLA2 (the paper's MATa1-SLA2 junction). MATalpha1 is single-exon, minus strand, 191 aa. The MATa1 model is partial (103 aa, two exons with a 217 bp intron of "
          "198 reads that joins it toward SLA2); the uninterrupted frame through the intron scores lower (a1 60) and runs through a low-complexity repeat. "
          "SLA2 is a diverged copy: 256 aa (472,708-473,478, minus), 5-prime partial (no start codon); a second SLA2 piece (473,680-474,606, 309 aa) is not recorded. GENETIC CODE: NCBI table 12 (CUG = Ser).")
NOTE_R = ("Right (palindromic) cassette of the active MAT locus on chromosome IV of S. schoenii GDD12 (CM143129.1; the paper's span 528,044-537,777; an SLA2 anchor between mirrored MAT arrays, "
          "Kriti et al. G3 2026, PMID 41849659). HOMOTHALLIC (Dost et al. 2024, PMID 38492364). The record covers the arm that RNA-seq resolves: SLA2 533,981-535,828 (plus, 615 aa, 5-prime partial: no start codon in the ORF), "
          "MATa1 at 536,046-536,944 (minus, 211 aa; a1 profile 85 bits, alpha2 27; introns of 205 bp, 221 reads, and 58 bp, 218 reads, the conserved 58 bp junction the paper reports at 536,635) and MATa2 at 537,379-538,170 "
          "(plus, 263 aa, ATG to stop, single exon, head-to-head with MATa1; a2 profile hit, results/2026-10-30_mat_hmms). The a2 profile finds only three a2-like ORFs in the whole genome, all on chromosome IV: this one, its mirror "
          "at 527,654-528,454 (minus) and one at 502,064-502,864 at the central cluster. The mirrored arm (SLA2 529,992-532,490; MATa1 529,186-529,746) is not recorded: short reads cannot separate the two inverted copies, "
          "and the paper's mirrored junction at 529,128 is not seen. GENE MODELS ARE OURS (tier 3), see results/2026-10-29_ss_records/README.md. The library is not stranded. GENETIC CODE: NCBI table 12 (CUG = Ser).")
ra = record("44080_gdd12_MATss_MATa-MATalpha", ["MATa", "MATalpha"], LEFT, "MATalpha2->MATalpha1->MATa1->sla2", NOTE_L, None)
rb = record("44080_gdd12_MATss_MATa", "MATa", RIGHT, "sla2->MATa1->MATa2", NOTE_R, None)
for r in (ra, rb):
    open(f"{OUT}{r['record_id']}.yaml", "w").write(yaml.safe_dump(r, sort_keys=False))
    print(r["record_id"], [(g["name"], g["start"], g["end"], g["strand"], len(g["exons"])) for g in r["genes"]])
