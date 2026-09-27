"""Aggregate the full Basidiomycota run into one per-genome table.

Reads every chunk's runs/*/detection_report.yaml, evidence_diagnostics.jsonl
and wall_seconds; joins BFD samples.csv (taxonomy) and tables/asm_stats.parquet
(assembly size, contigs, N50). Writes genomes.tsv (one row per genome) and
loci.tsv (one row per reported locus). Run with /usr/bin/python3.12.
"""
import csv, glob, json, os, sys
import yaml
import pyarrow.parquet as pq

HERE = os.path.dirname(os.path.abspath(__file__))
SAMPLES = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"
ASM = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/tables/asm_stats.parquet"

samp = {r["ASMID"]: r for r in csv.DictReader(open(SAMPLES))}
# samples.csv leaves SUBPHYLUM blank for most Basidiomycota; derive it from CLASS.
SUBPHYLUM_OF_CLASS = {
    "Agaricomycetes": "Agaricomycotina", "Tremellomycetes": "Agaricomycotina", "Dacrymycetes": "Agaricomycotina",
    "Wallemiomycetes": "Wallemiomycotina", "Bartheletiomycetes": "Agaricomycotina",
    "Ustilaginomycetes": "Ustilaginomycotina", "Exobasidiomycetes": "Ustilaginomycotina",
    "Malasseziomycetes": "Ustilaginomycotina", "Moniliellomycetes": "Ustilaginomycotina",
    "Pucciniomycetes": "Pucciniomycotina", "Microbotryomycetes": "Pucciniomycotina",
    "Cystobasidiomycetes": "Pucciniomycotina", "Agaricostilbomycetes": "Pucciniomycotina",
    "Atractiellomycetes": "Pucciniomycotina", "Classiculomycetes": "Pucciniomycotina",
    "Mixiomycetes": "Pucciniomycotina", "Tritirachiomycetes": "Pucciniomycotina",
    "Spiculogloeomycetes": "Pucciniomycotina", "Cryptomycocolacomycetes": "Pucciniomycotina",
    "Entorrhizomycetes": "Entorrhizomycota", "Peribolosporomycetes": "Ustilaginomycotina",
    "Geminibasidiomycetes": "Agaricomycotina",
}
asm = {r["ASMID"]: r for r in pq.read_table(
    ASM, columns=["ASMID", "contig_count", "total_length_bp", "N50_bp", "total_n_bases"]).to_pylist()}

# genome -> report dir; a later chunk (e.g. Pucciniales_long) overrides an earlier one
dirs = {}
for chunk in sorted(glob.glob(os.path.join(HERE, "*", "runs"))):
    for d in glob.glob(os.path.join(chunk, "*")):
        a = os.path.basename(d)
        if os.path.exists(os.path.join(d, "detection_report.yaml")) and os.path.getsize(os.path.join(d, "detection_report.yaml")):
            if a not in dirs or dirs[a][1] is None or "Pucciniales_long" in chunk:
                dirs[a] = (os.path.basename(os.path.dirname(chunk)), d)
        elif a not in dirs:
            dirs[a] = (os.path.basename(os.path.dirname(chunk)), None)

grows, lrows = [], []
for a, (chunk, d) in sorted(dirs.items()):
    s = samp.get(a, {})
    st = asm.get(a, {})
    g = dict(genome=a, chunk=chunk, species=s.get("SPECIES", ""),
             subphylum=s.get("SUBPHYLUM", "") or SUBPHYLUM_OF_CLASS.get(s.get("CLASS", ""), "?"),
             class_=s.get("CLASS", ""), order=s.get("ORDER", ""), family=s.get("FAMILY", ""),
             size_bp=st.get("total_length_bp", ""), contigs=st.get("contig_count", ""),
             n50=st.get("N50_bp", ""), n_bases=st.get("total_n_bases", ""))
    if d is None:
        g.update(status="no_report")
        grows.append(g)
        continue
    try:
        g["wall_s"] = int(open(os.path.join(d, "wall_seconds")).read().strip())
    except Exception:
        g["wall_s"] = ""
    r = yaml.load(open(os.path.join(d, "detection_report.yaml")), Loader=yaml.CSafeLoader) or {}
    det = r.get("detected") or []
    g.update(status="called" if det else "uncalled", routing=r.get("routing_mode", ""),
             not_searched_reason=r.get("not_searched_reason") or "",
             routing_error=bool(r.get("routing_error")),
             n_loci=len(det), families_called="|".join(sorted({x["family"] for x in det})),
             n_suppressed=len(r.get("suppressed_loci") or []),
             suppressed_flank_carried=r.get("suppressed_flank_carried") or 0,
             gap_at_locus=len(r.get("assembly_gap_at_locus") or []),
             zygosity=(r.get("zygosity") or {}).get("status", "") if isinstance(r.get("zygosity"), dict) else (r.get("zygosity") or ""),
             homothallic=sum(1 for x in det if x.get("locus_class") == "homothallic_candidate"),
             unverified=sum(1 for x in det if x.get("verification") == "unverified"),
             idiomorph_unmodelled=sum(1 for x in det if x.get("idiomorph_unmodelled")))
    # withheld-locus shape: best suppressed locus's genes and polished count
    sup = r.get("suppressed_loci") or []
    g["sup_max_polished"] = max((x.get("polished_genes") or 0 for x in sup), default="")
    g["sup_max_genes"] = max((len(x.get("genes_found") or []) for x in sup), default="")
    nd = r.get("not_detected") or []
    g["nd_best_fraction"] = max((x.get("best_fraction_found") or 0 for x in nd), default="")
    g["nd_reasons"] = "|".join(sorted({(x.get("reason") or "")[:60] for x in nd}))
    # diagnostics: capped clusters, admitted clusters
    capped = admitted = 0
    ev = os.path.join(d, "evidence_diagnostics.jsonl")
    if os.path.exists(ev):
        for line in open(ev):
            try:
                e = json.loads(line)
            except Exception:
                continue
            if e.get("kind") == "evidence":
                admitted += bool(e.get("admitted"))
                capped += bool(e.get("polish_capped"))
    g["admitted_clusters"], g["capped_clusters"] = admitted, capped
    grows.append(g)
    for x in det:
        lrows.append(dict(genome=a, order=g["order"], subphylum=g["subphylum"], family_called=x["family"],
                          contig=x["contig"], start=x["start"], end=x["end"], idiomorph=x.get("idiomorph"),
                          confidence=x.get("confidence"), locus_class=x.get("locus_class"),
                          detection_pass=x.get("detection_pass"), polished_genes=x.get("polished_genes"),
                          genes_found="|".join(x.get("genes_found") or []),
                          verification=x.get("verification") or "",
                          idiomorph_unmodelled=bool(x.get("idiomorph_unmodelled"))))

def write(path, rows):
    keys = []
    for row in rows:
        for k in row:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=keys, delimiter="\t", restval="")
        w.writeheader(); w.writerows(rows)

write(os.path.join(HERE, "genomes.tsv"), grows)
write(os.path.join(HERE, "loci.tsv"), lrows)
print(len(grows), "genomes;", len(lrows), "loci")
