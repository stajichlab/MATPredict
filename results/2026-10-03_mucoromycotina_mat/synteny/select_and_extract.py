#!/usr/bin/env python3
"""Pick 28 loci (one Plus + one Minus for 14 genera) and write one GenBank file each.

Selection (calls.tsv): status called, locus_class mat_locus, flank genes on both
sides, confidence high, no taxon override. Per genus, prefer one species with both
idiomorphs (largest smaller-margin); else the largest-margin locus per idiomorph.
Genera are ordered by the nf_phyling species tree (rooted on Umbelopsis).

Region: locus span (flank to flank) +- PAD bp. Gene models: the genome annotation
(LCG/Jena funannotate .gbk; BFD .gff3 + the same assembly). Each annotation CDS that
overlaps a MATPredict gene model on the same strand by >= 50% takes the MATPredict
name; a MATPredict gene with no overlapping annotation CDS is added from its model
(small MAT genes are often missed by annotation). Every CDS gets a unique locus_tag
<sample>_<n> so clinker can group genes by function.
"""
import csv, gzip, os, re, sys
from Bio import SeqIO, Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
import yaml

C = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-03_mucoromycotina_mat"
O = f"{C}/synteny"
PAD = 3000
DROP = {"Helicostylum", "Parasitella"}
LCG_ANN = "/bigdata/stajichlab/shared/projects/ZyGoLife/LCG/Annotation/annotate/{g}/annotate_results/{g}.gbk"
JENA_ANN = "/bigdata/stajichlab/shared/projects/ZyGoLife/Mucor_Jena/annotation/{g}/annotate_results"
BFD_GFF = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/input/gff3/{n}.gff3"
BFD_DNA = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes/{a}.fa.gz"
SPTREE = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-01_lcg_name_check/nf_out/protein/buildtree/mucoromycota_odb12/fasttree/protein-lcg_namecheck_v1-taxa_941.mucoromycota_odb12.fasttree.support.treefile"

def fnum(x):
    try: return float(x)
    except Exception: return -1e9

rows = [r for r in csv.DictReader(open(f"{C}/calls.tsv"), delimiter="\t")
        if r["status"] == "called" and r["locus_class"] == "mat_locus" and r["flank_status"] == "both_sides"
        and r["confidence"] == "high" and not r["override"]]
SAMPLES = {x["ASMID"]: x for x in csv.DictReader(open("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"))}

def bfd_gff(asmid):
    samp = SAMPLES[asmid]
    base = samp["SPECIES_IN"]
    if samp["STRAIN"] and samp["STRAIN"] not in base:
        base += " " + samp["STRAIN"]
    p = BFD_GFF.format(n=re.sub(r"[ /]+", "_", base))
    return p if os.path.exists(p) else None

def ann_path(r):
    if r["source"] == "LCG":
        p = LCG_ANN.format(g=r["genome"])
        return p if os.path.exists(p) else None
    if r["source"] == "JENA":
        d = JENA_ANN.format(g=r["genome"])
        if os.path.isdir(d):
            for x in os.listdir(d):
                if x.endswith(".gbk"):
                    return os.path.join(d, x)
        return None
    return bfd_gff(r["genome"])

rows = [r for r in rows if ann_path(r)]
for r in rows:
    r["annotated_rank"] = 0 if r["source"] in ("LCG", "JENA") else 1
    r["genus"] = r["name"].split()[0]
    r["species"] = " ".join(r["name"].split()[:2])
by_genus = {}
for r in rows:
    by_genus.setdefault(r["genus"], []).append(r)
picks = []
for g, rs in by_genus.items():
    if g in DROP or {x["idiomorph"] for x in rs} < {"Plus", "Minus"}:
        continue
    best = None
    for sp in {x["species"] for x in rs}:
        P = [x for x in rs if x["species"] == sp and x["idiomorph"] == "Plus"]
        M = [x for x in rs if x["species"] == sp and x["idiomorph"] == "Minus"]
        if P and M and "sp." not in sp:
            key = lambda x: (-x["annotated_rank"], fnum(x["margin"]))
            p, m = max(P, key=key), max(M, key=key)
            score = min(fnum(p["margin"]), fnum(m["margin"]))
            if best is None or score > best[0]:
                best = (score, p, m)
    if best:
        picks += [best[1], best[2]]
    else:
        for idio in ("Plus", "Minus"):
            picks.append(max([x for x in rs if x["idiomorph"] == idio], key=lambda x: (-x["annotated_rank"], fnum(x["margin"]))))

# order genera by the species tree
tree = Phylo.read(SPTREE, "newick")
umb = [t for t in tree.get_terminals() if "Umbelopsis" in t.name]
tree.root_with_outgroup(tree.common_ancestor(umb))
tree.ladderize()
order = {}
for i, t in enumerate(tree.get_terminals()):
    g = re.sub(r"^(LCG|JENA|BFD)__", "", t.name).split("_")[0]
    order.setdefault(g, i)
picks.sort(key=lambda r: (order.get(r["genus"], 10**6), r["idiomorph"] != "Plus"))

def mp_genes(r):
    d = yaml.safe_load(open(f"{C}/runs/{r['source']}__{r['genome']}/detection_report.yaml"))
    loc = [x for x in d["detected"] if x["contig"] == r["contig"] and str(x["start"]) == r["locus_start"]][0]
    out = []
    for e in loc["gene_evidence"]:
        if e["contig"] != r["contig"] or not str(e.get("status", "")).startswith("polished"):
            continue
        # sexM and sexP share an HMG box, so both references model the same core
        # gene; keep only the called idiomorph's core gene.
        if e["gene"] in ("sexP", "sexM") and e["gene"] != {"Plus": "sexP", "Minus": "sexM"}[r["idiomorph"]]:
            continue
        out.append(e)
    return out

def bfd_record(r, start, end):
    gff = bfd_gff(r["genome"])
    seq = None
    with gzip.open(BFD_DNA.format(a=r["genome"]), "rt") as fh:
        for rec in SeqIO.parse(fh, "fasta"):
            if rec.id == r["contig"]:
                seq = rec.seq; break
    cds = {}
    for l in open(gff):
        f = l.rstrip("\n").split("\t")
        if len(f) < 9 or f[0] != r["contig"] or f[2] != "CDS":
            continue
        a = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
        cds.setdefault(a.get("Parent", a.get("ID")), []).append((int(f[3]), int(f[4]), f[6]))
    feats = []
    for pid, ex in cds.items():
        s0, e0 = min(x[0] for x in ex), max(x[1] for x in ex)
        if e0 < start or s0 > end:
            continue
        strand = 1 if ex[0][2] == "+" else -1
        parts = [FeatureLocation(a - 1, b, strand=strand) for a, b, _ in sorted(ex, reverse=strand < 0)]
        loc = parts[0] if len(parts) == 1 else CompoundLocation(parts)
        feats.append(SeqFeature(loc, type="CDS", qualifiers={"note": [pid]}))
    return SeqRecord(seq, id=r["contig"]), feats, gff

def ann_record(r):
    path = ann_path(r)
    for rec in SeqIO.parse(path, "genbank"):
        if rec.id == r["contig"] or rec.name == r["contig"]:
            return rec, [f for f in rec.features if f.type == "CDS"], path
    raise KeyError(f"{r['contig']} not in {path}")

def overlap(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0) + 1)

summary = []
for i, r in enumerate(picks):
    start, end = max(1, int(r["locus_start"]) - PAD), int(r["locus_end"]) + PAD
    rec, feats, src = bfd_record(r, start, end) if r["source"] == "BFD" else ann_record(r)
    end = min(end, len(rec.seq))
    sample = f"{i+1:02d}_{r['genus']}_{r['idiomorph']}"
    genes = mp_genes(r)
    new = []
    used = set()
    for f in feats:
        f0, f1 = int(f.location.start) + 1, int(f.location.end)
        if f1 < start or f0 > end:
            continue
        # An annotation model that covers >= 30% of two MATPredict genes is a fused
        # model (e.g. sexP + rnhA): drop it and use the MATPredict models instead.
        hits = [j for j, e in enumerate(genes)
                if (1 if e["strand"] == "+" else -1) == f.location.strand
                and overlap(f0, f1, e["start"], e["end"]) >= 0.3 * (e["end"] - e["start"] + 1)]
        if len(hits) >= 2:
            continue
        name, best = "other", 0.0
        for j, e in enumerate(genes):
            es = 1 if e["strand"] == "+" else -1
            if es != f.location.strand:
                continue
            ov = overlap(f0, f1, e["start"], e["end"])
            if ov > 0:
                used.add(j)   # an annotation CDS covers this model: never add it again
            frac = ov / min(f1 - f0 + 1, e["end"] - e["start"] + 1)
            if frac >= 0.3 and frac > best:
                name, best = e["gene"], frac
        new.append((f.location, name))
    for j, e in enumerate(genes):
        if j in used:
            continue
        strand = 1 if e["strand"] == "+" else -1
        ex = e.get("exons") or [{"start": e["start"], "end": e["end"]}]
        parts = [FeatureLocation(x["start"] - 1, x["end"], strand=strand) for x in sorted(ex, key=lambda x: x["start"], reverse=strand < 0)]
        new.append((parts[0] if len(parts) == 1 else CompoundLocation(parts), e["gene"] + "*"))
    sub = rec.seq[start - 1:end]
    out = SeqRecord(sub, id=sample[:16], name=sample[:16],
                    description=f"{r['name']} {r['source']} {r['genome']} {r['contig']}:{start}-{end} {r['idiomorph']}",
                    annotations={"molecule_type": "DNA"})
    k = 0
    for loc, name in sorted(new, key=lambda x: int(x[0].start)):
        shifted = loc._shift(-(start - 1)) if hasattr(loc, "_shift") else loc
        if int(shifted.start) < 0 or int(shifted.end) > len(sub):
            continue
        k += 1
        nt = shifted.extract(sub)
        nt = nt[: len(nt) - len(nt) % 3]
        prot = str(nt.translate(table=1)).rstrip("*")
        out.features.append(SeqFeature(shifted, type="CDS", qualifiers={
            "locus_tag": [f"{sample}_{k}"], "gene": [name.rstrip("*")], "product": [name],
            "translation": [prot.replace("*", "X")]}))
    SeqIO.write(out, f"{O}/gbk/{sample}.gbk", "genbank")
    summary.append(dict(sample=sample, name=r["name"], source=r["source"], genome=r["genome"], idiomorph=r["idiomorph"],
                        contig=r["contig"], region=f"{start}-{end}", locus_size_inner=r["locus_size_inner"],
                        flanks=f"{r['left_flank']}|{r['right_flank']}", margin=r["margin"], n_cds=k,
                        mat_genes=",".join(sorted({n for _, n in new if n != 'other'})), annotation=src))
with open(f"{O}/selection.tsv", "w") as fo:
    w = csv.DictWriter(fo, fieldnames=list(summary[0]), delimiter="\t", lineterminator="\n")
    w.writeheader(); w.writerows(summary)
for s in summary:
    print(s["sample"], s["name"], s["source"], s["flanks"], s["locus_size_inner"], s["n_cds"], s["mat_genes"])
