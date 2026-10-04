#!/usr/bin/env python3
"""Collect the Mucoromycotina MAT campaign into one table and one protein FASTA.

Inputs: runs/<SOURCE>__<id>/detection_report.yaml and detected_loci.gff3 (run with
--emit-cds-fasta), inputs.tsv, curator names (LCG, Jena), BFD samples.csv,
db/taxon_overrides.tsv.

Outputs:
  calls.tsv         one row per detected Mucoromycota:MAT locus (plus one row per
                    genome with no call, locus fields empty)
  core_proteins.faa translation of the called core gene (sexP for Plus, sexM for
                    Minus) when it has a polished model (exonerate or miniprot)
  flank_proteins.faa translations of polished flank genes (for synteny checks)

Locus size (curator ruling 2026-10-03): inner ends of the flanking genes, i.e. from
the end of the nearest flank gene on one side of the core gene to the start of the
nearest flank gene on the other side. Only flank genes with a polished gene model
(exonerate/miniprot) count. Only when flank genes lie on both sides of the
core gene on the same contig; otherwise locus_size_inner is empty and
flank_status says why.
"""
import csv, os, re, sys, glob
import yaml

C = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-03_mucoromycotina_mat"
R = "/bigdata/stajichlab/jstajich/projects/MATPredict/results"
WT = "/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-7c7ed99"
CORE = {"Plus": "sexP", "Minus": "sexM"}

def names():
    out = {}
    for r in csv.DictReader(open(f"{R}/2026-09-28_lcg_holdout/curator_table.tsv"), delimiter="\t"):
        out[("LCG", r["org"])] = dict(name=r.get("curator_taxonomy") or r["org"].replace("_", " "),
                                     name_in_source=r["org"].replace("_", " "),
                                     label=r.get("curator_mating_type", ""))
    for r in csv.DictReader(open(f"{R}/2026-09-28_mucor_jena_holdout/curator_table.tsv"), delimiter="\t"):
        out[("JENA", r["strain"])] = dict(name=r.get("curator_taxonomy") or r.get("curator_name_as_given") or r["strain"],
                                        name_in_source=r.get("curator_name_as_given") or r["strain"],
                                        label=r.get("curator_known_mating_type", ""))
    for r in csv.DictReader(open("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv")):
        out[("BFD", r["ASMID"])] = dict(name=(r["SPECIES"] or r["SPECIES_IN"]) + (f" {r['STRAIN']}" if r["STRAIN"] else ""),
                                       name_in_source=r["SPECIES_IN"], label="", order=r["ORDER"])
    return out

def overrides():
    o = {}
    for l in open(f"{WT}/db/taxon_overrides.tsv"):
        if l.startswith("#") or l.startswith("genome_id"):
            continue
        f = l.rstrip("\n").split("\t")
        if len(f) >= 9:
            o[f[0]] = f
    return o

def gff_translations(path):
    """{(contig, start, end, gene): translation} from the --emit-cds-fasta GFF3."""
    gene_by_id, tr = {}, {}
    if not os.path.exists(path):
        return tr
    for l in open(path):
        if l.startswith("#"):
            continue
        f = l.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        a = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
        if f[2] == "gene":
            gene_by_id[a["ID"]] = (f[0], int(f[3]), int(f[4]), a.get("Name"))
        elif f[2] == "CDS" and "translation" in a:
            g = gene_by_id.get(a.get("Parent"))
            if g:
                tr[g] = a["translation"]
    return tr

def main():
    nm, ov = names(), overrides()
    inputs = list(csv.DictReader(open(f"{C}/inputs.tsv"), delimiter="\t"))
    rows, core_fa, flank_fa = [], [], []
    for inp in inputs:
        src, gid = inp["source"], inp["id"]
        d = f"{C}/runs/{src}__{gid}"
        rep = f"{d}/detection_report.yaml"
        n = nm.get((src, gid), {})
        base = dict(source=src, genome=gid, name=n.get("name", ""), name_in_source=n.get("name_in_source", ""),
                    file_label=n.get("label", ""), bfd_order=inp.get("order", ""),
                    override=("" if gid not in ov else ov[gid][5] + ": " + ov[gid][2]),
                    override_use=("" if gid not in ov else ov[gid][6]))
        if not os.path.exists(rep) or os.path.getsize(rep) == 0:
            rows.append(dict(base, status="no_report")); continue
        doc = yaml.safe_load(open(rep)) or {}
        det = [x for x in (doc.get("detected") or []) if x.get("family") == "Mucoromycota:MAT"]
        if not det:
            rows.append(dict(base, status="uncalled", n_calls=0)); continue
        tr = gff_translations(f"{d}/detected_loci.gff3")
        for i, x in enumerate(det):
            idio = x.get("idiomorph")
            core_name = CORE.get(idio)
            ev = x.get("gene_evidence") or []
            core = [e for e in ev if e.get("gene") == core_name and e.get("contig") == x["contig"]]
            core = max(core, key=lambda e: e.get("bitscore") or 0) if core else None
            polished = bool(core and str(core.get("status", "")).startswith("polished"))
            prot = tr.get((core["contig"], core["start"], core["end"], core_name)) if core else None
            # Only flank genes with a polished model bound the locus. An unpolished
            # tblastn hit can be a weak, distant match (e.g. a 27%-identity algA hit
            # 48 kb from the core gene gave a false 49 kb locus in S. racemosum).
            flanks = [e for e in ev if str(e.get("role", "")).startswith("flanking") and e.get("contig") == x["contig"]
                      and str(e.get("status", "")).startswith("polished")]
            size = left = right = ""
            if core is None:
                fstat = "no_core_gene_model"
            else:
                L = [e for e in flanks if e["end"] < core["start"]]
                Rr = [e for e in flanks if e["start"] > core["end"]]
                if L and Rr:
                    lf, rf = max(L, key=lambda e: e["end"]), min(Rr, key=lambda e: e["start"])
                    size, left, right = rf["start"] - lf["end"] - 1, lf["gene"], rf["gene"]
                    fstat = "both_sides"
                else:
                    fstat = "flank_one_side" if (L or Rr) else "no_flank"
            loc_id = f"{src}__{gid}__L{i}"
            row = dict(base, status="called", n_calls=len(det), locus_id=loc_id, idiomorph=idio,
                       confidence=x.get("confidence"), locus_class=x.get("locus_class"),
                       detection_pass=x.get("detection_pass"), margin=x.get("idiomorph_margin"),
                       contig=x["contig"], locus_start=x["start"], locus_end=x["end"],
                       genes_found=",".join(x.get("genes_found") or []),
                       core_gene=core_name or "", core_status=(core or {}).get("status", ""),
                       core_start=(core or {}).get("start", ""), core_end=(core or {}).get("end", ""),
                       core_strand=(core or {}).get("strand", ""), core_identity=(core or {}).get("identity", ""),
                       core_aa=len(prot) if prot else "",
                       flank_status=fstat, left_flank=left, right_flank=right, locus_size_inner=size,
                       two_idiomorphs=bool(doc.get("two_idiomorphs")))
            rows.append(row)
            if prot and polished:
                core_fa.append(f">{loc_id}|{idio}|{core_name}\n{prot.replace('*', '')}\n")
            for e in flanks:
                p = tr.get((e["contig"], e["start"], e["end"], e["gene"]))
                if p and str(e.get("status", "")).startswith("polished"):
                    flank_fa.append(f">{loc_id}|{e['gene']}\n{p.replace('*', '')}\n")
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(f"{C}/calls.tsv", "w") as fo:
        w = csv.DictWriter(fo, fieldnames=keys, delimiter="\t", lineterminator="\n", restval="")
        w.writeheader(); w.writerows(rows)
    open(f"{C}/core_proteins.faa", "w").write("".join(core_fa))
    open(f"{C}/flank_proteins.faa", "w").write("".join(flank_fa))
    from collections import Counter
    print("genomes", len(inputs), Counter(r["status"] for r in rows if r.get("status") != "called"),
          "called loci", sum(r["status"] == "called" for r in rows),
          "genomes called", len({r["genome"] for r in rows if r["status"] == "called"}))
    print("idiomorph", Counter(r.get("idiomorph") for r in rows if r["status"] == "called"))
    print("flank_status", Counter(r.get("flank_status") for r in rows if r["status"] == "called"))
    print("core proteins", len(core_fa), "flank proteins", len(flank_fa))

if __name__ == "__main__":
    main()
