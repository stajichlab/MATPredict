#!/usr/bin/env python3
"""CAAX precursor scan (commit c8412b7, basidio-anchors f017ecb): isolate its effect.

Arms, all per genome:
  on   = run-f017ecb           (scan on for Basidiomycota:PR)
  off  = run-f017ecb-nocaax    (same code and db, scan removed from the roster)
  full = results/2026-09-26_basidiomycota_full (run-ad1f865), Agaricales only

Per family: PR/B-locus calls gained/lost/changed (on vs off, on vs full),
precursor ORFs per call, gained receptors' nearest neighbour in the step-1
STE3 tree (inside a small subclade around a curated Agaricales mating
receptor?), identity to curated mating receptors, receptor on the same contig
as HD, HD calls unchanged, runtime.

Usage: compare.py   (writes calls_*.tsv, gained.tsv, per_family.tsv, summary.txt)
"""
import collections
import csv
import glob
import gzip
import os
import subprocess
import sys
import tarfile
import io

import yaml
from Bio import Phylo
from Bio.Seq import Seq

HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.dirname(HERE)
LIB = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes"
BIN = "/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin"
STEP1 = f"{R}/2026-09-27_receptor_explore"
FULL = f"{R}/2026-09-26_basidiomycota_full"
PR_LIKE = ("Basidiomycota:PR", "Basidiomycota:Balpha", "Basidiomycota:Bbeta")
HD_LIKE = ("Basidiomycota:HD", "Basidiomycota:Aalpha", "Basidiomycota:Abeta",
           "Basidiomycota:aLocus", "Basidiomycota:bLocus", "Basidiomycota:MAT")


def fam(d):
    return f"{d['family']['phylum']}:{d['family']['locus_name']}" if isinstance(d["family"], dict) else str(d["family"])


def load_arm(outdir):
    rep = {}
    for p in glob.glob(f"{outdir}/runs/*/detection_report.yaml"):
        g = p.split("/")[-2]
        rep[g] = yaml.safe_load(open(p))
        w = os.path.join(os.path.dirname(p), "wall_seconds")
        rep[g]["_wall"] = float(open(w).read().strip()) if os.path.exists(w) else None
    return rep


def load_full(genomes):
    rep = {}
    for p in glob.glob(f"{FULL}/*/runs/*/detection_report.yaml"):
        g = p.split("/")[-2]
        if g in genomes:
            rep[g] = yaml.safe_load(open(p))
            w = os.path.join(os.path.dirname(p), "wall_seconds")
            rep[g]["_wall"] = float(open(w).read().strip()) if os.path.exists(w) else None
    return rep


def calls(r, families):
    out = []
    for d in (r or {}).get("detected") or []:
        if fam(d) in families:
            out.append(d)
    return out


def key(d):
    return (fam(d), d["contig"], d["start"], d["end"])


def overlap(a, b):
    return fam(a) == fam(b) and a["contig"] == b["contig"] and a["start"] <= b["end"] and b["start"] <= a["end"]


def diff(on, off):
    """gained/lost/changed PR-like calls between two arms of one genome."""
    g, l, c = [], [], []
    for d in on:
        m = [e for e in off if overlap(d, e)]
        if not m:
            g.append(d)
        elif any((e["confidence"], e.get("idiomorph"), e.get("locus_class")) != (d["confidence"], d.get("idiomorph"), d.get("locus_class")) for e in m):
            c.append((d, m[0]))
    for e in off:
        if not any(overlap(d, e) for d in on):
            l.append(e)
    return g, l, c


def read_contigs(asm, wanted):
    seqs, name = {}, None
    with gzip.open(f"{LIB}/{asm}.fa.gz", "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                if name in wanted:
                    seqs[name] = []
            elif name in seqs:
                seqs[name].append(line.strip())
    return {k: "".join(v).upper() for k, v in seqs.items()}


def translate_evidence(e, seqs):
    s = seqs.get(e["contig"])
    if not s:
        return ""
    exons = sorted((x["start"], x["end"]) for x in (e.get("exons") or [])) or [(e["start"], e["end"])]
    parts = [s[a - 1:b] for a, b in exons]
    cds = "".join(parts)
    if e["strand"] == "-":
        cds = str(Seq(cds).reverse_complement())
    best = ""
    for f in range(3):
        sub = cds[f:]
        sub = sub[: len(sub) - len(sub) % 3]
        for piece in str(Seq(sub).translate()).split("*"):
            if len(piece) > len(best):
                best = piece
    return best


def mating_subclade_tips():
    """Union of the smallest clades (>=12 tips) around each curated Agaricales
    mating receptor (Coprinopsis 5346_, Schizophyllum 5334_) in the step-1 tree."""
    tree = Phylo.read(f"{STEP1}/tree_ft.treefile", "newick")
    tips = set()
    refs = [t for t in tree.get_terminals() if t.name.startswith(("REF|5346_", "REF|5334_"))]
    for t in refs:
        path = tree.get_path(t)
        for node in reversed(path[:-1]):
            n = node.count_terminals()
            if n >= 12:
                if n <= 60:
                    tips |= {x.name for x in node.get_terminals()}
                break
    return tips, len(refs)


def main():
    lists = {"agaricales_panel": f"{HERE}/lists/agaricales_panel.tsv",
             "polyrussu_controls": f"{HERE}/lists/polyrussu_controls.tsv"}
    meta = {r["genome"]: r for r in csv.DictReader(open(f"{FULL}/genomes.tsv"), delimiter="\t")}
    samp = {r["ASMID"]: r for r in csv.DictReader(open("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"))}
    on, off = {}, {}
    for name in lists:
        on.update(load_arm(f"{HERE}/on_{name}"))
        off.update(load_arm(f"{HERE}/off_{name}"))
    agar = [l.split("\t")[0] for l in open(lists["agaricales_panel"])]
    full = load_full(set(agar))
    genomes = sorted(set(on) | set(off))
    tips, nrefs = mating_subclade_tips()

    fam_of = lambda g: (meta.get(g, {}).get("family") or samp.get(g, {}).get("FAMILY") or "?")
    order_of = lambda g: (meta.get(g, {}).get("order") or samp.get(g, {}).get("ORDER") or "?")

    per = collections.defaultdict(collections.Counter)
    gained_rows, all_rows = [], []
    hd_changed = []
    need = collections.defaultdict(set)
    for g in genomes:
        o, f_ = on.get(g), off.get(g)
        if o is None or f_ is None:
            per[("?", "missing")]["genomes"] += 1
            continue
        k = (order_of(g), fam_of(g))
        per[k]["genomes"] += 1
        pon, poff = calls(o, PR_LIKE), calls(f_, PR_LIKE)
        per[k]["PR_on"] += len([d for d in pon if fam(d) == "Basidiomycota:PR"])
        per[k]["PR_off"] += len([d for d in poff if fam(d) == "Basidiomycota:PR"])
        per[k]["Bab_on"] += len([d for d in pon if fam(d) != "Basidiomycota:PR"])
        per[k]["Bab_off"] += len([d for d in poff if fam(d) != "Basidiomycota:PR"])
        gn, ls, ch = diff(pon, poff)
        per[k]["gained"] += len(gn); per[k]["lost"] += len(ls); per[k]["changed"] += len(ch)
        per[k]["genomes_gaining"] += 1 if gn else 0
        hon, hoff = calls(o, HD_LIKE), calls(f_, HD_LIKE)
        hg, hl, hc = diff(hon, hoff)
        if hg or hl or hc:
            hd_changed.append(g)
        per[k]["HD_on"] += len(hon)
        if g in full:
            pfull = calls(full[g], PR_LIKE)
            fg, fl, fc = diff(pon, pfull)
            per[k]["vsfull_gained"] += len(fg); per[k]["vsfull_lost"] += len(fl); per[k]["vsfull_changed"] += len(fc)
            hfg, hfl, hfc = diff(hon, calls(full[g], HD_LIKE))
            per[k]["HD_vsfull_diff"] += 1 if (hfg or hfl or hfc) else 0
        if o.get("_wall") and f_.get("_wall"):
            per[k]["wall_on"] += o["_wall"]; per[k]["wall_off"] += f_["_wall"]
        hd_contigs = {d["contig"] for d in hon}
        for d in pon:
            caax = [e for e in d["gene_evidence"] if e.get("method") == "caax_scan"]
            rec = [e for e in d["gene_evidence"] if e["gene"] == "pheromone_receptor"]
            row = dict(genome=g, order=k[0], family=k[1], species=samp.get(g, {}).get("SPECIES_IN", ""),
                       locus=fam(d), contig=d["contig"], start=d["start"], end=d["end"],
                       confidence=d["confidence"], locus_class=d.get("locus_class"),
                       genes=",".join(d.get("genes_found", [])),
                       caax_orfs=sum(e.get("orf_count") or 0 for e in caax) if caax else 0,
                       caax_motifs=",".join(sorted({e.get("caax_motif") or "" for e in caax})),
                       receptor_identity=max((e["identity"] for e in rec), default=None),
                       receptor_status=",".join(sorted({e["status"] for e in rec})),
                       receptor_ref=",".join(sorted({e["reference_record"] for e in rec})),
                       same_contig_as_HD=d["contig"] in hd_contigs,
                       gained=any(d is x for x in gn))
            all_rows.append(row)
            if row["gained"]:
                gained_rows.append((row, rec))
                for e in rec:
                    need[g].add(e["contig"])

    # nearest neighbour of each gained receptor in the step-1 STE3 set
    tmp = os.path.join(os.environ.get("SCRATCH", "/tmp"), "caax_cmp")
    os.makedirs(tmp, exist_ok=True)
    q = os.path.join(tmp, "gained_receptors.faa")
    with open(q, "w") as fo:
        for g, ctgs in need.items():
            seqs = read_contigs(g, ctgs)
            for row, rec in gained_rows:
                if row["genome"] != g:
                    continue
                for i, e in enumerate(rec):
                    p = translate_evidence(e, seqs)
                    if len(p) >= 100:
                        fo.write(f">{g}|{row['contig']}|{row['start']}|{i}\n{p}\n")
    db = os.path.join(tmp, "step1")
    subprocess.run([f"{BIN}/makeblastdb", "-in", f"{STEP1}/ste3_all.faa", "-dbtype", "prot", "-out", db],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    out = subprocess.run([f"{BIN}/blastp", "-query", q, "-db", db, "-outfmt", "6 qseqid sseqid pident bitscore",
                          "-max_target_seqs", "5", "-evalue", "1e-5", "-num_threads", "8"],
                         check=True, capture_output=True, text=True).stdout
    nn = {}
    for line in out.splitlines():
        qq, ss, pid, bs = line.split("\t")
        if qq not in nn or float(bs) > nn[qq][2]:
            nn[qq] = (ss, float(pid), float(bs))
    by_call = collections.defaultdict(list)
    for qq, v in nn.items():
        g, c, s, i = qq.split("|")
        by_call[(g, c, int(s))].append(v)
    for row, _ in gained_rows:
        v = by_call.get((row["genome"], row["contig"], row["start"]), [])
        best = max(v, key=lambda x: x[2]) if v else None
        row["step1_nn"] = best[0] if best else ""
        row["step1_nn_ident"] = best[1] if best else None
        row["nn_in_mating_subclade"] = bool(best) and best[0] in tips
        row["nn_is_curated_mating"] = bool(best) and best[0].startswith("REF|")

    cols = list(all_rows[0].keys()) if all_rows else []
    with open(f"{HERE}/calls_on.tsv", "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=cols, delimiter="\t"); w.writeheader(); w.writerows(all_rows)
    gcols = cols + ["step1_nn", "step1_nn_ident", "nn_in_mating_subclade", "nn_is_curated_mating"]
    with open(f"{HERE}/gained.tsv", "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=gcols, delimiter="\t", extrasaction="ignore"); w.writeheader()
        w.writerows([r for r, _ in gained_rows])
    fields = ["genomes", "PR_off", "PR_on", "gained", "lost", "changed", "genomes_gaining", "Bab_off", "Bab_on",
              "HD_on", "vsfull_gained", "vsfull_lost", "vsfull_changed", "HD_vsfull_diff", "wall_off", "wall_on"]
    with open(f"{HERE}/per_family.tsv", "w") as fo:
        fo.write("order\tfamily\t" + "\t".join(fields) + "\n")
        for k in sorted(per):
            fo.write(f"{k[0]}\t{k[1]}\t" + "\t".join(str(round(per[k][f], 1)) for f in fields) + "\n")

    tot = collections.Counter()
    for k, c in per.items():
        tot.update(c)
    gr = [r for r, _ in gained_rows]
    with open(f"{HERE}/summary.txt", "w") as fo:
        pr = lambda *a: print(*a, file=fo)
        pr(f"genomes with both arms: {sum(1 for g in genomes if g in on and g in off)} of {len(genomes)}")
        pr(f"mating subclade tips: {len(tips)} around {nrefs} curated Agaricales receptors")
        pr(f"PR calls off -> on: {tot['PR_off']} -> {tot['PR_on']}; gained {tot['gained']}, lost {tot['lost']}, changed {tot['changed']}")
        pr(f"Balpha/Bbeta calls off -> on: {tot['Bab_off']} -> {tot['Bab_on']}")
        pr(f"HD-type calls changed on vs off: {len(hd_changed)} genomes {hd_changed[:10]}")
        pr(f"vs full run (Agaricales): gained {tot['vsfull_gained']}, lost {tot['vsfull_lost']}, changed {tot['vsfull_changed']}; HD differ in {tot['HD_vsfull_diff']} genomes")
        if tot["wall_off"]:
            pr(f"wall on/off total: {tot['wall_on']/3600:.2f} h / {tot['wall_off']/3600:.2f} h (x{tot['wall_on']/tot['wall_off']:.2f})")
        pr(f"gained calls: {len(gr)}; with a caax precursor: {sum(1 for r in gr if r['caax_orfs'])}")
        orfs = sorted(r["caax_orfs"] for r in gr)
        if orfs:
            pr(f"  caax ORFs per gained call: median {orfs[len(orfs)//2]}, max {orfs[-1]}")
        pr(f"  confidence: {collections.Counter(r['confidence'] for r in gr)}")
        pr(f"  receptor status: {collections.Counter(r['receptor_status'] for r in gr)}")
        pr(f"  same contig as HD: {sum(1 for r in gr if r['same_contig_as_HD'])}")
        pr(f"  step-1 NN inside a curated-mating subclade: {sum(1 for r in gr if r.get('nn_in_mating_subclade'))}; NN is a curated mating receptor: {sum(1 for r in gr if r.get('nn_is_curated_mating'))}; no NN: {sum(1 for r in gr if not r.get('step1_nn'))}")
        ids = sorted(r["receptor_identity"] for r in gr if r["receptor_identity"] is not None)
        if ids:
            pr(f"  receptor identity to curated mating receptors: min {ids[0]:.1f} median {ids[len(ids)//2]:.1f} max {ids[-1]:.1f}")
        pr("per family (gained / genomes gaining / NN in subclade):")
        byf = collections.defaultdict(list)
        for r in gr:
            byf[(r["order"], r["family"])].append(r)
        for k2, lst in sorted(byf.items()):
            pr(f"  {k2[0]}/{k2[1]}: {len(lst)} calls in {len({r['genome'] for r in lst})} genomes; NN in subclade {sum(1 for r in lst if r.get('nn_in_mating_subclade'))}; median receptor id {sorted(r['receptor_identity'] or 0 for r in lst)[len(lst)//2]:.1f}")
    print(open(f"{HERE}/summary.txt").read())


if __name__ == "__main__":
    main()
