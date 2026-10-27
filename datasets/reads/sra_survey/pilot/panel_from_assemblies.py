#!/usr/bin/env python3
"""Cut a multi-lineage reads-type panel for one species from BFD assemblies and the `detect` calls on them.

Why: one reference locus per idiomorph misses lineages whose MAT locus is too divergent for 31-mers (a Pyricularia oryzae panel from a rice-blast strain shared only
1 to 4% of its k-mers with the loci of other *P. oryzae* lineages and left their reads `low_depth`). Here every `detect` locus of the species is cut from its
assembly, loci are clustered by shared k-mers (a FracMinHash sketch, scale 1/8), and one representative per cluster and idiomorph goes into the panel FASTA
(several records in one file; `reads-type` pools the k-mers of all records of an idiomorph). Shared k-mers (in both idiomorphs) stay the single-copy control.

With --reports (the campaign's reports_all.tar.zst) the panel is built from CORES and FLANKS, not whole loci. The core of a locus is the extent of its own
core_MAT genes (MAT1-1-* for a MAT1-1 locus, MAT1-2-* for MAT1-2) plus a margin; everything else is flank. Loci are clustered into LINEAGES by their flanks (so a
MAT1-1 and a MAT1-2 locus of one lineage fall together); each idiomorph's file holds its own cores and the flanks of EVERY representative locus of both idiomorphs.
Flank k-mers are then in both idiomorph sets (the single-copy depth control, never the call), and only core k-mers are idiomorph-specific. Without this a lineage
whose flank occurs in one idiomorph's file only would put flank k-mers among that idiomorph's "unique" k-mers and miscall strains of that lineage.

    panel_from_assemblies.py --species "Pyricularia oryzae" --genomes genomes.tsv --loci loci.tsv --reports reports_all.tar.zst --genome-dir DIR --out OUT --prefix pyricularia
"""
import argparse, collections, csv, gzip, hashlib, multiprocessing as mp, os, re, statistics, subprocess, sys, tarfile

K = 31; SCALE = 8
COMP = str.maketrans("ACGT", "TGCA")


def canon_sketch(seq):
    seq = seq.upper(); out = set()
    for i in range(len(seq) - K + 1):
        km = seq[i:i + K]
        if set(km) <= set("ACGT"):
            rc = km.translate(COMP)[::-1]; c = km if km <= rc else rc
            if int(hashlib.md5(c.encode()).hexdigest()[:8], 16) % SCALE == 0:
                out.add(c)
    return out


def contig(path, name):
    seq = []; on = False
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if on: break
                on = line[1:].split()[0] == name; continue
            if on: seq.append(line.strip())
    return "".join(seq)


def core_extents(tar_path, wanted):
    """{(genome, contig, start, end): (core_start, core_end)} from the reports of `wanted` genomes: the extent of the core_MAT genes of the locus's own idiomorph."""
    import yaml
    L = getattr(yaml, "CSafeLoader", yaml.SafeLoader); out = {}
    proc = subprocess.Popen(["zstd", "-dc", tar_path], stdout=subprocess.PIPE)
    with tarfile.open(fileobj=proc.stdout, mode="r|") as tf:
        for m in tf:
            if not (m.isfile() and m.name.endswith("/detection_report.yaml")): continue
            gn = os.path.basename(os.path.dirname(m.name))
            if gn not in wanted: continue
            d = yaml.load(tf.extractfile(m).read(), Loader=L) or {}
            for loc in d.get("detected") or []:
                idio = loc.get("idiomorph") or ""
                pref = {"MAT1-1": "MAT1-1-", "MAT1-2": "MAT1-2-"}.get(idio)
                if not pref: continue
                gs = [g for g in loc.get("gene_evidence") or [] if g.get("role") == "core_MAT" and str(g.get("gene", "")).startswith(pref)
                      and g.get("start") is not None and g.get("status", "").startswith("polished")]
                if not gs:
                    gs = [g for g in loc.get("gene_evidence") or [] if g.get("role") == "core_MAT" and str(g.get("gene", "")).startswith(pref) and g.get("start") is not None]
                if gs:
                    out[(gn, loc["contig"], int(loc["start"]), int(loc["end"]))] = (min(g["start"] for g in gs), max(g["end"] for g in gs))
    proc.wait()
    return out


def work(a):
    genome, c, s, e, path = a
    try:
        sq = contig(path, c)[s - 1:e]
    except Exception as ex:
        return genome, None, None
    return genome, sq, canon_sketch(sq)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--species", required=True); ap.add_argument("--genomes", required=True); ap.add_argument("--loci", required=True)
    ap.add_argument("--genome-dir", required=True); ap.add_argument("--out", required=True); ap.add_argument("--prefix", required=True)
    ap.add_argument("--max-span", type=int, default=60000); ap.add_argument("--min-span", type=int, default=9000)
    ap.add_argument("--min-cluster-size", type=int, default=2); ap.add_argument("--contain", type=float, default=0.4, help="shared sketch fraction (of the smaller) that joins a cluster")
    ap.add_argument("--max-reps", type=int, default=6); ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--reports", help="reports_all.tar.zst of the campaign: build the panel from cores and flanks (recommended)"); ap.add_argument("--margin", type=int, default=1000)
    ap.add_argument("--max-flank", type=int, default=6000, help="bases of flank kept on each side of a core (nearest the core): every lineage then adds about the same number of shared k-mers, so the shared depth of a minor lineage is not diluted by a major one")
    ap.add_argument("--min-flank", type=int, default=1500, help="a locus with less flank on both sides together is not used as a representative")
    a = ap.parse_args()
    g = {r["genome"]: r for r in csv.DictReader(open(a.genomes), delimiter="\t") if r["species"].startswith(a.species)}
    L = [r for r in csv.DictReader(open(a.loci), delimiter="\t") if r["genome"] in g and r["locus_class"] in ("mat_locus", "partial_locus")
         and r["idiomorph"] in ("MAT1-1", "MAT1-2") and a.min_span <= int(r["end"]) - int(r["start"]) <= a.max_span]
    print(f"{len(g)} genomes, {len(L)} loci within the span limits", flush=True)
    jobs = [(r["genome"], r["contig"], int(r["start"]), int(r["end"]), os.path.join(a.genome_dir, r["genome"] + ".fa.gz")) for r in L]
    with mp.Pool(a.jobs) as pool:
        res = pool.map(work, jobs, chunksize=4)
    items = []
    for r, (gn, sq, sk) in zip(L, res):
        if sq: items.append(dict(genome=gn, contig=r["contig"], start=int(r["start"]), end=int(r["end"]), idio=r["idiomorph"], conf=r["confidence"], cls=r["locus_class"], seq=sq, sketch=sk))
    print(len(items), "loci cut", flush=True)
    if a.reports:
        ext = core_extents(a.reports, {r["genome"] for r in L})
        ci = []
        for x in items:
            e = ext.get((x["genome"], x["contig"], x["start"], x["end"]))
            if not e: continue
            c0 = max(0, e[0] - a.margin - x["start"]); c1 = min(len(x["seq"]), e[1] + a.margin - x["start"] + 1)
            left, core, right = x["seq"][:c0], x["seq"][c0:c1], x["seq"][c1:]
            left, right = left[-a.max_flank:] if a.max_flank else left, right[:a.max_flank] if a.max_flank else right
            if len(left) + len(right) < a.min_flank or len(core) < 300: continue
            x.update(left=left, right=right, core=core, fsketch=canon_sketch(left) | canon_sketch(right), csketch=canon_sketch(core)); ci.append(x)
        print(len(ci), "loci with a core and enough flank (of", len(items), ")", flush=True)
        os.makedirs(a.out, exist_ok=True)
        lin = []     # lineage clusters by flank
        for x in sorted(ci, key=lambda x: (x["cls"] != "mat_locus", x["conf"] != "high", -(len(x["left"]) + len(x["right"])))):
            for cl in lin:
                sh = len(x["fsketch"] & cl["fsk"]) / max(1, min(len(x["fsketch"]), len(cl["fsk"])))
                if sh >= a.contain:
                    cl["members"].append(x); cl["fsk"] |= x["fsketch"]; break
            else:
                lin.append(dict(members=[x], fsk=set(x["fsketch"])))
        lin.sort(key=lambda c: -len(c["members"]))
        reps = {i: [] for i in ("MAT1-1", "MAT1-2")}; flanks = []
        with open(os.path.join(a.out, f"{a.prefix}_lineages.tsv"), "w") as fh:
            fh.write("lineage\tloci\tMAT1-1\tMAT1-2\trep_MAT1-1\trep_MAT1-2\tkept\n")
            for n, cl in enumerate(lin, 1):
                kept = n <= a.max_reps and len(cl["members"]) >= a.min_cluster_size
                row = []
                for i in ("MAT1-1", "MAT1-2"):
                    ms = [m for m in cl["members"] if m["idio"] == i]
                    ms.sort(key=lambda m: (m["cls"] != "mat_locus", m["conf"] != "high", -(len(m["left"]) + len(m["right"]))))
                    row.append(ms[0] if ms else None)
                    if ms and kept: reps[i].append((n, ms[0]))
                if kept:
                    for r in row:
                        if r: flanks += [(n, r, "L", r["left"]), (n, r, "R", r["right"])]
                f = lambda r: f"{r['genome']}:{r['contig']}:{r['start']}-{r['end']}" if r else ""
                fh.write(f"{n}\t{len(cl['members'])}\t{sum(1 for m in cl['members'] if m['idio'] == 'MAT1-1')}\t{sum(1 for m in cl['members'] if m['idio'] == 'MAT1-2')}\t{f(row[0])}\t{f(row[1])}\t{'yes' if kept else 'no'}\n")
        def wr(fh, name, desc, seq):
            if len(seq) >= 200: fh.write(f">{name} {desc}\n" + "\n".join(seq[j:j + 80] for j in range(0, len(seq), 80)) + "\n")
        for i in ("MAT1-1", "MAT1-2"):
            seen = []
            with open(os.path.join(a.out, f"{a.prefix}_{i}.fasta"), "w") as fh:
                for n, r in reps[i]:                        # cores of this idiomorph (a core that repeats an earlier one adds no k-mer)
                    if any(len(r["csketch"] & s) / max(1, min(len(r["csketch"]), len(s))) >= 0.9 for s in seen): continue
                    seen.append(r["csketch"]); wr(fh, f"{a.prefix}_{i}_core_L{n}", f"{r['genome']}:{r['contig']}:{r['start']}-{r['end']} core", r["core"])
                for n, r, side, seq in flanks:              # flanks of every representative, in BOTH idiomorph files: shared k-mers
                    wr(fh, f"{a.prefix}_flank_L{n}_{r['idio']}_{side}", f"{r['genome']}:{r['contig']}", seq)
        print("lineages (loci):", [len(c["members"]) for c in lin]); return
    med = {i: statistics.median(len(x["seq"]) for x in items if x["idio"] == i) for i in ("MAT1-1", "MAT1-2") if any(x["idio"] == i for x in items)}
    # greedy clustering per idiomorph: best loci first (mat_locus, high confidence, span near the median)
    rank = lambda x: (x["cls"] != "mat_locus", x["conf"] != "high", abs(len(x["seq"]) - med[x["idio"]]))
    clusters = {i: [] for i in med}
    for x in sorted(items, key=rank):
        for cl in clusters[x["idio"]]:
            sh = len(x["sketch"] & cl["rep"]["sketch"]) / max(1, min(len(x["sketch"]), len(cl["rep"]["sketch"])))
            if sh >= a.contain:
                cl["members"].append(x); break
        else:
            clusters[x["idio"]].append(dict(rep=x, members=[x]))
    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, f"{a.prefix}_clusters.tsv"), "w") as fh:
        fh.write("idiomorph\tcluster\tmembers\trep_genome\trep_contig\trep_start\trep_end\trep_len\tkept\n")
        for i, cls in clusters.items():
            cls.sort(key=lambda c: -len(c["members"]))
            for n, cl in enumerate(cls, 1):
                r = cl["rep"]; kept = len(cl["members"]) >= a.min_cluster_size and n <= a.max_reps
                fh.write(f"{i}\t{n}\t{len(cl['members'])}\t{r['genome']}\t{r['contig']}\t{r['start']}\t{r['end']}\t{len(r['seq'])}\t{'yes' if kept else 'no'}\n")
                cl["kept"] = kept; cl["n"] = n
    for i, cls in clusters.items():
        with open(os.path.join(a.out, f"{a.prefix}_{i}.fasta"), "w") as fh:
            for cl in cls:
                if not cl["kept"]: continue
                r = cl["rep"]
                fh.write(f">{a.prefix}_{i}_c{cl['n']} {r['genome']}:{r['contig']}:{r['start']}-{r['end']} members={len(cl['members'])}\n")
                fh.write("\n".join(r["seq"][j:j + 80] for j in range(0, len(r["seq"]), 80)) + "\n")
    print({i: [len(c["members"]) for c in cls] for i, cls in clusters.items()})


if __name__ == "__main__":
    main()
