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
from collections import Counter
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


PAIRS = [("MAT1-1", "MAT1-2"), ("a", "alpha"), ("P", "M"), ("A", "B")]      # idiomorph name pairs, in order of preference


def norm_idio(name):
    """One name per idiomorph across the reports' spellings: MATa / MTLa / A (when paired with alpha) -> a; MATalpha / MTLalpha / alpha -> alpha; MAT1-1, MAT1-2, P, M, A, B kept. '' for undetermined."""
    n = (name or "").strip()
    if not n or n.lower() in ("undetermined", "none", "na"): return ""
    if n in ("MAT1-1", "MAT1-2"): return n
    t = re.sub(r"^(MAT|MTL)", "", n, flags=re.I)
    return "alpha" if t.lower() in ("alpha", "α") else "a" if n == "a" or n.lower() in ("mata", "mtla") else t if t in ("P", "M", "A", "B") else n


def choose_names(counts, partner_counts=None):
    """The two idiomorph names of a species from the idiomorph counts of its loci (`partner_counts`: those of its genus, to complete a species that shows only one). The first preferred pair that is present in both
    names; a species with one name takes its pair's partner when the genus has it; otherwise the two commonest."""
    c = {k: v for k, v in counts.items() if k}
    pc = {k: v for k, v in (partner_counts or {}).items() if k}
    for a_, b_ in PAIRS:
        if a_ in c and b_ in c: return [a_, b_]
    if c and pc:
        for a_, b_ in PAIRS + [("a", "A")]:
            for x, y in ((a_, b_), (b_, a_)):
                if x in c and y in pc: return [a_, b_]
    top = [k for k, _ in sorted({**pc, **c}.items(), key=lambda kv: -kv[1])]
    for a_, b_ in PAIRS:
        if any(n in c for n in (a_, b_)): return [a_, b_]
    return (top + ["", ""])[:2]


def core_extents(tar_path, wanted, names=("MAT1-1", "MAT1-2"), norm=None):
    """{(genome, contig, start, end): (core_start, core_end)} from the reports of `wanted` genomes: the extent of the core_MAT genes of the locus's own idiomorph."""
    import yaml
    L = getattr(yaml, "CSafeLoader", yaml.SafeLoader); out = {}
    import shutil
    tar_path = os.path.realpath(tar_path)       # zstd refuses to read a symbolic link
    zstd = shutil.which("zstd") or "/opt/linux/rocky/8.x/x86_64/pkgs/miniconda3/py39_4.12.0/bin/zstd"      # not on PATH in a non-login Slurm shell
    proc = subprocess.Popen([zstd, "-dc", tar_path], stdout=subprocess.PIPE)
    with tarfile.open(fileobj=proc.stdout, mode="r|") as tf:
        for m in tf:
            if not (m.isfile() and m.name.endswith("/detection_report.yaml")): continue
            gn = os.path.basename(os.path.dirname(m.name))
            if gn not in wanted: continue
            d = yaml.load(tf.extractfile(m).read(), Loader=L) or {}
            for loc in d.get("detected") or []:
                idio = (norm or norm_idio)(loc.get("idiomorph"))
                if idio not in names: continue
                pref = {"MAT1-1": "MAT1-1-", "MAT1-2": "MAT1-2-"}.get(idio, "")      # other names (yeast MATa/MATalpha, MTLA1...): every core_MAT gene of the locus
                gs = [g for g in loc.get("gene_evidence") or [] if g.get("role") == "core_MAT" and str(g.get("gene", "")).startswith(pref)
                      and g.get("start") is not None and g.get("status", "").startswith("polished")]
                if not gs:
                    gs = [g for g in loc.get("gene_evidence") or [] if g.get("role") == "core_MAT" and str(g.get("gene", "")).startswith(pref) and g.get("start") is not None]
                if gs:
                    out[(gn, loc["contig"], int(loc["start"]), int(loc["end"]))] = (min(g["start"] for g in gs), max(g["end"] for g in gs))
    if proc.wait(): raise RuntimeError(f"zstd failed on {tar_path}")
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
    ap.add_argument("--diverse", action="store_true", help="rank lineages by the number of distinct species in them before their size (for a genus panel: a few rare lineages of other species beat many loci of the commonest one)")
    ap.add_argument("--per-lineage", action="store_true", help="also write one complete panel per lineage in OUT/L<n>/ (cores and flanks of that lineage only; a lineage lacking one idiomorph borrows the core of the first lineage that has it). A multi-lineage panel dilutes breadth and depth over cores a strain does not carry; typing against each lineage panel and keeping the best-matching one does not.")
    ap.add_argument("--classes", nargs="+", default=["mat_locus", "partial_locus"], help="locus classes used (add idiomorph_gene_only for lineages whose loci are only ever reported as a gene, such as the Taphrinomycotina)")
    ap.add_argument("--names", nargs=2, help="the two idiomorph names (default: chosen from the loci: MAT1-1/MAT1-2, a/alpha, P/M or A/B)")
    ap.add_argument("--min-own", type=int, default=2, help="with --borrow: an idiomorph with fewer than this many loci in the species also takes loci from the genus")
    ap.add_argument("--borrow", help="genus prefix (for example 'Aspergillus '): an idiomorph with no locus in the species takes its loci from the other species of this genus (flank clustering then pairs them with the species' lineages)")
    ap.add_argument("--min-flank", type=int, default=1500, help="a locus with less flank on both sides together is not used as a representative")
    a = ap.parse_args()
    gall = {r["genome"]: r for r in csv.DictReader(open(a.genomes), delimiter="\t")}
    g = {k: r for k, r in gall.items() if r["species"].startswith(a.species)}
    ok_loc = lambda r: r["locus_class"] in a.classes and a.min_span <= int(r["end"]) - int(r["start"]) <= a.max_span
    loci = list(csv.DictReader(open(a.loci), delimiter="\t"))
    own = [r for r in loci if r["genome"] in g and ok_loc(r)]
    genus_loci = [r for r in loci if a.borrow and r["genome"] in gall and gall[r["genome"]]["species"].startswith(a.borrow) and ok_loc(r)]
    has_alpha = any(norm_idio(r["idiomorph"]) == "alpha" for r in own + genus_loci)
    nm = lambda x: "a" if has_alpha and norm_idio(x) == "A" else norm_idio(x)       # a yeast that has alpha: its A is MATa
    names = a.names or choose_names(Counter(nm(r["idiomorph"]) for r in own), Counter(nm(r["idiomorph"]) for r in genus_loci))
    L = [dict(r, idiomorph=nm(r["idiomorph"])) for r in own if nm(r["idiomorph"]) in names]
    own_set = set(g)
    if a.borrow:
        cnt = Counter(r["idiomorph"] for r in L); missing = [n for n in names if n and cnt[n] < a.min_own]      # no locus, or fewer than --min-own: the genus tops it up
        extra = [dict(r, idiomorph=nm(r["idiomorph"])) for r in genus_loci if r["genome"] not in g and nm(r["idiomorph"]) in missing]
        print("borrowed", len(extra), "loci of", missing, "from", len({r["genome"] for r in extra}), "genomes of", a.borrow, flush=True)
        L += extra; g = {**g, **{k: gall[k] for k in {r["genome"] for r in extra}}}
    open(os.path.join(a.out, f"{a.prefix}_idiomorphs.txt"), "w").write(" ".join(names) + "\n") if os.makedirs(a.out, exist_ok=True) is None else None
    print(f"{len(g)} genomes, {len(L)} loci within the span limits", flush=True)
    jobs = [(r["genome"], r["contig"], int(r["start"]), int(r["end"]), os.path.join(a.genome_dir, r["genome"] + ".fa.gz")) for r in L]
    with mp.Pool(a.jobs) as pool:
        res = pool.map(work, jobs, chunksize=4)
    items = []
    for r, (gn, sq, sk) in zip(L, res):
        if sq: items.append(dict(genome=gn, contig=r["contig"], start=int(r["start"]), end=int(r["end"]), idio=r["idiomorph"], conf=r["confidence"], cls=r["locus_class"], seq=sq, sketch=sk))
    print(len(items), "loci cut", flush=True)
    if a.reports:
        ext = core_extents(a.reports, {r["genome"] for r in L}, names, nm)
        ci = []
        for x in items:
            e = ext.get((x["genome"], x["contig"], x["start"], x["end"]))
            if not e: continue
            c0 = max(0, e[0] - a.margin - x["start"]); c1 = min(len(x["seq"]), e[1] + a.margin - x["start"] + 1)
            left, core, right = x["seq"][:c0], x["seq"][c0:c1], x["seq"][c1:]
            left, right = left[-a.max_flank:] if a.max_flank else left, right[:a.max_flank] if a.max_flank else right
            if len(left) + len(right) < a.min_flank or len(core) < 300: continue
            x["sp"] = g.get(x["genome"], {}).get("species", "")
            x.update(left=left, right=right, core=core, fsketch=canon_sketch(left) | canon_sketch(right), csketch=canon_sketch(core)); ci.append(x)
        print(len(ci), "loci with a core and enough flank (of", len(items), ")", flush=True)
        os.makedirs(a.out, exist_ok=True)
        lin = []     # lineage clusters by flank. The species' own loci first; a borrowed locus (another species of the genus) only joins a lineage the species already has (its flank matches); the borrowed loci that
        #              match nothing are kept aside and used only when a lineage of the species still lacks an idiomorph altogether, so a big genus cannot bury the species' own lineages.
        for x in ci: x["own"] = x["genome"] in own_set
        aside = []
        for x in sorted(ci, key=lambda x: (not x["own"], x["cls"] != "mat_locus", x["conf"] != "high", -(len(x["left"]) + len(x["right"])))):
            for cl in lin:
                sh = len(x["fsketch"] & cl["fsk"]) / max(1, min(len(x["fsketch"]), len(cl["fsk"])))
                if sh >= a.contain:
                    cl["members"].append(x); cl["fsk"] |= x["fsketch"]; break
            else:
                if x["own"]: lin.append(dict(members=[x], fsk=set(x["fsketch"])))
                else: aside.append(x)
        have_idio = {m["idio"] for cl in lin for m in cl["members"]}
        if any(n not in have_idio for n in names) and aside:      # nothing borrowed shares a flank with the species: cluster the borrowed loci among themselves and add the largest lineages
            extra = []
            for x in aside:
                if x["idio"] in have_idio: continue
                for cl in extra:
                    if len(x["fsketch"] & cl["fsk"]) / max(1, min(len(x["fsketch"]), len(cl["fsk"]))) >= a.contain: cl["members"].append(x); cl["fsk"] |= x["fsketch"]; break
                else: extra.append(dict(members=[x], fsk=set(x["fsketch"])))
            extra.sort(key=lambda c: -len(c["members"])); lin += extra[:max(1, a.max_reps // 2)]
            print("no borrowed locus shares a flank with the species: added", min(len(extra), max(1, a.max_reps // 2)), "borrowed lineages", flush=True)
        else: print("borrowed loci that joined a lineage of the species:", sum(1 for cl in lin for m in cl["members"] if not m["own"]), "; set aside:", len(aside), flush=True)
        for cl in lin: cl["has_own"] = any(m["own"] for m in cl["members"])
        lin.sort(key=(lambda c: (not c["has_own"], -len({m["sp"] for m in c["members"]}), -len(c["members"]))) if a.diverse else (lambda c: (not c["has_own"], -len(c["members"]))))
        reps = {i: [] for i in names}; flanks = []
        with open(os.path.join(a.out, f"{a.prefix}_lineages.tsv"), "w") as fh:
            fh.write("lineage\tloci\t" + "\t".join(names) + "\t" + "\t".join("rep_" + n for n in names) + "\tkept\n")
            for n, cl in enumerate(lin, 1):
                kept = n <= a.max_reps and len(cl["members"]) >= a.min_cluster_size
                row = []
                for i in names:
                    ms = [m for m in cl["members"] if m["idio"] == i]
                    ms.sort(key=lambda m: (m["cls"] != "mat_locus", m["conf"] != "high", -(len(m["left"]) + len(m["right"]))))
                    row.append(ms[0] if ms else None)
                    if ms and kept: reps[i].append((n, ms[0]))
                if kept:
                    for r in row:
                        if r: flanks += [(n, r, "L", r["left"]), (n, r, "R", r["right"])]
                f = lambda r: f"{r['genome']}:{r['contig']}:{r['start']}-{r['end']}" if r else ""
                fh.write(f"{n}\t{len(cl['members'])}\t{sum(1 for m in cl['members'] if m['idio'] == names[0])}\t{sum(1 for m in cl['members'] if m['idio'] == names[1])}\t{f(row[0])}\t{f(row[1])}\t{'yes' if kept else 'no'}\n")
        def wr(fh, name, desc, seq):
            if len(seq) >= 200: fh.write(f">{name} {desc}\n" + "\n".join(seq[j:j + 80] for j in range(0, len(seq), 80)) + "\n")
        for i in names:
            seen = []
            with open(os.path.join(a.out, f"{a.prefix}_{i}.fasta"), "w") as fh:
                for n, r in reps[i]:                        # cores of this idiomorph (a core that repeats an earlier one adds no k-mer)
                    if any(len(r["csketch"] & s) / max(1, min(len(r["csketch"]), len(s))) >= 0.9 for s in seen): continue
                    seen.append(r["csketch"]); wr(fh, f"{a.prefix}_{i}_core_L{n}", f"{r['genome']}:{r['contig']}:{r['start']}-{r['end']} core", r["core"])
                for n, r, side, seq in flanks:              # flanks of every representative, in BOTH idiomorph files: shared k-mers
                    wr(fh, f"{a.prefix}_flank_L{n}_{r['idio']}_{side}", f"{r['genome']}:{r['contig']}", seq)
        if a.per_lineage:
            first = {i: next((r for n, r in reps[i]), None) for i in reps}
            by_n = {i: dict(reps[i]) for i in reps}
            for n, cl in enumerate(lin, 1):
                if n > a.max_reps: break
                cores = {i: by_n[i].get(n) or first[i] for i in names}
                if not all(cores.values()): continue
                ld = os.path.join(a.out, f"L{n}"); os.makedirs(ld, exist_ok=True)
                src = [by_n[j].get(n) for j in names if by_n[j].get(n)] or [cores[names[0]], cores[names[1]]]
                for i in names:
                    with open(os.path.join(ld, f"{a.prefix}_{i}.fasta"), "w") as fh:
                        r = cores[i]; wr(fh, f"{a.prefix}_{i}_core_L{n}", f"{r['genome']}:{r['contig']}:{r['start']}-{r['end']} core", r["core"])
                        for r2 in src:
                            wr(fh, f"{a.prefix}_flank_L{n}_{r2['idio']}_L", f"{r2['genome']}:{r2['contig']}", r2["left"]); wr(fh, f"{a.prefix}_flank_L{n}_{r2['idio']}_R", f"{r2['genome']}:{r2['contig']}", r2["right"])
        print("lineages (loci):", [len(c["members"]) for c in lin]); return

    med = {i: statistics.median(len(x["seq"]) for x in items if x["idio"] == i) for i in names if any(x["idio"] == i for x in items)}
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
