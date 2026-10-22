"""Standalone model-extension / chooser prototype (read-only analysis). Not part of MATPredict.
usage:
  python modelkit.py survey RUNS_DIR GENOME_DIR TAXONOMY_CSV OUT.tsv [NPROC]
  python modelkit.py test   RUNS_DIR GENOME_DIR GENOME_LIST REFS.faa OUT.tsv   (needs mafft on PATH)
"""
import gzip, sys, csv, subprocess, tempfile, os
from pathlib import Path
from multiprocessing import Pool
import yaml

L = getattr(yaml, "CSafeLoader", yaml.SafeLoader)
B = "TCAG"; AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {B[i // 16] + B[(i // 4) % 4] + B[i % 4]: AA[i] for i in range(64)}
comp = str.maketrans("ACGTacgtN", "TGCAtgcaN")
rc = lambda s: s.translate(comp)[::-1]
EXT = 3000   # bp scanned past an end


def tr(s, code=1):
    out = []
    for i in range(0, len(s) - 2, 3):
        a = CODON.get(s[i:i + 3].upper(), "X")
        if code == 4 and s[i:i + 3].upper() == "TGA":
            a = "W"
        out.append(a)
    return "".join(out)


def contigs(path, names):
    got, name, buf = {}, None, []
    with gzip.open(path, "rt") as f:
        for l in f:
            if l[0] == ">":
                if name in names:
                    got[name] = "".join(buf)
                    if len(got) == len(names):
                        return got
                name = l[1:].split()[0]; buf = []
            elif name in names:
                buf.append(l.strip())
    if name in names:
        got[name] = "".join(buf)
    return got


def build(seq, exons, strand, code):
    """-> dict(cds, prot, down(3' flank in transcript orientation), up(5' flank))"""
    ex = sorted((a, b) for a, b in exons)
    lo, hi = ex[0][0], ex[-1][1]
    cds = "".join(seq[a - 1:b] for a, b in ex)
    if strand == "+":
        down = seq[hi:hi + EXT]; up = seq[max(0, lo - 1 - EXT):lo - 1]
    else:
        cds = rc(cds)
        down = rc(seq[max(0, lo - 1 - EXT):lo - 1]); up = rc(seq[hi:hi + EXT])
    return cds, down, up


def analyse(cds, down, up, code):
    n3 = len(cds) // 3
    prot = tr(cds, code)
    term = prot.endswith("*")
    body = prot[:-1] if term else prot
    internal = body.count("*")
    longest = max((len(x) for x in body.split("*")), default=0)
    # 3' extension, only meaningful when no terminal stop
    ext_aa, ext_seq, ext_state = 0, "", "terminal_stop"
    if not term:
        t = tr(cds[n3 * 3:] + down, code)
        k = t.find("*")
        if k >= 0:
            ext_seq, ext_aa, ext_state = t[:k], k, "stop_found"
        else:
            ext_seq, ext_aa, ext_state = t, len(t), "no_stop_in_%dbp" % EXT
    # 5' look-up: in-frame upstream codons back to previous stop (frame set by the model start)
    u = up[len(up) % 3:] if len(up) % 3 else up
    ut = tr(u, code)
    k = ut.rfind("*")
    seg = ut[k + 1:]
    up_to_stop = len(seg) if k >= 0 else None       # aa between preceding stop and model start
    m = seg.find("M")
    up_M = (len(seg) - m) if m >= 0 else 0           # aa from furthest upstream M to model start
    return dict(prot=body, term=term, internal=internal, longest=longest, ext_aa=ext_aa, ext_seq=ext_seq,
                ext_state=ext_state, starts_M=body[:1] == "M", up_to_stop=up_to_stop, up_M=up_M,
                up_seq=seg[m:] if m >= 0 else "")


def gene_models(ev, seq, code):
    out = {}
    ex = [(a["start"], a["end"]) for a in ev["exons"]]
    cds, down, up = build(seq, ex, ev["strand"], code)
    out["exo"] = analyse(cds, down, up, code)
    alt = ev.get("alternate_model")
    if alt and alt.get("exons"):
        c2, d2, u2 = build(seq, [(a["start"], a["end"]) for a in alt["exons"]], alt["strand"], code)
        out["mini"] = analyse(c2, d2, u2, code)
    return out


def load_tax(p):
    t = {}
    for r in csv.DictReader(open(p)):
        t[r["ASM_ACCESSION"]] = (r["ORDER"], r["GENUS"])
    return t


def survey_one(args):
    runs, gdir, g, tax = args
    rows = []
    try:
        d = yaml.load(open(Path(runs) / g / "detection_report.yaml"), Loader=L)
        code = d.get("genetic_code") or 1
        evs = [e for l in d["detected"] for e in l["gene_evidence"]
               if e.get("role") == "core_MAT" and e.get("exons") and e.get("method") in ("exonerate_refine", "miniprot_refine")]
        if not evs:
            return rows
        sq = contigs(Path(gdir) / f"{g}.fa.gz", {e["contig"] for e in evs})
        for e in evs:
            if e["contig"] not in sq:
                continue
            m = gene_models(e, sq[e["contig"]], code)
            x = m["exo"]
            r = dict(genome=g, order=tax.get(g, ("", ""))[0], genus=tax.get(g, ("", ""))[1], gene=e["gene"], method=e["method"],
                     status=e["status"], identity=e.get("identity"), coverage=e.get("coverage"), code=code,
                     exo_aa=len(x["prot"]), exo_term=x["term"], exo_internal=x["internal"], exo_longest=x["longest"],
                     exo_ext_aa=x["ext_aa"], exo_ext_state=x["ext_state"], exo_startM=x["starts_M"], exo_up_to_stop=x["up_to_stop"], exo_up_M=x["up_M"],
                     has_alt=False)
            if "mini" in m:
                y = m["mini"]
                r.update(has_alt=True, mini_aa=len(y["prot"]), mini_term=y["term"], mini_internal=y["internal"], mini_longest=y["longest"],
                         mini_ext_aa=y["ext_aa"], mini_ext_state=y["ext_state"], mini_identity=(e.get("alternate_model") or {}).get("identity"))
            rows.append(r)
    except Exception as ex:
        rows.append(dict(genome=g, gene="ERROR", status=repr(ex)[:100]))
    return rows


def survey(runs, gdir, taxp, outp, nproc=16):
    tax = load_tax(taxp)
    gl = sorted(p.name for p in Path(runs).iterdir() if (p / "detection_report.yaml").exists())
    allrows = []
    with Pool(nproc) as pool:
        for rs in pool.imap_unordered(survey_one, [(runs, gdir, g, tax) for g in gl], chunksize=4):
            allrows += rs
    cols = []
    for r in allrows:
        for k in r:
            if k not in cols: cols.append(k)
    with open(outp, "w") as f:
        w = csv.DictWriter(f, cols, delimiter="\t"); w.writeheader(); w.writerows(allrows)


# ---------------- chooser (test mode) ----------------
def read_fa(p):
    d, n = {}, None
    for l in open(p):
        l = l.rstrip()
        if l.startswith(">"): n = l[1:]; d[n] = []
        elif n: d[n].append(l)
    return {k: "".join(v) for k, v in d.items()}


def align_score(cand, ref):
    """pairwise mafft; returns (ref_cov, ident_over_aligned, overhang_unaligned_cand_aa, aligned_cols)"""
    with tempfile.NamedTemporaryFile("w", suffix=".faa", delete=False) as f:
        f.write(f">ref\n{ref}\n>cand\n{cand}\n"); fn = f.name
    o = subprocess.run(["mafft", "--quiet", "--localpair", "--maxiterate", "100", fn], capture_output=True, text=True).stdout
    os.unlink(fn)
    s = {}; n = None
    for l in o.splitlines():
        if l.startswith(">"): n = l[1:]; s[n] = []
        else: s[n].append(l.strip())
    a, b = "".join(s["ref"]).upper(), "".join(s["cand"]).upper()
    both = [(x, y) for x, y in zip(a, b) if x != "-" and y != "-"]
    ident = sum(x == y for x, y in both)
    over = sum(1 for x, y in zip(a, b) if x == "-" and y != "-")
    return len(both) / len(ref), (ident / len(both) if both else 0.0), over, len(both), ident


def support_all(cand, refs):
    """score of one candidate against EVERY ref: {ref: dict}"""
    out = {}
    c = cand.replace("*", "X")
    for rn, rs in refs.items():
        cov, idn, over, nb, ni = align_score(c, rs)
        out[rn] = dict(score=ni - 0.5 * over, ref_len=len(rs), ref_cov=cov, ident=idn, overhang=over)
    return out


def after_stop_support(prot, ref):
    """ref columns aligned to candidate residues AFTER its first internal stop (0 when no stop). Readthrough/frameshift/intron evidence."""
    k = prot.find("*")
    if k < 0:
        return 0
    with tempfile.NamedTemporaryFile("w", suffix=".faa", delete=False) as f:
        f.write(f">ref\n{ref}\n>cand\n{prot.replace('*', 'X')}\n"); fn = f.name
    o = subprocess.run(["mafft", "--quiet", "--localpair", "--maxiterate", "100", fn], capture_output=True, text=True).stdout
    os.unlink(fn)
    s = {}; n = None
    for l in o.splitlines():
        if l.startswith(">"): n = l[1:]; s[n] = []
        else: s[n].append(l.strip())
    a, b = "".join(s["ref"]), "".join(s["cand"])
    pos = 0; n_after = 0
    for x, y in zip(a, b):
        if y != "-":
            pos += 1
        if y != "-" and pos > k + 1 and x != "-":
            n_after += 1
    return n_after


def choose(models, refs):
    cands = {}
    for k, m in models.items():
        cands[k] = m["prot"]
        if m["ext_state"] == "stop_found" and m["ext_aa"] >= 20:
            cands[k + "+ext3"] = m["prot"] + m["ext_seq"]
    internal = {k: models[k.split("+")[0]]["internal"] for k in cands}
    mat = {k: support_all(p, refs) for k, p in cands.items()}
    clean = [k for k in cands if internal[k] == 0] or list(cands)
    # one reference for all candidates: the one that best supports the best clean candidate
    ref = max(refs, key=lambda r: max(mat[k][r]["score"] for k in clean))
    rows = {k: dict(aa=len(p), internal=internal[k], ref=ref, **mat[k][ref]) for k, p in cands.items()}
    cl = {k: v for k, v in rows.items() if v["internal"] == 0}
    if cl:
        best = max(cl, key=lambda k: (round(cl[k]["score"] / 3), cl[k]["aa"]))   # scores within ~3 are ties -> longer
        flag = "clean"
    else:
        best = max(rows, key=lambda k: rows[k]["score"]); flag = "interrupted"
    v = rows[best]
    if v["ref_cov"] < 0.8 and v["ident"] >= 0.85:
        flag += "+REMNANT(ident>=85,cov<80)"
    elif v["ref_cov"] < 0.8:
        flag += "+partial(cov<80)"
    for k in rows:
        rows[k]["after_stop_ref_cols"] = after_stop_support(cands[k], refs[ref]) if internal[k] else 0
    return best, flag, rows


def one_genome(a):
    runs, gdir, grp, g, allref = a
    d = yaml.load(open(Path(runs) / g / "detection_report.yaml"), Loader=L)
    code = d.get("genetic_code") or 1
    out = []
    for l in d["detected"]:
        for e in l["gene_evidence"]:
            if e.get("role") != "core_MAT" or e["gene"] != "MAT1-1-1" or not e.get("exons"):
                continue
            sq = contigs(Path(gdir) / f"{g}.fa.gz", {e["contig"]})[e["contig"]]
            models = gene_models(e, sq, code)
            best, flag, rows = choose(models, allref)
            for k, v in rows.items():
                out.append(dict(group=grp, genome=g[:22], cand=k, chosen=(k == best), flag=flag if k == best else "",
                                aa=v["aa"], internal=v["internal"], ref=v["ref"].split("|")[1][:28], ref_len=v["ref_len"],
                                ref_cov=round(v["ref_cov"], 3), ident=round(v["ident"], 3),
                                overhang=v["overhang"], score=round(v["score"], 1), after_stop_ref_cols=v["after_stop_ref_cols"],
                                exo_aa=len(models["exo"]["prot"]), exo_longest=models["exo"]["longest"], exo_startM=models["exo"]["starts_M"],
                                exo_ext3=models["exo"]["ext_aa"], up_to_stop=models["exo"]["up_to_stop"], up_M=models["exo"]["up_M"],
                                orig_status=e["status"], orig_cov=e.get("coverage")))
    return out


def test(runs, gdir, glist, refp, outp):
    allref = read_fa(refp)
    jobs = [(runs, gdir, *line.split(), allref) for line in open(glist)]
    out = []
    with Pool(16) as pool:
        for r in pool.imap(one_genome, jobs):
            out += r
    with open(outp, "w") as f:
        w = csv.DictWriter(f, list(out[0]), delimiter="\t"); w.writeheader(); w.writerows(out)


if __name__ == "__main__":
    c = sys.argv[1]
    if c == "survey":
        survey(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], int(sys.argv[6]) if len(sys.argv) > 6 else 16)
    else:
        test(*sys.argv[2:7])
