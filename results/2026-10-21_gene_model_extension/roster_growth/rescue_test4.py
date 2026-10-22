"""Stage 0b rescue test. python rescue_test.py TRUTH_TSV WINDOWS_FA WINDOWS_TSV REFS_FA OUT_TSV [MAXINTRON]
For every curated MAT1-1-1/MAT1-2-1 gene with a deposited exon structure: withhold its record (radius `record`) or its genus (`genus`), truncate every
remaining reference of that gene to its first F of residues (a simulated fast-evolving tail: the alignment can only reach F of the reference), run exonerate as
production does (`protein2genome --refine region`, all references in one query, highest-scoring model kept), then extend the model in frame to the next stop and compare
with the deposited CDS end. Exonerate and the extension are the real tools; nothing here is scored by homology of the tail."""
import csv, subprocess, sys, tempfile, os
truth = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
meta = {r["key"]: r for r in csv.DictReader(open(sys.argv[3]), delimiter="\t")}
MAXI = sys.argv[6] if len(sys.argv) > 6 else "5000"   # search.py MAX_INTRON_BP
def rd(f):
    d, k = {}, None
    for l in open(f):
        l = l.strip()
        if l.startswith(">"): k = l[1:]; d[k] = ""
        elif k: d[k] += l
    return d
W, R = rd(sys.argv[2]), rd(sys.argv[4])
aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"; B = "TCAG"
CODON = {B[i//16]+B[(i//4)%4]+B[i%4]: aa[i] for i in range(64)}
comp = str.maketrans("ACGTacgt", "TGCAtgca"); rc = lambda s: s.translate(comp)[::-1]
tr = lambda s: "".join(CODON.get(s[i:i+3], "X") for i in range(0, len(s)-2, 3))
genus = lambda sp: sp.split("_")[0]

def extend(seq, exons, strand, maxbp=3000):
    """exons: genomic (1-based closed) spans of a CDS in the tool frame. Returns (protein, ext_aa, state, new_end_coordinate_of_last_codon_before_stop, next_bases)."""
    ex = sorted(exons)
    order = ex[::-1] if strand == "-" else ex
    cds = "".join(rc(seq[a-1:b]) if strand == "-" else seq[a-1:b] for a, b in order)
    rest = len(cds) % 3
    cds_in = cds[:len(cds)-rest]
    if strand == "+":
        end = max(b for _, b in ex); tail = seq[end - rest:end - rest + maxbp] if rest else seq[end:end+maxbp]
        # leftover bases are the last `rest` bases of the CDS; walk codons starting at the first leftover base
        tail = seq[end-rest:end+maxbp]
    else:
        start = min(a for a, _ in ex); tail = rc(seq[max(0, start-1-maxbp):start-1+rest])
    ext = tr(tail); stop = ext.find("*")
    prot = tr(cds_in)
    if stop < 0: return prot, None, "no_stop", rest
    return prot, stop, "stop", rest

def last_sense_end(exons, strand):  # 3' genomic coordinate of the last coding base
    return min(a for a, _ in exons) if strand == "-" else max(b for _, b in exons)

out = csv.writer(open(sys.argv[5], "w"), delimiter="\t")
out.writerow("key gene order species radius frac n_refs model_found strand model_aa model_internal_stops closed_by_alignment ext_aa ext_state truth_aa truth_end_nostop model_end gap_before_aa final_end delta_aa exact donor_GT identity score best_ref q_start q_end truth_introns_after_edge".split())
for t in truth:
    key = f"{t['record']}|{t['gene']}"
    if key not in W: continue
    seq = W[key]; off = int(meta[key]["seq_start"]) - 1
    ex_t = sorted((int(a)-off, int(b)-off) for a, b in (e.split("-") for e in t["exons"].split(";")))
    strand_t = t["strand"]
    tsum = sum(b-a+1 for a, b in ex_t)
    tprot = int(t["prot_aa"]); includes_stop = (tsum == 3*tprot+3)
    tend = last_sense_end(ex_t, strand_t)
    tend_nostop = (tend + 3 if strand_t == "-" else tend - 3) if includes_stop else tend
    sp_t = [h.split("|")[2] for h in R if h.split("|")[0] == t["record"]]
    sp_t = sp_t[0] if sp_t else "NA"
    for radius in ("genus",):
        for frac in (1.0, 0.7, 0.5, 0.3):
            refs = {h: s for h, s in R.items() if h.split("|")[1] == t["gene"] and h.split("|")[0] != t["record"]
                    and (radius == "record" or genus(h.split("|")[2]) != genus(sp_t))}
            if not refs: continue
            with tempfile.TemporaryDirectory() as td:
                with open(td+"/q.fa", "w") as q:
                    for h, s in refs.items(): q.write(f">{h.split('|')[0]}\n{s[:max(30, int(round(frac*len(s))))]}\n")
                open(td+"/t.fa", "w").write(f">w\n{seq}\n")
                try:
                    o = subprocess.run(["exonerate", "--model", "protein2genome", "--query", td+"/q.fa", "--target", td+"/t.fa", "--maxintron", MAXI,
                                        "--refine", "region", "--showtargetgff", "yes", "--showalignment", "no", "--showvulgar", "yes"], capture_output=True, text=True, timeout=600)
                except subprocess.TimeoutExpired:
                    out.writerow([key, t["gene"], t["order"], sp_t, radius, frac, len(refs), "timeout"] + [""]*20); continue
            recs = []
            vul = [l.split() for l in o.stdout.splitlines() if l.startswith("vulgar:")]
            for l in o.stdout.splitlines():
                f = l.split("\t")
                if len(f) < 9: continue
                if f[2] == "gene":
                    try: sc = float(f[5])
                    except ValueError: sc = 0
                    recs.append([sc, f, [], vul[len(recs)] if len(vul) > len(recs) else None])
                elif f[2] == "exon" and recs: recs[-1][2].append((int(f[3]), int(f[4])))
            # keep the highest-scoring model that overlaps the truth locus (production selects by raw score within a window around the locus)
            recs = [r for r in recs if r[2] and min(a for a, _ in r[2]) <= max(b for _, b in ex_t) and max(b for _, b in r[2]) >= min(a for a, _ in ex_t)]
            if not recs:
                out.writerow([key, t["gene"], t["order"], sp_t, radius, frac, len(refs), "none"] + [""]*20); continue
            sc, gl, exs, vg = max(recs, key=lambda r: r[0]); strand = gl[6]
            qid = vg[1] if vg else ''; qend = int(vg[3]) if vg else ''; qstart = int(vg[2]) if vg else ''
            ident = [p.split(" ")[1] for p in gl[8].split(" ; ") if p.startswith("identity ")]
            prot, ext, state, rest = extend(seq, exs, strand)
            mend = last_sense_end(exs, strand)
            closed = (ext == 0)
            fin = (mend - rest + 3*(ext or 0)) if strand == "+" else (mend + rest - 3*(ext or 0))
            sign = 1 if strand_t == "+" else -1
            gap = sign*(tend_nostop - mend)/3 if strand == strand_t else None
            delta = sign*(fin - tend_nostop)/3 if (strand == strand_t and ext is not None) else None
            # truth introns lying beyond the model end, in transcript direction
            tin = 0
            for (a1, b1), (a2, b2) in zip(ex_t[:-1], ex_t[1:]):
                if strand_t == "+" and a2 - 1 > mend and b1 >= mend - 3*3: tin += 1
                if strand_t == "-" and b1 + 1 < mend and a2 <= mend + 3*3: tin += 1
            donor = ""
            if ext not in (None, 0):
                nb = seq[mend:mend+2] if strand == "+" else rc(seq[mend-3:mend-1])
                donor = "GT" if nb == "GT" else "other"
            out.writerow([key, t["gene"], t["order"], sp_t, radius, frac, len(refs), "yes", strand, len(prot), prot[:-1].count("*"), closed, ext, state, tprot, tend_nostop, mend,
                          None if gap is None else round(gap, 1), fin, None if delta is None else round(delta, 1), None if delta is None else (abs(delta) < 0.5), donor, ident[0] if ident else "", sc, qid, qstart, qend, tin])
