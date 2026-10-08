#!/usr/bin/env python3
"""Hit tally in the stretches outside a locus's own genes, for the called loci with >= MIN_BEYOND bp beyond their core.

For each such locus: tblastn (detect's settings: -seg no, default e-value 10, the genome's genetic code) of the run's own
_reference.faa over the contig window; hits that overlap the stretch left of the core (start..core_start-1) or right of it
(core_end+1..end) are counted as strong (e < 1e-3) or weak (e >= 1e-3), only if the hit lies WHOLLY inside the stretch. The hit that sets each span edge is reported.
usage: tally.py RUNS_DIR OUT_TSV [MIN_BEYOND]"""
import glob, gzip, os, subprocess, sys, tempfile
import yaml

runs, out = sys.argv[1], sys.argv[2]
MIN = int(sys.argv[3]) if len(sys.argv) > 3 else 10000
LIB = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes"
STRONG = 1e-3

def contig_seq(path, name):
    seq, on = [], False
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line[0] == ">":
                if on: break
                on = line[1:].split()[0] == name
            elif on:
                seq.append(line.strip())
    return "".join(seq)

rows = []
for f in sorted(glob.glob(f"{runs}/*/detection_report.yaml")):
    g = f.split("/")[-2]; y = yaml.safe_load(open(f)); code = y.get("genetic_code") or 1
    for x in y.get("detected") or []:
        cs = x.get("core_span")
        if not cs or cs["beyond_core_bp"] < MIN: continue
        seq = contig_seq(f"{LIB}/{g}.fa.gz", x["contig"])
        lo, hi = max(1, x["start"] - 2000), min(len(seq), x["end"] + 2000)
        with tempfile.TemporaryDirectory() as d:
            open(f"{d}/c.fa", "w").write(f">c\n{seq[lo-1:hi]}\n")
            subprocess.run(["makeblastdb", "-in", f"{d}/c.fa", "-dbtype", "nucl", "-out", f"{d}/db"], check=True, capture_output=True)
            res = subprocess.run(["tblastn", "-query", f"{runs}/{g}/_reference.faa", "-db", f"{d}/db", "-evalue", "10", "-seg", "no",
                                  "-db_gencode", str(code), "-outfmt", "6 qseqid pident length sstart send evalue bitscore"],
                                 check=True, capture_output=True, text=True).stdout
        hits = []
        for l in res.splitlines():
            q, pid, ln, ss, se, ev, bs = l.split("\t")
            a, b = sorted((int(ss), int(se))); hits.append((a + lo - 1, b + lo - 1, q, float(pid), int(ln), float(ev), float(bs)))
        def tally(a0, b0):
            sel = [h for h in hits if h[0] >= a0 and h[1] <= b0]   # wholly inside the stretch (straddlers belong to the core genes)
            strong = [h for h in sel if h[5] < STRONG]; weak = [h for h in sel if h[5] >= STRONG]
            best = min(sel, key=lambda h: h[5]) if sel else None
            return len(sel), len(strong), len(weak), (best[5] if best else None), best
        L = (x["start"], cs["start"] - 1) if cs["start"] > x["start"] else None
        R = (cs["end"] + 1, x["end"]) if x["end"] > cs["end"] else None
        edge = {}
        for side, pos, key in (("left", x["start"], 0), ("right", x["end"], 1)):
            at = [h for h in hits if h[key] == pos]
            edge[side] = min(at, key=lambda h: h[5]) if at else None
        rows.append(dict(genome=g[:34], family=x["family"].split(":")[1], contig=x["contig"], start=x["start"], end=x["end"], core=(cs["start"], cs["end"]),
                         beyond=cs["beyond_core_bp"], left=(L, tally(*L) if L else None), right=(R, tally(*R) if R else None), edge=edge))
with open(out, "w") as o:
    o.write("genome\tfamily\tcontig\tspan\tcore\tbeyond_bp\tside\tstretch_bp\thits\tstrong\tweak\tbest_e\tbest_hit_in_stretch(query;pos;e;aa;pid)\tedge_hit(query;e;bits;aa;pid)\n")
    for r in rows:
        for side in ("left", "right"):
            st, t = r[side]
            if st is None: continue
            e = r["edge"][side]
            ed = f"{e[2].split('|')[0][:28]}|{e[2].split('|')[-1]};{e[5]:.1e};{e[6]:.0f};{e[4]};{e[3]:.0f}" if e else "no tblastn hit at the edge"
            b = t[4]
            bh = f"{b[2].split('|')[0][:24]}|{b[2].split('|')[-1]};{b[0]}-{b[1]};{b[5]:.1e};{b[4]};{b[3]:.0f}" if b else ""
            o.write(f"{r['genome']}\t{r['family']}\t{r['contig']}\t{r['start']}-{r['end']}\t{r['core'][0]}-{r['core'][1]}\t{r['beyond']}\t{side}\t{st[1]-st[0]+1}\t{t[0]}\t{t[1]}\t{t[2]}\t{'' if t[3] is None else format(t[3], '.1e')}\t{bh}\t{ed}\n")
print("loci tallied:", len(rows))
