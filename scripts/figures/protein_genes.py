"""Draw-ready gene boxes on contigs by tblastn of curated A. fumigatus record proteins (same method for every track).

genes_on(contigs: {name: seq}, db_dir) -> {contig: [(start, end, strand, label)]}
Queries: COX13, APN2, SLA2, MAT1-1-1 and the MAT1-2-1 remnant from the A1163 record; MAT1-2-1 and MAT1-2-4 from the Af293 record.
A gene is drawn when its tblastn HSPs (same contig and strand, merged when < 3 kb apart) cover >= 40% of the protein at >= 70% identity.
The remnant is drawn only where MAT1-1-1 is within 3 kb (otherwise its hit is the 3' part of the real MAT1-2-1).
"""
import subprocess
import tempfile
from pathlib import Path

ORDER = [("746128_a1163_MAT_MAT1-1", {"COX13": "COX13", "APN2": "APN2", "SLA2": "SLA2", "MAT1-1-1": "MAT1-1-1", "MAT1-2-1": "remnant"}),
         ("746128_af293_MAT_MAT1-2", {"MAT1-2-1": "MAT1-2-1", "MAT1-2-4": "MAT1-2-4"})]


def query_proteins(db_dir):
    q = {}
    for rec, names in ORDER:
        name = None
        for line in open(Path(db_dir) / rec / "proteins.faa"):
            if line.startswith(">"):
                nm = dict(kv.split("=", 1) for kv in line[1:].strip().split("|")[1:])["name"]
                name = names.get(nm)
                if name and not any(k == name for k in q):
                    q[name] = []
                else:
                    name = None if (name and name in q and q[name]) else name
            elif name:
                q[name].append(line.strip())
    return {k: "".join(v) for k, v in q.items()}


def genes_on(contigs, db_dir, min_cov=0.40, min_pid=70.0, merge=3000):
    prot = query_proteins(db_dir)
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        (d / "q.faa").write_text("".join(f">{k}\n{v}\n" for k, v in prot.items()))
        (d / "s.fa").write_text("".join(f">{k}\n{v}\n" for k, v in contigs.items()))
        subprocess.run(["makeblastdb", "-in", str(d / "s.fa"), "-dbtype", "nucl", "-out", str(d / "db")], check=True, capture_output=True)
        out = subprocess.run(["tblastn", "-query", str(d / "q.faa"), "-db", str(d / "db"), "-evalue", "1e-8", "-outfmt",
                              "6 qseqid sseqid pident length qstart qend sstart send"], check=True, capture_output=True, text=True).stdout
    hsps = {}
    for line in out.splitlines():
        q, s, pid, ln, qs, qe, ss, se = line.split("\t")
        if float(pid) < min_pid:
            continue
        strand = "+" if int(ss) < int(se) else "-"
        hsps.setdefault((q, s, strand), []).append((min(int(ss), int(se)), max(int(ss), int(se)), int(qs), int(qe)))
    found = {}
    for (q, s, strand), hs in hsps.items():
        hs.sort(); groups, cur = [], [hs[0]]
        for h in hs[1:]:
            if h[0] - cur[-1][1] <= merge:
                cur.append(h)
            else:
                groups.append(cur); cur = [h]
        groups.append(cur)
        for g in groups:
            cov = len({p for h in g for p in range(h[2], h[3] + 1)}) / len(prot[q])
            if cov >= min_cov:
                found.setdefault(s, []).append((min(h[0] for h in g), max(h[1] for h in g), strand, q, cov))
    res = {}
    for s, gs in found.items():
        keep = []
        for g in gs:
            if g[3] == "remnant" and not any(o[3] == "MAT1-1-1" and abs(o[0] - g[1]) < 3000 or abs(g[0] - o[1]) < 3000 and o[3] == "MAT1-1-1" for o in gs):
                continue
            keep.append(g)
        # one best (largest coverage) per gene name per contig
        best = {}
        for g in keep:
            if g[3] not in best or g[4] > best[g[3]][4]:
                best[g[3]] = g
        res[s] = [(g[0], g[1], g[2], g[3]) for g in best.values()]
    return res
