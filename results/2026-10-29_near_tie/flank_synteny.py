"""Flank check for the family_ambiguous loci of the v8 replay: which family's flank genes sit beside the locus?

For each ambiguous locus cut the contig window locus +/- WINDOW bp, search the flank proteins of both families
with miniprot (SLA2, CHA1, APN2, DIC1 from the MATsc, MATyl and MATog records; PAP1, OBP1, PIK1 from the MTL records) and report the nearest hit
of each gene with its distance to the locus edge and its identity. MATsc-type neighbourhood: SLA2 (or CHA1) next to the locus; MTL-type: PIK1/OBP1/PAP1.
usage: python3 flank_synteny.py <yeast_routing dir> <db dir> <library dir> <out.tsv>   (run on the cluster with miniprot on PATH)"""
import collections, glob, gzip, json, os, subprocess, sys, tempfile
import yaml

R, DB, LIB, OUT = sys.argv[1:5]
WINDOW = 40000
FLANK = {"sla2": "SLA2", "cha1": "CHA1", "apn2": "APN2", "dic1": "DIC1", "PAP1": "PAP1", "OBP1": "OBP1", "PIK1": "PIK1"}

def flank_refs(path):
    """One protein per (gene, record) from the curated records, named <gene>|<record>."""
    n = 0
    with open(path, "w") as out:
        for faa in glob.glob(f"{DB}/Ascomycota/*/*/proteins.faa"):
            rec = os.path.basename(os.path.dirname(faa))
            if not any(f in rec for f in ("_MATsc_", "_MTL_", "_MATyl_", "_MATog_")):
                continue
            name = seq = None
            for line in open(faa):
                if line.startswith(">"):
                    if name in FLANK and seq: out.write(f">{FLANK[name]}|{rec}\n{seq}\n"); n += 1
                    name = next((x[5:] for x in line.strip().split("|") if x.startswith("name=")), None); seq = ""
                else:
                    seq += line.strip()
            if name in FLANK and seq: out.write(f">{FLANK[name]}|{rec}\n{seq}\n"); n += 1
    return n

def genome_seq(asm, contig):
    want, seq, on = contig, [], False
    with gzip.open(f"{LIB}/{asm}.fa.gz", "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if on: break
                on = line[1:].split()[0] == want
            elif on:
                seq.append(line.strip())
    return "".join(seq)

work = tempfile.mkdtemp(dir=R)
refs = os.path.join(work, "flank_refs.faa"); print("flank reference proteins:", flank_refs(refs), flush=True)
rows = []
for d in sorted(glob.glob(f"{R}/out_v8/runs/*")):
    asm = os.path.basename(d); f = d + "/detection_report.yaml"
    if not os.path.exists(f): continue
    y = yaml.safe_load(open(f))
    for l in y.get("detected") or []:
        if not l.get("family_ambiguous"): continue
        seq = genome_seq(asm, l["contig"])
        lo, hi = max(1, l["start"] - WINDOW), min(len(seq), l["end"] + WINDOW)
        fa = os.path.join(work, "win.fa"); open(fa, "w").write(f">w\n{seq[lo - 1:hi]}\n")
        gff = subprocess.run(["miniprot", "--gff", "-t", "2", fa, refs], capture_output=True, text=True).stdout
        best = {}
        for g in gff.splitlines():
            c = g.split("\t")
            if len(c) < 9 or c[2] != "mRNA": continue
            a = dict(x.split("=", 1) for x in c[8].split(";") if "=" in x)
            gene = a.get("Target", "").split()[0].split("|")[0]
            ident = float(a.get("Identity", 0)) * 100; s0, e0 = lo - 1 + int(c[3]), lo - 1 + int(c[4])
            dist = 0 if (s0 <= l["end"] and e0 >= l["start"]) else min(abs(s0 - l["end"]), abs(e0 - l["start"]))
            if gene not in best or ident > best[gene][0]: best[gene] = (round(ident, 1), dist, s0, e0, c[6])
        rows.append(dict(genome=asm, family=l["family"].split(":")[1], idiomorph=l.get("idiomorph"), contig=l["contig"], start=l["start"], end=l["end"],
                         identity=l["family_ambiguous"]["identities"], **{k: (f"{v[0]}@{v[1]}" if (v := best.get(k)) else "") for k in sorted(set(FLANK.values()))}))
cols = ["genome", "family", "idiomorph", "contig", "start", "end"] + sorted(set(FLANK.values()))
with open(OUT, "w") as out:
    out.write("\t".join(cols) + "\n")
    for r in rows: out.write("\t".join(str(r[c]) for c in cols) + "\n")
print(len(rows), "loci")
