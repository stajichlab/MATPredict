"""S. schoenii chromosome IV (CM143129.1): StringTie transcripts and TransDecoder ORFs in the MAT region, identified by phmmer
against the curated Saccharomycotina MAT proteins (refs.faa). Usage: chr4_models.py data/ refs.faa A B"""
import re, sys, collections, pyhmmer
from pyhmmer.easel import SequenceFile, Alphabet, TextSequence, DigitalSequenceBlock
d, refs, A, B = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
tx = {}
for l in open(f"{d}/stringtie.gtf"):
    f = l.split("\t")
    if f[0] != "CM143129.1" or f[2] not in ("transcript", "exon"): continue
    tid = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
    t = tx.setdefault(tid, {"strand": f[6], "exons": [], "span": None, "cov": None})
    if f[2] == "transcript":
        t["span"] = (int(f[3]), int(f[4])); t["cov"] = float(re.search(r'cov "([^"]+)"', f[8]).group(1))
    else: t["exons"].append((int(f[3]), int(f[4])))
keep = {k: v for k, v in tx.items() if v["span"] and A <= v["span"][0] and v["span"][1] <= B}
pep = {}
name = None
for l in open(f"{d}/longest_orfs.pep"):
    if l.startswith(">"):
        name = l[1:].split()[0]; hdr = l; pep[name] = [hdr.strip(), ""]
    elif name: pep[name][1] += l.strip()
alpha = Alphabet.amino()
q = []
for tid in keep:
    for n, (h, s) in pep.items():
        if n.startswith(tid + ".p") and s:
            q.append(TextSequence(name=n.encode(), sequence=s.rstrip("*")).digitize(alpha))
with SequenceFile(refs, digital=True, alphabet=alpha) as sf: db = sf.read_block()
out = {}
for hits in pyhmmer.phmmer(q, db, E=1e-2):
    qn = hits.query.name if isinstance(hits.query.name, str) else hits.query.name.decode()
    if len(hits): 
        h = hits[0]; out[qn] = (h.name if isinstance(h.name, str) else h.name.decode(), h.score, h.evalue)
for tid, t in sorted(keep.items(), key=lambda kv: kv[1]["span"]):
    print(tid, t["strand"], t["span"], f"cov={t['cov']:.0f}", "exons=" + ",".join(f"{a}-{b}" for a, b in sorted(t["exons"])))
    for n, (h, s) in pep.items():
        if n.startswith(tid + ".p") and s:
            hit = out.get(n)
            print("   ", h[:110].split(" ")[0], h.split()[-2] if False else "", f"len={len(s.rstrip('*'))}", ("-> %s %.0f bits E=%.1g" % (hit[0].split("|")[0][-45:] + ("|" + hit[0].split("name=")[1].split("|")[0] if "name=" in hit[0] else ""), hit[1], hit[2])) if hit else "")
