import pyhmmer, glob
from pyhmmer.easel import SequenceFile, Alphabet
from pyhmmer.plan7 import HMMFile
R="/bigdata/stajichlab/jstajich/projects/MATPredict/results"
V=f"{R}/2026-10-24_ctg_mtla1_verify"; T=f"{R}/2026-10-24_ctg_mtla1_records/testdb/Ascomycota/Serinales"
W="/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-56e11cd/db/Ascomycota/Serinales"
a=Alphabet.amino()
seqs=[]
def add(path, tag):
    with SequenceFile(path, digital=True, alphabet=a) as sf:
        for s in sf:
            nm=(s.name.decode() if isinstance(s.name,bytes) else s.name).split("|"); n=(tag+":"+(nm[2].replace("name=","") if len(nm)>2 else nm[1] if len(nm)>1 else nm[0]))
            s.name=n.encode(); seqs.append(s)
for rec in sorted(glob.glob(T+"/*_MTL_A")):
    r=rec.split("/")[-1]
    if r.split("_")[0] in ("717740","43968","4924","1303940","45512","45502","45501"): add(rec+"/proteins.faa", r.split("_")[0]+"_"+r.split("_")[1][:6])
for r in ("36911_cl143_MTL_A","498019_b11243_MTL_A","5476_sc5314_MTL_A","498019_b11205_MTL_A"): add(f"{W}/{r}/proteins.faa", r.split("_")[0]+"_"+r.split("_")[1][:6]+"(ref)")
with HMMFile("/srv/projects/db/pfam/current/Pfam-A.hmm") as hf:
    for hits in pyhmmer.hmmscan(seqs, hf, E=1e-3, cpus=4):
        q=hits.query.name if isinstance(hits.query.name,str) else hits.query.name.decode()
        out=[]
        for h in hits:
            if h.included: out.append("%s(%.1e)"%(h.name if isinstance(h.name,str) else h.name.decode(), h.evalue))
        print(q, len(hits.query.sequence) if hasattr(hits.query,"sequence") else "", " ".join(out) or "none")
