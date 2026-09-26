import sys,gzip,re
# usage: region.py fa.gz contig start end  -> print N runs and base composition of window
fa,ctg,s,e=sys.argv[1],sys.argv[2],int(sys.argv[3]),int(sys.argv[4])
seq=[];keep=False
for line in gzip.open(fa,'rt'):
    if line.startswith('>'):
        if keep: break
        keep=line[1:].split()[0]==ctg; continue
    if keep: seq.append(line.strip())
seq=''.join(seq); w=seq[s-1:e]
print(f"{ctg} len={len(seq)} window {s}-{e}: N={w.upper().count('N')} lower={sum(c.islower() for c in w)}")
for m in re.finditer(r'[Nn]{10,}',w): print(f"  Nrun {s+m.start()}-{s+m.end()-1} ({m.end()-m.start()} bp)")
if len(sys.argv)>5: open(sys.argv[5],'w').write(f">{ctg}_{s}_{e}\n{w}\n")
