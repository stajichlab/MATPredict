import sys
L=[];n=0;cur=0;name=None;ends={}
seqs={}
for line in open(sys.argv[1]):
    if line.startswith('>'):
        if name: L.append(cur)
        name=line[1:].split()[0]; cur=0
    else:
        s=line.strip(); cur+=len(s); n+=s.upper().count('N')
if name: L.append(cur)
L.sort(reverse=True); t=sum(L); c=0
for x in L:
    c+=x
    if c>=t/2: n50=x; break
print(f"contigs={len(L)} total={t} N50={n50} largest={L[0]} Ns={n}")
