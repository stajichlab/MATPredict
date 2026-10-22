import re, csv
S='/private/tmp/claude-501/-Users-jstajich-projects-MATPredict/519c0d9b-c67c-487f-bb14-2f7c2115d61e/scratchpad'
sel=list(csv.DictReader(open(f'{S}/heldout/selected.tsv'),delimiter='\t'))
want={r['protein']:r for r in sel}
txt=open(f'{S}/heldout/gp_all.txt').read(); recs=re.split(r'^//\s*$',txt,flags=re.M)
seqs={}
for r in recs:
    r=r.strip('\n'); v=re.search(r'^VERSION\s+(\S+)',r,re.M)
    if not v or v.group(1) not in want or 'ORIGIN' not in r: continue
    seqs[v.group(1)]=''.join(re.findall(r'\d+ ([a-z ]+)',r[r.index('ORIGIN'):])).replace(' ','').upper()
print('sequences found',len(seqs),'of',len(want))
with open('refs_aug.faa','w') as out:
    for l in open(f'{S}/s0b/refs.faa'): out.write(l)
    for p,r in want.items():
        sp=r['organism'].replace(' ','_').replace('|','_')
        out.write(f">HO_{p}|{r['which']}|{sp}\n{seqs[p]}\n")
# intron positions (aa) for the deposits, from their exon structures
tr=list(csv.DictReader(open(f'{S}/heldout/truth_heldout.tsv'),delimiter='\t'))
rows=list(csv.DictReader(open(f'{S}/s0b/ref_introns_all.tsv'),delimiter='\t'))
for t in tr:
    ex=[tuple(map(int,e.split('-'))) for e in t['exons'].split(';')]; cum=0; pos=[]
    for a,b in ex[:-1]: cum+=b-a+1; pos.append(round(cum/3,1))
    rows.append(dict(record=t['record'],gene=t['gene'],strand=t['strand'],nexons=len(ex),cds_aa=round(sum(b-a+1 for a,b in ex)/3,1),intron_aa=';'.join(map(str,pos))))
w=csv.DictWriter(open('ref_introns_aug.tsv','w',newline=''),rows[0].keys(),delimiter='\t'); w.writeheader(); w.writerows(rows); print(len(rows),'reference structures')
# combined truth + windows (34 curated + 106 deposits)
def cat(files,out):
    with open(out,'w') as o:
        for f in files: o.write(open(f).read())
cat([f'{S}/s0b/truth_windows.fa',f'{S}/heldout/heldout_windows.fa'],'all_windows.fa')
h=open(f'{S}/s0b/truth.tsv').read().split('\n')[0]
t1=list(csv.DictReader(open(f'{S}/s0b/truth.tsv'),delimiter='\t')); t2=list(csv.DictReader(open(f'{S}/heldout/truth_heldout.tsv'),delimiter='\t'))
cols=list(t1[0].keys()); allt=[{c:r.get(c,'') for c in cols} for r in t1+t2]
m1=list(csv.DictReader(open(f'{S}/s0b/truth_windows.tsv'),delimiter='\t')); m2=list(csv.DictReader(open(f'{S}/heldout/heldout_windows.tsv'),delimiter='\t'))
w=csv.DictWriter(open('all_windows.tsv','w',newline=''),m1[0].keys(),delimiter='\t'); w.writeheader(); w.writerows(m1+m2)
n=5
for i in range(n):
    w=csv.DictWriter(open(f'truth_chunk{i}.tsv','w',newline=''),cols,delimiter='\t'); w.writeheader(); w.writerows(allt[i::n])
print(len(allt),'genes in',n,'chunks')
