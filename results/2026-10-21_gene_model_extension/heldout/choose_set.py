import re, csv, json, time, urllib.request, urllib.parse, collections, random
rows=list(csv.DictReader(open('candidates_direct.tsv'),delimiter='\t'))
txt=open('gp_all.txt').read(); recs=re.split(r'^//\s*$',txt,flags=re.M)
lin={}
for r in recs:
    r=r.strip('\n'); v=re.search(r'^VERSION\s+(\S+)',r,re.M)
    if not v: continue
    m=re.search(r'ORGANISM\s+.+\n((?:\s{12}.+\n?)+)',r)
    if m:
        l=' '.join(m.group(1).split()); lin[v.group(1)]=l
def order(l):
    for t in re.split(r';\s*',l.replace('.','')): 
        if t.endswith('ales'): return t
    return 'NA'
for r in rows: r['order']=order(lin.get(r['protein'],''))
# nucleotide lengths
accs=sorted({r['nuc'] for r in rows}); slen={}
for i in range(0,len(accs),150):
    b=accs[i:i+150]
    u="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"; data=urllib.parse.urlencode(dict(db='nuccore',id=','.join(b),retmode='json')).encode()
    j=json.load(urllib.request.urlopen(urllib.request.Request(u,data=data),timeout=120)); time.sleep(0.5)
    for k,v in j['result'].items():
        if k=='uids': continue
        slen[v['accessionversion']]=int(v['slen'])
rows=[r for r in rows if r['nuc'] in slen]
for r in rows:
    spans=[(int(a),int(b)) for a,b in re.findall(r'[<>]?(\d+)\.\.[<>]?(\d+)',r['coded_by'])]
    r['lo']=min(a for a,b in spans); r['hi']=max(b for a,b in spans); r['slen']=slen[r['nuc']]; r['down']=r['slen']-r['hi'] if not r['complement'] else r['lo']-1
    r['up']=r['lo']-1 if not r['complement'] else r['slen']-r['hi']
print('candidates with a known sequence length:',len(rows))
for thr in (0,300,1000,2000,3000): print(f' downstream flank (3 prime side) >= {thr} bp: {sum(r["down"]>=thr for r in rows)}')
good=[r for r in rows if r['down']>=1000]
print('with >= 1000 bp downstream:',len(good),' genera:',len({r['genus'] for r in good}),' pubmed:',sum(str(r['pubmed'])=='True' for r in good))
print(collections.Counter(r['nexons'] for r in good), collections.Counter(r['which'] for r in good))
random.seed(1)
# select: max 2 per genus (one per gene where possible), prefer PubMed-cited, spread exon counts
by=collections.defaultdict(list)
for r in good: by[r['genus']].append(r)
sel=[]
for g,l in sorted(by.items()):
    l.sort(key=lambda r:(str(r['pubmed'])!='True',-int(r['nexons']),random.random()))
    seen=set(); 
    for r in l:
        k=(r['which'],r['nexons'])
        if k in seen: continue
        seen.add(k); sel.append(r)
        if len([x for x in sel if x['genus']==g])>=2: break
print('selected',len(sel),' genera',len({r['genus'] for r in sel}),collections.Counter(r['nexons'] for r in sel),collections.Counter(r['which'] for r in sel),'pubmed',sum(str(r['pubmed'])=='True' for r in sel))
print(collections.Counter(r['order'] for r in sel).most_common())
w=csv.DictWriter(open('selected.tsv','w',newline=''),list(sel[0].keys()),delimiter='\t'); w.writeheader(); w.writerows(sel)
