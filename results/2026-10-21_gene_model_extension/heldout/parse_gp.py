import re, csv, collections, glob, yaml
txt=open('gp_all.txt').read()
recs=re.split(r'^//\s*$',txt,flags=re.M)
refs_g=set(); ref_sp=set()
for md in glob.glob('/Users/jstajich/projects/MATPredict-roster/db/*/*/*/metadata.yaml'):
    if '/candidates/' in md or '/classifiers/' in md: continue
    d=yaml.safe_load(open(md)); sp=(d.get('organism') or {}).get('species') or ''
    if sp: refs_g.add(sp.split()[0]); ref_sp.add(sp)
rows=[]
for r in recs:
    r=r.strip('\n')
    if not r.startswith('LOCUS'): continue
    acc=re.search(r'^ACCESSION\s+(\S+)',r,re.M); ver=re.search(r'^VERSION\s+(\S+)',r,re.M)
    org=re.search(r'^\s+ORGANISM\s+(.+)$',r,re.M); 
    lin=re.search(r'ORGANISM.*?\n((?:\s{12}.+\n)+)',r)
    prod=re.search(r'/product="([^"]+)',r.replace('\n                     ',' ')); gene=re.search(r'/gene="([^"]+)',r)
    flat=r.replace('\n                     ',' ')
    cb=re.search(r'/coded_by="([^"]+)"',flat)
    dbs=re.search(r'^DBSOURCE\s+(.+)$',r,re.M)
    ln=int(re.search(r'LOCUS\s+\S+\s+(\d+) aa',r).group(1)) if re.search(r'LOCUS\s+\S+\s+(\d+) aa',r) else 0
    cbs=cb.group(1) if cb else ''
    nacc=re.findall(r'([A-Z]{1,6}_?\d+\.\d+):',cbs)
    spans=re.findall(r'[<>]?(\d+)\.\.[<>]?(\d+)',cbs)
    rows.append(dict(protein=ver.group(1) if ver else '',organism=org.group(1).strip() if org else '',gene=gene.group(1) if gene else '',product=prod.group(1) if prod else '',aa=ln,
        coded_by=cbs,nuc=nacc[0] if nacc else '',nexons=len(spans),partial=('<' in cbs or '>' in cbs),complement='complement' in cbs,
        refseq=ver.group(1).startswith(('NP_','XP_','YP_')) if ver else False,nuc_refseq=(nacc[0].startswith(('NC_','NW_','NZ_','NT_','XM_','NM_')) if nacc else False)))
print('records',len(rows))
c=collections.Counter()
for r in rows:
    g=(r['gene'] or r['product']).upper()
    r['which']='MAT1-1-1' if '1-1-1' in g or '1-1-1' in r['product'] else ('MAT1-2-1' if '1-2-1' in r['product'] or '1-2-1' in g else '')
    r['genus']=r['organism'].split()[0] if r['organism'] else ''
    r['in_ref_genus']=r['genus'] in refs_g
rows=[r for r in rows if r['which']]
print('with a MAT1-1-1 / MAT1-2-1 name:',len(rows))
step=[('has coded_by',lambda r:r['coded_by']),('not a RefSeq protein',lambda r:not r['refseq']),('nucleotide not RefSeq/genome-annotation (direct deposit)',lambda r:r['nuc'] and not r['nuc_refseq']),
      ('CDS has >= 2 exons',lambda r:r['nexons']>=2),('complete CDS (no < or >)',lambda r:not r['partial']),('genus not in the 55-reference set',lambda r:not r['in_ref_genus'])]
cur=rows
for n,f in step:
    cur=[r for r in cur if f(r)]; print(f'{n:60s} {len(cur)}')
# one per species
seen={};
for r in cur: seen.setdefault((r['organism'],r['which']),r)
print('distinct species x gene:',len(seen),' genera:',len({r["genus"] for r in seen.values()}))
w=csv.DictWriter(open('candidates.tsv','w',newline=''),list(next(iter(seen.values())).keys()),delimiter='\t'); w.writeheader(); w.writerows(seen.values())
print(collections.Counter(r['genus'] for r in seen.values()).most_common(40))
