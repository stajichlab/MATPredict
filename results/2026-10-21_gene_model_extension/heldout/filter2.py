import re, csv, collections
txt=open('gp_all.txt').read(); recs=re.split(r'^//\s*$',txt,flags=re.M)
pub={}; sub={}
for r in recs:
    r=r.strip('\n'); v=re.search(r'^VERSION\s+(\S+)',r,re.M)
    if not v: continue
    pub[v.group(1)]=bool(re.search(r'PUBMED\s+\d+',r)); sub[v.group(1)]=('Direct Submission' in r) 
rows=list(csv.DictReader(open('candidates.tsv'),delimiter='\t'))
direct=re.compile(r'^[A-Z]{1,2}\d{5,6}\.\d+$')
kept=[r for r in rows if direct.match(r['nuc'])]
print('candidates (species x gene) from the previous filters:',len(rows))
print('nucleotide is an old-style direct deposit (1-2 letters + 5-6 digits):',len(kept),' genera:',len({r['genus'] for r in kept}))
withpub=[r for r in kept if pub.get(r['protein'])]
print('  of which the protein record cites a PubMed article:',len(withpub),' genera:',len({r['genus'] for r in withpub}))
print('  accession patterns of the rest (genome annotation):',collections.Counter(re.sub(r'\d','N',r['nuc'])[:8] for r in rows if not direct.match(r['nuc'])).most_common(5))
for r in kept: r['pubmed']=pub.get(r['protein'],False)
w=csv.DictWriter(open('candidates_direct.tsv','w',newline=''),list(kept[0].keys()),delimiter='\t'); w.writeheader(); w.writerows(kept)
print(collections.Counter(r['nexons'] for r in kept))
print(collections.Counter(r['which'] for r in kept))
print(sorted(collections.Counter(r['genus'] for r in kept).items(),key=lambda x:-x[1])[:45])
