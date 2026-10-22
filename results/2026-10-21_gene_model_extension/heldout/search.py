import urllib.request, urllib.parse, json, time, re, sys
def esearch(term, retmax=5000):
    u="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"+urllib.parse.urlencode(dict(db='protein',term=term,retmax=retmax,retmode='json'))
    j=json.load(urllib.request.urlopen(u,timeout=60)); time.sleep(0.4); return int(j['esearchresult']['count']), j['esearchresult']['idlist']
terms={
 'MAT1-1-1 gene':'("MAT1-1-1"[Gene Name]) AND Fungi[Organism]',
 'MAT1-2-1 gene':'("MAT1-2-1"[Gene Name]) AND Fungi[Organism]',
 'title MAT1-1-1':'("MAT1-1-1"[Title] OR "mating type protein 1-1-1"[Title] OR "MAT1-1-1 protein"[Title]) AND Fungi[Organism]',
 'title MAT1-2-1':'("MAT1-2-1"[Title] OR "mating type protein 1-2-1"[Title] OR "MAT1-2-1 protein"[Title]) AND Fungi[Organism]',
 'refseq excluded combined':'("MAT1-1-1"[Title] OR "MAT1-2-1"[Title] OR "MAT1-1-1"[Gene Name] OR "MAT1-2-1"[Gene Name]) AND Fungi[Organism] NOT refseq[filter]',
}
ids=set()
for k,t in terms.items():
    n,l=esearch(t); print(k,n); ids|=set(l)
print('union',len(ids)); open('ids.txt','w').write('\n'.join(sorted(ids)))
