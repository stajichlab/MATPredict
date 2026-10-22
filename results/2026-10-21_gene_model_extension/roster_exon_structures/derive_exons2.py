"""RefSeq proteins: take the CDS feature (protein_id match) from the genomic record slice. python derive_exons2.py"""
import csv, re, time, urllib.request
rows={ (r['record'],r['gene']):r for r in csv.DictReader(open('derived_exons.tsv'),delimiter='\t')}
miss={(r['record'],r['gene']):r for r in csv.DictReader(open('missing_refs.tsv'),delimiter='\t')}
P={};k=None
for l in open('refs.faa'):
    if l[0]=='>': k='|'.join(l[1:].strip().split('|')[:2]); P[k]=''
    elif k: P[k]+=l.strip()
for key,r in list(rows.items()):
    if r['status']!='length_mismatch' or not r['coded_by'].startswith(('XM_','NM_')) : continue
    m=miss[key]; acc=m['accession']; gs,ge=int(m['gstart']),int(m['gend'])
    s0=max(1,gs-200)
    url=f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={acc}&rettype=gb&retmode=text&seq_start={s0}&seq_stop={ge+200}"
    txt=urllib.request.urlopen(url,timeout=90).read().decode(); time.sleep(0.5)
    feat=re.split(r'\n     (?=\S)',txt)
    got=None
    for f in feat:
        if f.startswith('CDS') and f'/protein_id="{r["protein"]}"' in f.replace('\n                     ',''):
            loc=re.match(r'CDS\s+(\S+)',f.replace('\n                     ','')).group(1); got=loc; break
    if not got: rows[key]['status']='refseq_cds_not_found'; continue
    comp='complement' in got
    spans=[(int(a)+s0-1,int(b)+s0-1) for a,b in re.findall(r'(\d+)\.\.(\d+)',got)]
    plen=len(P.get(f"{key[0]}|{key[1]}",'')); cds=sum(b-a+1 for a,b in spans)
    tr=spans[::-1] if comp else spans
    rows[key].update(status='ok' if cds in (3*plen,3*plen+3) else 'length_mismatch',coded_by='genomic CDS '+got,coded_by_accession=acc,strand_cb='-' if comp else '+',nexons=len(spans),
        exons_transcript_order=';'.join(f'{a}-{b}' for a,b in tr),cds_nt=cds,db_protein_aa=plen,includes_stop=cds==3*plen+3)
# partial CDS deposits
for key,r in rows.items():
    if r['status']=='length_mismatch' and '<' in r['coded_by'] or r['status']=='length_mismatch' and '>' in r['coded_by']: r['status']='partial_cds'
cols=list(next(iter(rows.values())).keys())
w=csv.DictWriter(open('derived_exons.tsv','w',newline=''),cols,delimiter='\t'); w.writeheader(); w.writerows(rows.values())
import collections; print(collections.Counter(r['status'] for r in rows.values()))
for r in rows.values():
    if r['status']!='ok': print(r['record'][:34].ljust(34),r['gene'],r['status'].ljust(18),r['nexons'],r['cds_nt'],r['db_protein_aa'],r['coded_by'][:60])
