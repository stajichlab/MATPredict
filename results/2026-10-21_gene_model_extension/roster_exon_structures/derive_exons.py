"""Derive exon structures for curated references lacking them, from the NCBI protein record's coded_by. python derive_exons.py missing_refs.tsv refs.faa OUT_TSV"""
import csv, re, sys, time, urllib.request
rows=list(csv.DictReader(open(sys.argv[1]),delimiter='\t'))
P={};k=None
for l in open(sys.argv[2]):
    if l[0]=='>': k=l[1:].strip().split('|')[:2]; k='|'.join(k); P[k]=''
    elif k: P[k]+=l.strip()
out=[]
for r in rows:
    pid=(r['protein'] or '').replace('ncbi_protein:','')
    url=f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id={pid}&rettype=gp&retmode=text"
    txt=None
    for a in range(3):
        try: txt=urllib.request.urlopen(url,timeout=60).read().decode(); break
        except Exception: time.sleep(2)
    time.sleep(0.4)
    rec=dict(record=r['record'],gene=r['gene'],protein=pid,accession=r['accession'],strand=r['strand'],gstart=r['gstart'],gend=r['gend'])
    if not txt: out.append({**rec,'status':'fetch_failed'}); continue
    m=re.search(r'/coded_by="([^"]+)"',txt.replace('\n                     ',''))
    if not m: out.append({**rec,'status':'no_coded_by'}); continue
    cb=m.group(1); comp='complement' in cb
    acc=re.findall(r'([A-Z]{1,2}_?\d+\.\d+):',cb)
    spans=[(int(a),int(b)) for a,b in re.findall(r'(\d+)\.\.(\d+)',cb)]
    plen=len(P.get(f"{r['record']}|{r['gene']}",''))
    cds=sum(b-a+1 for a,b in spans)
    # transcript order: complement(join(a..b,c..d)) lists ascending coordinates; transcript order is descending
    tr=spans[::-1] if comp else spans
    same_acc = (len(set(acc))<=1) and (not acc or acc[0]==r['accession'])
    status='ok' if (plen and cds in (3*plen,3*plen+3) and same_acc) else ('length_mismatch' if plen else 'no_protein_in_db')
    if status=='ok' and not same_acc: status='accession_differs'
    out.append({**rec,'status':status,'coded_by':cb,'coded_by_accession':';'.join(sorted(set(acc))),'strand_cb':'-' if comp else '+','nexons':len(spans),
                'exons_transcript_order':';'.join(f'{a}-{b}' for a,b in tr),'cds_nt':cds,'db_protein_aa':plen,'includes_stop':cds==3*plen+3})
cols=['record','gene','protein','accession','strand','gstart','gend','status','coded_by','coded_by_accession','strand_cb','nexons','exons_transcript_order','cds_nt','db_protein_aa','includes_stop']
w=csv.DictWriter(open(sys.argv[3],'w',newline=''),cols,delimiter='\t',extrasaction='ignore',restval=''); w.writeheader(); w.writerows(out)
import collections; print(collections.Counter(o['status'] for o in out))
for o in out: print(o['record'][:34].ljust(34),o['gene'],o['status'].ljust(16),o.get('nexons',''),o.get('strand_cb',''),o['strand'],o.get('coded_by','')[:70])
