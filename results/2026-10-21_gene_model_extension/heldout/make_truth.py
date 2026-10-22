import csv, re
rows=list(csv.DictReader(open('selected.tsv'),delimiter='\t'))
out=[]
for r in rows:
    spans=[(int(a),int(b)) for a,b in re.findall(r'[<>]?(\d+)\.\.[<>]?(\d+)',r['coded_by'])]
    tr=spans[::-1] if r['complement'] in ('True',True) else spans
    out.append(dict(record='HO_'+r['protein'],order=r['order'],species=r['organism'],taxid='',gene=r['which'],accession=r['nuc'],strand='-' if r['complement']=='True' else '+',
                    gstart=min(a for a,b in spans),gend=max(b for a,b in spans),exons=';'.join(f'{a}-{b}' for a,b in tr),nexons=len(spans),prot_aa=r['aa'],status='deposit',tier=''))
w=csv.DictWriter(open('truth_heldout.tsv','w',newline=''),out[0].keys(),delimiter='\t'); w.writeheader(); w.writerows(out); print(len(out),'held-out genes')
