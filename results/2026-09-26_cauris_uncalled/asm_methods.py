import csv,json,urllib.request,time,sys
rows=list(csv.reader(open('cauris_samples.csv')))
accs=[r[0].split('_')[0]+'_'+r[0].split('_')[1] for r in rows]
out=open('asm_methods.tsv','w'); out.write('asmid\taccession\tlevel\tmethod\ttech\tbioproject\tcontigs\n')
m={}
for i in range(0,len(accs),100):
    batch=accs[i:i+100]
    url='https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/'+','.join(batch)+'/dataset_report?page_size=100'
    for t in range(4):
        try: d=json.load(urllib.request.urlopen(url,timeout=120)); break
        except Exception as e: time.sleep(5)
    for r in d.get('reports',[]):
        a=r['assembly_info']; s=r['assembly_stats']
        m[r['accession']]=(a.get('assembly_level'),a.get('assembly_method'),a.get('sequencing_tech'),a.get('bioproject_accession'),s.get('number_of_contigs'))
    time.sleep(0.5)
for r,a in zip(rows,accs):
    v=m.get(a,('NA',)*5); out.write('\t'.join([r[0],a]+[str(x) for x in v])+'\n')
