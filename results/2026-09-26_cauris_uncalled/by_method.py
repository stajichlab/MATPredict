import csv,yaml,os,collections,re
S='../2026-09-26_serinales_all_882aa01/runs'
tab=collections.defaultdict(collections.Counter); ex=collections.defaultdict(list)
w=open('method_calls.tsv','w'); w.write('asmid\tmethod\tbioproject\tcontigs\tcall\n')
for r in csv.DictReader(open('asm_methods.tsv'),delimiter='\t'):
    meth=re.sub(r'[ ;,].*','',r['method'] or 'NA').lower()
    p=f"{S}/{r['asmid']}/detection_report.yaml"
    if not os.path.exists(p): call='no_report'
    else:
        y=yaml.safe_load(open(p)); det=y.get('detected') or []
        call='+'.join(sorted({d['idiomorph'] for d in det})) or 'uncalled'
    tab[meth][call]+=1
    w.write('\t'.join([r['asmid'],r['method'],r['bioproject'],r['contigs'],call])+'\n')
for m in sorted(tab,key=lambda k:-sum(tab[k].values())):
    print(f"{m:25s} {sum(tab[m].values()):4d} {dict(tab[m])}")
