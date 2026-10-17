"""Across the v0.6.0 Ascomycota campaign: loci whose cross-match resolution picked an accessory/HMG-named MAT1-1 gene over MAT1-2-1.
usage: campaign_pattern.py REPORT_DIR genomes.tsv out.tsv"""
import sys,glob,os,csv,collections,yaml
L=getattr(yaml,"CSafeLoader",yaml.SafeLoader)
gen={r['genome']:r for r in csv.DictReader(open(sys.argv[2]),delimiter='\t')}
rows=[]
for p in glob.glob(sys.argv[1]+'/**/detection_report.yaml',recursive=True):
    g=os.path.basename(os.path.dirname(p)); doc=yaml.load(open(p),Loader=L)
    for c in doc.get('detected') or []:
        if c.get('family')!='Ascomycota:MAT': continue
        ev={e['gene']:e for e in c.get('gene_evidence') or []}
        for r in c.get('idiomorph_resolutions') or []:
            if r['loser']=='MAT1-2-1' and r['winner'] in('MAT1-1-2','MAT1-1-3','MAT1-1-4','MAT1-1-5'):
                lb=ev.get('MAT1-2-1',{}).get('bitscore'); wb=ev.get(r['winner'],{}).get('bitscore')
                gi=gen.get(g,{})
                rows.append(dict(genome=g,species=gi.get('species'),order=gi.get('order'),label=c['idiomorph'],winner=r['winner'],loser_bits=lb,winner_bits=wb,core_alpha='MAT1-1-1' in ev,alpha_bits=ev.get('MAT1-1-1',{}).get('bitscore'),basis=r['basis'],wref=ev.get(r['winner'],{}).get('reference_record')))
w=csv.DictWriter(open(sys.argv[3],'w'),fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows);print(len(rows))
