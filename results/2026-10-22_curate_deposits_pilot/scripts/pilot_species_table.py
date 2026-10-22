import yaml, csv, sys
from pathlib import Path
P=Path(sys.argv[1])
for r in csv.DictReader(open(P/'replay_genomes.tsv'),delimiter='\t'):
    if r['group']!='pilot_species': continue
    out=[]
    for s in ('base','pilot'):
        d=yaml.safe_load(open(P/'runs'/s/r['asm']/'detection_report.yaml'))
        calls=[]
        for l in d.get('detected') or []:
            core=';'.join(f"{g['gene']}:{g.get('reference_record','').split('_')[0]}:{g.get('identity')}/{g.get('coverage')}:{len(g.get('exons') or [])}ex" for g in l['gene_evidence'] if g['gene'] in('MAT1-1-1','MAT1-2-1') and g.get('status','').startswith('polished'))
            calls.append(f"{l['idiomorph']}/{l['confidence']}/{l['locus_class']} [{core}]")
        out.append(' | '.join(calls) or 'no call')
    print(f"{r['species']}\t{r['asm']}\tBASE: {out[0]}\tPILOT: {out[1]}")
