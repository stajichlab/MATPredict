import sys,subprocess,os,tempfile
sys.path.insert(0,"/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-26_v070_submit")
import taph_mat_candidates as t
P="/bigdata/stajichlab/jstajich/projects/MATPredict"; I=P+"/results/sra_reads_typing/wave_inputs/genomes"; O=P+"/results/2026-10-28_taphrinomycotina_mat_candidates"
tests=[("GCA_051599425.1_Sai_col_YB2330_v1","AP039747.1",177119,198616),("GCF_001661265.1_Saico1","NW_017567017.1",185699,210006),("GCA_001929475.1_NeoIirr1.0","LXFE01000411.1",5949,25065)]
rlen={}
for h,s in t.read_fasta(O+"/references.faa"): rlen[h.split()[0]]=len(s)
for g,c,s,e in tests:
    seqs=t.load_genome(f"{I}/{g}.fa.gz"); w=seqs[c][max(0,s-1):e]
    d=tempfile.mkdtemp(dir="/tmp"); open(d+"/w.fa","w").write(f">{c}\n{w}\n")
    for tag,args in (("default",[]),("relaxed",["--outs","0.1","--outc","0.1","-k","5","-n","1"])):
        p=subprocess.run([P+"/.pixi/envs/default/bin/miniprot","--gff","-N","30"]+args+[d+"/w.fa",O+"/references.faa"],capture_output=True,text=True)
        hits=t.parse_gff(p.stdout); best={}
        for h in hits:
            ql=rlen.get(h["target"],1); cov=(h["tend"]-h["tstart"]+1)/ql
            k=(h["start"]//500)
            if k not in best or h["ident"]*cov>best[k][0]: best[k]=(h["ident"]*cov,h,cov)
        print(g,c,tag,len(hits),"alignments;",len(best),"loci")
        for _,(sc,h,cov) in sorted(best.items(),key=lambda kv:-kv[1][0])[:6]:
            print("   ",h["target"].split("|")[1],"id %.2f cov %.2f"%(h["ident"],cov),"win",h["start"]+s-1,h["end"]+s-1,"CDS",len(h["cds"]))
