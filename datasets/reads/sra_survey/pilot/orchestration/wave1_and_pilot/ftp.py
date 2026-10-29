import csv
W="/bigdata/stajichlab/jstajich/projects/MATPredict/results/sra_reads_typing/wave_inputs/asm_summary_fungi.txt"
hdr=None; want=("GCF_001477535","GCF_001477545","GCF_000002945","GCF_001661265","GCA_005281805","GCF_000","GCA_051599425","GCA_001929475")
for l in open(W):
    if l.startswith("#assembly_accession"): hdr=l[1:].rstrip("\n").split("\t"); continue
    if l.startswith("#"): continue
    p=l.rstrip("\n").split("\t"); d=dict(zip(hdr,p))
    o=d["organism_name"]
    if any(d["assembly_accession"].startswith(w) for w in want[:5]+want[6:]) or o.startswith(("Protomyces","Saitoella","Neolecta","Taphrina deformans","Pneumocystis jirovecii","Pneumocystis carinii","Schizosaccharomyces pombe")):
        if o.startswith(("Schizosaccharomyces pombe")) and not d["assembly_accession"].startswith("GCF_000002945"): continue
        print(d["assembly_accession"],"|",o,"|",d["assembly_level"],"|",d["refseq_category"],"|",d["ftp_path"])
