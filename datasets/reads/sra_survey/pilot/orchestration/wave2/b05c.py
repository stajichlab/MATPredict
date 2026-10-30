import csv
R="/bigdata/stajichlab/jstajich/projects/MATPredict/results/sra_reads_typing/wave_inputs"
rd=lambda p:list(csv.DictReader(open(p,newline=""),delimiter="\t"))
M=[r for r in rd(R+"/fungi_meta.tsv") if r["run_accession"] in ("SRR29713114","SRR29713115","SRR29713116")]
for r in M: print({k:v[:120] for k,v in r.items() if v})
S=[r for r in rd(R+"/studies.tsv") if r["study"]=="PRJNA1132156"]; print(S)
W=[r for r in rd(R+"/fungi_wgs_illumina.tsv") if r["run_accession"] in ("SRR29713114","SRR29713115","SRR29713116")]
for r in W: print(r["run_accession"],r["sample_accession"],r["library_layout"],r["base_count"],r["read_count"],r["center_name"],r["collection_date"])
# other B05.10 runs in the survey
B=[r for r in rd(R+"/fungi_meta.tsv") if r["strain"].upper().replace(" ","") in ("B05.10","B0510") or "B05.10" in r.get("sample_title","")]
print(len(B),"B05.10 runs in the survey:",[(r["run_accession"],r["study_title"][:60]) for r in B][:12])
