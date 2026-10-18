import sys
import pandas as pd
pd.set_option('display.width', 280); pd.set_option('display.max_columns', 60); pd.set_option('display.max_rows', 300); pd.set_option('display.max_colwidth', 70)
D = sys.argv[1]
L = pd.read_csv(D + '/loci_audit_flagged.tsv', sep='\t', low_memory=False)
G = pd.read_csv(D + '/genome_audit_verdict.tsv', sep='\t', low_memory=False)
c = L[L.called]
what = sys.argv[2]
if what == 'disagree':
    x = c[c.vote_vs_profile == 'disagree']
    print(x[['genome', 'species', 'idiomorph', 'confidence', 'locus_class', 'win_gene', 'own_id', 'own_cov', 'prof_set', 'prof_call', 'prof_margin', 'prof_order_call', 'genes_found', 'pfam_MAT1-1-1', 'pfam_MAT1-2-1']].to_string())
elif what == 'genus':
    g = sys.argv[3]
    x = c[c.genus == g]
    print(x.groupby(['species', 'idiomorph', 'locus_class', 'confidence']).size())
    print(x[['own_id', 'own_cov', 'win_bits', 'a_bits', 'a_id', 'a_cov', 'h_bits', 'h_id', 'h_cov', 'prof_alpha_on_alpha', 'prof_hmg_on_hmg']].describe().T)
    print(x.genes_found.value_counts().head(8))
    print(x.groupby(['idiomorph', 'vote_vs_profile', 'domain_support']).size())
    print(x[['pfam_MAT1-1-1', 'pfam_MAT1-2-1']].apply(lambda s: s.fillna('-').str.slice(0, 40)).value_counts().head(8))
elif what == 'lowid':
    x = c[(c.own_id < 40) | (c.locus_class == 'idiomorph_gene_only')]
    print(x.groupby(['order', 'genus', 'idiomorph', 'locus_class']).size().sort_values(ascending=False).head(40))
elif what == 'sample':
    g = sys.argv[3]
    x = c[c.genus == g].head(int(sys.argv[4]))
    print(x[['genome', 'species', 'idiomorph', 'confidence', 'locus_class', 'win_gene', 'win_bits', 'own_id', 'own_cov', 'a_bits', 'a_id', 'a_cov', 'h_bits', 'h_id', 'h_cov', 'prof_call', 'prof_margin', 'res_pairs', 'genes_found', 'n_loci_genome']].to_string())
