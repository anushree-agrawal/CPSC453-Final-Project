#!/usr/bin/python3

"""
Holly Zhou, CPSC 453
Purpose: reduces the original framingham dataset from ~10000 lines to ~80 lines;
removes the linkage disequilibrium confounding factor by removing repeated SNPs 
with the same miRNA, keeping the SNP with the lowest P-value 
"""

import pandas as pd

snp_data = pd.read_table('data_files/SNP_eQTLs_Raw.txt', sep = '\t', index_col = 0)
df = pd.DataFrame(snp_data)
df = df.sort_values('miRNA_FHS', ascending=False)
out_file = 'data_files/SNP_eQTLs_cut.txt'

g = df.groupby('miRNA_FHS')
for name, group in g:
	# get list of miRNAs with similar distances 
	similar_dist = group[group['SNP.pos'] - group['SNP.pos'].median() < 1000000]
	# save miRNA with smallest Pval 
	min_pval = similar_dist[similar_dist.Pval == min(similar_dist.Pval)]
	# drop similar miRNAs, add back in the min Pval row 
	new = group.drop(group['SNP.pos'][group['SNP.pos'] - group['SNP.pos'].median() < 1000000].index)
	new = new.append(min_pval)
	new.to_csv(out_file, sep='\t', index=True, header=False, mode='a')