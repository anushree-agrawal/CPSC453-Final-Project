#!/usr/bin/python3

"""
Holly Zhou, CPSC 453
Purpose: given file with snpID and corresponding SNP locations, check if any 
SNPs are in enhancer regions; pipe the output into an outfile
"""

import pandas as pd
import sys

raw_data = pd.read_table('data_files/enhancers_blood.txt', sep='\t', index_col=None)
# c1 = pd.read_table('data_files/c1_analysis_norm_maf_tstat_beta_dist.txt', sep='\t', index_col=0)
# c2 = pd.read_table('data_files/c2_analysis_norm_maf_tstat_beta_dist.txt', sep='\t', index_col=0)

c1 = pd.read_table('data_files/ch_6_high_snps.txt', sep='\t', index_col=0)
c2 = pd.read_table('data_files/ch_6_high_snps_2.txt', sep='\t', index_col=0)

df = pd.DataFrame(raw_data)
df_c1 = pd.DataFrame(c1)
df_c2 = pd.DataFrame(c2)

df = df.groupby('chromosome')
df_c1 = df_c1.groupby('chr.SNP')
df_c2 = df_c2.groupby('chr.SNP')

def check_range(name, cluster, enhancer):
	enhancers = []
	for ix, row in enhancer.iterrows():
		# print (row)
		if row['start'] < cluster < row['end']:
			# print snpID if found in enhancer region 
			print(name)
			enhancers.append(name)
	return enhancers

print("c1")
c1_enhancers = {}
for name, group in df_c1:
	# check if the snp in df_c1 falls in one of the regions of df.get_group(name) 
	print(name)
	g = df.get_group(name)
	c1_enhancers[name] = group.apply(lambda row: check_range(row.name, row['SNP.pos'], g), axis=1)

print("c2")
c2_enhancers = {}
for name, group in df_c2:
	g = df.get_group(name)
	c2_enhancers[name] = group.apply(lambda row: check_range(row.name, row['SNP.pos'], g), axis=1)

