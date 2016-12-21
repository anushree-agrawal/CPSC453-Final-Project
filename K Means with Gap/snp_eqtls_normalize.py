#!/usr/bin/python3

'''
Holly Zhou, CPSC 453
Purpose: Normalizes the reduced framingham dataset by dividing specified rows by 
their maximum value; prevents overweighting of large values; in the command line
specify whether want to build file for first (get_partition) or second part of 
the k_means (k_means_gap)
'''

import pandas as pd
import sys
import random 

if len(sys.argv) != 2:
	print("incorrect number of arguments")
	sys.exit()

if sys.argv[1] == 'clusters':
	snp_data = pd.read_table('data_files/raw_reduct.txt', sep = '\t', index_col = 'miRNA_FHS')

elif sys.argv[1] == 'kmeans':
	snp_data = pd.read_table('data_files/raw_reduct.txt', sep = '\t', index_col = 0)

df = pd.DataFrame(snp_data)

def normalize(names):
	raw = df[names]
	return raw / raw.max()

# norm_distances = normalize(['abs_dist_btwn_SNP_and_miRNA(kb)'])
# norm_dist_pval = normalize(['abs_dist_btwn_SNP_and_miRNA(kb)', 'Pval'])
# norm_dist_pval_h2q = normalize(['abs_dist_btwn_SNP_and_miRNA(kb)', 'Pval', 'h2q'])
# norm_tval_dist = normalize(['Tvalue', 'abs_dist_btwn_SNP_and_miRNA(kb)'])
norm_tstat_beta = normalize(['Tvalue', 'Estimate'])
norm_maf_beta = normalize(['MAF', 'Estimate'])
norm_maf_beta_tstat = normalize(['MAF', 'Estimate', 'Tvalue'])
norm_maf_tstat = normalize(['MAF', 'Tvalue'])
norm_maf_tstat_beta_dist = normalize(['MAF', 'Tvalue', 'Estimate', 'abs_dist_btwn_SNP_and_miRNA(kb)'])
norm_maf_tstat_beta_dist_h2q = normalize(['MAF', 'Tvalue', 'Estimate', 'abs_dist_btwn_SNP_and_miRNA(kb)', 'h2q'])
norm_h2q_beta_maf = normalize(['h2q', 'Estimate', 'MAF'])
norm_h2q_beta_tstat = normalize(['h2q', 'Estimate', 'Tvalue'])
# rand = random.sample(range(0,len(dist_pval_h2q) - 1),1000)
# trunc_norm_dist_pval_h2q = norm_dist_pval_h2q.iloc[rand]

# for k_means_w_gap
if sys.argv[1] == 'clusters':
	print("building cluster files")
	# norm_dist_pval.to_csv('data_files/norm_dist_pval_h2q.txt', sep='\t', index=True, header=False)
	norm_tstat_beta.to_csv('data_files/new_cluster_norm_tstat_beta.txt', sep='\t', index=True, header=False) 
	norm_maf_beta.to_csv('data_files/new_cluster_norm_maf_beta.txt', sep='\t', index=True, header=False)
	norm_maf_beta_tstat.to_csv('data_files/new_cluster_norm_maf_beta_tstat.txt', sep='\t', index=True, header=False)
	norm_maf_tstat.to_csv('data_files/new_cluster_norm_maf_tstat.txt', sep='\t', index=True, header=False)
	norm_maf_tstat_beta_dist.to_csv('data_files/new_cluster_norm_maf_tstat_beta_dist.txt', sep='\t', index=True, header=False)
	norm_maf_tstat_beta_dist_h2q.to_csv('data_files/new_cluster_norm_maf_tstat_beta_dist_h2q.txt', sep='\t', index=True, header=False)
	norm_h2q_beta_maf.to_csv('data_files/new_cluster_norm_h2q_beta_maf.txt', sep='\t', index=True, header=False)
	norm_h2q_beta_tstat.to_csv('data_files/new_cluster_norm_h2q_beta_tstat.txt', sep='\t', index=True, header=False)

# for get_partition
elif sys.argv[1] == 'kmeans':
	print("building kmeans files")
	# norm_distances.to_csv('norm_distances_test.txt', sep='\t', index=False, header=False)
	# norm_dist_pval.to_csv('data_files/norm_dist_pval.txt', sep='\t', index=False, header=False)
	# norm_dist_pval_h2q.to_csv('data_files/norm_dist_pval_h2q.txt', sep='\t', index=False, header=False)
	# norm_tval_dist.to_csv('data_files/norm_tval_dist.txt', sep='\t', index=False, header=False)
	norm_tstat_beta.to_csv('data_files/norm_tstat_beta.txt', sep='\t', index=False, header=False) 
	norm_maf_beta.to_csv('data_files/norm_maf_beta.txt', sep='\t', index=False, header=False)
	norm_maf_beta_tstat.to_csv('data_files/norm_maf_beta_tstat.txt', sep='\t', index=False, header=False)
	norm_maf_tstat.to_csv('data_files/norm_maf_tstat.txt', sep='\t', index=False, header=False)
	norm_maf_tstat_beta_dist.to_csv('data_files/norm_maf_tstat_beta_dist.txt', sep='\t', index=False, header=False)
	norm_maf_tstat_beta_dist_h2q.to_csv('data_files/norm_maf_tstat_beta_dist_h2q.txt', sep='\t', index=False, header=False)
	norm_h2q_beta_maf.to_csv('data_files/norm_h2q_beta_maf.txt', sep='\t', index=False, header=False)
	norm_h2q_beta_tstat.to_csv('data_files/norm_h2q_beta_tstat.txt', sep='\t', index=False, header=False)

else:
	print("incorrect argument")

