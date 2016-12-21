#!/usr/bin/python2

'''
Holly Zhou, CPSC 453
Written with assistance of Declan Clarke, Gerstein Lab at Yale University
'''

import numpy as np
import random
import matplotlib.pyplot as plt
import sys
import os
from scipy.stats import mode

'''
Given a an output file from get_partition.py (example below) and the full Framingham list of mirQTLs, this script 
uses the "represeantitve" eQTLs in the get_partition.py to process data for all corresponding eQTLs.

Script usage:
python extract_all_eQTLs_from_seed_lists.py output__cluster_norm_tstat_beta__2.txt /Users/admin/Desktop/rsch/eQTLs/h_zhou/course_materials/Supp_data_4___cis_miR-eQTLs_at_FDR_less_than_0p1__append.txt data_for_clust_1.txt data_for_clust_2.txt


'''

'''
##  If the user did not enter the correct number of args, notify the user and terminate the program
if (len(sys.argv) != 5):
   sys.stderr.write("input expected:\n")
   sys.stderr.write("   <input_file_of_representative_eQTLs> \n")
   sys.stderr.write("   <full_eqtls file> \n")
   sys.stderr.write("   <data_for_clust_1> \n")
   sys.stderr.write("   <data_for_clust_2> \n")
   sys.exit()
'''

##  Assign input variables
input_file_of_representative_eQTLs = open(sys.argv[1], "r")
full_eqtls = open(sys.argv[2], "r")
data_for_clust_1 = open(sys.argv[3], "w")
data_for_clust_2 = open(sys.argv[4], "w")



list_of_all_miRNAs = list()
cluster___2___repr_eQTLs = {}
for line in input_file_of_representative_eQTLs:
    ln_elms = list()
    ln_elms = line.split()
    if len(ln_elms) > 0:
        if ln_elms[0] == "Cluster":
            #clust_num = int(ln_elms[1][:-1])

            stripped_string = ""
            for char in line:
                if (char != ":") and (char != "[") and (char != "'") and (char != ",") and (char != "]"):
                    stripped_string = stripped_string + char
            stripped_string_elms = list()
            stripped_string_elms = stripped_string.split()

            repr_eQTLs = list()
            i = 2
            while i < len(stripped_string_elms):
                eQTL = stripped_string_elms[i]
                if eQTL not in repr_eQTLs:
                    repr_eQTLs.append(eQTL)
                    if eQTL not in list_of_all_miRNAs:
                        list_of_all_miRNAs.append(eQTL)
                i += 1
            clust_num = int(stripped_string_elms[1])

            cluster___2___repr_eQTLs[clust_num] = repr_eQTLs
input_file_of_representative_eQTLs.close()

## fill dictionaries w/data
for line in full_eqtls:
    if line[0:2] == "rs":
        ln_elms = list()
        ln_elms = line.split()

        rs_ID_val = ln_elms[0]
        maf = float(ln_elms[3])
        beta = float(ln_elms[4])
        t_val = float(ln_elms[6])
        p_val = float(ln_elms[7])
        h2q = float(ln_elms[8])
        fdr = float(ln_elms[9])
        snp_funct = ln_elms[13]
        mir_type = ln_elms[21]
        abs_dist_btwn_SNP_and_miRNA = float(ln_elms[26])
        neg_log_p_val = float(ln_elms[28])

        miRNA = ln_elms[1]
        if miRNA in cluster___2___repr_eQTLs[1]:
            data_for_clust_1.write(rs_ID_val + "\t" + str(maf) + "\t" + str(beta) + "\t" + str(t_val) + "\t" + str(p_val) + "\t" + str(h2q) + "\t" + str(fdr) + "\t" + snp_funct + "\t" + mir_type + "\t" + str(abs_dist_btwn_SNP_and_miRNA) + "\t" + str(neg_log_p_val) + "\n")
        elif miRNA in cluster___2___repr_eQTLs[2]:
            data_for_clust_2.write(rs_ID_val + "\t" + str(maf) + "\t" + str(beta) + "\t" + str(t_val) + "\t" + str(p_val) + "\t" + str(h2q) + "\t" + str(fdr) + "\t" + snp_funct + "\t" + mir_type + "\t" + str(abs_dist_btwn_SNP_and_miRNA) + "\t" + str(neg_log_p_val) + "\n")
        else:
            print "ERROR: cannot get data for mirna " + miRNA
            sys.exit("ERROR: cannot get data for mirna " + miRNA)

full_eqtls.close()
data_for_clust_1.close()
data_for_clust_2.close()


print "\n\n fin \n\n"



