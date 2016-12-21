#!/usr/bin/python2

'''
Holly Zhou, CPSC 453
Written with assistance of Declan Clarke, Gerstein Lab at Yale University
'''

from __future__ import division
import numpy as np
import sys
import os
import random
from collections import Counter
# import matplotlib.pyplot as plt
# import os.path

###   Script to identify optimal K value in k-means clustering, as described in:
###       http://datasciencelab.wordpress.com/2013/12/12/clustering-with-k-means-in-python/
###       http://datasciencelab.wordpress.com/2013/12/27/finding-the-k-in-k-means-clustering/

###   USAGE: 
###   python k_means_w_gap.py ./sample_input_files/ errors_log.txt 10 10000 0.50 10 > output.txt


#  If the user did not enter the correct number of args, notify the user and terminate the program
if (len(sys.argv) != 7):
   sys.stderr.write("Usage: " + sys.argv[0] + " <inputFile> <outputDir>\n")
   sys.stderr.write("where:\n")
   sys.stderr.write("   <inputFileDir> is a directory containing the rmsd_vals* (or qh_vals*) files generated from STAMP runs\n")
   sys.stderr.write("   <errors> is a file to which potential errors/failures are printed \n")
   sys.stderr.write("   <num_runs> is an int specifying how many times the algorithm should be run on a matrix. Note that the k-means algorithm w/the gap statistic has stochastic elements, so multiple runs should be performed to search for a 'consensus k value'\n")
   sys.stderr.write("   <num_attempts_to_try_gap_calc_on_small_clusts> is an int specifying the number of attempts to perform -- sometimes (especially for small matices), the algorithm fails to converge proprely. This int specifies the number of attempts to perform before giving up (unlikely to be a problem in large matrices)\n")
   sys.stderr.write("   <fraction_for_mode_thresh> is the threshold used to determine the 'consensus k value' -- if a particlar k value (for example, 3) is acheived in at least fraction_for_mode_thresh of the specified num_runs, then this k value is taken as the consensus value. \n")
   sys.stderr.write("   <B> Create B reference datasets -- for context, see the line that has: BWkbs = np.zeros(B) within the WordPress site (links above)\n")
   sys.exit()


##  Assign input variables
input_dir = sys.argv[1]
errors_file = sys.argv[2]
num_runs = int(sys.argv[3])
num_attempts_to_try_gap_calc_on_small_clusts = int(sys.argv[4])
fraction_for_mode_thresh = float(sys.argv[5])
B = int(sys.argv[6])


#####  Basic K-means clustering algorithm implemented in the next 4 functions
## Put each point in X into a cluster based on its proximity to one of the mu values, 
## and return the respective clusters (ie, cluster indeces and members in each cluster)
def cluster_points(X, mu):
    clusters  = {}
    for x in X:
        bestmukey = min([(i[0], np.linalg.norm(x-mu[i[0]])) for i in enumerate(mu)], key=lambda t:t[1])[0]
        try:
            clusters[bestmukey].append(x)
        except KeyError:
            clusters[bestmukey] = [x]
    return clusters

def reevaluate_centers(clusters):
    newmu = []
    keys = sorted(clusters.keys())
    for k in keys:
        #print "rand keys: " + str(k)
        newmu.append(np.mean(clusters[k], axis = 0))
    return newmu

def has_converged(mu, oldmu):
    return (set([tuple(a) for a in mu]) == set([tuple(a) for a in oldmu]))

def find_centers(X, K):
    ## Assign the K oldmu's and the K mu's to random points within X
    oldmu = random.sample(X, K)  # Initialize to K random centers
    mu = random.sample(X, K)
    while not has_converged(mu, oldmu):  ##  Based on mu and oldmu (the seeds) -- first assign points to centers (and thus clusters). Then, find the centers of THOSE clusters, and keep iterating (ie, standard k-means clustering)
        oldmu = mu  ##  oldmu is the previous cluster centers
        clusters = cluster_points(X, mu)  ##  Assign all points in X to clusters
        mu = reevaluate_centers(clusters) ##  Based on cluster assignments, find new centers
    return(mu, clusters)



#######  The purpose of the next several functions is to find the optimal K in k-means clustering by using the gap statistic 
###  Given a set of mu's and the clusters (ie, points assinged to each mu) -- calculate the Wk statistic 
###  (qualitatively, this is a measure of the intra-cluster variance, a metric we'd like to minimize for small k)
def Wk(mu, clusters):
    K = len(mu)
    return sum([np.linalg.norm(mu[i]-c)**2/(2*len(c)) for i in range(K) for c in clusters[i]])

##  To calc Wk* (ie, the null), we need to set up points that are randomply distributed in such a way that 
##  the randomly distributed points lie within the same overall interval as the oritinal dataset
def bounding_box(X):
    dimensions = X.shape
    num_points = dimensions[0]
    num_dimensions = dimensions[1]
    min_bounds = list() ## the lenght of this list will be equal to the number of dimensions
    max_bounds = list() ## '''
    d = 0
    while (d < num_dimensions):
        min_bounds.append(min(X,key=lambda a:a[d])[d])       
        max_bounds.append(max(X,key=lambda a:a[d])[d])
        d += 1
    return(min_bounds, max_bounds)

def gap_statistic(X):
    dimensions = X.shape
    num_points = dimensions[0]
    num_dimensions = dimensions[1]
    min_bounds, max_bounds = bounding_box(X)
    if (num_points >= 10):
        range_of_k_vals_to_test = 8
    else:
        range_of_k_vals_to_test = num_points - 1
    ks = range(1,range_of_k_vals_to_test)
    Wks = np.zeros(len(ks))
    Wkbs = np.zeros(len(ks))
    sk = np.zeros(len(ks))
    for indk, k in enumerate(ks):
        mu, clusters = find_centers(X,k)
        Wks[indk] = np.log(Wk(mu, clusters))
        BWkbs = np.zeros(B)  # Making B reference datasets
        for i in range(B):
            Xb = []
            for n in range(len(X)):
                rand_coords = list()
                rc = 0
                while(rc < num_dimensions):
                    rand_coords.append(random.uniform(min_bounds[rc],max_bounds[rc]))
                    rc += 1
                Xb.append(rand_coords)
            Xb = np.array(Xb)
            mu, clusters = find_centers(Xb,k)
            BWkbs[i] = np.log(Wk(mu, clusters))
        Wkbs[indk] = sum(BWkbs)/B
        sk[indk] = np.sqrt(sum((BWkbs-Wkbs[indk])**2)/B)
    sk = sk*np.sqrt(1+1/B)
    return(ks, Wks, Wkbs, sk)

def find_mode_in_list(data_list):
    data = Counter(data_list)
    data.most_common()   # Returns all unique items and their counts
    mode_and_frequency = data.most_common(1)  # Returns the highest occurring item
    return mode_and_frequency



##  Start running calculations at this step (ie, the code below calls and uses the functions defined above)
##  Make a list of the non-empty MDS files
tot_num_files = 0;
files_list = list() 
files_list.append(input_dir)
# for filename in os.listdir(input_dir):
#     full_file_name = input_dir + "/" + filename
#     file_size = os.path.getsize(full_file_name)
#     if file_size > 0:
#         files_list.append(filename)

##  Iterate through all files in the directory
errors_file_hndle = open(errors_file, "w")
tot_num_corrupt_files = 0;
for filename in files_list:
    X = np.empty([])   ## make sure X is empty before starting
    many_attempts_to_run_gap_stat = "NO"  ## starting out, there are NO attempts made yet
    # full_file_name = input_dir + "/" + filename
    full_file_name = files_list[0]

    ##  Extract the coordinates, and store data in a list
    matrix_file = open(full_file_name, "r")
    file_lines = list()
    for line in matrix_file:
        file_lines.append(line)
    num_data_points = len(file_lines)
    print "\n\nFile:  " + str(filename)

    ##  Perform singletons k-means clustering search
    buckets = []
    for line in file_lines:
        line_elems = list()
        line_elems = line.split()
        vals = list()
        for line_element in line_elems:
            val = float(line_element)
            vals.append(val)
        buckets.append(vals)
    X = np.array(buckets)
    dimensions = X.shape

    ## Run many iterations to get vector of k values
    k_values = list()
    confidently_assigned_k = "CANNOT_CONFIDENTLY_ASSIGN_K_VAL" ## assume this by default -- ie, the matrix is crappy until proven otherwise
    b = 0  ## b is just an index to keep track of the number of runs performed on a given matrix
    while b < num_runs:
        num_attempts = 0
        while True:
            try:
                ks, logWks, logWkbs, sk = gap_statistic(X)
                ######   Determine the optimal k value
                k_ind = 0
                gaps = list()
                while k_ind < len(ks):
                    gap_val = logWkbs[k_ind] - logWks[k_ind]
                    gaps.append(gap_val)
                    k_ind += 1
                k_ind = 0
                gap_minus_sks = list()
                while k_ind < len(gaps):
                    gap_minus_sk = gaps[k_ind] - sk[k_ind]
                    gap_minus_sks.append(gap_minus_sk)
                    k_ind += 1
                ####   determine the optimal k here
                optimal_k = 0
                k_ind = 0
                while k_ind < (len(gaps)-1):
                    actual_k = k_ind + 1
                    value_for_k = gaps[k_ind] - gap_minus_sks[k_ind+1]
                    if value_for_k > 0:
                        optimal_k = actual_k
                        break # just to break the loop
                    k_ind += 1
                k_values.append(optimal_k)
                '%s %d %s %d' % ('attempts', num_attempts, 'k val', optimal_k)
                break
            except:
                num_attempts += 1
                if num_attempts >= num_attempts_to_try_gap_calc_on_small_clusts:
                    many_attempts_to_run_gap_stat = "YES"
                    break
        b += 1
    mode_and_frequency = find_mode_in_list(k_values)
    if len(mode_and_frequency) > 0:  ## if len==0, then there was a convergence issue
        mode = mode_and_frequency[0][0]
        mode_frequency = mode_and_frequency[0][1]
        mode_fraction = mode_frequency/len(k_values)
        print "    List of k values from simulations:  " + str(k_values)
        print "    num_data_points:  " + str(num_data_points)
        print "    Mode:  " + str(mode)
        print "    Frequency:  " + str(mode_frequency)
        print "    Fraction: " + str(mode_fraction)
        if mode_fraction >= fraction_for_mode_thresh:
            confidently_assigned_k = mode
    print "\n\t\tFinal K value determined for this matrix (using the gap statistic):  " + str(confidently_assigned_k) + "\n\n\n"

    ## If something went wrong, print full description to output error log file:
    if many_attempts_to_run_gap_stat == "YES":
        errors_file_hndle.write("Clust w/" + str(num_data_points) + " members seems to be too small -- many attempts to run gap calc\n")
        tot_num_corrupt_files += 1

    matrix_file.close()
    print "\n\n\n"
    tot_num_files += 1


print "tot num files:  " + str(tot_num_files)
print "tot_num_corrupt_files:  " + str(tot_num_corrupt_files)

print "\n\n\n\n   ---  Run Complete  ---   \n\n\n"


