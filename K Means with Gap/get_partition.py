#!/usr/bin/python2

"""
Holly Zhou, CPSC 453
Written with assistance of Declan Clarke, Gerstein Lab at Yale University
"""

import numpy as np
import random
import matplotlib.pyplot as plt
import sys
import os
from scipy.stats import mode

'''
Given a K value and a list of labeled points (in n-dimensional space), this script runs the standard k-means
algorithm on the points, thereby clustering them into K clusters. The assignment of points into distinct clusters
is called a "partiton". For example, the pointes {A,B,C,D,E} may be partitioned into the 2 clusters {A,C} & {B,D,E}.
Note that the k-means algorithm is is stochastic in nature. Thus, each time the algorithm is run, it may generate
a different partition. As such, this script runs the k-means algorithm many times, and it ultimately reports to the
user the most commonly-generated partition. Further details and info are provided in the DataScienceLab site:
   https://datasciencelab.wordpress.com/2013/12/12/clustering-with-k-means-in-python/

Script usage:
python get_partition.py input_points_file.txt 3 1000 10000 > output_file_from_k_means.txt

'''


##  If the user did not enter the correct number of args, notify the user and terminate the program
if (len(sys.argv) != 5):
   sys.stderr.write("Usage: " + sys.argv[0] + " <inputFile> <outputDir>\n")
   sys.stderr.write("where:\n")
   sys.stderr.write("   <input_points_file> is a file containing the list of labeled points (in n-dimensional space) to be clustered \n")
   sys.stderr.write("   <K> is the assigned k value -- this clusters data into k clusters, and the value should be obtained from the k-means w/ gap statistic algorithm \n")
   sys.stderr.write("   <num_iterations> is the number of times that this script runs the k-means algorithm -- as a result_from_k_means, a total of num_iterations partitions (not necessarily distinct) are generated \n")
   sys.stderr.write("   <num_attempted_loops_to_build_partition> -- sometimes, the k-means algorithm fails to converge properly. Thus, num_attempted_loops_to_build_partition attempts are made to generate a partition. \n")
   sys.exit()

##  Assign input variables
input_points_file = open(sys.argv[1], "r")
K = int(sys.argv[2])
num_iterations = int(sys.argv[3])
num_attempted_loops_to_build_partition = int(sys.argv[4])



##  The functions below were taken & adapted from the WordPress DataScienceLab site:
##  https://datasciencelab.wordpress.com/2013/12/12/clustering-with-k-means-in-python/
def reevaluate_centers(mu, clusters):
    newmu = []
    keys = sorted(clusters.keys())
    for k in keys:
        newmu.append(np.mean(clusters[k], axis = 0))
    return newmu


def cluster_points(X, mu):
    clusters  = {}
    for x in X:
        bestmukey = min([(i[0], np.linalg.norm(x-mu[i[0]])) \
                    for i in enumerate(mu)], key=lambda t:t[1])[0]
        try:
            clusters[bestmukey].append(x)
        except KeyError:
            clusters[bestmukey] = [x]
    return clusters


def has_converged(mu, oldmu):
    return (set([tuple(a) for a in mu]) == set([tuple(a) for a in oldmu]))


def find_centers(X, K):
    # Initialize to K random centers
    oldmu = random.sample(X, K)
    mu = random.sample(X, K)
    while not has_converged(mu, oldmu):
        oldmu = mu
        # Assign all points in X to clusters
        clusters = cluster_points(X, mu)
        # Reevaluate centers
        mu = reevaluate_centers(oldmu, clusters)
    return(mu, clusters)



##  Read in the list of points from the input file, and store all of their coordinates in "X"
coordinates_list = list()
input_file_lines = list()
point_label___2___coordinates = {}
for point in input_points_file:
    input_file_lines.append(point)
    point_elements = list()
    point_elements = point.split()
    label = point_elements[0]
    num_dimensions = len(point_elements) - 1  ## Note that the number of dimensions = "len(point_elements) - 1" because the first "element" within point_elements is not a numerical data point, but rather a label (for exapmle, "eQTL_4")

    ## Generate a list of all of the coordiante values for the current point
    coordinates_for_point = list()
    dimension = 1
    while dimension <= num_dimensions:
        coordinates_for_point.append(float(point_elements[dimension]))
        dimension += 1
    coordinates_list.append(coordinates_for_point)

    point_label___2___coordinates[label] = coordinates_list

tot_num_points = len(point_label___2___coordinates)
X = np.array([(coordinates_list[i]) for i in range(tot_num_points)])


## Run K-means clustering num_iterations times -- this will produce num_iterations partitions (not necessarily distinct)
iteration = 0
list_of_partition_dictionaries = list()
iteration___2___partition = {}
while iteration < num_iterations:
    ##  First perform K means to build the partition. Note that we 
    ##  may need to attempt this many times (hence the loops here)
    t1 = 0
    while t1 < num_attempted_loops_to_build_partition:
        ##  For reasons that are not 100% clear, the call to find_centers 
        ##  does not always work -- so try it multiple attempts
        t = 0
        while t < num_attempted_loops_to_build_partition:
            try:
                result_from_k_means = find_centers(X, K)
                break
            except:
                t += 1
        if t > num_attempted_loops_to_build_partition - 1:
            print "\n\nERROR  --  UNABLE TO GET 'result_from_k_means' FROM CALL TO function find_centers! -- may need to increase the parameter <num_attempted_loops_to_build_partition> \n\n"
            sys.exit()
        if len(result_from_k_means[1]) == K:
            break
        t1 += 1
    if t1 > num_attempted_loops_to_build_partition - 1:
        sys.exit()

    ##  Determine the list of eQTLs within each cluster. This is stored 
    ##  as the dictionary cluster_index___2___list_of_eQTLs_in_clust
    cluster_index___2___list_of_eQTLs_in_clust = {}
    k = 0
    while k < K:
        ##  In this loop, we're looping through each of the K clusters, and for each 
        ##  cluster, we're determing the associated eQTLs within that cluster.
        cluster_index = k + 1
        list_of_eQTLs_in_clust = list()
        j = 0 #  Each j represents a different set of coordinates within cluster k
        while j < len(result_from_k_means[1][k]):
            ##  Now go through each eQTL in the original input file, and see which one(s) 
            ##  has/have the coordinates that match result_from_k_means[1][k][j]
            for point in input_file_lines:
                num_matching_coords = 0
                point_elements = list()
                point_elements = point.split()
                dimension = 1
                while dimension <= num_dimensions:
                    coordinate_for_point = float(point_elements[dimension])
                    if coordinate_for_point in result_from_k_means[1][k][j]:
                        num_matching_coords += 1
                    dimension += 1
                if num_matching_coords == num_dimensions:  ## Then this eQTL is within cluster k!
                    eQTL_label = point_elements[0]
                    list_of_eQTLs_in_clust.append(eQTL_label)
            j += 1
        cluster_index___2___list_of_eQTLs_in_clust[cluster_index] = list_of_eQTLs_in_clust
        k += 1

    ##  Make sure we have ALL eQTLs assigned to a cluster! --- Note that 
    ##  an error may (potentially) occur because of a float estimation
    cl = 1
    tot_num_points_assigned_to_a_clust = 0
    while cl <= K:
        tot_num_points_assigned_to_a_clust = tot_num_points_assigned_to_a_clust + len(cluster_index___2___list_of_eQTLs_in_clust[cl])
        cl += 1
    if tot_num_points_assigned_to_a_clust != tot_num_points:
        print "\n\nERROR for: " + str(root_molecule) + "  --  NOT EVERY MOLECULE HAS BEEN ASSIGNED TO A CLUSTER! \n\n"
        sys.exit()

    ##  Now that we have a partition, let's "regularze" the ordering of everything in this partition such 
    ##  that an identical partition w/different ordering (the ordering in the dict is arbitrary) is recognized 
    ##  as being the same ordering when calculating the most common partitoin (below). The ordering is performed 
    ##  alphabetically, so alphebetize each element within a cluster, and then alphebitze the clusters
    ##  (with a "cluster name" being assinged based on the first alphebical element in that cluster). For example, 
    ##  the partition:
    ##      {1: ['eQTL_3', 'eQTL_5', 'eQTL_2'], 2: ['eQTL_4', 'eQTL_9', 'eQTL_1'], 3: ['eQTL_7', 'eQTL_6', 'eQTL_8']}
    ##  is first treated to alphatize each element within the 3 clusters, giving:
    ##      {1: ['eQTL_2', 'eQTL_3', 'eQTL_5'], 2: ['eQTL_1', 'eQTL_4', 'eQTL_9'], 3: ['eQTL_6', 'eQTL_7', 'eQTL_8']}
    ##  and then the clusters are alphabetically assigned based on the 1st element within each cluster, giving:
    ##      {1: ['eQTL_1', 'eQTL_4', 'eQTL_9'], 2: ['eQTL_2', 'eQTL_3', 'eQTL_5'], 3: ['eQTL_6', 'eQTL_7', 'eQTL_8']}
    cluster_index___2___list_of_eQTLs_in_clust_alph = {}
    cluster_index___2___cluster_name = {} ## again, the "cluster_name" is really just the first alphabetical element within a cluster
    list_of_cluster_names = list()
    cluster_name___2___list_of_eQTLs_in_clust_alph = {}
    clst = 1
    while clst <= K:
        cluster_index___2___list_of_eQTLs_in_clust_alph[clst] = sorted(cluster_index___2___list_of_eQTLs_in_clust[clst])
        cl_name = cluster_index___2___list_of_eQTLs_in_clust_alph[clst][0]
        cluster_index___2___cluster_name[clst] = cl_name
        cluster_name___2___list_of_eQTLs_in_clust_alph[cl_name] = cluster_index___2___list_of_eQTLs_in_clust_alph[clst]
        list_of_cluster_names.append(cluster_index___2___cluster_name[clst])
        clst += 1
    list_of_cluster_names.sort()
    cluster_index___2___list_of_eQTLs_in_clust_alph_final = {}
    clst = 1
    r = 0
    while clst <= K:
        cluster_index___2___list_of_eQTLs_in_clust_alph_final[clst] = cluster_name___2___list_of_eQTLs_in_clust_alph[list_of_cluster_names[r]]
        clst += 1
        r += 1
    list_of_partition_dictionaries.append(cluster_index___2___list_of_eQTLs_in_clust_alph_final)
    #iteration___2___partition[iteration] = cluster_index___2___list_of_eQTLs_in_clust
    iteration += 1

##  Determine the "consensus partition" (most_common_partition) based 
##  on the num_iterations times we performed k-means clustering
most_common_partition = mode(list_of_partition_dictionaries)
fract = float(most_common_partition[1][0]) / float(num_iterations)
print "Most common partition from running the k-means algorithm:"
k = 1
while k <= K:
    print "    Cluster " + str(k) + ":  " + str(most_common_partition[0][0][k])
    k += 1


print "\n\nFraction of times that this particular partition was obtained = " + str(fract)

