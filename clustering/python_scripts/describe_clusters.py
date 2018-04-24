#! /usr/bin/python3
###############################################################
###		EXTRACT INFORMATION FOR CLUSTERS	    ###
###				        		    ###
### Autor: Diego Granados				    ###
### Date:23/04/18					    ###	
### Requieres several libraries listed below                ###
###############################################################

# Reading the libraries

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tic
import seaborn as sns
import hdbscan
import pandas as pd

# Function for reading files #

def read_data(filename):
        with open(filename, 'r') as f:
                x = []
                y = []
                d=f.readlines()
                for line in d:
                        p = line.split()
                        x.append(float(p[0]))
                        y.append(float(p[1]))
        data = np.column_stack((x, y))
        return data

# Actually reading the files

L117K = read_data("../data/at-angl_L117K_open_all-dis.dat")
L117Q = read_data('../data/at-angl_L117Q_open_all-dis.dat')
L117R = read_data('../data/at-angl_L117R_open_all-dis.dat')
LAO = read_data('../data/at-angl_2LAO_open_all-dis.dat')

# Convert data to radians for using harversine distance

factor=3.141516/180
LAO=factor*LAO; L117K=factor*L117K; L117R=factor*L117R; L117Q=factor*L117Q

# Find the clusters 

clusters_LAO = hdbscan.HDBSCAN( min_cluster_size=50, metric='haversine').fit(LAO)
clusters_L117K = hdbscan.HDBSCAN( min_cluster_size=50, metric='haversine').fit(L117K)
clusters_L117R = hdbscan.HDBSCAN( min_cluster_size=50, metric='haversine').fit(L117R)
clusters_L117Q = hdbscan.HDBSCAN( min_cluster_size=50, metric='haversine').fit(L117Q)

# List that contains the fitting of all the systems
L = [ clusters_LAO, clusters_L117K, clusters_L117R, clusters_L117Q ]

# This list is for naming the files
models=["LAO","L117K","L117R","L117Q"]

# Build summary files for every system

for i in range(len(L)) :
	cluster,counts = np.unique(L[i].labels_, return_counts=True)
	# Percentage of the simulation time occupied by every cluster
	percentage=[(x/5000)*100 for x in counts]
	summary=pd.DataFrame(list(zip(cluster,counts,percentage)))
	summary.to_csv("../results/summary_{0}.csv".format(models[i]),index=False,header=False)

# Obtain the frames that compose every cluster for all the systems. Then save them in a text file

index=list(range(1,(len(L[0].labels_)+1)))

for i in (range(len(L))):
	bigger_cluster=L[i].labels_.max()+1
	indexed_labels=list(zip(index,list(L[i].labels_)))
	for f in range(-1,bigger_cluster):
		tmp=[x[:][0] for x in indexed_labels if x[:][1]==f]
		frames=pd.DataFrame(tmp)
		name="../results/frames_in_cluster_{0}_for_{1}.dat".format(f,models[i])
		frames.to_csv(name,index=False,header=False)

# Distribution of frames by simulation for every system. Then save a text file with a table
for i in range(len(L)):
	bigger_cluster=L[i].labels_.max()+1
        indexed_labels=list(zip(index,list(L[i].labels_)))
	dist_mat=np.zeros((bigger_cluster+1,10))
	for f in range(-1,bigger_cluster):
		dist_mat[f+1,0]=len([x[:][1] for x in indexed_labels if x[:][0] < 500 and x[:][1]==f])
		for g in range(1,10):
			dist_mat[f+1,g]=len([x[:][1] for x in indexed_labels if x[:][0]>=g*500 and x[:][0] < (g*500)+500 and x[:][1]==f])
	summary_by_replica=pd.DataFrame(dist_mat)
	name="../results/summary_by_replica_in_{0}.dat".format(models[i])
	summary_by_replica.to_csv(name,sep=' ', header=False, index=False)
