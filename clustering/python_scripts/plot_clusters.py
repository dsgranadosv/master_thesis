#! /usr/bin/python3

# Loading libraries #

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tic
import seaborn as sns
import hdbscan

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

# Reading the files #

L117K = read_data("../data/at-angl_L117K_open_all-dis.dat")
L117Q = read_data('../data/at-angl_L117Q_open_all-dis.dat')
L117R = read_data('../data/at-angl_L117R_open_all-dis.dat')
LAO = read_data('../data/at-angl_2LAO_open_all-dis.dat')

# Function for plotting clusters, 
# an adaptation of the function defined
# in the HDBSCAN documentatation: http://hdbscan.readthedocs.io/en/latest/comparing_clustering_algorithms.html

def plot_clusters(dataset,args,kwds):
	labels = hdbscan.HDBSCAN(*args, **kwds).fit_predict(dataset)
	palette = sns.color_palette('deep', np.unique(labels).max() + 1)
	colors = [palette[x] if x >= 0 else (0, 0, 0) for x in labels]
	plt.scatter(dataset.T[0], dataset.T[1],  c=colors)
	plt.xlim(1.3,2.4)
	plt.ylim(0.3,2.2)

# Convert data to radians for using harversine distance

factor=3.141516/180
LAO=factor*LAO; L117K=factor*L117K; L117R=factor*L117R; L117Q=factor*L117Q

###### Plotting the clusters #######
sns.set()
fig,[ax1,ax2,ax3,ax4]=plt.subplots(1,4,figsize=(18, 6),dpi=300)
# LAO
plt.subplot(1,4,1)
plot_clusters(LAO,(),{'min_cluster_size':50,'metric':'haversine'})
plt.tick_params(axis='both',labelsize=15)
plt.ylabel('Ángulo de torsión (radianes)', fontsize = 17)
plt.title("LAO",fontsize=17)
# L117K
plt.subplot(1,4,2)
clusters_L117K=plot_clusters(L117K,(),{'min_cluster_size':50,'metric':'haversine'})
plt.tick_params(axis='y',labelleft=False)
plt.tick_params(axis='both',labelsize=15)
plt.title("L117K",fontsize=17)
# L117R
plt.subplot(1,4,3)
clusters_L117R=plot_clusters(L117R,(),{'min_cluster_size':50,'metric':'haversine'})
plt.tick_params(axis='y',labelleft=False)
plt.tick_params(axis='both',labelsize=15)
plt.title("L117R",fontsize=17)
# L117Q
plt.subplot(1,4,4)
clusters_L117Q=plot_clusters(L117Q,(),{'min_cluster_size':50,'metric':'haversine'})
plt.tick_params(axis='y',labelleft=False)
plt.tick_params(axis='both',labelsize=15)
plt.title("L117Q",fontsize=17)
# Adjusting white space and additional information
plt.subplots_adjust(wspace=0.1)
fig.suptitle('Agrupamiento de conformaciones basado en densidad', fontsize=20)
fig.text(0.5, 0.018, 'Ángulo de apertura (radianes)', ha='center', fontsize=17)
# Changing size and saving the figure
fig.savefig('../results/clustering.png', dpi = 300, bbox_inches='tight')
plt.close()
