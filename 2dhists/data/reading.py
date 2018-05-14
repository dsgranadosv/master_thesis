# Loading libraries #
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import hdbscan

# Definition of function for reading files #
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

# Reading of the files #
L117K = read_data('at-angl_L117K_open_all-dis.dat')
L117Q = read_data('at-angl_L117Q_open_all-dis.dat')
L117R = read_data('at-angl_L117R_open_all-dis.dat')
LAO = read_data('at-angl_2LAO_open_all-dis.dat')
# Function for plotting clusters, 
# an adaptation of the function defined
# in the HDBSCAN documentatation: http://hdbscan.readthedocs.io/en/latest/comparing_clustering_algorithms.html
def plot_clusters(dataset,args,kwds):
	device=hdbscan.HDBSCAN(*args, **kwds)
	labels = hdbscan.HDBSCAN(*args, **kwds).fit_predict(dataset)
	palette = sns.color_palette('deep', np.unique(labels).max() + 1)
	colors = [palette[x] if x >= 0 else (0, 0, 0) for x in labels]
	plt.scatter(dataset.T[0], dataset.T[1],  c=colors)
	#plt.xlim(80,130)
	#plt.ylim(10,130)
	return device

#Conver data to radians to use harversine distance
factor=3.141516/180
LAO=factor*LAO; L117K=factor*L117K; L117R=factor*L117R; L117Q=factor*L117Q
#Plotting the clusters and saving them for posterior ussage
fig, ax = plt.subplots(nrows=1, ncols=4)
plt.subplot(1,4,1)
clusters_LAO=plot_clusters(LAO,(),{'min_cluster_size':50,'metric':'haversine'})
plt.subplot(1,4,2)
clusters_L117K=plot_clusters(L117K,(),{'min_cluster_size':50,'metric':'haversine'})
plt.subplot(1,4,3)
clusters_L117R=plot_clusters(L117R,(),{'min_cluster_size':50,'metric':'haversine'})
plt.subplot(1,4,4)
clusters_L117Q=plot_clusters(L117Q,(),{'min_cluster_size':50,'metric':'haversine'})
plt.show()
