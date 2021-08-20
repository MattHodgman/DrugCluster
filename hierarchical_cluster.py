from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys

print('Loading distance matrix...')
distance_matrix = pd.read_csv(sys.argv[1], delimiter='\t', index_col='Drug') # get distance matrix
print('Converting to condensed matrix...')
X = squareform(distance_matrix) # convert to condensed form
print('Performing hierarchical clustering...')
Z = linkage(X, 'ward') # cluster
clusters = fcluster(Z, 100,'maxclust') # get cluster labels


# data = np.empty((0, 2)) # init empty array to store mean silhouette scores and numebr of clusters
# # loop through values of k
# for k in range(10,210,10):
#     print(f'Extracting clusters (k={k})...')
#     clusters = fcluster(Z, k,'maxclust') # get cluster labels
#     print('\tCalculating mean silhouette score...')
#     mean_s = silhouette_score(distance_matrix, clusters, metric="precomputed") # calculate mean silhouette score
#     data = np.append(data, np.array([[k, mean_s]]), axis=0)

# # plot
# print('Plotting...')
# plt.scatter(x=data[:,0], y=data[:,1])
# plt.title('ChEBML Drug Hierarchical Clustering Solutions')
# plt.xlabel('k')
# plt.ylabel('mean silhouette score')
# plt.show()

# OLD
results = distance_matrix.index.to_frame().drop(columns='Drug') # create df with drug names as index
results['Cluster'] = clusters # match cluster labels to Drug names
results.to_csv('all_drugs_hierarchical_clusters.tsv', sep='\t')