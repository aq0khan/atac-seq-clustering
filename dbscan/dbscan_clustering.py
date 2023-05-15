import sys, os
from pathlib import Path

root_path = Path(__file__).parents[2]
sys.path.append(str(root_path))

import common
from sklearn.cluster import DBSCAN
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics import *

def scores(adata, clustering):
	''' compute scores '''
	# ground truth needed
	ari = adjusted_rand_score(adata.obs['cell_type'], clustering.labels_)
	nmi = normalized_mutual_info_score(adata.obs['cell_type'], clustering.labels_)
	homogeneity = homogeneity_score(adata.obs['cell_type'], clustering.labels_)
	# no ground truth needed
	X = adata.obsm['X_lsi']
	silhouette = silhouette_score(X, clustering.labels_)

	return (ari, silhouette, nmi, homogeneity)

def dbscan_euclidian(adata):
	# compute clustering
	X = adata.obsm['X_lsi']
	clustering = DBSCAN(eps=10, min_samples=2, metric = 'euclidean', n_jobs=-1).fit(X)
	return scores(adata, clustering)

def dbscan_manhattan(adata):
	# compute clustering
	X = adata.obsm['X_lsi']
	clustering = DBSCAN(eps=0.5, min_samples=5, metric = 'manhattan', n_jobs=-1).fit(X)
	return scores(adata, clustering)

lsi_path = ''
file_name = ''
common.run(Path(lsi_path) / file_name, dbscan_euclidian)
