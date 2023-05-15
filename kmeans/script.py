from pathlib import Path
import sys

from sklearn.metrics import normalized_mutual_info_score

root_path = Path(__file__).parents[2]
sys.path.append(str(root_path))

import common
import episcanpy as epi

n_clusters=10
def kmeans_euclidian(data):
	epi.tl.kmeans(data, num_clusters = 4)
	return (epi.tl.ARI(data, "kmeans", "cell_type"),
		epi.tl.silhouette(data, "cell_type"),
		0, # normalized_mutual_info_score(data.obs['cell_type'], data.obs['kmeans'].labels_), #errors because labels_ doesn't exist
		epi.tl.homogeneity(data, "kmeans", "cell_type"),
	)

from sklearn.cluster import KMeans
def kmeans_lsi(adata):
	# compute clustering
	X = adata.obsm['X_lsi']
	clustering = KMeans(n_clusters=n_clusters, random_state=0).fit(X)
	return clustering

def kmeans_pca(adata):
	X = adata.obsm['X_pca']
	clustering = KMeans(n_clusters=n_clusters, random_state=0).fit(X)
	return clustering

common.run_all(kmeans_lsi, directory='', pca_prep=False)
common.run_all(kmeans_pca, directory="", pca_prep=True)
