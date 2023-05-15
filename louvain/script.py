import episcanpy.api as epi
import anndata as ad
import numpy as np
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics import *
import time
import tracemalloc
from sklearn.metrics import normalized_mutual_info_score as NMI

file_names = ["adrenal","cerebellum","cerebrum","eye","heart","intestine", "kidney", "liver","lung", "muscle","pancreas","placenta", "spleen", "stomach", "thymus"]

LSI_PATH = "/work/lect0077/ATAC/h5ad/LSI/"
PCA_PATH = "/hpcwork/lect0077/ATAC/h5ad/PCA/"


def louvain_cluster(name, path):
    adata = ad.read(path + name + ".h5ad")
    start_time = time.time()
    tracemalloc.start()
    epi.pp.neighbors(adata, n_neighbors=15, n_pcs=None, use_rep=None, knn=True, random_state=0, method='gauss', metric='euclidean',  metric_kwds={}, copy=False)
    epi.tl.louvain(adata)
    
def calculate_score(adata, label)
    if X_lsi in adata.obsm:
        sil = silhouette_score(adata.obsm['X_lsi'], adata.obs['louvain'])
    else :
        sil = silhouette_score(adata.obsm['X_pca'], adata.obs['louvain'])
    
    ari = epi.tl.ARI(adata, 'louvain', 'cell_type')
    homogenity = epi.tl.homogeneity(adata, 'louvain', 'cell_type')
    nmi = NMI(adata.obs['cell_type'], adata.obs['louvain'])

    curr_mem_size, peak_mem_size = tracemalloc.get_traced_memory()
    end_time = time.time()
    total_time =  end_time-start_time
    memory_used = str((int)(peak_mem_size/(10**6))) + 'Mb'
    with open("louv.txt", 'a+') as f:
        f.write(path + name + ".h5ad")
        f.write("\n" + "ari   " + str(ari) )
        f.write("\n" + "silhouette     :   "  + str(sil))
        f.write("\n" + "nmi  :   " + str(nmi))
        f.write('\n'+ "homogenity   :  " + str(homogenity))
        f.write('\n' + "total_time   :   " + str(total_time))
        f.write('\n' + "memory_used   :   " + memory_used)
        f.write('\n')


for name in file_names:
    agglomerative_cluster(name , LSI_PATH)

for name in file_names:
    agglomerative_cluster(name , PCA_PATH)
