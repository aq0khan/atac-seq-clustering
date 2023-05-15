import episcanpy.api as epi
import anndata as ad
import numpy as np
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics import *
from sklearn.metrics import normalized_mutual_info_score as NMI
import time
import tracemalloc
file_names = [""]

LSI_PATH = ""
PCA_PATH = ""

def agglomerative_cluster(name, path):
    adata = ad.read(path + name + ".h5ad")
    start_time = time.time()
    tracemalloc.start()
    #use these parameters as it worked for thesse datasets
    cluster = AgglomerativeClustering(n_clusters=10, affinity='euclidean', linkage='ward')
    adata.obs['hclust_10'] = cluster.fit_predict(adata.obsm['X_lsi']).astype(str)
    calculate_score(adata)

def calculate_score(adata)
    ari = epi.tl.ARI(adata, 'hclust_10', 'cell_type')
    homogenity = epi.tl.homogeneity(adata, 'hclust_10', 'cell_type')
    

def calculate_score(adata, label)
    if X_lsi in adata.obsm:
        sil = silhouette_score(adata.obsm['X_lsi'], adata.obs['hclust_10'])
    else
        sil = silhouette_score(adata.obsm['X_pca'], adata.obs['hclust_10'])
    
    ari = epi.tl.ARI(adata, 'hclust_10', 'cell_type')
    homogenity = epi.tl.homogeneity(adata, 'hclust_10', 'cell_type')
    nmi = NMI(adata.obs['cell_type'], adata.obs['hclust_10'])
    curr_mem_size, peak_mem_size = tracemalloc.get_traced_memory()
    end_time = time.time()
    total_time =  end_time-start_time
    memory_used = str((int)(peak_mem_size/(10**6))) + 'Mb'
    with open("aggl.txt", 'a+') as f:
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

