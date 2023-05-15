import episcanpy.api as epi
import anndata as ad
from sklearn.cluster import MeanShift
import numpy as np
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics import *
from sklearn.metrics import normalized_mutual_info_score as NMI
import time
import tracemalloc
file_names = ["adrenal","cerebellum","cerebrum","eye","heart","intestine", "kidney", "liver","lung", "muscle","pancreas","placenta", "spleen", "stomach", "thymus"]

LSI_PATH = "/work/lect0077/ATAC/h5ad/LSI/"
PCA_PATH = "/hpcwork/lect0077/ATAC/h5ad/PCA/"

def mean_shift(name, path, data_type):
    adata = ad.read(path + name + ".h5ad")
    start_time = time.time()
    tracemalloc.start()
    #making bandwidht more than 10 or default makes only 1 cluster. This gives better result overall 
    ms = MeanShift(bandwidth=10).fit(adata.obsm[data_type])
    labels = ms.labels_

    print(len(np.unique(labels)))
    X = adata.obsm[data_type]
    
    sil = silhouette_score(X, labels)
    ari = adjusted_rand_score(adata.obs['cell_type'], labels)
    homogenity = homogeneity_score(adata.obs['cell_type'], labels)
    nmi =  NMI(adata.obs['cell_type'], labels)
    curr_mem_size, peak_mem_size = tracemalloc.get_traced_memory()
    end_time = time.time()
    total_time =  end_time-start_time
    memory_used = str((int)(peak_mem_size/(10**6))) + 'Mb'
    with open("MeanShiftdataPca.txt", 'a+') as f:
        f.write(path + name + ".h5ad")
        f.write("\n" + "ari   " + str(ari) )
        f.write("\n" + "silhouette     :   "  + str(sil))
        f.write("\n" + "nmi  :   " + str(nmi))
        f.write('\n'+ "homogenity   :  " + str(homogenity))
        f.write('\n' + "total_time   :   " + str(total_time))
                f.write('\n' + "total_time   :   " + str(total_time))
        f.write('\n' + "memory_used   :   " + memory_used)
        f.write('\n' + "Number of clusters  :     " + str(len(np.unique(labels))))
        f.write('\n')

for name in file_names:
    agglomerative_cluster(name , LSI_PATH,'X_pca')

for name in file_names:
    agglomerative_cluster(name , PCA_PATH, 'X_lsi')
