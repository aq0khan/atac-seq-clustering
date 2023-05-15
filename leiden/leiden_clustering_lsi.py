from inspect import trace
import sys, getopt, os, time
from tracemalloc import Snapshot, start

# parse command line options
plot = False
script_dir = os.path.dirname(os.path.abspath(__file__))
res_id = time.time()
input_file = None
input_dir = None
args = sys.argv[1:]
param_setting = '1'
try:
    opts, args = getopt.getopt(args, 'hi:c:p', [])
except getopt.GetoptError:
    print('gmm_clustering.py -i <inputfile-or-dir> -c <setting> -p')
    print('-p is an optional parameter for plotting the clustering results.')
    sys.exit()

for opt, arg in opts:
    if opt == '-h':
        print('gmm_clustering.py -i <inputfile-or-dir> -p')
        print('-p is an optional parameter for plotting the clustering results.')
        sys.exit()
    elif opt == '-p':
        plot = True
    elif opt == '-i':
        if os.path.isfile(arg):
            input_file = arg
        elif os.path.isdir(arg):
            input_dir = arg
        else:
            print('Illegal input provided.')
            sys.exit()
    elif opt == '-c':
        param_setting = arg

res_file = os.path.join(script_dir,'results','res_'+'LEIDEN_LSI'+'_'+str(res_id))

if input_file is None and input_dir is None:
    print('An input file or directory must be provided.')
    sys.exit()
else:
    print('Input file: ', input_file)
    print('Input directory: ', input_dir)
if plot:
    print('Clustering results will be plotted.')

import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import episcanpy.api as epi
from sklearn.mixture import GaussianMixture as GMM
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics import *
import tracemalloc

def gmm_routine(file, setting, res_file):
    # read h5ad data
    # e.g. 'BioInfo/ATAC/h5ad/spleen.h5ad'
    adata = ad.read(file)

    # find out the number of clusters
    cell_types = adata.obs['cell_type'].unique().tolist()
    # don't consider the cluster with "unknown" cells
    out = ['nknown' in e for e in cell_types]
    cell_types = [d for (d, remove) in zip(cell_types, out) if not remove]
    cluster_num = len(cell_types)

    sc.pp.neighbors(adata, use_rep="X_lsi", metric="cosine")

    # REPEAT 10 TIMES BECAUSE RESULTS ARE ALWAYS DIFFERENT
    ari_list = []
    nmi_list = []
    hom_list = []
    sil_list = []
    time_list = []
    mem_list = []
    for i in range(10):

        # cluster using GMM
        start_time = time.time()
        tracemalloc.start()

        X = adata.obsm['X_lsi']

        epi.tl.getNClusters(adata, n_cluster=cluster_num, method='leiden')

        end_time = time.time()
        curr_mem_size, peak_mem_size = tracemalloc.get_traced_memory()

        # calculate statistics
        # ground truth needed
        ari = adjusted_rand_score(adata.obs['cell_type'], adata.obs['leiden'])
        nmi = normalized_mutual_info_score(adata.obs['cell_type'], adata.obs['leiden'])
        homogeneity = homogeneity_score(adata.obs['cell_type'], adata.obs['leiden'])
        # no ground truth needed
        silhouette = silhouette_score(X, adata.obs['leiden'])

        ari_list.append(ari)
        nmi_list.append(nmi)
        hom_list.append(homogeneity)
        sil_list.append(silhouette)
        time_list.append(end_time-start_time)
        mem_list.append(peak_mem_size/(10**6))
        i += 1

    # save statistics
    with open(res_file, 'a+') as f:
        spacing = '\t'
        f.write('GMM clustering results for'+ spacing+ file.split('/')[-1]+ '\n')
        f.write(str(np.mean(ari_list))+ '\n')
        f.write(str(np.mean(sil_list))+ '\n')
        f.write(str(np.mean(nmi_list))+ '\n')
        f.write(str(np.mean(hom_list))+ '\n')
        # Runtime in seconds:
        f.write(str(np.mean(time_list))+ '\n')
        # Peak memory consumption in Mb
        f.write(str(np.mean(mem_list))+ 'Mb\n')
        f.write('\n')

    #epi.pl.umap(adata, color=['leiden', 'cell_type'], wspace=0.3, show=False)
    if plot:
        from matplotlib import pyplot as plt
        sc.settings.set_figure_params(dpi=80, color_map='gist_earth')
        with plt.rc_context():
            epi.pl.umap(adata, color=['leiden', 'cell_type'], wspace=0.3, show=False)
            plt.savefig('BioInfo/ATAC/test1.png', bbox_inches="tight")
    
    print(file.split('/')[-1], 'done')

if not input_file is None:
    gmm_routine(input_file, param_setting, res_file)
elif not input_dir is None:
    for file in os.listdir(input_dir):
        filename = os.fsdecode(file)
        if filename.endswith('.h5ad'):
            gmm_routine(os.path.join(input_dir,file), param_setting, res_file)

#gmm_routine('~/spp-data/h5ad/spleen.h5ad', 'PCA30', 'spleen.h5ad')
#gmm_routine('~/spp-data/h5ad/heart.h5ad', 'PCA30', 'heart.h5ad')
#gmm_routine('~/spp-data/h5ad/kidney.h5ad', 'PCA30', 'kidney.h5ad')
