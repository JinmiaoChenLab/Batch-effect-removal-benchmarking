
# coding: utf-8

import bbknn
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib import rcParams
import scanpy as sc
import h5py
import os
import time
from datetime import timedelta
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()


dirname = os.getcwd()
print(dirname)
save_dir = dirname
data_dir = '/batch_norm/datasets/dataset8_Mouse_brain/filtered_genes_and_cells/'


# f = h5py.File(os.path.join(data_dir,'dropviz_and_nuclei_combined_filtered_UMI.h5'), 'r')
# keys = list(f.keys())
# k2 = [x for x in keys if x not in ['gene_names', 'cell_names']]

# myData = np.array(f[k2[0]])
# print(myData.shape)
# gene_names = f['gene_names']
# cell_names = f['cell_names']
# print(gene_names.shape)
# print(cell_names.shape)

# # Decoding utf-8
# # take ~ 10 mins
# gene_names = [x.decode() for x in gene_names]  
# print(gene_names[1:3])

# sample_adata = pd.read_csv(os.path.join(data_dir,'dropviz_and_nuclei_combined_filtered_cell_info.txt'),header=0, index_col=0, sep='\t')
# print(sample_adata.values.shape)
# print(sample_adata.keys())


# adata = sc.AnnData(np.transpose(myData))  # if transposed, just use myData
# #adata = sc.AnnData(myData) 
# adata.obs_names = sample_adata.index
# adata.var_names = gene_names

# adata.obs['cell_type'] = sample_adata.loc[adata.obs_names,['cell_type']]
# adata.obs['batch'] = sample_adata.loc[adata.obs_names,['batch']]
# # adata.obs['tissue'] = sample_adata.loc[adata.obs_names,['tissue']]
# # adata.obs['organ'] = sample_adata.loc[adata.obs_names,['organ']]
# # adata.obs['batchlb'] = sample_adata.loc[adata.obs_names,['batchlb']]
# print(adata)
# adata.write(os.path.join(save_dir,'dataset8.h5ad'))

# Save output into h5ad, easy to access 
adata1 = sc.read_h5ad(os.path.join(save_dir,'dataset8.h5ad'))


sc.pp.normalize_per_cell(adata1, counts_per_cell_after=1e4)
filter_result = sc.pp.filter_genes_dispersion(adata1.X, min_mean=0.0125, max_mean=2, min_disp=0.5)
sc.pl.filter_genes_dispersion(filter_result)
print([sum([i[0] for i in filter_result]),len(filter_result)])


adata1 = adata1[:, filter_result.gene_subset]
print(adata1)


sc.pp.log1p(adata1)
#sc.pp.scale(adata1, max_value=10)
sc.tl.pca(adata1, svd_solver='arpack')
adata1.obsm['X_pca'] *= -1  # multiply by -1 to match Seurat
sc.pl.pca_variance_ratio(adata1, log=True)


num_pcs = 20
sc.pp.neighbors(adata1,n_pcs=num_pcs, n_neighbors=20)
# sc.tl.umap(adata1)
# color_group = ['cell_type','batch']
# sc.pl.umap(adata1, color=color_group,save='raw.png')


# # BBKNN first input option
# # input: batch vector in ann data
# t1 = time.time()
# adata_bbknn = bbknn.bbknn(adata1, copy=True, neighbors_within_batch=5, trim=0, n_pcs=num_pcs, batch_key='batch') #approx=False,
# t2 = time.time()
# print('Took '+str(timedelta(seconds=t2-t1)))
# sc.tl.umap(adata_bbknn)
# color_group = ['cell_type','batch']
# sc.pl.umap(adata_bbknn, color=color_group,save='bbknn.png')
# adata_bbknn.write_h5ad(os.path.join(save_dir,'bbknn.h5ad'))
# save_runtime(t1,t2,'bbknn')

# # BBKNN uses adata.obs['batch'] by default to split the batches
# # BBKNN also offers a trim parameter, which filters each cell's connections 
# # to a user-specified number of highest connectivities. 
# # This can help increase cell type separation at a potential slight cost to batch mixing
# t3 = time.time()
# adata_bbknn_trim = bbknn.bbknn(adata1,copy=True,neighbors_within_batch=5,trim=50,n_pcs=num_pcs,batch_key='batch') #approx=False,
# t4 = time.time()
# print('Took '+str(timedelta(seconds=t4-t3)))
# sc.tl.umap(adata_bbknn_trim)
# color_group = ['cell_type','batch']
# sc.pl.umap(adata_bbknn_trim, color=color_group,save='bbknn_trim.png')
# adata_bbknn_trim.write_h5ad(os.path.join(save_dir,'bbknn_trim.h5ad'))
# save_runtime(t3,t4,'bbknn_trim')

# Test different options and choose the best 
# By default, BBKNN uses annoy to compute approximate neighbours. 
# This has the potential to improve batch mixing, particularly in larger data, as connections will be present in the graph between more representatives of different batches. BBKNN does offer the option to compute exact neighbours, and does so by default via faiss. 
# The default metric changes from angular to Euclidean.
t5 = time.time()
adata_bbknn_faiss = bbknn.bbknn(adata1,copy=True,neighbors_within_batch=5,approx=False,n_pcs=num_pcs,batch_key='batch') #approx=False,
t6 = time.time()
print('Took '+str(timedelta(seconds=t6-t5)))
sc.tl.umap(adata_bbknn_faiss)
#sc.tl.louvain(adata_bbknn_faiss)
color_group = ["cell_type","batch"]
sc.pl.umap(adata_bbknn_faiss, color=color_group,save='bbknn_faiss.png')
adata_bbknn_faiss.write_h5ad(os.path.join(save_dir,'bbknn_faiss.h5ad'))
save_runtime(t5,t6,'bbknn_faiss')

# Test different options and choose the best 
# ckdtree
# However, faiss is quite difficult to install, and temperamental in execution 
# (e.g. it doesn't work via reticulate on my machine). Setting use_faiss to False, 
# or just failing to have faiss installed, makes BBKNN use scipy.spatial.cKDTree at a performance loss.
t7 = time.time()
adata_bbknn_ckdtree = bbknn.bbknn(adata1, copy=True, neighbors_within_batch=5, approx=False, n_pcs=30,use_faiss=False,batch_key='batch') #approx=False,
t8 = time.time()
print('Took '+str(timedelta(seconds=t8-t7)))
sc.tl.umap(adata_bbknn_ckdtree)
#sc.tl.louvain(adata_bbknn_ckdtree)
color_group = ["cell_type","batch"]
sc.pl.umap(adata_bbknn_ckdtree, color=color_group,save='bbknn_ckdtree.png')
adata_bbknn_ckdtree.write_h5ad(os.path.join(save_dir,'bbknn_ckdtree.h5ad'))
save_runtime(t7,t8,'bbknn_ckdtree')

