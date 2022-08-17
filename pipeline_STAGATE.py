import argparse
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import os
import sys
import STAGATE
import sklearn

parser=argparse.ArgumentParser()
parser.add_argument('-cd',"--counts_df",required=True)
parser.add_argument("-ld","--location_df",required=True)
args =  parser.parse_args()

counts = pd.read_csv(args.counts_df,index_col = 0)
coor_df = pd.read_csv(args.location_df,index_col = 0)
print("The file has been successfully read in")

adata = sc.AnnData(counts.T) 
adata.var_names_make_unique() 
coor_df = coor_df.loc[adata.obs_names, ['y', 'x']] 
adata.obsm["spatial"] = coor_df.to_numpy() 
sc.pp.calculate_qc_metrics(adata, inplace=True) 
sc.pp.filter_genes(adata, min_cells=50)
print('After flitering: ', adata.shape) 
#Normalization
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000) 
sc.pp.normalize_total(adata, target_sum=1e4) 
sc.pp.log1p(adata)
STAGATE.Cal_Spatial_Net(adata, rad_cutoff=50) 
adata = STAGATE.train_STAGATE(adata, alpha=0) 
sc.pp.neighbors(adata, use_rep='STAGATE') 
sc.tl.umap(adata)  
sc.tl.louvain(adata, resolution=0.8) 

sc.pl.embedding(adata, basis="spatial", color="louvain",s=6, show=False, title='STAGATE')
plt.savefig("./spatial_domain.pdf",dpi=300,bbox_inches = 'tight')

