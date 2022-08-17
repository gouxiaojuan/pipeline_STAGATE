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

adata = sc.AnnData(counts.T) #把counts变为带注释的数据矩阵。
adata.var_names_make_unique() #var_names_make_unique
coor_df = coor_df.loc[adata.obs_names, ['y', 'x']] #将adata.obs_names作为行名，并按y,x作为列名
adata.obsm["spatial"] = coor_df.to_numpy() #to_numpy 方法将 DataFrame 转换为 NumPy 数组
sc.pp.calculate_qc_metrics(adata, inplace=True) #计算质量控制指标。
sc.pp.filter_genes(adata, min_cells=50)
print('After flitering: ', adata.shape) 
#Normalization
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000) #选取高变基因
sc.pp.normalize_total(adata, target_sum=1e4) #标准化每个细胞的计数。
sc.pp.log1p(adata)#对 数据矩阵进行对数。
STAGATE.Cal_Spatial_Net(adata, rad_cutoff=50) #构建空间邻居网络。
adata = STAGATE.train_STAGATE(adata, alpha=0) #训练图注意力自动编码器。
sc.pp.neighbors(adata, use_rep='STAGATE') #计算观测值的邻域图
sc.tl.umap(adata) #使用 UMAP [McInnes18]_ 嵌入邻域图。
sc.tl.louvain(adata, resolution=0.8) #将细胞聚集成不同的空间域

sc.pl.embedding(adata, basis="spatial", color="louvain",s=6, show=False, title='STAGATE')
plt.savefig("./spatial_domain.pdf",dpi=300,bbox_inches = 'tight')

