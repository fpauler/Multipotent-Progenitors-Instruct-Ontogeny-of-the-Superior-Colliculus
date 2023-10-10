# bash

# used via conda environment
# source /mnt/nfs/clustersw/Debian/bullseye/anaconda3/2022.05/activate_anaconda3_2022.05.txt
# conda activate velocyto

# Note that velocyto has to be called to produce the loom file
# details on calling velocyto in the readme file

# script to analyse data with scVelo
# from: https://github.com/basilkhuder/Seurat-to-RNA-Velocity

# python 3.10.6
import anndata
import scvelo as scv # v0.2.4
import pandas as pd
import numpy as np
import matplotlib as plt

# read meta data
# files are produced during analysis with R
sample_obs = pd.read_csv("./files/cellID_obs.csv")
umap_cord = pd.read_csv("./files/cell_embeddings.csv")
meta_data = pd.read_csv("./files/meta_data.csv")

# define loom files:

# read the loom files
# give absolute path here
ann_data = anndata.read_loom("possorted_genome_bam_IWTQ9.loom")

# change rownames of pandas data frame to match the Seurat labels
# and do the filtering
new_rownames = list()
for i in range(0, len(ann_data.obs.index)):
   tmp = ann_data.obs.index[i].split(":")[1]
   tmp = tmp.replace('x', '-1')
   new_rownames.append(tmp)

#use these new rownames
ann_data.obs.index = new_rownames
#filter dataframe
ann_data = ann_data[np.isin(ann_data.obs.index, sample_obs['x'])]

#this block is necessary to match the UMAP coordinates to the Velocity dataframe
all_samples_index = pd.DataFrame(ann_data.obs.index)
all_samples_index = all_samples_index.rename(columns = {0:'Cell ID'})

umap_cord = umap_cord.rename(columns = {'Unnamed: 0':'Cell ID'})
umap_ordered = all_samples_index.merge(umap_cord,on="Cell ID")

umap_ordered = umap_ordered.iloc[:,1:]

ann_data.obsm['X_umap'] = umap_ordered.values

# do the same with coloring - to match to the UMAp colors of Seurat
meta_data = meta_data.rename(columns = {'cellID':'Cell ID'})
meta_data = meta_data[['Cell ID', 'seurat_clusters']]

meta_ordered = all_samples_index.merge(meta_data,on="Cell ID")

meta_ordered = meta_ordered.iloc[:,1:]

ann_data.obs['cluster'] = meta_ordered.values

#do the Velocity analysis
scv.pp.filter_and_normalize( ann_data )
scv.pp.moments( ann_data )
scv.tl.velocity( ann_data, mode = "stochastic")
scv.tl.velocity_graph( ann_data )

#do the plotting
scv.pl.velocity_embedding_stream(ann_data, basis='umap', save = './files/scVelo_stochastic.png', dpi=600, title = None, legend_loc = 'right', legend_fontoutline = 0, legend_fontsize = 0)
scv.pl.velocity_embedding_stream(ann_data, basis='umap', save = './files/scVelo_stochastic_noDots.png', dpi=600, title = None, legend_loc = 'right', legend_fontoutline = 0, legend_fontsize = 0, size=0)
