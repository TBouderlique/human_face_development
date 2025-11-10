#!/usr/bin/env python
# coding: utf-8

# In[1]:


import cellrank
import scvelo as scv
import os
os.chdir("/home/yakov/face_data")
path = os.getcwd()
print(path)
import platform
print(platform.python_version())


# In[2]:


import seaborn as sns
import pandas as pd
import random
import numpy as np
import loompy as lp
import scanpy as scp
import scvelo as scv
import matplotlib.pyplot as plt
import scanpy.external as sce
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = ['sans-serif']
plt.rcParams['font.sans-serif'] = ['Arial']
random.seed(101)
get_ipython().run_line_magic('matplotlib', 'inline')


# In[3]:


plt.rcParams['figure.figsize'] = [12, 8]
plt.rcParams['figure.dpi'] = 100 # 200 e.g. is really fine, but slower
scp.set_figure_params(scanpy=True, fontsize=14 )


# In[ ]:


# General Embedding
adata = scp.read_h5ad('DATA/all_samples-no_cell_cycle_260325.h5ad')
#adata_old = scp.read_h5ad('DATA/all_old-no_cell_cycle.h5ad')
#adata_new = scp.read_h5ad('DATA/all_new-no_cell_cycle.h5ad')
#adata = adata_old.concatenate(adata_new, join='outer')
#scp.pp.scale(adata, max_value=10)
#adata.raw=adata
scp.tl.pca(adata,n_comps = 30, svd_solver='arpack',use_highly_variable=True)
if adata.raw is not None:
    adata = adata.raw.to_adata()
    del adata.raw
#sce.pp.harmony_integrate(adata, 'Sample_type')
#scp.pp.neighbors(adata, n_neighbors = 15, n_pcs = 25, use_rep = 'X_pca_harmony' )
scp.pp.neighbors(adata, n_neighbors = 30, n_pcs = 25, use_rep = 'X_pca' )
scp.tl.leiden(adata, resolution = 0.4)
scp.tl.umap(adata,n_components = 2)
scp.pl.umap(adata, color=['leiden'], save = '-all.png')
#for obsp_key in list(adata.obsp.keys()):
#    del adata.obsp[obsp_key]
#for varm_key in list(adata.varm.keys()):
#    del adata.varm[varm_key]
#for uns_key in list(adata.uns.keys()):
#    del adata.uns[uns_key]
#if adata.raw is not None:
#    adata = adata.raw.to_adata()
#    del adata.raw
adata.write('DATA/all_samples.h5ad')


# In[ ]:


# Mesenchymal Embedding
adata = scp.read_h5ad('DATA/all_samples.h5ad')
scp.pl.umap(adata, color=['leiden'], legend_loc="on data", legend_fontsize = 10)
adata_mes = adata[adata.obs['leiden'].isin(['0','1','2','3','4','5','6','7','8','9','10','11','13','19']),:].copy()
adata_mes.obs.leiden=adata_mes.obs.leiden.cat.rename_categories(range(len(adata_mes.obs.leiden.cat.categories))).astype(str)
scp.pl.umap(adata_mes, color=['leiden'], legend_loc="on data")

cell_cycle_genes = pd.read_csv('../cell_cycle_genes/all_human_genes.csv')
s_genes = cell_cycle_genes.all_genes[:43]
g2m_genes = cell_cycle_genes.all_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes.all_genes if x in adata_mes.var_names]
adata_mes.raw = adata_mes.copy()
scp.pp.scale(adata_mes)
scp.tl.score_genes_cell_cycle(adata_mes, s_genes=s_genes, g2m_genes=g2m_genes)
scp.pp.regress_out(adata_mes, ['S_score', 'G2M_score'])
scp.tl.pca(adata_mes,n_comps = 30, svd_solver='arpack',use_highly_variable=True)
if adata_mes.raw is not None:
    adata_mes = adata_mes.raw.to_adata()
    del adata_mes.raw
sce.pp.harmony_integrate(adata_mes, 'Sample')
scp.pp.neighbors(adata_mes, n_neighbors = 30, n_pcs = 25, use_rep = 'X_pca_harmony' )
#scp.pp.neighbors(adata_mes, n_neighbors = 30, n_pcs = 25, use_rep = 'X_pca' )
scp.tl.leiden(adata_mes, resolution = 1.5)
scp.tl.umap(adata_mes,n_components = 2)
scp.pl.umap(adata_mes, color=['leiden'], save = 'all-mes-harmony-sample.pdf')
adata_mes.write('DATA/all-mes-harmony-sample.h5ad')


# In[67]:


# DotPlot
adata = scp.read_h5ad('DATA/all_samples.h5ad')
adata.obs['annotation'] = "Mesenchyme"
for cell in adata.obs.loc[adata.obs.leiden == '0'].index:
    adata.obs.at[cell, 'annotation'] = "Fibroblast"
for cell in adata.obs.loc[adata.obs.leiden == '1'].index:
    adata.obs.at[cell, 'annotation'] = "Osteogenic"
for cell in adata.obs.loc[adata.obs.leiden == '2'].index:
    adata.obs.at[cell, 'annotation'] = "Mesenchyme"
for cell in adata.obs.loc[adata.obs.leiden == '3'].index:
    adata.obs.at[cell, 'annotation'] = "Dermal"    
for cell in adata.obs.loc[adata.obs.leiden == '4'].index:
    adata.obs.at[cell, 'annotation'] = "Progenitor"
for cell in adata.obs.loc[adata.obs.leiden == '5'].index:
    adata.obs.at[cell, 'annotation'] = "Stroma"
for cell in adata.obs.loc[adata.obs.leiden == '6'].index:
    adata.obs.at[cell, 'annotation'] = "Chondrogenic"
for cell in adata.obs.loc[adata.obs.leiden == '7'].index:
    adata.obs.at[cell, 'annotation'] = "Myoblast"
for cell in adata.obs.loc[adata.obs.leiden == '8'].index:
    adata.obs.at[cell, 'annotation'] = "Meningeal"
for cell in adata.obs.loc[adata.obs.leiden == '9'].index:
    adata.obs.at[cell, 'annotation'] = "Mesenchyme"
for cell in adata.obs.loc[adata.obs.leiden == '10'].index:
    adata.obs.at[cell, 'annotation'] = "Subdermal"
for cell in adata.obs.loc[adata.obs.leiden == '11'].index:
    adata.obs.at[cell, 'annotation'] = "Mesenchyme"
for cell in adata.obs.loc[adata.obs.leiden == '12'].index:
    adata.obs.at[cell, 'annotation'] = "Neuronal"
for cell in adata.obs.loc[adata.obs.leiden == '13'].index:
    adata.obs.at[cell, 'annotation'] = "Perivascular"
for cell in adata.obs.loc[adata.obs.leiden == '14'].index:
    adata.obs.at[cell, 'annotation'] = "Glia"
for cell in adata.obs.loc[adata.obs.leiden == '15'].index:
    adata.obs.at[cell, 'annotation'] = "Endothelia"
for cell in adata.obs.loc[adata.obs.leiden == '16'].index:
    adata.obs.at[cell, 'annotation'] = "Neuronal"    
for cell in adata.obs.loc[adata.obs.leiden == '17'].index:
    adata.obs.at[cell, 'annotation'] = "Immune"
for cell in adata.obs.loc[adata.obs.leiden == '18'].index:
    adata.obs.at[cell, 'annotation'] = "Blood"
for cell in adata.obs.loc[adata.obs.leiden == '19'].index:
    adata.obs.at[cell, 'annotation'] = "Myoblast"
for cell in adata.obs.loc[adata.obs.leiden == '20'].index:
    adata.obs.at[cell, 'annotation'] = "Neuronal"
for cell in adata.obs.loc[adata.obs.leiden == '21'].index:
    adata.obs.at[cell, 'annotation'] = "Mesenchyme"
for cell in adata.obs.loc[adata.obs.leiden == '21'].index:
    adata.obs.at[cell, 'annotation'] = "Epithelia"
for cell in adata.obs.loc[adata.obs.leiden == '22'].index:
    adata.obs.at[cell, 'annotation'] = "Neuronal" 
for cell in adata.obs.loc[adata.obs.leiden == '23'].index:
    adata.obs.at[cell, 'annotation'] = "Immune"     

#scp.pl.umap(adata, color=['annotation'],use_raw=False) 

#chondro is always orange '#ff7f0e', 
#dermal is always red '#d62728', 
#osteo always blue '#1f77b4', 
#teno/connective always brown '#8c564b', 
#progenitors purple '#aa40fc',
#dental and fox that don't always show up

#adata.uns['annotation_colors'][0] = '#ff7f0e'
#adata.uns['annotation_colors'][1] = '#8c564b'
#adata.uns['annotation_colors'][2] = '#279e68'
#adata.uns['annotation_colors'][3] = '#1f77b4'
#adata.uns['annotation_colors'][4] = '#d62728'
#adata.uns['annotation_colors'][5] = '#aa40fc'


#sc.pl.umap(adata, color=['leiden'],use_raw=False)

scp.set_figure_params(dpi_save=300, dpi=300)
scp.pl.umap(adata, color=['annotation'],use_raw=False, size=5, save = '_full.pdf')

markers = ["ALAS2","HBB","NFKBIA","PECAM1","EPCAM","TUBB3","CKB","SOX10","S100B","MYF5","TNNT1","COL2A1","COL9A1","SOX9","PTN","RUNX2","ACTA2","SULT1E1","TWIST2","VIM","PDGFRA","OGN","RARB","GJA1","WNT5A","CXCL12","DCN","POSTN","COL5A1","H19","HMGA2","CACNA1C","ROBO1","RERE","COL14A1"]
scp.pl.dotplot(adata, markers, groupby='annotation', dendrogram=False, swap_axes = False, save = 'Fig1E-DotPlot.pdf', standard_scale="var", categories_order = ["Blood","Immune","Endothelia","Epithelia","Neuronal","Glia","Myoblast","Chondrogenic","Osteogenic","Perivascular","Dermal","Subdermal","Meningeal","Stroma","Fibroblast","Progenitor","Mesenchyme"])
scp.pl.umap(adata, color=['leiden'], save = 'Fig1D-General.pdf')
adata.write('DATA/all_samples.h5ad')


# In[66]:


adata = scp.read_h5ad('DATA/all_samples.h5ad')
markers = ["ALAS2","HBB","NFKBIA","PECAM1","EPCAM","TUBB3","CKB","SOX10","S100B","MYF5","TNNT1","COL2A1","COL9A1","SOX9","PTN","RUNX2","ACTA2","SULT1E1","TWIST2","VIM","PDGFRA","OGN","RARB","GJA1","WNT5A","CXCL12","DCN","POSTN","COL5A1","H19","HMGA2","CACNA1C","ROBO1","RERE","COL14A1"]
scp.pl.dotplot(adata, markers, groupby='annotation', dendrogram=False, swap_axes = False, save = 'Fig1E-DotPlot.pdf', standard_scale="var", categories_order = ["Blood","Immune","Endothelia","Epithelia","Neuronal","Glia","Myoblast","Chondrogenic","Osteogenic","Perivascular","Dermal","Subdermal","Meningeal","Stroma","Fibroblast","Progenitor","Mesenchyme"])
adata.write('DATA/all_samples.h5ad')


# In[40]:


adata = scp.read_h5ad('DATA/all_samples.h5ad')
scp.pp.log1p(adata)
scp.tl.rank_genes_groups(adata, 'leiden')
pd.set_option("display.max_columns", 60)
pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(10)


# In[7]:


adata = scp.read_h5ad('DATA/all-mes-harmony-sample.h5ad')
scp.pl.umap(adata, color=['leiden'], save = 'Fig1D-Mesenchymal.pdf')


# In[5]:


# CytoTRACE
adata = scp.read_h5ad('DATA/all_samples.h5ad')
scp.pl.umap(adata, color=['leiden'], legend_loc="on data", legend_fontsize = 10)
adata_mes = adata[adata.obs['leiden'].isin(['0','1','2','3','4','5','6','8','9','10','11','13']),:].copy()
adata_mes.obs.leiden=adata_mes.obs.leiden.cat.rename_categories(range(len(adata_mes.obs.leiden.cat.categories))).astype(str)
scp.pl.umap(adata_mes, color=['leiden'], legend_loc="on data")
adata = adata_mes.copy()

scp.pp.neighbors(adata)
adata.layers["spliced"] = adata.X
adata.layers["unspliced"] = adata.X
#scv.pp.moments(adata, n_neighbors = 15, n_pcs = 25 )
scv.pp.moments(adata)
ctk = cellrank.kernels.CytoTRACEKernel(adata)
ckt = ctk.compute_cytotrace().compute_transition_matrix()
scp.set_figure_params(dpi=150)
scp.pl.umap(adata,color="ct_score", save = 'Fig1H-cytoTRACE.pdf', size = 10)
adata.write('DATA/all-samples-Fig1H.h5ad')


# In[ ]:


adata


# In[4]:


# Age Origin
adata = scp.read_h5ad('DATA/all_samples.h5ad')
scp.pl.umap(adata, color=['leiden'], legend_loc="on data")
k = adata.obs.groupby(['leiden', 'Age'])['Age'].count()
k.index
for cell in adata.obs_names:
    cells_count = int(k[str(adata.obs['leiden'][cell]), '6']) + int(k[str(adata.obs['leiden'][cell]), '6.5']) + int(k[str(adata.obs['leiden'][cell]), '7']) + int(k[str(adata.obs['leiden'][cell]), '7.5']) + int(k[str(adata.obs['leiden'][cell]), '9']) + int(k[str(adata.obs['leiden'][cell]), '9.5']) + int(k[str(adata.obs['leiden'][cell]), '10']) + int(k[str(adata.obs['leiden'][cell]), '10.5']) + int(k[str(adata.obs['leiden'][cell]), '11']) + int(k[str(adata.obs['leiden'][cell]), '11.5'])
    sum_age = 6 * int(k[str(adata.obs['leiden'][cell]), '6']) + 6.5 * int(k[str(adata.obs['leiden'][cell]), '6.5']) + 7 * int(k[str(adata.obs['leiden'][cell]), '7']) + 7.5 * int(k[str(adata.obs['leiden'][cell]), '7.5']) + 9 * int(k[str(adata.obs['leiden'][cell]), '9']) + 9.5 * int(k[str(adata.obs['leiden'][cell]), '9.5']) + 10 * int(k[str(adata.obs['leiden'][cell]), '10']) + 10.5 * int(k[str(adata.obs['leiden'][cell]), '10.5']) + 11 * int(k[str(adata.obs['leiden'][cell]), '11']) + 11.5 * int(k[str(adata.obs['leiden'][cell]), '11.5'])
    adata.obs.loc[cell, 'ClusterAge' ] = sum_age/cells_count
scp.pl.umap(adata, color=['ClusterAge'],cmap='turbo', size=1, title = 'Average age of cluster (week)', save = 'Fig1G-Age.pdf')
scp.pl.umap(adata, color=['Origin'], save = 'Fig1F-Origin.pdf', title = 'Positional origin of cell')


# In[6]:


scp.pl.umap(adata, color=['ClusterAge'],cmap='turbo', size=1, title = 'Average age of cluster (week)', save = 'Fig1G-Age.pdf')


# In[7]:


import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# Define custom colormap
colors = [
    (0.6, 0.8, 1.0),   # light blue
    (0.7, 1.0, 0.7),   # light green
    (1.0, 1.0, 0.6),   # yellow
    (1.0, 0.8, 0.6),   # light orange
    (1.0, 0.6, 0.6)    # light red
]
cmap_custom = LinearSegmentedColormap.from_list("soft_heat", colors, N=256)

# Plot with Scanpy
scp.pl.umap(
    adata,
    color=['ClusterAge'],
    cmap=cmap_custom,
    size=1,
    title='Average age of cluster (week)',
    save='Fig1G-Age.pdf'
)


# In[14]:


# Mesenchymal CytoTRACE Origin
adata = scp.read_h5ad('DATA/all-mes-harmony-sample.h5ad')
scp.pp.neighbors(adata)
adata.layers["spliced"] = adata.X
adata.layers["unspliced"] = adata.X
#scv.pp.moments(adata, n_neighbors = 15, n_pcs = 25 )
scv.pp.moments(adata)
ctk = cellrank.kernels.CytoTRACEKernel(adata)
ckt = ctk.compute_cytotrace().compute_transition_matrix()
scp.set_figure_params(dpi=150)
scp.pl.umap(adata,color="ct_score", cmap="Spectral_r",  save = 'Fig1H-Mesenchymal-cytoTRACE.pdf', size = 10)
scp.pl.umap(adata, color=['Origin'], save = 'Fig1F-Mesenchymal-Origin.pdf', title = 'Positional origin of cell')
adata.write('DATA/all-mes-harmony-sample-Fig1H-Fig1F.h5ad')


# In[14]:


# Mesenchymal Age
adata = scp.read_h5ad('DATA/all-mes-harmony-sample-Fig1H-Fig1F.h5ad')
scp.pl.umap(adata, color=['leiden'], legend_loc="on data")
k = adata.obs.groupby(['leiden', 'Age'])['Age'].count()
k.index
for cell in adata.obs_names:
    cells_count = int(k[str(adata.obs['leiden'][cell]), '6']) + int(k[str(adata.obs['leiden'][cell]), '6.5']) + int(k[str(adata.obs['leiden'][cell]), '7']) + int(k[str(adata.obs['leiden'][cell]), '7.5']) + int(k[str(adata.obs['leiden'][cell]), '9']) + int(k[str(adata.obs['leiden'][cell]), '9.5']) + int(k[str(adata.obs['leiden'][cell]), '10']) + int(k[str(adata.obs['leiden'][cell]), '10.5']) + int(k[str(adata.obs['leiden'][cell]), '11']) + int(k[str(adata.obs['leiden'][cell]), '11.5'])
    sum_age = 6 * int(k[str(adata.obs['leiden'][cell]), '6']) + 6.5 * int(k[str(adata.obs['leiden'][cell]), '6.5']) + 7 * int(k[str(adata.obs['leiden'][cell]), '7']) + 7.5 * int(k[str(adata.obs['leiden'][cell]), '7.5']) + 9 * int(k[str(adata.obs['leiden'][cell]), '9']) + 9.5 * int(k[str(adata.obs['leiden'][cell]), '9.5']) + 10 * int(k[str(adata.obs['leiden'][cell]), '10']) + 10.5 * int(k[str(adata.obs['leiden'][cell]), '10.5']) + 11 * int(k[str(adata.obs['leiden'][cell]), '11']) + 11.5 * int(k[str(adata.obs['leiden'][cell]), '11.5'])
    adata.obs.loc[cell, 'ClusterAge' ] = sum_age/cells_count
scp.pl.umap(adata, color=['ClusterAge'], title = 'Average age of cluster (week)', save = 'Fig1G-Mesenchymal-Age.pdf')


# In[ ]:


adata = scp.read_h5ad('DATA/all-mes-harmony-sample.h5ad')
scp.tl.embedding_density(adata, basis='umap', groupby = 'Age', key_added = 'kdeAge')
scp.pl.embedding_density(adata, basis='umap', key = 'kdeAge', save = 'kde-mes-Age.png')
scp.tl.embedding_density(adata, basis='umap', groupby = 'Origin', key_added = 'kdeOrigin')
scp.pl.embedding_density(adata, basis='umap', key = 'kdeOrigin', save = 'kde-mes-Origin.png')
scp.tl.embedding_density(adata, basis='umap', groupby = 'Sample_type', key_added = 'kdeSample_type')
scp.pl.embedding_density(adata, basis='umap', key = 'kdeSample_type', save = 'kde-mes-Sample_type.png')


# In[ ]:


adata = scp.read_h5ad('DATA/all_samples.h5ad')
scp.tl.embedding_density(adata, basis='umap', groupby = 'Age', key_added = 'kdeAge')
scp.pl.embedding_density(adata, basis='umap', key = 'kdeAge', save = 'kde-all-Age.png')
scp.tl.embedding_density(adata, basis='umap', groupby = 'Origin', key_added = 'kdeOrigin')
scp.pl.embedding_density(adata, basis='umap', key = 'kdeOrigin', save = 'kde-all-Origin.png')
scp.tl.embedding_density(adata, basis='umap', groupby = 'Sample_type', key_added = 'kdeSample_type')
scp.pl.embedding_density(adata, basis='umap', key = 'kdeSample_type', save = 'kde-all-Sample_type.png')


# In[ ]:


scp.pl.umap(adata, color=['leiden'], legend_loc="on data", legend_fontsize = 10)
scp.pl.umap(adata, color=['Sample'], legend_loc="on data")
scp.pl.umap(adata, color=['Age'], legend_loc="on data")
scp.pl.umap(adata, color=['Origin'], legend_loc="on data")
scp.pl.umap(adata, color=['leiden', 'Age', 'Origin'])
scp.pl.umap(adata, color=['leiden'])
scp.pl.umap(adata, color=['Sample'])
scp.pl.umap(adata, color=['Age'])
scp.pl.umap(adata, color=['Origin'])
scp.pl.umap(adata, color=['Sample', 'Age', 'Origin','leiden'], legend_loc="on data")
scp.pl.umap(adata, color=['leiden'], legend_loc="on data")
scp.pl.umap(adata, color=['Sample'])
scp.pl.umap(adata, color=['leiden', 'VCAN', 'SOX6','MSX1','SCX','RGS5', 'COL3A1', 'COL5A1', 'RUNX2','COL1A1','TNMD', 'COL2A1'], legend_loc="on data")


# In[ ]:


adata = scp.read_h5ad('DATA/all_samples-harmony5.h5ad')
scp.pl.umap(adata, color=['leiden'], legend_loc="on data", legend_fontsize = 10)
adata_mes = adata[adata.obs['leiden'].isin(['0','1','2','3','4','5','6','7','8','9','10','12','14','19']),:].copy()
adata_mes.obs.leiden=adata_mes.obs.leiden.cat.rename_categories(range(len(adata_mes.obs.leiden.cat.categories))).astype(str)
scp.pl.umap(adata_mes, color=['leiden'], legend_loc="on data")
scp.tl.pca(adata_mes,n_comps = 30, svd_solver='arpack',use_highly_variable=True)
sce.pp.harmony_integrate(adata_mes, 'Sample_type')
scp.pp.neighbors(adata_mes, n_neighbors = 15, n_pcs = 25, use_rep = 'X_pca_harmony' )
scp.tl.leiden(adata_mes)
scp.tl.umap(adata_mes,n_components = 2)
scp.pl.umap(adata_mes, color=['leiden'], save = 'all-mes-harmony-sample_type.png')
adata_mes.write('DATA/all-mes-harmony-sample_type.h5ad')


# In[ ]:


scp.pl.umap(adata, color=['leiden'], legend_loc="on data", legend_fontsize = 10, size=15)
scp.pl.umap(adata, color=['Sample'], legend_loc="on data", size=15)
scp.pl.umap(adata, color=['Age'], legend_loc="on data", size=15)
scp.pl.umap(adata, color=['Origin'], legend_loc="on data", size=15)
scp.pl.umap(adata, color=['leiden', 'Age', 'Origin'],size=15)
scp.pl.umap(adata, color=['leiden'],size=15)
scp.pl.umap(adata, color=['Sample'],size=15)
scp.pl.umap(adata, color=['Age'],size=15)
scp.pl.umap(adata, color=['Origin'],size=15)
scp.pl.umap(adata, color=['Sample', 'Age', 'Origin','leiden'], legend_loc="on data",size=15)
scp.pl.umap(adata, color=['leiden'], legend_loc="on data",size=15)
scp.pl.umap(adata, color=['Sample'],size=15)
scp.pl.umap(adata, color=['leiden', 'VCAN', 'SOX6','MSX1','SCX','RGS5', 'COL3A1', 'COL5A1', 'RUNX2','COL1A1','TNMD', 'COL2A1'], legend_loc="on data",size=15)


# In[ ]:


adata = scp.read_h5ad('DATA/all_samples-harmony5.h5ad')
scp.pl.umap(adata, color=['leiden'], legend_loc="on data", legend_fontsize = 10)
adata_mes = adata[adata.obs['leiden'].isin(['0','1','2','3','4','5','6','7','8','9','10','12','14','19']),:].copy()
adata_mes.obs.leiden=adata_mes.obs.leiden.cat.rename_categories(range(len(adata_mes.obs.leiden.cat.categories))).astype(str)
scp.pl.umap(adata_mes, color=['leiden'], legend_loc="on data")
scp.tl.pca(adata_mes,n_comps = 30, svd_solver='arpack',use_highly_variable=True)
#sce.pp.harmony_integrate(adata_mes, 'Sample_type')
scp.pp.neighbors(adata_mes, n_neighbors = 15, n_pcs = 25, use_rep = 'X_pca' )
scp.tl.leiden(adata_mes)
scp.tl.umap(adata_mes,n_components = 2)
scp.pl.umap(adata_mes, color=['leiden'], save = 'all-mes.png')
adata_mes.write('DATA/all-mes.h5ad')


# In[ ]:


scp.pp.neighbors(adata_mes, n_neighbors = 30, n_pcs = 25, use_rep = 'X_pca' )
scp.tl.leiden(adata_mes)
scp.tl.umap(adata_mes,n_components = 2)
scp.pl.umap(adata_mes, color=['leiden'], save = 'all-mes30.png')
adata_mes.write('DATA/all-mes30.h5ad')

