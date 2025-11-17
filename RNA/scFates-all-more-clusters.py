#!/usr/bin/env python
# coding: utf-8

# In[1]:


import cellrank
import scvelo as scv
import os
os.chdir("/home/yakov/face_data/scFates")
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


# In[4]:


import os, sys
os.environ['R_HOME'] = sys.exec_prefix+"/lib/R/"
import scFates as scf
import warnings
warnings.filterwarnings("ignore")
import os, sys
# to avoid any possible jupyter crashes due to rpy2 not finding the R install on conda
os.environ['R_HOME'] = sys.exec_prefix+"/lib/R/"
from anndata import AnnData
import numpy as np
import pandas as pd
import scanpy as sc
import scFates as scf
import palantir
import matplotlib.pyplot as plt
sc.settings.verbosity = 3
sc.settings.logfile = sys.stdout

## fix palantir breaking down some plots
import seaborn
seaborn.reset_orig()
get_ipython().run_line_magic('matplotlib', 'inline')
sc.set_figure_params()
scf.set_figure_pubready()


# In[11]:


adata = scp.read_h5ad('all_samples_55clusters.h5ad')


# In[12]:


#adata.obs['general_annotation'] = "Mesenchyme"
for cell in adata.obs.loc[adata.obs.leiden == '0'].index:
    adata.obs.at[cell, 'general_annotation'] = "Myoblast1"
    adata.obs.at[cell, 'specific_annotation'] = "Myoblast1"
for cell in adata.obs.loc[adata.obs.leiden == '1'].index:
    adata.obs.at[cell, 'general_annotation'] = "Osteoblast"
    adata.obs.at[cell, 'specific_annotation'] = "Osteoblasts and pre-skeletal mesenchyme"
for cell in adata.obs.loc[adata.obs.leiden == '2'].index:
    adata.obs.at[cell, 'general_annotation'] = "CXCL14 DerM"
    adata.obs.at[cell, 'specific_annotation'] = "CXCL14+ dermal mesenchyme"
for cell in adata.obs.loc[adata.obs.leiden == '3'].index:
    adata.obs.at[cell, 'general_annotation'] = "FOXP EOM"
    adata.obs.at[cell, 'specific_annotation'] = "FOXP+ extraocular mesenchyme"  
for cell in adata.obs.loc[adata.obs.leiden == '4'].index:
    adata.obs.at[cell, 'general_annotation'] = "PTN FB"
    adata.obs.at[cell, 'specific_annotation'] = "PTN+ w9 ear skeletal fibroblasts"
for cell in adata.obs.loc[adata.obs.leiden == '5'].index:
    adata.obs.at[cell, 'general_annotation'] = "SalivGlaM"
    adata.obs.at[cell, 'specific_annotation'] = "Salivary gland mesenchyme"
for cell in adata.obs.loc[adata.obs.leiden == '6'].index:
    adata.obs.at[cell, 'general_annotation'] = "LinguaM"
    adata.obs.at[cell, 'specific_annotation'] = "Lingual mesenchyme"
for cell in adata.obs.loc[adata.obs.leiden == '7'].index:
    adata.obs.at[cell, 'general_annotation'] = "Early MPC"
    adata.obs.at[cell, 'specific_annotation'] = "Early mesenchymal progenitor cell"
for cell in adata.obs.loc[adata.obs.leiden == '8'].index:
    adata.obs.at[cell, 'general_annotation'] = "Nasal StroM1"
    adata.obs.at[cell, 'specific_annotation'] = "PAX3+ nasal stromal mesenchyme"
for cell in adata.obs.loc[adata.obs.leiden == '9'].index:
    adata.obs.at[cell, 'general_annotation'] = "Nasal StroM2"
    adata.obs.at[cell, 'specific_annotation'] = "PAX7/PAX9+ nasal stromal mesenchyme"
for cell in adata.obs.loc[adata.obs.leiden == '10'].index:
    adata.obs.at[cell, 'general_annotation'] = "Chondrocyte1"
    adata.obs.at[cell, 'specific_annotation'] = "Chondrocyte 1"
for cell in adata.obs.loc[adata.obs.leiden == '11'].index:
    adata.obs.at[cell, 'general_annotation'] = "Pericyte"
    adata.obs.at[cell, 'specific_annotation'] = "Perivascular cells and pericytes"
for cell in adata.obs.loc[adata.obs.leiden == '12'].index:
    adata.obs.at[cell, 'general_annotation'] = "DentaM"
    adata.obs.at[cell, 'specific_annotation'] = "Dental mesenchyme"
for cell in adata.obs.loc[adata.obs.leiden == '13'].index:
    adata.obs.at[cell, 'general_annotation'] = "DLK1/CXCL stroM"
    adata.obs.at[cell, 'specific_annotation'] = "DLK1+ CXCL12+ stromal mesenchyme"
for cell in adata.obs.loc[adata.obs.leiden == '14'].index:
    adata.obs.at[cell, 'general_annotation'] = "Chondrocyte2"
    adata.obs.at[cell, 'specific_annotation'] = "Chondrocyte 2"
for cell in adata.obs.loc[adata.obs.leiden == '15'].index:
    adata.obs.at[cell, 'general_annotation'] = "ConnectiveMPC"
    adata.obs.at[cell, 'specific_annotation'] = "Connective mesenchymal progenitor cell"
for cell in adata.obs.loc[adata.obs.leiden == '16'].index:
    adata.obs.at[cell, 'general_annotation'] = "SOX11 MPC"
    adata.obs.at[cell, 'specific_annotation'] = "SOX11+ mesenchyme progenitor cell"    
for cell in adata.obs.loc[adata.obs.leiden == '17'].index:
    adata.obs.at[cell, 'general_annotation'] = "Ant FB"
    adata.obs.at[cell, 'specific_annotation'] = "Anterior fibroblast"
for cell in adata.obs.loc[adata.obs.leiden == '18'].index:
    adata.obs.at[cell, 'general_annotation'] = "Post MPC"
    adata.obs.at[cell, 'specific_annotation'] = "Posterior and mandibular mesenchyme"
for cell in adata.obs.loc[adata.obs.leiden == '19'].index:
    adata.obs.at[cell, 'general_annotation'] = "DerMPC"
    adata.obs.at[cell, 'specific_annotation'] = "Dermal progenitor mesenchyme"
for cell in adata.obs.loc[adata.obs.leiden == '20'].index:
    adata.obs.at[cell, 'general_annotation'] = "Mature ant DerFB"
    adata.obs.at[cell, 'specific_annotation'] = "Mature anterior dermal fibroblast"
for cell in adata.obs.loc[adata.obs.leiden == '21'].index:
    adata.obs.at[cell, 'general_annotation'] = "Fibroblast"
    adata.obs.at[cell, 'specific_annotation'] = "Fibroblast"
for cell in adata.obs.loc[adata.obs.leiden == '21'].index:
    adata.obs.at[cell, 'general_annotation'] = "NPC1"
    adata.obs.at[cell, 'specific_annotation'] = "Neural progenitor cell 1"
for cell in adata.obs.loc[adata.obs.leiden == '22'].index:
    adata.obs.at[cell, 'general_annotation'] = "Glia"
    adata.obs.at[cell, 'specific_annotation'] = "Glia"
for cell in adata.obs.loc[adata.obs.leiden == '23'].index:
    adata.obs.at[cell, 'general_annotation'] = "DerFB"
    adata.obs.at[cell, 'specific_annotation'] = "Dermal fibroblast"
for cell in adata.obs.loc[adata.obs.leiden == '24'].index:
    adata.obs.at[cell, 'general_annotation'] = "Meninge"
    adata.obs.at[cell, 'specific_annotation'] = "Leptomeninges-like cells"
for cell in adata.obs.loc[adata.obs.leiden == '25'].index:
    adata.obs.at[cell, 'general_annotation'] = "Mature POSTN FB"
    adata.obs.at[cell, 'specific_annotation'] = "Mature POSTN+ fibroblast"
for cell in adata.obs.loc[adata.obs.leiden == '26'].index:
    adata.obs.at[cell, 'general_annotation'] = "Tenocyte"
    adata.obs.at[cell, 'specific_annotation'] = "Tenocyte"
for cell in adata.obs.loc[adata.obs.leiden == '27'].index:
    adata.obs.at[cell, 'general_annotation'] = "POM_ScleraM"
    adata.obs.at[cell, 'specific_annotation'] = "Periocular and sclera mesenchyme"
for cell in adata.obs.loc[adata.obs.leiden == '28'].index:
    adata.obs.at[cell, 'general_annotation'] = "Stress FB1"
    adata.obs.at[cell, 'specific_annotation'] = "Stress signature fibroblast"
for cell in adata.obs.loc[adata.obs.leiden == '29'].index:
    adata.obs.at[cell, 'general_annotation'] = "Stress FB2"
    adata.obs.at[cell, 'specific_annotation'] = "Stress signature fibroblast2"
for cell in adata.obs.loc[adata.obs.leiden == '30'].index:
    adata.obs.at[cell, 'general_annotation'] = "PrePerivaM"
    adata.obs.at[cell, 'specific_annotation'] = "Perivascular mesenchyme"
for cell in adata.obs.loc[adata.obs.leiden == '31'].index:
    adata.obs.at[cell, 'general_annotation'] = "Endothelia"
    adata.obs.at[cell, 'specific_annotation'] = "Endothelial cell"
for cell in adata.obs.loc[adata.obs.leiden == '32'].index:
    adata.obs.at[cell, 'general_annotation'] = "Myoblast2"
    adata.obs.at[cell, 'specific_annotation'] = "Myoblast 2"
for cell in adata.obs.loc[adata.obs.leiden == '33'].index:
    adata.obs.at[cell, 'general_annotation'] = "Ant_MC_2"
    adata.obs.at[cell, 'specific_annotation'] = "Anterior (eye and maxillary region) mesenchyme"
for cell in adata.obs.loc[adata.obs.leiden == '34'].index:
    adata.obs.at[cell, 'general_annotation'] = "PCM"
    adata.obs.at[cell, 'specific_annotation'] = "Perichondrial mesenchyme"
for cell in adata.obs.loc[adata.obs.leiden == '35'].index:
    adata.obs.at[cell, 'general_annotation'] = "Subdermal"
    adata.obs.at[cell, 'specific_annotation'] = "Subdermal"
for cell in adata.obs.loc[adata.obs.leiden == '36'].index:
    adata.obs.at[cell, 'general_annotation'] = "Skeletal MPC"
    adata.obs.at[cell, 'specific_annotation'] = "Skeletal progenitor mesenchyme"
for cell in adata.obs.loc[adata.obs.leiden == '37'].index:
    adata.obs.at[cell, 'general_annotation'] = "Myeloid"
    adata.obs.at[cell, 'specific_annotation'] = "Immune"
for cell in adata.obs.loc[adata.obs.leiden == '38'].index:
    adata.obs.at[cell, 'general_annotation'] = "Erythrocytes"
    adata.obs.at[cell, 'specific_annotation'] = "Erythrocytes"
for cell in adata.obs.loc[adata.obs.leiden == '39'].index:
    adata.obs.at[cell, 'general_annotation'] = "Chondroblast"
    adata.obs.at[cell, 'specific_annotation'] = "Chondroblast"
for cell in adata.obs.loc[adata.obs.leiden == '40'].index:
    adata.obs.at[cell, 'general_annotation'] = "Neurons 1"
    adata.obs.at[cell, 'specific_annotation'] = "Neurons 1" 
for cell in adata.obs.loc[adata.obs.leiden == '41'].index:
    adata.obs.at[cell, 'general_annotation'] = "Myoblast3"
    adata.obs.at[cell, 'specific_annotation'] = "Myoblast3"
for cell in adata.obs.loc[adata.obs.leiden == '42'].index:
    adata.obs.at[cell, 'general_annotation'] = "Neurons 2"
    adata.obs.at[cell, 'specific_annotation'] = "Neurons 2"
for cell in adata.obs.loc[adata.obs.leiden == '43'].index:
    adata.obs.at[cell, 'general_annotation'] = "WNT5A stroM"
    adata.obs.at[cell, 'specific_annotation'] = "WNT5A+ stroma"
for cell in adata.obs.loc[adata.obs.leiden == '44'].index:
    adata.obs.at[cell, 'general_annotation'] = "Neurons 3"
    adata.obs.at[cell, 'specific_annotation'] = "Neurons 3"
for cell in adata.obs.loc[adata.obs.leiden == '45'].index:
    adata.obs.at[cell, 'general_annotation'] = "Peripheral neurons"
    adata.obs.at[cell, 'specific_annotation'] = "Peripheral neurons"
for cell in adata.obs.loc[adata.obs.leiden == '46'].index:
    adata.obs.at[cell, 'general_annotation'] = "Epithelia"
    adata.obs.at[cell, 'specific_annotation'] = "Epithelia"
for cell in adata.obs.loc[adata.obs.leiden == '47'].index:
    adata.obs.at[cell, 'general_annotation'] = "Ant MPC"
    adata.obs.at[cell, 'specific_annotation'] = "Anterior mesenchymal progenitor"
for cell in adata.obs.loc[adata.obs.leiden == '48'].index:
    adata.obs.at[cell, 'general_annotation'] = "NPC2"
    adata.obs.at[cell, 'specific_annotation'] = "Neural progenitor cell 2"
for cell in adata.obs.loc[adata.obs.leiden == '49'].index:
    adata.obs.at[cell, 'general_annotation'] = "Chondrocyte3"
    adata.obs.at[cell, 'specific_annotation'] = "Nasal chondrocyte"
for cell in adata.obs.loc[adata.obs.leiden == '50'].index:
    adata.obs.at[cell, 'general_annotation'] = "Keratocytes"
    adata.obs.at[cell, 'specific_annotation'] = "Keratocytes"
for cell in adata.obs.loc[adata.obs.leiden == '51'].index:
    adata.obs.at[cell, 'general_annotation'] = "T cells"
    adata.obs.at[cell, 'specific_annotation'] = "T cells" 
for cell in adata.obs.loc[adata.obs.leiden == '52'].index:
    adata.obs.at[cell, 'general_annotation'] = "Doublet"
    adata.obs.at[cell, 'specific_annotation'] = "Doublet" 
for cell in adata.obs.loc[adata.obs.leiden == '53'].index:
    adata.obs.at[cell, 'general_annotation'] = "Erythrocytes2"
    adata.obs.at[cell, 'specific_annotation'] = "Erythrocytes2" 
for cell in adata.obs.loc[adata.obs.leiden == '54'].index:
    adata.obs.at[cell, 'general_annotation'] = "B cells"
    adata.obs.at[cell, 'specific_annotation'] = "B cells" 
for cell in adata.obs.loc[adata.obs.leiden == '55'].index:
    adata.obs.at[cell, 'general_annotation'] = "Melanocytes"
    adata.obs.at[cell, 'specific_annotation'] = "Melanocytes"     


# In[13]:


adata.obs


# In[14]:


import pandas as pd

meta = adata.obs.copy()
meta['UMAP_1'] = adata.obsm['X_umap'][:,0]
meta['UMAP_2'] = adata.obsm['X_umap'][:,1]
meta.to_csv('cell_metadata.csv')


# In[7]:


adata = scp.read_h5ad('all_samples_55clusters.h5ad')


# In[7]:


scp.pp.neighbors(adata)
adata.layers["spliced"] = adata.X
adata.layers["unspliced"] = adata.X
#scv.pp.moments(adata, n_neighbors = 15, n_pcs = 25 )
scv.pp.moments(adata)
ctk = cellrank.kernels.CytoTRACEKernel(adata)
ckt = ctk.compute_cytotrace().compute_transition_matrix()
scp.set_figure_params(dpi=150)
scp.pl.umap(adata,color="ct_score",)


# In[29]:


scp.pl.umap(adata,color="ct_score", save = 'Fig1H-cytoTRACE.pdf', size = 10)


# In[12]:


import pandas as pd

# Define column names
cluster_col = 'leiden'
score_col = 'ct_score'

# Calculate average ct_score per Leiden cluster
avg_scores = adata.obs.groupby(cluster_col)[score_col].mean()

# Convert to DataFrame
avg_scores_df = avg_scores.reset_index()
avg_scores_df.columns = [cluster_col, 'avg_ct_score']

# Save to CSV
avg_scores_df.to_csv('avg_ct_score_per_cluster.csv', index=False)

print("✅ Saved average ct_score per Leiden cluster to avg_ct_score_per_cluster.csv")
avg_scores_df.head()


# In[40]:


# Ensure ct_score is numeric
adata.obs['ct_score'] = pd.to_numeric(adata.obs['ct_score'], errors='coerce')

# Calculate average ct_score per Leiden cluster
cluster_avg_ct = adata.obs.groupby('leiden')['ct_score'].mean()

# Map the cluster average to all cells
adata.obs['Cluster_ct_score'] = adata.obs['leiden'].map(cluster_avg_ct)

print("✅ Added adata.obs['Cluster_ct_score'] with average ct_score per Leiden cluster")


# In[17]:


import pandas as pd

cluster_col = 'leiden'
age_col = 'Age'

# Convert Age to numeric (ignore any non-numeric entries)
adata.obs[age_col] = pd.to_numeric(adata.obs[age_col], errors='coerce')

# Drop missing or invalid values
adata.obs = adata.obs.dropna(subset=[age_col])

# Calculate average Age per Leiden cluster
avg_age = adata.obs.groupby(cluster_col)[age_col].mean()

# Convert to DataFrame
avg_age_df = avg_age.reset_index()
avg_age_df.columns = [cluster_col, 'avg_Age']

# Save to CSV
avg_age_df.to_csv('avg_Age_per_cluster.csv', index=False)

print("✅ Saved average Age per Leiden cluster to avg_Age_per_cluster.csv")
avg_age_df.head()


# In[19]:


# Age Origin
import pandas as pd
import numpy as np

# Make sure Age is numeric
adata.obs['Age'] = adata.obs['Age'].astype(float)

# Compute average age per cluster
cluster_age = adata.obs.groupby('leiden')['Age'].mean()

# Assign average cluster age back to all cells
adata.obs['ClusterAge'] = adata.obs['leiden'].map(cluster_age)

print("✅ Added adata.obs['ClusterAge'] (average Age per Leiden cluster)")


# In[28]:


import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import scanpy as sc

# Make sure ClusterAge is numeric
adata.obs['ClusterAge'] = adata.obs['ClusterAge'].astype(float)

# Define custom colormap
colors = [
    (0.5, 0.8, 1.0),  # light blue
    (0.6, 1.0, 0.6),  # light green
    (1.0, 1.0, 0.6),  # yellow
    (1.0, 0.8, 0.6),  # light orange
    (1.0, 0.6, 0.6)   # light red
]
my_cmap = LinearSegmentedColormap.from_list("blue_green_yellow_red", colors, N=256)

# Plot UMAP using scanpy (not scp) and force the colorbar
sc.pl.umap(
    adata,
    color='ClusterAge',
    cmap=my_cmap,
    colorbar_loc='right',   # explicitly request colorbar
    size=5,
    title='Average age of cluster (week)',
    save='_Age_custom.pdf',
    show=True
)


# In[33]:


scp.pl.umap(adata, color=['ClusterAge'],cmap='turbo', size=1, title = 'Average age of cluster (week)', save = 'Age-55clusters.pdf')


# In[ ]:


import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import scanpy as sc

# Make sure ClusterAge is numeric
adata.obs['ClusterAge'] = adata.obs['ClusterAge'].astype(float)

# Define custom colormap
colors = [
    (0.5, 0.8, 1.0),  # light blue
    (0.6, 1.0, 0.6),  # light green
    (1.0, 1.0, 0.6),  # yellow
    (1.0, 0.8, 0.6),  # light orange
    (1.0, 0.6, 0.6)   # light red
]
my_cmap = LinearSegmentedColormap.from_list("blue_green_yellow_red", colors, N=256)

# Plot UMAP using scanpy (not scp) and force the colorbar
sc.pl.umap(
    adata,
    color='ClusterAge',
    cmap=my_cmap,
    colorbar_loc='right',   # explicitly request colorbar
    size=5,
    title='Average age of cluster (week)',
    save='_Age_custom.pdf',
    show=True
)


# In[39]:


#scp.pl.umap(adata, color=['ClusterAge'],cmap='turbo', size=20, title = 'Average age of cluster (week)', save = 'Age-55clusters.pdf')
scp.pl.umap(
    adata,
    color=['ClusterAge'],
    cmap='turbo',
    size=10,
    title='Average age of cluster (week)',
    save='Age-55clusters.pdf',
    show=True
)
# Then increase resolution when saving
plt.savefig('Age-55clusters_highdpi.pdf', dpi=300)


# In[44]:


import scanpy as sc
import matplotlib.pyplot as plt

# Ensure the new column is numeric
adata.obs['Cluster_ct_score'] = pd.to_numeric(adata.obs['Cluster_ct_score'], errors='coerce')

# Plot UMAP colored by cluster-average ct_score using Turbo
sc.pl.umap(
    adata,
    color='Cluster_ct_score',
    cmap='turbo',       # <--- use only cmap
    size=20,            # adjust marker size
    title='Average ct_score per cluster',
    colorbar_loc='right',  # show the scale bar
    save='_Cluster_ct_score_turbo.pdf',
    show=True
)

# Optional: save a high-resolution version manually
plt.savefig('Cluster_ct_score_umap_turbo_highdpi.pdf', dpi=300)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


ct_score


# In[ ]:





# In[ ]:


adata = scp.read_h5ad('all_samples_55clusters.h5ad')


# In[8]:


#adata.obs['general_annotation'] = "Mesenchyme"
for cell in adata.obs.loc[adata.obs.leiden == '0'].index:
    adata.obs.at[cell, 'general_annotation'] = "Myoblast"
    adata.obs.at[cell, 'specific_annotation'] = "Myoblast"
for cell in adata.obs.loc[adata.obs.leiden == '1'].index:
    adata.obs.at[cell, 'general_annotation'] = "Osteoblast"
    adata.obs.at[cell, 'specific_annotation'] = "Osteoblast"
for cell in adata.obs.loc[adata.obs.leiden == '2'].index:
    adata.obs.at[cell, 'general_annotation'] = "Dermal"
    adata.obs.at[cell, 'specific_annotation'] = "CXCL14+ dermal mesenchme"
for cell in adata.obs.loc[adata.obs.leiden == '3'].index:
    adata.obs.at[cell, 'general_annotation'] = "Periocular"
    adata.obs.at[cell, 'specific_annotation'] = "FOXP+ extraocular mesenchyme"  
for cell in adata.obs.loc[adata.obs.leiden == '4'].index:
    adata.obs.at[cell, 'general_annotation'] = "Osteogenic"
    adata.obs.at[cell, 'specific_annotation'] = "PTN+ w9 ear skeletal fibroblasts"
for cell in adata.obs.loc[adata.obs.leiden == '5'].index:
    adata.obs.at[cell, 'general_annotation'] = "Gland"
    adata.obs.at[cell, 'specific_annotation'] = "salivary or submandibular gland mesenchyme"
for cell in adata.obs.loc[adata.obs.leiden == '6'].index:
    adata.obs.at[cell, 'general_annotation'] = "lingual"
    adata.obs.at[cell, 'specific_annotation'] = "MEIS+ lingual mesenchyme"
for cell in adata.obs.loc[adata.obs.leiden == '7'].index:
    adata.obs.at[cell, 'general_annotation'] = "Progenitor"
    adata.obs.at[cell, 'specific_annotation'] = "Progenitors"
for cell in adata.obs.loc[adata.obs.leiden == '8'].index:
    adata.obs.at[cell, 'general_annotation'] = "Stroma"
    adata.obs.at[cell, 'specific_annotation'] = "PAX3+ nasal dermal mesenchyme"
for cell in adata.obs.loc[adata.obs.leiden == '9'].index:
    adata.obs.at[cell, 'general_annotation'] = "Stroma"
    adata.obs.at[cell, 'specific_annotation'] = "PAX9/PAX7+ nasal perichondrium"
for cell in adata.obs.loc[adata.obs.leiden == '10'].index:
    adata.obs.at[cell, 'general_annotation'] = "Chondrocyte"
    adata.obs.at[cell, 'specific_annotation'] = "Chondrocyte 1"
for cell in adata.obs.loc[adata.obs.leiden == '11'].index:
    adata.obs.at[cell, 'general_annotation'] = "Perivascular"
    adata.obs.at[cell, 'specific_annotation'] = "Pericyte"
for cell in adata.obs.loc[adata.obs.leiden == '12'].index:
    adata.obs.at[cell, 'general_annotation'] = "Dental"
    adata.obs.at[cell, 'specific_annotation'] = "Dental 1"
for cell in adata.obs.loc[adata.obs.leiden == '13'].index:
    adata.obs.at[cell, 'general_annotation'] = "Stroma"
    adata.obs.at[cell, 'specific_annotation'] = "DLK1+ CXCL12+ stroma"
for cell in adata.obs.loc[adata.obs.leiden == '14'].index:
    adata.obs.at[cell, 'general_annotation'] = "Chondrocyte"
    adata.obs.at[cell, 'specific_annotation'] = "Chondrocyte 2"
for cell in adata.obs.loc[adata.obs.leiden == '15'].index:
    adata.obs.at[cell, 'general_annotation'] = "Connective"
    adata.obs.at[cell, 'specific_annotation'] = "connective progenitor mesenchyme"
for cell in adata.obs.loc[adata.obs.leiden == '16'].index:
    adata.obs.at[cell, 'general_annotation'] = "Mesenchyme"
    adata.obs.at[cell, 'specific_annotation'] = "SOX11+ progenitor mesenchyme"    
for cell in adata.obs.loc[adata.obs.leiden == '17'].index:
    adata.obs.at[cell, 'general_annotation'] = "Fibroblast"
    adata.obs.at[cell, 'specific_annotation'] = "Fibroblast 1"
for cell in adata.obs.loc[adata.obs.leiden == '18'].index:
    adata.obs.at[cell, 'general_annotation'] = "Mesenchyme"
    adata.obs.at[cell, 'specific_annotation'] = "posterior mandibular mesenchyme"
for cell in adata.obs.loc[adata.obs.leiden == '19'].index:
    adata.obs.at[cell, 'general_annotation'] = "Mesenchyme"
    adata.obs.at[cell, 'specific_annotation'] = "Dermal progenitor mesenchyme"
for cell in adata.obs.loc[adata.obs.leiden == '20'].index:
    adata.obs.at[cell, 'general_annotation'] = "Dermal"
    adata.obs.at[cell, 'specific_annotation'] = "Mature anterior dermal fibroblast"
for cell in adata.obs.loc[adata.obs.leiden == '21'].index:
    adata.obs.at[cell, 'general_annotation'] = "Fibroblast"
    adata.obs.at[cell, 'specific_annotation'] = "Fibroblast"
for cell in adata.obs.loc[adata.obs.leiden == '21'].index:
    adata.obs.at[cell, 'general_annotation'] = "Neural"
    adata.obs.at[cell, 'specific_annotation'] = "Neural progenitors"
for cell in adata.obs.loc[adata.obs.leiden == '22'].index:
    adata.obs.at[cell, 'general_annotation'] = "Glia"
    adata.obs.at[cell, 'specific_annotation'] = "Glia"
for cell in adata.obs.loc[adata.obs.leiden == '23'].index:
    adata.obs.at[cell, 'general_annotation'] = "Dermal"
    adata.obs.at[cell, 'specific_annotation'] = "Dermal fibroblast"
for cell in adata.obs.loc[adata.obs.leiden == '24'].index:
    adata.obs.at[cell, 'general_annotation'] = "Meninges"
    adata.obs.at[cell, 'specific_annotation'] = "Leptomeninges"
for cell in adata.obs.loc[adata.obs.leiden == '25'].index:
    adata.obs.at[cell, 'general_annotation'] = "Mesenchyme"
    adata.obs.at[cell, 'specific_annotation'] = "POSTN+ late mesenchyme"
for cell in adata.obs.loc[adata.obs.leiden == '26'].index:
    adata.obs.at[cell, 'general_annotation'] = "Connective"
    adata.obs.at[cell, 'specific_annotation'] = "Tenocyte"
for cell in adata.obs.loc[adata.obs.leiden == '27'].index:
    adata.obs.at[cell, 'general_annotation'] = "Meninges-like"
    adata.obs.at[cell, 'specific_annotation'] = "Eye meningeal-like"
for cell in adata.obs.loc[adata.obs.leiden == '28'].index:
    adata.obs.at[cell, 'general_annotation'] = "Dermal"
    adata.obs.at[cell, 'specific_annotation'] = "Stress dermis"
for cell in adata.obs.loc[adata.obs.leiden == '29'].index:
    adata.obs.at[cell, 'general_annotation'] = "Dermal"
    adata.obs.at[cell, 'specific_annotation'] = "Dermal mesenchyme"
for cell in adata.obs.loc[adata.obs.leiden == '30'].index:
    adata.obs.at[cell, 'general_annotation'] = "Perivascular"
    adata.obs.at[cell, 'specific_annotation'] = "Pre-perivascular"
for cell in adata.obs.loc[adata.obs.leiden == '31'].index:
    adata.obs.at[cell, 'general_annotation'] = "Endothelia"
    adata.obs.at[cell, 'specific_annotation'] = "Endothelia"
for cell in adata.obs.loc[adata.obs.leiden == '32'].index:
    adata.obs.at[cell, 'general_annotation'] = "Myoblast"
    adata.obs.at[cell, 'specific_annotation'] = "Myoblast"
for cell in adata.obs.loc[adata.obs.leiden == '33'].index:
    adata.obs.at[cell, 'general_annotation'] = "Mesenchyme"
    adata.obs.at[cell, 'specific_annotation'] = "Eye dermal mesenchyme"
for cell in adata.obs.loc[adata.obs.leiden == '34'].index:
    adata.obs.at[cell, 'general_annotation'] = "Meninges"
    adata.obs.at[cell, 'specific_annotation'] = "Outer dura"
for cell in adata.obs.loc[adata.obs.leiden == '35'].index:
    adata.obs.at[cell, 'general_annotation'] = "Subdermal"
    adata.obs.at[cell, 'specific_annotation'] = "Subdermal"
for cell in adata.obs.loc[adata.obs.leiden == '36'].index:
    adata.obs.at[cell, 'general_annotation'] = "Skeletal"
    adata.obs.at[cell, 'specific_annotation'] = "Skeletal progenitor mesenchyme"
for cell in adata.obs.loc[adata.obs.leiden == '37'].index:
    adata.obs.at[cell, 'general_annotation'] = "Immune"
    adata.obs.at[cell, 'specific_annotation'] = "Myeloid"
for cell in adata.obs.loc[adata.obs.leiden == '38'].index:
    adata.obs.at[cell, 'general_annotation'] = "Blood"
    adata.obs.at[cell, 'specific_annotation'] = "Erythrocytes"
for cell in adata.obs.loc[adata.obs.leiden == '39'].index:
    adata.obs.at[cell, 'general_annotation'] = "Preskeletal"
    adata.obs.at[cell, 'specific_annotation'] = "Chondroblast"
for cell in adata.obs.loc[adata.obs.leiden == '40'].index:
    adata.obs.at[cell, 'general_annotation'] = "Neural"
    adata.obs.at[cell, 'specific_annotation'] = "Neurons 1" 
for cell in adata.obs.loc[adata.obs.leiden == '41'].index:
    adata.obs.at[cell, 'general_annotation'] = "Myoblast"
    adata.obs.at[cell, 'specific_annotation'] = "Myoblast"
for cell in adata.obs.loc[adata.obs.leiden == '42'].index:
    adata.obs.at[cell, 'general_annotation'] = "Neural"
    adata.obs.at[cell, 'specific_annotation'] = "Neurons 2"
for cell in adata.obs.loc[adata.obs.leiden == '43'].index:
    adata.obs.at[cell, 'general_annotation'] = "Stroma"
    adata.obs.at[cell, 'specific_annotation'] = "WNT5A+ stroma"
for cell in adata.obs.loc[adata.obs.leiden == '44'].index:
    adata.obs.at[cell, 'general_annotation'] = "Neural"
    adata.obs.at[cell, 'specific_annotation'] = "Neurons 3"
for cell in adata.obs.loc[adata.obs.leiden == '45'].index:
    adata.obs.at[cell, 'general_annotation'] = "Neural"
    adata.obs.at[cell, 'specific_annotation'] = "Peripheral neurons"
for cell in adata.obs.loc[adata.obs.leiden == '46'].index:
    adata.obs.at[cell, 'general_annotation'] = "Epithelial"
    adata.obs.at[cell, 'specific_annotation'] = "Epithelia"
for cell in adata.obs.loc[adata.obs.leiden == '47'].index:
    adata.obs.at[cell, 'general_annotation'] = "Progenitors"
    adata.obs.at[cell, 'specific_annotation'] = "Anterior mesenchymal progenitor"
for cell in adata.obs.loc[adata.obs.leiden == '48'].index:
    adata.obs.at[cell, 'general_annotation'] = "Neural"
    adata.obs.at[cell, 'specific_annotation'] = "Neural progenitor"
for cell in adata.obs.loc[adata.obs.leiden == '49'].index:
    adata.obs.at[cell, 'general_annotation'] = "Chondrocyte"
    adata.obs.at[cell, 'specific_annotation'] = "Nasal chondrocyte"
for cell in adata.obs.loc[adata.obs.leiden == '50'].index:
    adata.obs.at[cell, 'general_annotation'] = "Cornea"
    adata.obs.at[cell, 'specific_annotation'] = "Keratocytes"
for cell in adata.obs.loc[adata.obs.leiden == '51'].index:
    adata.obs.at[cell, 'general_annotation'] = "Immune"
    adata.obs.at[cell, 'specific_annotation'] = "T cells" 
for cell in adata.obs.loc[adata.obs.leiden == '52'].index:
    adata.obs.at[cell, 'general_annotation'] = "Eye"
    adata.obs.at[cell, 'specific_annotation'] = "Eye neural mesenchymal doublet" 
for cell in adata.obs.loc[adata.obs.leiden == '53'].index:
    adata.obs.at[cell, 'general_annotation'] = "Blood"
    adata.obs.at[cell, 'specific_annotation'] = "Blood" 
for cell in adata.obs.loc[adata.obs.leiden == '54'].index:
    adata.obs.at[cell, 'general_annotation'] = "Immune"
    adata.obs.at[cell, 'specific_annotation'] = "B cells" 
for cell in adata.obs.loc[adata.obs.leiden == '55'].index:
    adata.obs.at[cell, 'general_annotation'] = "Melanocytes"
    adata.obs.at[cell, 'specific_annotation'] = "Melanocytes"     


# In[8]:


import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Full gene list
genes = [ "DIRAS3","RNU5D-1","TP73","LRRTM1","GPR1","ZDBF2","NAP1L5","RHOBTB3","ERAP2","VTRNA2-1","NR3C1","HMMR",
"ADTRP","FAM50B","PXDC1","LIN28B","AIM1","PLAGL1","HYMAI","SLC22A2*","SLC22A3*","GRB10","DDC","GLI3",
"HECW1","HOXA4","PEG10","MAGI2","SGCE","PPP1R9A","TFPI2","CCDC71L","COPG2IT1","CPA4","MEST","MESTIT1",
"KLF14","SVOPL","PTPRN2","DLGAP2","KCNK9","PEG13","MYC","ZFAT","ZFAT-AS1","GLIS3","VIM","WT1-AS","WT1",
"KCNQ1DN","OSBPL5","PHLDA2","INS","H19","SLC22A18AS","CDKN1C","IGF2AS","KCNQ1","KCNQ1OT1","IGF2","SLC22A18",
"ANO1","ZC3H12C","NTM","ST8SIA1","RBP5","ATP5F1EP2","RB1","SMOC1","DIO3","MEG3","DLK1","MEG8","SNORD114-1",
"SNORD113-1","DIO3OS","RTL1","MAGEL2","UBE3A","MKRN3","SNORD115-48","SNORD115@","SNORD116","ATP10A",
"SNORD109A","SNORD107","SNORD108","SNRPN","SNORD109B","PWAR6","PWCR1","NPAP1","NDN","SNORD64","SNURF",
"RASGRF1","IRAIN","ZNF597","PRR25","NAA60","CMTM1","ZFP90","ZNF396","TCEB3C","PARD6G","PARD6G-AS1","ADNP2",
"DNMT1","AXL","ZIM2","PEG3","MIR371A","NLRP2","PEG3-AS1","MIMT1","BLCAP","NNAT","MCTS2","GDAP1L1","SGK2",
"GNAS","L3MBTL1","MIR298","GNAS-AS1","MIR296","MIR125B2","DSCAM","DGCR6L","DGCR6" ]

# Filter genes present in adata
existing_genes = [g for g in genes if g in adata.var_names]

# Compute average per cluster (specific_annotation)
cluster_avg = adata.to_df()[existing_genes].groupby(adata.obs['specific_annotation']).mean().T

# Plot clustered heatmap of averages
with PdfPages("imprinted_genes_cluster_avg2.pdf") as pdf:
    g = sns.clustermap(
        cluster_avg,
        cmap='viridis',
        standard_scale=0,
        figsize=(25, 40),
        row_cluster=True,
        col_cluster=True,
        xticklabels=True,
        yticklabels=True
    )
    pdf.savefig(g.fig)
    plt.close(g.fig)


# In[31]:


import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Full gene list
genes = [ "NDUFA4P1", "GFI1", "DIRAS3", "FUCA1", "RNU5D-1","BMP8B","DVL1","WDR8","TP73","RPL22",
          "PRDM16","PEX10","TMEM52","HSPA6","PTPN14","OBSCN","HIST3H2BB","OR11L1","LRRTM1","OTX1",
          "VAX2","COMMD1","CCDC85A","ABCG8","CYP1B1","ZFP36L2","GPR1","ZDBF2","TIGD1","MYEOV2",
          "ALDH1L1","ZIC1","HES1","FGFRL1","SPON2","KIAA1530","NAP1L5","NFKB1","GAB1","ADAMTS16",
          "CDH18","ERAP2","RHOBTB3","CSF2","VTRNA2-1","PRIM2","BTNL2","ADTRP","PXDC1","FAM50B",
          "MRAP2","C6orf117","AIM1","LIN28B","PLAGL1","HYMAI","SLC22A2*","IGF2R","SLC22A3*","BRP44L",
          "GRB10","DDC","GLI3","HECW1","EVX1","HOXA2","HOXA3","HOXA4","HOXA5","HOXA11","RAC1",
          "TMEM60","MAGI2","PEG10","SGCE","ASB4","PPP1R9A","PON2","CALCR","PON1","PON3","DLX5",
          "TFPI2","CCDC71L","COPG2","CPA4","COPG2IT1","MEST","MESTIT1","KLF14","SVOPL","FASTK",
          "SLC4A2","PURG","DLGAP2","NKAIN3","LY6D","KCNK9","GPT","PEG13","ZFAT-AS1","ZFAT","GLIS3",
          "APBA1","C9orf85","FLJ46321","LMX1B","EGFL7","PHPT1","C9orf116","SFMBT2","GATA3","CTNNA3",
          "LDB1","NKX6-2","PAOX","C10orf93","C10orf91","VENTX","INPP5F V2","WT1-AS","WT1","PKP3",
          "AMPD3","KCNQ1OT1","ZNF215","OSBPL5","KCNQ1DN","ASCL2","NAP1L4","TSSC4","H19","TRPM5",
          "SLC22A18AS","TH","IGF2AS","SLC22A18","B4GALNT4","IGF2","CDKN1C","CD81","TSPAN32","PHLDA2",
          "IFITM1","KCNQ1","INS","RAB1B","ANO1","DHCR7","ZC3H12C","KBTBD3","SDHD","NTM","ABCC9",
          "ST8SIA1","RBP5","SLC38A4","SLC26A10","HOXC9","HOXC4","CDK4","E2F7","DCN","HNF1A","KIAA1545",
          "FBRSL1","ATP5F1EP2","HTR2A","RB1","FLJ40296","FAM70B","FOXG1","FERMT2","ESR2","SMOC1",
          "DIO3","MEG3","DLK1","MEG8","RTL1","DIO3OS","SNORD113-1","SNORD114-1","MAGEL2","UBE3A",
          "MKRN3","ZNF127AS","SNORD107","SNORD108","SNORD115@","SNORD109B","SNORD109A","ATP10A","SNRPN",
          "SNORD115-48","SNORD116","PWAR6","PWCR1","NPAP1","NDN","GABRA5","GABRB3","SNURF","GABRG3",
          "SNORD64","GATM","RASGRF1","MIR184","IRAIN","PRR25","ZNF597","SOX8","NAA60","SALL1",
          "C16orf57","CMTM1","ZFP90","ACD","FOXF1","CDH15","TP53","TMEM88","PYY2","HOXB2","HOXB3",
          "LOC100131170","IMPACT","BRUNOL4","FAM59A","ZNF396","TCEB3C","PARD6G","PPAP2C","DNMT1",
          "CHMP2A","TSHZ3","CHST8","ZNF225","ZNF264","MIMT1","PEG3","ZIM2","MZF1","LILRB4","ZNF229",
          "NLRP2","MIR371A","USP29","ZIM3","PEG3-AS1","C20orf82","ISM1","PSIMCT-1","BLCAP","NNAT",
          "MCTS2","HM13","GDAP1L1","SGK2","COL9A3","GNAS","L3MBTL1","GNASAS","SANG","MIR296","MIR298",
          "C20orf20","SIM2","DSCAM","DGCR6L","DGCR6","FLJ20464","TSIX","XIST" ]

# Filter genes present in adata
existing_genes = [g for g in genes if g in adata.var_names]

# Compute average per cluster (specific_annotation)
cluster_avg = adata.to_df()[existing_genes].groupby(adata.obs['specific_annotation']).mean().T

# Plot clustered heatmap of averages
with PdfPages("imprinted_genes_cluster_avg.pdf") as pdf:
    g = sns.clustermap(
        cluster_avg,
        cmap='viridis',
        standard_scale=0,
        figsize=(25, 40),
        row_cluster=True,
        col_cluster=True,
        xticklabels=True,
        yticklabels=True
    )
    pdf.savefig(g.fig)
    plt.close(g.fig)


# In[ ]:





# In[31]:


import scanpy as sc

sc.pp.log1p(adata) 
# Run differential expression analysis
sc.tl.rank_genes_groups(adata, groupby='leiden', method='t-test')  # or 'wilcoxon' for a more robust test


# In[32]:


import pandas as pd

# Extract ranking results
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names

# Number of top genes you want
n_genes = 50  # or 20, etc.

# Create a dictionary of top genes per cluster
top_genes = {group: result['names'][group][:n_genes] for group in groups}

# Convert to DataFrame
df_wide = pd.DataFrame(top_genes)

# Optionally add rank as index
df_wide.index.name = 'rank'
df_wide.index += 1

# Save to CSV
df_wide.to_csv('top_genes_per_cluster55.csv')


# In[ ]:





# In[ ]:





# In[ ]:





# In[7]:


adata.obs.to_csv('all_samples_55clusters.csv', index=True)


# In[5]:


adata = scp.read_h5ad('../DATA/all_samples.h5ad')


# In[11]:


import pandas as pd

# Load the CSV file (comma-separated)
csv_filename = '../face_genes_df_hg38_with_windows_v2.csv'
gene_df = pd.read_csv(csv_filename)

# Extract gene list from SYMBOL column
gene_list = gene_df['SYMBOL'].dropna().unique().tolist()

# Extract genes from your adata object
adata_genes = adata.var_names.tolist()

# Find genes in adata not in the CSV SYMBOL list
genes_not_in_csv = [gene for gene in adata_genes if gene not in gene_list]

print(f"Number of genes in adata not found in CSV SYMBOL list: {len(genes_not_in_csv)}")

# Save the output to a text file
output_filename = 'genes_not_in_csv.txt'
with open(output_filename, 'w') as f:
    for gene in genes_not_in_csv:
        f.write(gene + '\n')

print(f"List of genes not in CSV saved to '{output_filename}'")


# In[13]:


adata


# In[15]:


adata = scp.read_h5ad('all_samples_55clusters.h5ad')

# Save the index and a specific column to a new CSV
output_filename = 'cellIDs_with_clusters.csv'
adata.obs[['leiden']].to_csv(output_filename, index=True)

print(f"Saved index and 'SYMBOL' column to '{output_filename}'")


# In[ ]:





# In[8]:


data = adata.copy()


# In[11]:


scp.tl.leiden(data, resolution=2)
scp.tl.umap(data, n_components=2)
scp.pl.umap(data, color=['leiden'])
#data.write('all_samples_clusters.h5ad')


# In[24]:


data.write('all_samples_55clusters.h5ad')


# In[52]:


adata = scp.read_h5ad('all_samples_55clusters.h5ad')


# In[18]:


#adata.obs['general_annotation'] = "Mesenchyme"
for cell in adata.obs.loc[adata.obs.leiden == '0'].index:
    adata.obs.at[cell, 'general_annotation'] = "Myoblast"
    adata.obs.at[cell, 'specific_annotation'] = "Myoblast1"
for cell in adata.obs.loc[adata.obs.leiden == '1'].index:
    adata.obs.at[cell, 'general_annotation'] = "Osteoblast"
    adata.obs.at[cell, 'specific_annotation'] = "Osteoblast"
for cell in adata.obs.loc[adata.obs.leiden == '2'].index:
    adata.obs.at[cell, 'general_annotation'] = "Dermal"
    adata.obs.at[cell, 'specific_annotation'] = "CXCL14+ dermal mesenchme"
for cell in adata.obs.loc[adata.obs.leiden == '3'].index:
    adata.obs.at[cell, 'general_annotation'] = "Periocular"
    adata.obs.at[cell, 'specific_annotation'] = "FOXP+ extraocular mesenchyme"  
for cell in adata.obs.loc[adata.obs.leiden == '4'].index:
    adata.obs.at[cell, 'general_annotation'] = "Osteogenic"
    adata.obs.at[cell, 'specific_annotation'] = "PTN+ w9 ear skeletal fibroblasts"
for cell in adata.obs.loc[adata.obs.leiden == '5'].index:
    adata.obs.at[cell, 'general_annotation'] = "Gland"
    adata.obs.at[cell, 'specific_annotation'] = "salivary or submandibular gland mesenchyme"
for cell in adata.obs.loc[adata.obs.leiden == '6'].index:
    adata.obs.at[cell, 'general_annotation'] = "lingual"
    adata.obs.at[cell, 'specific_annotation'] = "MEIS+ lingual mesenchyme"
for cell in adata.obs.loc[adata.obs.leiden == '7'].index:
    adata.obs.at[cell, 'general_annotation'] = "Progenitor"
    adata.obs.at[cell, 'specific_annotation'] = "Progenitors"
for cell in adata.obs.loc[adata.obs.leiden == '8'].index:
    adata.obs.at[cell, 'general_annotation'] = "Stroma"
    adata.obs.at[cell, 'specific_annotation'] = "PAX3+ nasal dermal mesenchyme"
for cell in adata.obs.loc[adata.obs.leiden == '9'].index:
    adata.obs.at[cell, 'general_annotation'] = "Stroma"
    adata.obs.at[cell, 'specific_annotation'] = "PAX9/PAX7+ nasal perichondrium"
for cell in adata.obs.loc[adata.obs.leiden == '10'].index:
    adata.obs.at[cell, 'general_annotation'] = "Chondrocyte"
    adata.obs.at[cell, 'specific_annotation'] = "Chondrocyte 1"
for cell in adata.obs.loc[adata.obs.leiden == '11'].index:
    adata.obs.at[cell, 'general_annotation'] = "Perivascular"
    adata.obs.at[cell, 'specific_annotation'] = "Pericyte"
for cell in adata.obs.loc[adata.obs.leiden == '12'].index:
    adata.obs.at[cell, 'general_annotation'] = "Dental"
    adata.obs.at[cell, 'specific_annotation'] = "Dental 1"
for cell in adata.obs.loc[adata.obs.leiden == '13'].index:
    adata.obs.at[cell, 'general_annotation'] = "Stroma"
    adata.obs.at[cell, 'specific_annotation'] = "DLK1+ CXCL12+ stroma"
for cell in adata.obs.loc[adata.obs.leiden == '14'].index:
    adata.obs.at[cell, 'general_annotation'] = "Chondrocyte"
    adata.obs.at[cell, 'specific_annotation'] = "Chondrocyte 2"
for cell in adata.obs.loc[adata.obs.leiden == '15'].index:
    adata.obs.at[cell, 'general_annotation'] = "Connective"
    adata.obs.at[cell, 'specific_annotation'] = "connective progenitor mesenchyme"
for cell in adata.obs.loc[adata.obs.leiden == '16'].index:
    adata.obs.at[cell, 'general_annotation'] = "Mesenchyme"
    adata.obs.at[cell, 'specific_annotation'] = "SOX11+ progenitor mesenchyme"    
for cell in adata.obs.loc[adata.obs.leiden == '17'].index:
    adata.obs.at[cell, 'general_annotation'] = "Fibroblast"
    adata.obs.at[cell, 'specific_annotation'] = "Fibroblast 1"
for cell in adata.obs.loc[adata.obs.leiden == '18'].index:
    adata.obs.at[cell, 'general_annotation'] = "Mesenchyme"
    adata.obs.at[cell, 'specific_annotation'] = "posterior mandibular mesenchyme"
for cell in adata.obs.loc[adata.obs.leiden == '19'].index:
    adata.obs.at[cell, 'general_annotation'] = "Mesenchyme"
    adata.obs.at[cell, 'specific_annotation'] = "Dermal progenitor mesenchyme"
for cell in adata.obs.loc[adata.obs.leiden == '20'].index:
    adata.obs.at[cell, 'general_annotation'] = "Dermal"
    adata.obs.at[cell, 'specific_annotation'] = "Mature anterior dermal fibroblast"
for cell in adata.obs.loc[adata.obs.leiden == '21'].index:
    adata.obs.at[cell, 'general_annotation'] = "Fibroblast"
    adata.obs.at[cell, 'specific_annotation'] = "Fibroblast"
for cell in adata.obs.loc[adata.obs.leiden == '21'].index:
    adata.obs.at[cell, 'general_annotation'] = "Neural"
    adata.obs.at[cell, 'specific_annotation'] = "Neural progenitors"
for cell in adata.obs.loc[adata.obs.leiden == '22'].index:
    adata.obs.at[cell, 'general_annotation'] = "Glia"
    adata.obs.at[cell, 'specific_annotation'] = "Glia"
for cell in adata.obs.loc[adata.obs.leiden == '23'].index:
    adata.obs.at[cell, 'general_annotation'] = "Dermal"
    adata.obs.at[cell, 'specific_annotation'] = "Dermal fibroblast"
for cell in adata.obs.loc[adata.obs.leiden == '24'].index:
    adata.obs.at[cell, 'general_annotation'] = "Meninges"
    adata.obs.at[cell, 'specific_annotation'] = "Leptomeninges"
for cell in adata.obs.loc[adata.obs.leiden == '25'].index:
    adata.obs.at[cell, 'general_annotation'] = "Mesenchyme"
    adata.obs.at[cell, 'specific_annotation'] = "POSTN+ late mesenchyme"
for cell in adata.obs.loc[adata.obs.leiden == '26'].index:
    adata.obs.at[cell, 'general_annotation'] = "Connective"
    adata.obs.at[cell, 'specific_annotation'] = "Tenocyte"
for cell in adata.obs.loc[adata.obs.leiden == '27'].index:
    adata.obs.at[cell, 'general_annotation'] = "Meninges-like"
    adata.obs.at[cell, 'specific_annotation'] = "Eye meningeal-like"
for cell in adata.obs.loc[adata.obs.leiden == '28'].index:
    adata.obs.at[cell, 'general_annotation'] = "Dermal"
    adata.obs.at[cell, 'specific_annotation'] = "Stress dermis"
for cell in adata.obs.loc[adata.obs.leiden == '29'].index:
    adata.obs.at[cell, 'general_annotation'] = "Dermal"
    adata.obs.at[cell, 'specific_annotation'] = "Dermal mesenchyme"
for cell in adata.obs.loc[adata.obs.leiden == '30'].index:
    adata.obs.at[cell, 'general_annotation'] = "Perivascular"
    adata.obs.at[cell, 'specific_annotation'] = "Pre-perivascular"
for cell in adata.obs.loc[adata.obs.leiden == '31'].index:
    adata.obs.at[cell, 'general_annotation'] = "Endothelia"
    adata.obs.at[cell, 'specific_annotation'] = "Endothelia"
for cell in adata.obs.loc[adata.obs.leiden == '32'].index:
    adata.obs.at[cell, 'general_annotation'] = "Myoblast"
    adata.obs.at[cell, 'specific_annotation'] = "Myoblast2"
for cell in adata.obs.loc[adata.obs.leiden == '33'].index:
    adata.obs.at[cell, 'general_annotation'] = "Mesenchyme"
    adata.obs.at[cell, 'specific_annotation'] = "Eye dermal mesenchyme"
for cell in adata.obs.loc[adata.obs.leiden == '34'].index:
    adata.obs.at[cell, 'general_annotation'] = "Meninges"
    adata.obs.at[cell, 'specific_annotation'] = "Outer dura"
for cell in adata.obs.loc[adata.obs.leiden == '35'].index:
    adata.obs.at[cell, 'general_annotation'] = "Subdermal"
    adata.obs.at[cell, 'specific_annotation'] = "Subdermal"
for cell in adata.obs.loc[adata.obs.leiden == '36'].index:
    adata.obs.at[cell, 'general_annotation'] = "Skeletal"
    adata.obs.at[cell, 'specific_annotation'] = "Skeletal progenitor mesenchyme"
for cell in adata.obs.loc[adata.obs.leiden == '37'].index:
    adata.obs.at[cell, 'general_annotation'] = "Immune"
    adata.obs.at[cell, 'specific_annotation'] = "Myeloid"
for cell in adata.obs.loc[adata.obs.leiden == '38'].index:
    adata.obs.at[cell, 'general_annotation'] = "Blood"
    adata.obs.at[cell, 'specific_annotation'] = "Erythrocytes"
for cell in adata.obs.loc[adata.obs.leiden == '39'].index:
    adata.obs.at[cell, 'general_annotation'] = "Preskeletal"
    adata.obs.at[cell, 'specific_annotation'] = "Chondroblast"
for cell in adata.obs.loc[adata.obs.leiden == '40'].index:
    adata.obs.at[cell, 'general_annotation'] = "Neural"
    adata.obs.at[cell, 'specific_annotation'] = "Neurons 1" 
for cell in adata.obs.loc[adata.obs.leiden == '41'].index:
    adata.obs.at[cell, 'general_annotation'] = "Myoblast"
    adata.obs.at[cell, 'specific_annotation'] = "Myoblast3"
for cell in adata.obs.loc[adata.obs.leiden == '42'].index:
    adata.obs.at[cell, 'general_annotation'] = "Neural"
    adata.obs.at[cell, 'specific_annotation'] = "Neurons 2"
for cell in adata.obs.loc[adata.obs.leiden == '43'].index:
    adata.obs.at[cell, 'general_annotation'] = "Stroma"
    adata.obs.at[cell, 'specific_annotation'] = "WNT5A+ stroma"
for cell in adata.obs.loc[adata.obs.leiden == '44'].index:
    adata.obs.at[cell, 'general_annotation'] = "Neural"
    adata.obs.at[cell, 'specific_annotation'] = "Neurons 3"
for cell in adata.obs.loc[adata.obs.leiden == '45'].index:
    adata.obs.at[cell, 'general_annotation'] = "Neural"
    adata.obs.at[cell, 'specific_annotation'] = "Peripheral neurons"
for cell in adata.obs.loc[adata.obs.leiden == '46'].index:
    adata.obs.at[cell, 'general_annotation'] = "Epithelial"
    adata.obs.at[cell, 'specific_annotation'] = "Epithelia"
for cell in adata.obs.loc[adata.obs.leiden == '47'].index:
    adata.obs.at[cell, 'general_annotation'] = "Progenitors"
    adata.obs.at[cell, 'specific_annotation'] = "Anterior mesenchymal progenitor"
for cell in adata.obs.loc[adata.obs.leiden == '48'].index:
    adata.obs.at[cell, 'general_annotation'] = "Neural"
    adata.obs.at[cell, 'specific_annotation'] = "Neural progenitor"
for cell in adata.obs.loc[adata.obs.leiden == '49'].index:
    adata.obs.at[cell, 'general_annotation'] = "Chondrocyte"
    adata.obs.at[cell, 'specific_annotation'] = "Nasal chondrocyte"
for cell in adata.obs.loc[adata.obs.leiden == '50'].index:
    adata.obs.at[cell, 'general_annotation'] = "Cornea"
    adata.obs.at[cell, 'specific_annotation'] = "Keratocytes"
for cell in adata.obs.loc[adata.obs.leiden == '51'].index:
    adata.obs.at[cell, 'general_annotation'] = "Immune"
    adata.obs.at[cell, 'specific_annotation'] = "T cells" 
for cell in adata.obs.loc[adata.obs.leiden == '52'].index:
    adata.obs.at[cell, 'general_annotation'] = "Eye"
    adata.obs.at[cell, 'specific_annotation'] = "Eye neural mesenchymal doublet" 
for cell in adata.obs.loc[adata.obs.leiden == '53'].index:
    adata.obs.at[cell, 'general_annotation'] = "Blood"
    adata.obs.at[cell, 'specific_annotation'] = "Blood" 
for cell in adata.obs.loc[adata.obs.leiden == '54'].index:
    adata.obs.at[cell, 'general_annotation'] = "Immune"
    adata.obs.at[cell, 'specific_annotation'] = "B cells" 
for cell in adata.obs.loc[adata.obs.leiden == '55'].index:
    adata.obs.at[cell, 'general_annotation'] = "Melanocytes"
    adata.obs.at[cell, 'specific_annotation'] = "Melanocytes"     


# In[9]:


adata.obs.to_csv('all_samples_55clusters.csv', index=True)


# In[ ]:





# In[ ]:





# In[ ]:





# In[19]:


# DotPlot
#adata = scp.read_h5ad('all_samples_55clusters.h5ad')

scp.set_figure_params(dpi_save=300, dpi=300)
#scp.pl.umap(adata, color=['general_annotation'],use_raw=False, size=5)

markers = ["ALAS2","HBB","NFKBIA","PECAM1","EPCAM","MEIS1","SOX11","TUBB3","CKB","GRID2","SOX10","S100B","MYF5","TNNT1"
           ,"COL2A1","COL9A1","SOX9","FMOD","SP7","PTN","RUNX2","EYA1","MKX","PRRX1","ACTA2","SULT1E1",
           "CXCL14","TWIST2","BGN","GJA1","EMX2","WNT5A","MSX1","TBX3",
           "MSX1","ALPL","MEIS2","MIR503HG","FOXL2", "TBX18","CXCL12","POSTN","SIX1","COL8A2","MGP","TNMD","SCX","PAX9","EBF2","COL5A1",
           "DKK2","KERA","S100A10", "MITF","PTGDS", "PAX1"
]


#scp.pl.dotplot(adata, markers, groupby = 'annotation', dendrogram = False, swap_axes = False, standard_scale = "var", save = '_55.pdf', categories_order = ["Blood","Immune","Endothelia","Epithelia","Neuronal","Glia","Myoblast","Chondrogenic","Osteogenic","Perivascular","Dermal","Subdermal","Meningeal","Stroma","Fibroblast","Progenitor","Mesenchyme"])
scp.pl.dotplot(adata, markers, groupby = 'specific_annotation', dendrogram = False, swap_axes = False, standard_scale = "var", save = '_55test.pdf', 

categories_order = [
    # Blood
    "Erythrocytes",
    "Blood",
    "Myeloid",
    "T cells",
    "B cells",

    # Endothelia
    "Endothelia",

    # Epithelia
    "Epithelia",

    # Neuronal
    "Neural progenitors",
    "Neurons 1",
    "Neurons 2",
    "Neurons 3",
    "Peripheral neurons",
    "Neural progenitor",

    # Glia
    "Glia",

    # Myoblast
    "Myoblast1",
    "Myoblast2",
    "Myoblast3",

    # Chondrogenic
    "Chondrocyte 1",
    "Chondrocyte 2",
    "Nasal chondrocyte",
    "Chondroblast",

    # Osteogenic
    "Osteoblast",
    "PTN+ w9 ear skeletal fibroblasts",
    "Skeletal progenitor mesenchyme",

    # Perivascular
    "Pericyte",
    "Pre-perivascular",

    # Dermal
    "CXCL14+ dermal mesenchme",
    "Mature anterior dermal fibroblast",
    "Dermal fibroblast",
    "Stress dermis",
    "Dermal mesenchyme",
    "FOXP+ extraocular mesenchyme",  # Periocular

    # Subdermal
    "Subdermal",

    # Meningeal
    "Leptomeninges",
    "Outer dura",
    "Eye meningeal-like",

    # Stroma
    "PAX3+ nasal dermal mesenchyme",
    "PAX9/PAX7+ nasal perichondrium",
    "DLK1+ CXCL12+ stroma",
    "WNT5A+ stroma",

    # Fibroblast
    "Fibroblast 1",

    # Progenitor
    "Progenitors",
    "Anterior mesenchymal progenitor",
    "SOX11+ progenitor mesenchyme",
    "Dermal progenitor mesenchyme",
    "MEIS+ lingual mesenchyme",  # Lingual

    # Mesenchyme
    "posterior mandibular mesenchyme",
    "POSTN+ late mesenchyme",
    "Eye dermal mesenchyme",
    "connective progenitor mesenchyme",
    "Tenocyte",
    "Dental 1",
    "salivary or submandibular gland mesenchyme",
    "Keratocytes",
    "Eye neural mesenchymal doublet",
    "Melanocytes"
]

)


#scp.pl.dotplot(adata, markers, groupby = 'general_annotation', dendrogram = False, swap_axes = False, standard_scale="var", categories_order = ["Blood","Immune","Endothelia","Epithelia","Neuronal","Glia","Myoblast","Chondrogenic","Osteogenic","Perivascular","Dermal","Subdermal","Meningeal","Stroma","Fibroblast","Progenitor","Mesenchyme"])

#scp.pl.umap(adata, color=['leiden'])
#adata.write('DATA/all_samples.h5ad')


# In[68]:


import pandas as pd
import scanpy as sc

# 1️⃣ Load adata
adata = sc.read_h5ad('all_samples_55clusters.h5ad')

# 2️⃣ Load scores CSV
df_scores = pd.read_csv("corr_Universe_x_genes_ROZYsnpVARIATIONreal_flat.csv")

# 3️⃣ Clean cluster column
df_scores = df_scores.dropna(subset=['cluster'])
df_scores['cluster'] = df_scores['cluster'].astype(float).astype(int).astype(str)

# 4️⃣ Ensure leiden is string
adata.obs['leiden'] = adata.obs['leiden'].astype(str)

# 5️⃣ Create mappings for r2 and cor
cluster_to_r2 = df_scores.set_index('cluster')['r2'].to_dict()
cluster_to_cor = df_scores.set_index('cluster')['cor'].to_dict()

# 6️⃣ Map values into adata.obs
adata.obs['r2'] = adata.obs['leiden'].map(cluster_to_r2)
adata.obs['cor'] = adata.obs['leiden'].map(cluster_to_cor)

# 7️⃣ Plot UMAPs
sc.pl.umap(
    adata,
    color='r2',
    cmap='viridis',
    size=50,
    title='UMAP colored by Universe R2 score',
    save = ' Universe R2 score.pdf'
)

sc.pl.umap(
    adata,
    color='cor',
    cmap='coolwarm',  # diverging cmap makes sense for correlations
    size=50,
    title='UMAP colored by Universe Correlation',
    save = ' Universe Correlation score.pdf'
)


# 2️⃣ Load scores CSV
df_scores = pd.read_csv("corr_spatial.csv")

# 3️⃣ Clean cluster column
df_scores = df_scores.dropna(subset=['cluster'])
df_scores['cluster'] = df_scores['cluster'].astype(float).astype(int).astype(str)

# 4️⃣ Ensure leiden is string
adata.obs['leiden'] = adata.obs['leiden'].astype(str)

# 5️⃣ Create mappings for r2 and cor
cluster_to_r2 = df_scores.set_index('cluster')['r2'].to_dict()
cluster_to_cor = df_scores.set_index('cluster')['cor'].to_dict()

# 6️⃣ Map values into adata.obs
adata.obs['r2'] = adata.obs['leiden'].map(cluster_to_r2)
adata.obs['cor'] = adata.obs['leiden'].map(cluster_to_cor)

# 7️⃣ Plot UMAPs
sc.pl.umap(
    adata,
    color='r2',
    cmap='viridis',
    size=50,
    title='UMAP colored by Spatial R2 score',
    save = ' Spatial R2 score.pdf'
)

sc.pl.umap(
    adata,
    color='cor',
    cmap='coolwarm',  # diverging cmap makes sense for correlations
    size=50,
    title='UMAP colored by Spatial Correlation',
    save = ' Spatial Correlation score.pdf'
)


# In[ ]:





# In[ ]:





# In[ ]:





# In[5]:


adata1 = scp.read_h5ad('../DATA/all-mes-harmony-sample.h5ad') 


# In[6]:


import scanpy as sc

sc.pp.log1p(adata1) 
# Run differential expression analysis
sc.tl.rank_genes_groups(adata1, groupby='leiden', method='t-test')  # or 'wilcoxon' for a more robust test


# In[7]:


import pandas as pd

# Extract ranking results
result = adata1.uns['rank_genes_groups']
groups = result['names'].dtype.names

# Number of top genes you want
n_genes = 50  # or 20, etc.

# Create a dictionary of top genes per cluster
top_genes = {group: result['names'][group][:n_genes] for group in groups}

# Convert to DataFrame
df_wide = pd.DataFrame(top_genes)

# Optionally add rank as index
df_wide.index.name = 'rank'
df_wide.index += 1

# Save to CSV
df_wide.to_csv('top_genes_mes.csv')


# In[ ]:




