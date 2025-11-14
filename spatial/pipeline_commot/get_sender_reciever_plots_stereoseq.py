# Running COMMOT on stereo-seq data
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from tqdm import tqdm
import squidpy as sq
from scipy.stats import spearmanr
import commot as ct
import time
import sys

sc.settings.verbosity = 1
sns.set(font_scale=1)
sc.settings.set_figure_params(dpi=150)
sns.set_style("ticks")

from matplotlib import cm
from matplotlib.colors import ListedColormap

cm_color = cm.get_cmap("Reds", 128)
cm_grey = cm.get_cmap("Greys", 128)

Reds = ListedColormap(np.vstack((
    cm_grey(np.linspace(0.2, 0.2, 1)),
    cm_color(np.linspace(0.1, 1, 128)),
)))

felix_data_path="/home/felix/projects/facial/felix/data/reprocessed_data/"
batch = str(sys.argv[1])

import os
save_dir = f"/home/felix/projects/facial/felix/plots/reprocessed/{batch}_secreted_signaling" # Switch here

adata_dis500 = sc.read_h5ad(f"{felix_data_path}{batch}_adata_ccc_dis500_secreted_signaling.h5ad") # Switch here
adata = sc.read_h5ad(f"{felix_data_path}{batch}/processed/50.h5ad")


if not os.path.exists(save_dir):
    os.makedirs(save_dir)
    
import anndata as ad

ad_sender=ad.AnnData(adata_dis500.obsm['commot-cellchat-sum-sender'].loc[:,[i for i in adata_dis500.obsm['commot-cellchat-sum-sender'].columns if i.count('-')==1]], 
           obsm=adata_dis500.obsm)
ad_sender.obs['sum_sender']=pd.DataFrame(ad_sender.X,columns=ad_sender.var_names, index=ad_sender.obs_names).sum(axis=1).to_list()

ad_receiver=ad.AnnData(adata_dis500.obsm['commot-cellchat-sum-receiver'].loc[:,[i for i in adata_dis500.obsm['commot-cellchat-sum-receiver'].columns if i.count('-')==1]], 
           obsm=adata_dis500.obsm)
ad_receiver.obs['sum_receiver']=pd.DataFrame(ad_receiver.X,columns=ad_receiver.var_names, index=ad_receiver.obs_names).sum(axis=1).to_list()

# Print plots all sender together
import math
ncol=6
paths=sorted(list(set(ad_sender.var_names)))
nrows=math.ceil(len(paths)/5)
plt.figure(figsize=(30, nrows*6.5))
plt.subplots_adjust(hspace=0.2, wspace=.2)

for n, p in enumerate(paths):
    
    ax = plt.subplot(nrows, ncol, n + 1)
    
    sc.pl.embedding(ad_sender,#frameon=False,
                  color=p, 
                  show=False,
                  ax=ax,  
                  legend_fontoutline=3,
                  basis='spatial',
                  size=80,#marker='H'
               )
    _=ax.set_facecolor('black')
    ax.set_title(p +' pathway',y=1.1, pad=0)   
   
    ax.set(xlabel=None, ylabel=None)
    
    ax.invert_yaxis()
    ax.legend_ = None
    
plt.suptitle(fontsize=20, t=f"Secreted signaling receptor communication SENDER plots for {batch}", y=.98) # Switch here   
plt.savefig(f"{save_dir}/{batch}_COMMOT_spatial_sender.pdf",dpi=150)


# Print plots all receiver together

import math
ncol=6
paths=sorted(list(set(ad_receiver.var_names)))
nrows=math.ceil(len(paths)/5)
plt.figure(figsize=(30, nrows*6.5))
plt.subplots_adjust(hspace=0.2, wspace=.2)
  

for n, p in enumerate(paths):
    
    ax = plt.subplot(nrows, ncol, n + 1)
    
    sc.pl.embedding(ad_receiver,#frameon=False,
                  color=p, 
                  show=False,
                  ax=ax,  
                  legend_fontoutline=3,
                  basis='spatial',
                  size=80,#marker='H'
               )
    _=ax.set_facecolor('black')
    ax.set_title(p +' pathway',y=1.1, pad=0)   
    
    ax.set(xlabel=None, ylabel=None)
    
    ax.invert_yaxis()
    ax.legend_ = None

plt.suptitle(fontsize=20, t=f"Secreted signaling recptor communication RECEIVER plots for {batch}", y=.98)   # Switch here 
plt.savefig(f"{save_dir}/{batch}_COMMOT_spatial_receiver.pdf",dpi=150)











