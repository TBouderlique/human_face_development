import commot as ct
import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
import matplotlib.pyplot as plt
import time
import sys
import os
os.environ['R_HOME'] = '/home/felix/miniconda3/envs/commot/lib/R'
os.environ['R_USER'] = '/home/felix/miniconda3/envs/commot/lib/python3.7/site-packages/rpy2'

plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["font.family"] = ["sans-serif"]
plt.rcParams["font.sans-serif"] = ["Arial"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["savefig.bbox"]="tight"
plt.rcParams["savefig.format"]='svg'


sc.set_figure_params(color_map = 'bwr',vector_friendly=True, format='svg')
felix_data_path="/home/felix/projects/facial/felix/data/"

take = sys.argv[1]
datasets = ["070623_A1_eye_socket_visium",
    "070623_B1_sagital-mandible",
    "070623_D1_W8_frontal_mandible_visium"]
sample = datasets[int(take)]

save_dir = f"/home/felix/projects/facial/felix/plots/{sample}_secreted_signaling" # Switch here

adata_dis500 = sc.read_h5ad(f"{felix_data_path}{sample}_adata_ccc_dis500_secreted_siganling.h5ad") # Switch here

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
#fig, ax = plt.subplots(figsize=(ncol*5, nrows*5))
plt.figure(figsize=(30, nrows*4.5))
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

    ax.legend_ = None
    
plt.suptitle(fontsize=20, t=f"Secreted Signaling receptor communication SENDER plots for the {sample}", y=.98)    # Switch here  
plt.savefig(f"{save_dir}/{sample}_COMMOT_spatial_sender.svg",dpi=150)

# Print plots all receiver together
import math
ncol=6
paths=sorted(list(set(ad_receiver.var_names)))
nrows=math.ceil(len(paths)/5)
#fig, ax = plt.subplots(figsize=(ncol*5, nrows*5))
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
    
    ax.legend_ = None

plt.suptitle(fontsize=20, t=f"Secrted Signaling recptor communication RECEIVER plots for the {sample}", y=.98)    # Switch here
plt.savefig(f"{save_dir}/{sample}_COMMOT_spatial_receiver.pdf",dpi=150)

