# Use commot env. 

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
data_path = "/home/felix/data/processed_face/visium/h5ad/"

take = sys.argv[1] # Which dataset to use.
s_type = sys.argv[2] # Which siganling type to use. We use "Secreted Signaling"
db = sys.argv[3] # Which database to use. We use CellChat



datasets = ["070623_A1_eye_socket_visium",
    "070623_B1_sagital-mandible",
    "070623_D1_W8_frontal_mandible_visium"]
sample = datasets[int(take)]

save_dir = f"/home/felix/projects/facial/felix/plots/"


adata = sc.read_h5ad(f"{data_path}{sample}.h5ad")

adata=adata.raw.to_adata()
adata.raw=adata
adata_dis500=adata.copy()
adata_dis500.obsm['spatial'][:,1]=adata.obsm['spatial'][:,0]
adata_dis500.obsm['spatial'][:,0]=adata.obsm['spatial'][:,1]

df_cellchat =ct.pp.ligand_receptor_database(signaling_type=s_type, # Select signaling type. Ex Secreted Signaling
                                            database=db, # Select signaling db. Ex CellChat
                                            species='human')

df_cellchat_filtered = ct.pp.filter_lr_database(df_cellchat, adata_dis500, min_cell_pct=0.05)

print("Starting commot")
t0=time.time()
ct.tl.spatial_communication(adata_dis500,
    database_name='cellchat', df_ligrec=df_cellchat_filtered,
                            dis_thr=500, heteromeric=True, pathway_sum=True)
t1=time.time()
print(f"Elapsed time {t1-t0} for {adata_dis500.obsm['spatial'].shape[0]} locations ")



adata_dis500.write(f"{felix_data_path}{sample}_adata_ccc_dis500_{'_'.join(s_type.split(' '))}.h5ad") 

