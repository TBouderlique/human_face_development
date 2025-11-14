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

adata = sc.read_h5ad(f"{felix_data_path}{batch}/processed/50.h5ad")

adata_dis500=adata.copy()

adata_dis500.obsm["spatial"] = adata_dis500.obsm["X_spatial"].copy()
del adata_dis500.obsm["X_spatial"]

df_cellchat =ct.pp.ligand_receptor_database(signaling_type="Secreted Signaling", # Select signaling type
                                            database="CellChat", # Select signaling db
                                            species='human')


df_cellchat_filtered = ct.pp.filter_lr_database(df_cellchat, adata_dis500, min_cell_pct=0.05)

    
t0=time.time()
ct.tl.spatial_communication(adata_dis500,
                            database_name='cellchat',
                            df_ligrec=df_cellchat_filtered,
                            dis_thr=500,
                            heteromeric=True,
                            pathway_sum=True)
t1=time.time()
print(f"Elapsed time {t1-t0} for {adata_dis500.obsm['spatial'].shape[0]} locations ")



adata_dis500.write(f"{felix_data_path}{batch}_adata_ccc_dis500_secreted_signaling.h5ad") # Switch here
