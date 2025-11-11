import commot as ct
import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
import matplotlib.pyplot as plt
import pickle
import os
import sys

os.environ['R_HOME'] = '/home/felix/miniconda3/envs/commot/lib/R'
os.environ['R_USER'] = '/home/felix/miniconda3/envs/commot/lib/python3.7/site-packages/rpy2'

batch = str(sys.argv[1])
pathway = str(sys.argv[2])
db=str(sys.argv[3])

# Set paths
felix_data_path="/home/felix/projects/facial/felix/data/reprocessed_data/"
save_dir = f"{felix_data_path}deg_{batch}_{db}" 

# Create the save directory if it does not exist already
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
    
    
adata_dis500 = sc.read_h5ad(f"{felix_data_path}{batch}_adata_ccc_dis500_{db}.h5ad") 

adata_dis500.layers['counts']=adata_dis500.X

import sys
sys.path.append("/home/felix/projects/facial/felix/pipeline_commot/")
from help_tools_commot import communication_deg_detection_modified

df_deg, df_yhat = communication_deg_detection_modified(adata_dis500,
    database_name = 'cellchat', pathway_name=pathway, summary = 'receiver', batch=batch, pathway=pathway)
    

deg_result = {"df_deg": df_deg, "df_yhat": df_yhat}
with open(f"{felix_data_path}deg_{batch}_{db}/deg_{pathway}.pkl", 'wb') as handle:
    pickle.dump(deg_result, handle, protocol=pickle.HIGHEST_PROTOCOL)







