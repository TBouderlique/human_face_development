# Celltypist label transfer
# Use celltypist conda environvment
# After setting the correct model for the correct dataset below the script can be run.


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# This needs to be modified before running script, selkect your model here. 
dataset_to_model = {"C01939B2":"model_6", 
                    "A02989E1":"model_6", 
                    "A02993E1":"model_6", 
                    "A02994D6":"model_6", 
                    "A02994E6":"model_6", 
                    "C01939A4":"model_6", 
                    "C01939A5":"model_6", 
                    "C01939A6":"model_6" 
                   }

# Imports
import scanpy as sc
import celltypist
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Paths
DATA_PATH="/home/felix/projects/facial/felix/data"
FELIX_DATA_PATH="/home/felix/projects/facial/felix/data/reprocessed_data"

def predict_and_tranfer(query_adata, model_path):
    # CellTypist prediction without over-clustering and majority-voting.
    t_start = time.time()
    predictions = celltypist.annotate(query_adata, model = model_path)
    t_end = time.time()
    print(f"Time elapsed: {t_end - t_start} seconds")
    del query_adata.uns["neighbors"] # Removing neighborhood graph
    sc.pp.pca(query_adata) # Since X_pca did not work
    # CellTypist prediction with over-clustering and majority-voting.
    t_start = time.time()
    predictions = celltypist.annotate(query_adata, model = model_path, majority_voting = True)
    t_end = time.time()
    print(f"Time elapsed: {t_end - t_start} seconds")
    # Get an `AnnData` with predicted labels embedded into the cell metadata columns.
    adata = predictions.to_adata() #This uses predicted_labels as conf
    return adata

datasets = ["C01939B2", "A02989E1", "A02993E1",  "A02994D6",  "A02994E6", "C01939A4",  "C01939A5",  "C01939A6"]

model_dict={}
for x in range(1,11):
    model_dict[f"model_{str(x)}"] =f"{FELIX_DATA_PATH}/models/model_from_all_harmony_{str(x)}.pkl"
print(model_dict)


import os
if not os.path.exists(f"{FELIX_DATA_PATH}/annotated_stereo"):
    os.makedirs(f"{FELIX_DATA_PATH}/annotated_stereo")

for dataset in datasets:
    model_path=model_dict[dataset_to_model[dataset]]
    query_adata  = sc.read_h5ad(f"{FELIX_DATA_PATH}/{dataset}/processed/50.h5ad")
    adata = predict_and_tranfer(query_adata, model_path)
    adata.write_h5ad(f"{FELIX_DATA_PATH}/annotated_stereo/{dataset}_50.h5ad")