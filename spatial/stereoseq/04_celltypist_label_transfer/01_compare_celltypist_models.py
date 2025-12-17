# Compare celltypist models
# Uses celltypist conda env

# Imports
import scanpy as sc
import celltypist
import time
import numpy as np
import pandas as pd

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
results_dict={}
for model in model_dict:
    model_path=model_dict[model]
    model_res = {}
    for dataset in datasets:
        query_adata  = sc.read_h5ad(f"{FELIX_DATA_PATH}/{dataset}/processed/50.h5ad")
        adata = predict_and_tranfer(query_adata, model_path)
        tmp_dataset = {}
        for prec in [0.50, 0.60, 0.70, 0.80 , 0.90]:
            tmp_dataset[str(prec)] = [len(adata[adata.obs.conf_score>prec].obs), len(adata[adata.obs.conf_score>prec].obs)/len(adata.obs)]
        model_res[dataset] = tmp_dataset
    results_dict[model] = model_res
    
import os
if not os.path.exists("celltypist_model_comparison_results"):
    os.makedirs("celltypist_model_comparison_results")

import datetime
now = datetime.datetime.now()
date=str(now.date())

with open(f"celltypist_model_comparison_results/comparison_{date}.txt", "a") as f:

    for dataset in datasets:
        f.write(dataset + "\n")
        for prec in [0.50, 0.60, 0.70, 0.80 , 0.90]:
            f.write(f"Precentage = {prec}" + "\n")
            best = 0
            order=[]
            for model in model_dict:
                if best==0:
                    best=results_dict[model][dataset][str(prec)][1]
                    mod=model
                    continue
                if results_dict[model][dataset][str(prec)][1] > best:
                    order.append([best, mod])
                    best=results_dict[model][dataset][str(prec)][1]
                    mod=model
                else:
                    order.append([results_dict[model][dataset][str(prec)][1], model])

            order.append([best, mod])
            f.write(f"{mod} is best" + "\n")
            for val in reversed(order):
                f.write(f"{str(val[0])} {str(val[1])}" + "\n")
            f.write("-------------------------------" + "\n")
        f.write("=================================" + "\n")
        f.write("\n")
