# Uses celltypist conda env

# Imports
import scanpy as sc
import celltypist
import time
import numpy as np
import sys

# Paths
DATA_PATH="/home/felix/projects/facial/felix/data"
FELIX_DATA_PATH="/home/felix/projects/facial/felix/data/reprocessed_data"

# Load data
ref_adata = sc.read_h5ad(f"{DATA_PATH}/all-harmony.h5ad") # Integrated full dataset from RNA-seq analysis. 

max_iter_1=int(sys.argv[1]) # Original 5
max_iter_2=int(sys.argv[2]) # Original 100
n_genes=int(sys.argv[3]) # Original 100
save_suffix=sys.argv[4]

"""
# Annotate the reference -- old-harmony
annotations={
    "0":"mesenchyme A(near-dermal)",
    "1":"osteogenic",
    "2":"dermal",
    "3":"mesenchyme B (middle)",
    "4":"mesenchyme C (ebf)",
    "5":"chondrocytes",
    "6":"mesenchyme D (young)",
    "7":"myoblasts",
    "8":"mesenchyme E (fox)",
    "9":"mesenchyme F (middle)",
    "10":"perivascular",
    "11":"neural progenitors",
    "12":"glia",
    "13":"neruons",
    "14":"endothelia",
    "15":"immune",
    "16":"erythrocytes",
    "17":"myoblasts",
    "18":"mesengyme G (wierd)",
    "19":"interneurons",
    "20":"epithelia",
    "21":"photoreceptors",
    "22":"immune",
    "23":"melanocytes"
}
new_name = "annotation"
ref_adata.obs[new_name] = ""
for take in annotations:
    cluster = ref_adata[ref_adata.obs["leiden"] == take]
    ref_adata.obs.loc[cluster.obs_names, "annotation"] = annotations[take]
"""
# Annotate the all-harmony-55
annotations={
    "0":"Myoblast",
    "1":"Osteoblast",
    "2":"MC",
    "3":"MC",
    "4":"Fibroblast",
    "5":"MC",
    "6":"MC",
    "7":"MC progenitors",
    "8":"MC",
    "9":"MC",
    "10":"Chondrocyte",
    "11":"Pericyte",
    "12":"Dental",
    "13":"Stroma",
    "14":"Chondrocyte",
    "15":"MC progenitors",
    "16":"MC progenitors",
    "17":"Fibroblast",
    "18":"MC",
    "19":"MC progenitors",
    "20":"Fibroblast",
    "21":"Neural progenitors",
    "22":"Glia",
    "23":"Fibroblast",
    "24":"Leptomeninges",
    "25":"MC" ,
    "26":"Tenocyte" ,
    "27":"Eye meningeal-like" ,
    "28":"Stess dermis" ,
    "29":"MC" ,
    "30":"Pre-perivascular" ,
    "31":"Endothelia" ,
    "32":"Myoblast" ,
    "33":"MC" ,
    "34":"Perichondrial" ,
    "35":"Subdermal" ,
    "36":"MC progenitors" ,
    "37":"Myeloid" ,
    "38":"Erythrocytes" ,
    "39":"Chondroblast" ,
    "40":"Neurons" ,
    "41":"Myoblast" ,
    "42":"Neurons" ,
    "43":"Strom" ,
    "44":"Neurons" ,
    "45":"Peripheral neurons" ,
    "46":"Epithelia" ,
    "47":"MC progenitors" ,
    "48":"Neural progenitors" ,
    "49":"Chondrocyte" ,
    "50":"Keratocytes" ,
    "51":"T cells" ,
    "52":"Eye neural MC doublet" ,
    "53":"Blood" ,
    "54":"B cells" ,
    "55":"Melanocytes" ,
}
new_name = "annotation"
reference.obs[new_name] = ""
for take in annotations:
    cluster = reference[reference.obs["leiden"] == take]
    reference.obs.loc[cluster.obs_names, "annotation"] = annotations[take]
    
# Downsampling cells from each cell type to a given number
# Given number here: 500
# All cells from a given cell type will be selected if the cell type size is < 500.
sampled_cell_index = celltypist.samples.downsample_adata(ref_adata, mode = 'each', n_cells = 5000, by = 'annotation', return_index = True)


# Use `celltypist.train` to quickly train a rough CellTypist model.
# You can also set `mini_batch = True` to enable mini-batch training.
max_iter=max_iter_1 # Can improve model, original was 5
t_start = time.time()
model_fs = celltypist.train(ref_adata[sampled_cell_index], 'annotation', n_jobs = 10, max_iter = max_iter, use_SGD = True)
t_end = time.time()
print(f"Time elapsed: {t_end - t_start} seconds")
    
nr_genes=n_genes # Can be changed, original 100
gene_index = np.argpartition(np.abs(model_fs.classifier.coef_), -1*nr_genes, axis = 1)[:, -1*nr_genes:]
gene_index = np.unique(gene_index)
print(f"Number of genes selected: {len(gene_index)}")

# Add `check_expression = False` to bypass expression check with only a subset of genes.
max_iterations=max_iter_2 # Can improve model, original 100
t_start = time.time()
model = celltypist.train(ref_adata[sampled_cell_index, gene_index], 'annotation', check_expression = False, n_jobs = 10, max_iter = max_iterations)
t_end = time.time()
print(f"Time elapsed: {(t_end - t_start)/60} minutes")

# Save the model.
model.write(f"{FELIX_DATA_PATH}/models/model_from_all_harmony_{str(save_suffix)}.pkl")   
    
    


