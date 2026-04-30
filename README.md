
# Human embryonic face atlas 

Code related to Cell types, positional codes, and enhancers contributing to human facial individuality and pathology (Alek Erickson, Yakov Gerstein, and Rozalina Galimullina et al. 2026). This repository contains supporting code, Jupyter notebooks, and yml/conda/env files to analyse the embryonic human face data.


## Find the data

Explore the dataset interactively or download the H5 object directly using our CELLXGENE: https://adameykolab.hifo.meduniwien.ac.at/cellxgene_public/ 

Count matrix data has been submitted to Dryad, and is accessible via this link: https://doi.org/10.5061/dryad.7m0cfxq8t 

FASTQ files have been deposited to the Karolinska Institute Data Repository via the DORIS system from the Swedish National Data Service (SND): https://doi.org/10.48723/pq77-v833b 



## Description of values in the cell-metadata in the H5 file
 
### obs (cell metadata):
- n_genes / n_genes_by_counts: number of detected genes per cell
- total_counts: total UMI counts per cell
- total_counts_mt / pct_counts_mt: mitochondrial counts and percentage (QC metric)
- doublet_score / predicted_doublet: doublet detection results
- Sample: description of the age and region of the sample
- batch: sample-level flag used to integrate gene lists across sample batches
- Age: Age of the donor
- Origin: anatomical region of the tissue sample for each cell
- S_score / G2M_score / phase: cell cycle scoring and assigned phase
- leiden: cluster ID from Leiden clustering
- coarse_annotation / general_annotation / specific_annotation: hierarchical cell type labels

### var (gene metadata):
- n_cells: number of cells expressing each gene
- highly_variable: whether gene is selected as HVG
- means / dispersions / dispersions_norm: stats used for HVG selection

### uns (unstructured):
- annotation_colors / leiden_colors: color mappings
- hvg: parameters/results of HVG selection
- log1p: log-normalization info
- neighbors: kNN graph parameters
- pca / umap: embedding settings

### obsm (multi-dimensional embeddings):
- X_pca: PCA coordinates
- X_umap: UMAP coordinates

### obsp (pairwise matrices):
- connectivities: kNN graph (weighted adjacency)
- distances: pairwise neighbor distances



## Code 
Code for making many of the figures is available as Jupyter notebooks.
The package versions used to generate these figures are in this environment file







