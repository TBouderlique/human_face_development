
# Human embryonic face atlas 

Code related to Cell types, positional codes, and enhancers contributing to human facial individuality and pathology (Alek Erickson, Yakov Gerstein, and Rozalina Galimullina et al. 2026). This repository contains supporting code, Jupyter notebooks, and yml/conda/env files to analyse the embryonic human face data.


## Find the data

Explore the dataset interactively or download the H5 object directly using our CELLXGENE: https://adameykolab.hifo.meduniwien.ac.at/cellxgene_public/ 

Count matrix data has been submitted to Dryad, and is accessible via this link: https://doi.org/10.5061/dryad.7m0cfxq8t 

All raw FASTQ files have been deposited to the European Genome-Phenome Archive (EGA) and are available upon request from the following accession number EGAS50000002001 (https://ega-archive.org/studies/EGAS50000002001).


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
Code for analyzing the data is available as Jupyter notebooks.

### HiC_analysis
Notebooks used for calling TADs and classifying links according to HiC data
### RNA
Notebooks used for QC and processing of each scRNAseq sample, as well as trajectory analysis.
### multiome
Notebooks used for analysis of multiome data including links, DARs, and motif analysis.
### plots
Notebooks used to plot GWAS enrichment heatmaps and GRN schematics.  
### spatial
Notebooks used for analyzing Stereoseq and Visium spatial transcriptomics datasets
### table
Notebook used to construct link-peak table shown in Supplementary Table 2



