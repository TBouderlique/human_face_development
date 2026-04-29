
# Human embryonic face atlas 

Code related to Cell types, positional codes, and enhancers contributing to human facial individuality and pathology (Alek Erickson, Yakov Gerstein, and Rozalina Galimullina et al. 2026). This repository contains supporting code, Jupyter notebooks, and yml/conda/env files to analyse the embryonic human face data.


## Find the data

Explore the dataset interactively or download the H5 object directly using our CELLXGENE: https://adameykolab.hifo.meduniwien.ac.at/cellxgene_public/ 

Count matrix data has been submitted to Dryad, and is accessible via this link: https://doi.org/10.5061/dryad.7m0cfxq8t 

FASTQ files have been deposited to the Karolinska Institute Data Repository via the DORIS system from the Swedish National Data Service (SND): https://doi.org/10.48723/pq77-v833b 



## Description of values in the cell-metadata in the H5 file
 
- Annotation_Coarse: Primary clusters derived from leiden clustering at resolution = 1
- Annotation_Specific: Primary clusters derived from leiden clustering at resolution = 2
- Donor: Donor IDs
- Doublet Score: Output from SCrublet
- Embedding: 2D embedding as used in the paper
- Age: Age of the donor
- Region: anatomical Origin of the tissue sample for each cell



## Code 
Code for making many of the figures is available as Jupyter notebooks.
The package versions used to generate these figures are in this environment file
