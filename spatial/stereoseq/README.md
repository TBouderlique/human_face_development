# Stereoseq workflow for spatial analysis.
These folders are include the code and files used in the stereoseq analysis.

## 00 Reading the file.
Reads and saves the raw files into objects used later.

## 01 Preprocessing
Preprocesses the files some more and saves into bins.

## 02 Cell UMAP
Saves data on cell level.

## 03 Celltypist label transfer
Transfers the annotations from the full-RNA-seq data onto stereo-seq data.

## 04 Commot signal overlaid on celltypes
This requires that commot has been run before. The pipeline for commot is in another folder. Commot requires that you have run 01_preprocessing. Afterwards this folder can be used to plot commot signals on celltype plots.

## 05 Dotplot commot degs all
This will create a dotplot of the top global gene activation from pathways was retrieved from the COMMOT results.

## 06 Morans I Stereoseq
This is used to find spatial genes. Use st env. 