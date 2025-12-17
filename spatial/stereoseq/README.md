# Stereoseq workflow for spatial analysis.
These folders are include the code and files used in the stereoseq analysis.

## 00 Reading the file.
Reads and saves the raw files into objects used later.

## 01 Preprocessing commot
Preprocesses the files some more and saves into bins for commot.

## 02 Preprocessing Morans I
Preprocesses the files some more and saves into bins for Morans I. 

## 03 Cell UMAP
Saves data on cell level.

## 04 Celltypist label transfer
Transfers the annotations from the full-RNA-seq data onto stereo-seq data.

## 05 Commot signal overlaid on celltypes
This requires that commot has been run before. The pipeline for commot is in another folder. Commot requires that you have run 01_preprocessing. Afterwards this folder can be used to plot commot signals on celltype plots.

## 06 Dotplot commot degs all
This will create a dotplot of the top global gene activation from pathways was retrieved from the COMMOT results.

## 07 Morans I Stereoseq
This is used to find spatial genes. Use st env. 