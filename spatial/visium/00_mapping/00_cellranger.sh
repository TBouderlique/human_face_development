#!/bin/bash

folders=$(ls -d *alltissue*)

for fold in $folders
    do
   d=$(echo ${fold})
    echo $d
    #############    Mapping   #########################
    if [[ ! -d ${d}/outs ]]; then
   
    cellranger count --id ${d}_transcriptome --transcriptome ~/data/ref_genome/refdata-gex-GRCh38-2020-A  --fastqs ${d}/  --localcores 44  --localmem 100 --include-introns true
    
    mv ${d}_transcriptome ${d}
    fi
    
    ############################### Do denoising   ########################
    echo After mapping
    if [[ ! -d ${d}/denoised_dca ]]; then
    #export bash variables
    export d
    #get datatable
        if [[ ! -f ${d}/${d}_matrix.csv ]]; then
        /home/felix/miniconda3/envs/tangram-env/bin/python -c 'import scanpy as sc; import pandas as pd; import os; file=os.environ["'"d"'"]; data=sc.read_10x_h5("'"/home/felix/projects/facial/felix/for_github/visium/00_mapping/"'"+file+"'"/outs/filtered_feature_bc_matrix.h5"'"); pd.DataFrame(data.X.toarray(),index=data.obs_names, columns=data.var_names).T.to_csv("'"/home/felix/projects/facial/felix/for_github/visium/00_mapping/"'"+file+"'"/"'"+file+"'"_matrix.csv"'",chunksize=5000)' 
        fi
    echo After tangram
    #run dca
    /home/felix/miniconda3/envs/signac-seurat/bin/dca /home/felix/projects/facial/felix/for_github/visium/00_mapping/${d}/${d}_matrix.csv --threads 30 /home/felix/projects/facial/felix/for_github/visium/00_mapping/${d}/outs/denoised_dca
    fi
    
    ###run velocity####################################
    echo After dca
    if [[ ! -f ${d}/velocyto/${d}.loom ]]; then
    
    /home/felix/miniconda3/envs/RNAvelo/bin/velocyto run10x --samtools-memory 100 --samtools-threads 50 -m ~/data/ref_genome/GRCh38_rmsk.gtf  ${d}/ ~/data/ref_genome/refdata-gex-GRCh38-2020-A/genes/genes.gtf
    
    fi
    
 done

