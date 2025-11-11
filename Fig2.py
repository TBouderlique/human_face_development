#!/usr/bin/env python
# coding: utf-8

# In[6]:


import cellrank
import scvelo as scv
import os
os.chdir("/home/yakov/face_data")
path = os.getcwd()
print(path)
import platform
print(platform.python_version())


# In[7]:


import seaborn as sns
import pandas as pd
import random
import numpy as np
import loompy as lp
import scanpy as scp
import scvelo as scv
import matplotlib.pyplot as plt
import scanpy.external as sce
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = ['sans-serif']
plt.rcParams['font.sans-serif'] = ['Arial']
random.seed(101)
get_ipython().run_line_magic('matplotlib', 'inline')


# In[8]:


plt.rcParams['figure.figsize'] = [12, 8]
plt.rcParams['figure.dpi'] = 100 # 200 e.g. is really fine, but slower
scp.set_figure_params(scanpy=True, fontsize=14 )


# In[9]:


import matplotlib.pyplot as plt

def plot_umap(phenotype, suffix):
    # Plot without immediately saving
    scp.pl.umap(
        adata,
        color=[phenotype],
        use_raw=False,
        sort_order=True,
        size=20,  # You can change this back to 30 if preferred
        legend_loc='right margin',
        show=False
    )

    # Rasterize the points
    ax = plt.gca()
    if ax.collections:
        ax.collections[0].set_rasterized(True)

    # Save directly to PDF with reduced file size
    output_path = f'figures/{Folder}/Fig2-SNPs-PV6{Dist}-{suffix}.pdf'
    plt.savefig(output_path, dpi=150)
    plt.close()

    print(f"✅ Saved: {output_path}") 


# In[10]:


import os
import re
import pandas as pd
import numpy as np
import scanpy as scp
from anndata import AnnData
from scipy.sparse import csr_matrix

# Preload AnnData
adata = scp.read_h5ad('DATA/all_samples.h5ad')

# Define distances and files
dist_values = [500000, 50000, 5000]
files = [
    'SchafferPV6.csv',
    'NormalVariation.csv',
    'Disease.csv',
    'White.csv'
]

# Loop over files and distances
for file in files:
    base_name = os.path.splitext(os.path.basename(file))[0].replace(" ", "_").replace("-", "_")
    
    for Dist in dist_values:
        print(f"🔄 Processing: {file}, Distance: {Dist}")
        Folder = f'{base_name}_Dist{Dist}'
        os.makedirs(f'figures/{Folder}', exist_ok=True)

        # Load and preprocess gene-phenotype data
        human = pd.read_csv(file)

        # Check if the 'Segment' column exists
        if 'Segment' in human.columns:
            # Remove the 'segment_' prefix and split if necessary
            human['Phenotype'] = human['Segment'].str.replace('segment_', '', regex=False)

            # Flatten the list if it's a list of values (explode the DataFrame if needed)
            human = human.explode('Phenotype')

            # Convert the 'Phenotype' column to integers after exploding
            human['Phenotype'] = human['Phenotype'].astype(str)
        else:
            print("⚠️ Column 'Segment' not found in the DataFrame.")


        # Step 1: Normalize Phenotype string: replace _ and , with comma, then split
        human['Phenotype'] = human['Phenotype'].apply(lambda x: re.split(r'[_\,]+', str(x)))

        # Step 2: Explode into separate rows
        human = human.explode('Phenotype')

        # Step 3: Clean up
        human['Phenotype'] = human['Phenotype'].str.strip()  # or .astype(int) if needed  


        # Pre-filter human data to reduce redundancy
        human = human[human['Distance'] < Dist]
        Genes = list(set(human.Gene.unique()) & set(adata.var_names))

        # Prepare the gene expression matrix
        A = pd.DataFrame(adata[:, Genes].X.A, columns=Genes, index=adata.obs.index)

        # Initialize new columns for phenotypes
        for Phenotype in human.Phenotype.unique():

            adata.obs[Phenotype] = 0
            adata.obs[Phenotype + ' Norm'] = 0

        # Process each phenotype
        for Phenotype in human.Phenotype.unique():
            # Filter genes associated with the current phenotype
            human_Phenotype = human[human.Phenotype == Phenotype]
            human_Phenotype_Genes = list(set(human_Phenotype['Gene'].unique()) & set(adata.var_names))

            # Add expression to the phenotype columns
            for gene in human_Phenotype_Genes:
                exp = A[gene]
                adata.obs[Phenotype] += exp
                adata.obs[Phenotype + ' Norm'] += exp / A[gene].sum()

        # Compute general phenotype expressions
        adata.obs['General Phenotype'] = adata.obs[[Phenotype for Phenotype in human.Phenotype.unique()]].sum(axis=1)
        adata.obs['General Phenotype Norm'] = adata.obs[[Phenotype + ' Norm' for Phenotype in human.Phenotype.unique()]].sum(axis=1)

        ## Plot UMAPs for each phenotype and general phenotype
        #def plot_umap(phenotype, suffix):
        #    scp.pl.umap(adata, color=[phenotype], use_raw=False, sort_order=True, size=30,
        #                legend_loc='on data', save=f'Fig2-SNPs-PV6{Dist}-{suffix}.pdf')
        #    os.rename(f'figures/umapFig2-SNPs-PV6{Dist}-{suffix}.pdf', f'figures/Fig2-SNPs-PV6{Dist}-{suffix}.pdf')
        #
        #import matplotlib.pyplot as plt



        for Phenotype in human.Phenotype.unique():
            plot_umap(Phenotype, Phenotype)

        plot_umap('General Phenotype', 'General')

        # Compute dot plot data
        gene_sets = human.groupby('Phenotype')['Gene'].apply(list).to_dict()
        gene_set_matrix = []
        gene_set_names = []

        # Compute average expression for each gene set
        for name, genes in gene_sets.items():
            valid_genes = [g for g in genes if g in adata.var_names]
            if valid_genes:
                expr = adata[:, valid_genes].X
                avg_expr = np.asarray(expr.mean(axis=1)).flatten() if hasattr(expr, 'mean') else expr.mean(axis=1)
                gene_set_matrix.append(avg_expr)
                gene_set_names.append(name)


        scp.tl.dendrogram(adata, groupby='annotation')

        # Create a new AnnData object for dot plot
        X = np.vstack(gene_set_matrix).T  # Transpose to (cells, gene_sets)
        pseudo_adata = AnnData(
            X=csr_matrix(X),
            obs=adata.obs.copy(),
            var=pd.DataFrame(index=gene_set_names)
        )

        # Remove phenotype columns from the pseudo_adata to avoid duplication
        pseudo_adata.obs = pseudo_adata.obs.drop(columns=[Phenotype for Phenotype in human.Phenotype.unique() if Phenotype in pseudo_adata.obs.columns])

        pseudo_adata.uns = adata.uns.copy()

        # sort by average expression (descending)
        avg_exprs = [expr.mean() for expr in gene_set_matrix]
        gene_set_names = [name for _, name in sorted(zip(avg_exprs, gene_set_names), reverse=True)]


        # Plot dot plot
        scp.pl.dotplot(
            pseudo_adata,
            var_names=gene_set_names,
            groupby='annotation',
            standard_scale='var',
            dendrogram=True,
            save=f'Fig2c-SNPs-PV6{Dist}-dotplot.pdf'
        )

        os.rename(f'figures/dotplot_Fig2c-SNPs-PV6{Dist}-dotplot.pdf', f'figures/{Folder}/Fig2c-SNPs-PV6{Dist}-dotplot.pdf')


# In[14]:


#White old
import re
import pandas as pd
import numpy as np
import os
import scanpy as scp
from anndata import AnnData
from scipy.sparse import csr_matrix

# File paths
file = 'GENES - 100kb - greater 10-8  SNPs per segment.csv'
Dist = 50000
adata = scp.read_h5ad('DATA/all_samples.h5ad')
human = pd.read_csv(file)

# Check if the 'Segment' column exists
if 'Segment' in human.columns:
    # Remove the 'segment_' prefix and split if necessary
    human['Phenotype'] = human['Segment'].str.replace('segment_', '', regex=False)
    
    # Flatten the list if it's a list of values (explode the DataFrame if needed)
    human = human.explode('Phenotype')
    
    # Convert the 'Phenotype' column to integers after exploding
    human['Phenotype'] = human['Phenotype'].astype(str)
else:
    print("⚠️ Column 'Segment' not found in the DataFrame.")

    
# Step 1: Normalize Phenotype string: replace _ and , with comma, then split
human['Phenotype'] = human['Phenotype'].apply(lambda x: re.split(r'[_\,]+', str(x)))

# Step 2: Explode into separate rows
human = human.explode('Phenotype')

# Step 3: Clean up
human['Phenotype'] = human['Phenotype'].str.strip()  # or .astype(int) if needed  
    

# Pre-filter human data to reduce redundancy
human = human[human['Distance'] < Dist]
Genes = list(set(human.Gene.unique()) & set(adata.var_names))

# Prepare the gene expression matrix
A = pd.DataFrame(adata[:, Genes].X.A, columns=Genes, index=adata.obs.index)

# Initialize new columns for phenotypes
for Phenotype in human.Phenotype.unique():
    
    adata.obs[Phenotype] = 0
    adata.obs[Phenotype + ' Norm'] = 0

# Process each phenotype
for Phenotype in human.Phenotype.unique():
    # Filter genes associated with the current phenotype
    human_Phenotype = human[human.Phenotype == Phenotype]
    human_Phenotype_Genes = list(set(human_Phenotype['Gene'].unique()) & set(adata.var_names))
    
    # Add expression to the phenotype columns
    for gene in human_Phenotype_Genes:
        exp = A[gene]
        adata.obs[Phenotype] += exp
        adata.obs[Phenotype + ' Norm'] += exp / A[gene].sum()

# Compute general phenotype expressions
adata.obs['General Phenotype'] = adata.obs[[Phenotype for Phenotype in human.Phenotype.unique()]].sum(axis=1)
adata.obs['General Phenotype Norm'] = adata.obs[[Phenotype + ' Norm' for Phenotype in human.Phenotype.unique()]].sum(axis=1)

# Plot UMAPs for each phenotype and general phenotype
def plot_umap(phenotype, suffix):
    scp.pl.umap(adata, color=[phenotype], use_raw=False, sort_order=True, size=30,
                legend_loc='on data', save=f'Fig2-SNPs{Dist}-{suffix}.pdf')
    os.rename(f'figures/umapFig2-SNPs{Dist}-{suffix}.pdf', f'figures/Fig2-SNPs{Dist}-{suffix}.pdf')

    
    
for Phenotype in human.Phenotype.unique():
    plot_umap(Phenotype, Phenotype)

plot_umap('General Phenotype', 'General')

# Compute dot plot data
gene_sets = human.groupby('Phenotype')['Gene'].apply(list).to_dict()
gene_set_matrix = []
gene_set_names = []

# Compute average expression for each gene set
for name, genes in gene_sets.items():
    valid_genes = [g for g in genes if g in adata.var_names]
    if valid_genes:
        expr = adata[:, valid_genes].X
        avg_expr = np.asarray(expr.mean(axis=1)).flatten() if hasattr(expr, 'mean') else expr.mean(axis=1)
        gene_set_matrix.append(avg_expr)
        gene_set_names.append(name)


scp.tl.dendrogram(adata, groupby='annotation')
        
# Create a new AnnData object for dot plot
X = np.vstack(gene_set_matrix).T  # Transpose to (cells, gene_sets)
pseudo_adata = AnnData(
    X=csr_matrix(X),
    obs=adata.obs.copy(),
    var=pd.DataFrame(index=gene_set_names)
)

# Remove phenotype columns from the pseudo_adata to avoid duplication
pseudo_adata.obs = pseudo_adata.obs.drop(columns=[Phenotype for Phenotype in human.Phenotype.unique() if Phenotype in pseudo_adata.obs.columns])

pseudo_adata.uns = adata.uns.copy()

# Plot dot plot
scp.pl.dotplot(
    pseudo_adata,
    var_names=sorted(gene_set_names, key=int),
    groupby='annotation',
    standard_scale='var',
    dendrogram=True,
    save=f'Fig2c-SNPs{Dist}-dotplot.pdf'
)

os.rename(f'figures/dotplot_Fig2c-SNPs{Dist}-dotplot.pdf', f'figures/Fig2c-SNPs{Dist}-dotplot.pdf')


# In[11]:


#Disease old
import re
import pandas as pd
import numpy as np
import os
import scanpy as scp
from anndata import AnnData
from scipy.sparse import csr_matrix

# File paths
file = 'Disease.csv'
Dist = 500000
adata = scp.read_h5ad('DATA/all_samples.h5ad')
human = pd.read_csv(file)

# Check if the 'Segment' column exists
if 'Segment' in human.columns:
    # Remove the 'segment_' prefix and split if necessary
    human['Phenotype'] = human['Segment'].str.replace('segment_', '', regex=False)
    
    # Flatten the list if it's a list of values (explode the DataFrame if needed)
    human = human.explode('Phenotype')
    
    # Convert the 'Phenotype' column to integers after exploding
    human['Phenotype'] = human['Phenotype'].astype(str)
else:
    print("⚠️ Column 'Segment' not found in the DataFrame.")

    
# Step 1: Normalize Phenotype string: replace _ and , with comma, then split
human['Phenotype'] = human['Phenotype'].apply(lambda x: re.split(r'[_\,]+', str(x)))

# Step 2: Explode into separate rows
human = human.explode('Phenotype')

# Step 3: Clean up
human['Phenotype'] = human['Phenotype'].str.strip()  # or .astype(int) if needed  
    

# Pre-filter human data to reduce redundancy
human = human[human['Distance'] < Dist]
Genes = list(set(human.Gene.unique()) & set(adata.var_names))

# Prepare the gene expression matrix
A = pd.DataFrame(adata[:, Genes].X.A, columns=Genes, index=adata.obs.index)

# Initialize new columns for phenotypes
for Phenotype in human.Phenotype.unique():
    
    adata.obs[Phenotype] = 0
    adata.obs[Phenotype + ' Norm'] = 0

# Process each phenotype
for Phenotype in human.Phenotype.unique():
    # Filter genes associated with the current phenotype
    human_Phenotype = human[human.Phenotype == Phenotype]
    human_Phenotype_Genes = list(set(human_Phenotype['Gene'].unique()) & set(adata.var_names))
    
    # Add expression to the phenotype columns
    for gene in human_Phenotype_Genes:
        exp = A[gene]
        adata.obs[Phenotype] += exp
        adata.obs[Phenotype + ' Norm'] += exp / A[gene].sum()

# Compute general phenotype expressions
adata.obs['General Phenotype'] = adata.obs[[Phenotype for Phenotype in human.Phenotype.unique()]].sum(axis=1)
adata.obs['General Phenotype Norm'] = adata.obs[[Phenotype + ' Norm' for Phenotype in human.Phenotype.unique()]].sum(axis=1)

# Plot UMAPs for each phenotype and general phenotype
def plot_umap(phenotype, suffix):
    scp.pl.umap(adata, color=[phenotype], use_raw=False, sort_order=True, size=30,
                legend_loc='on data', save=f'Fig2-SNPsDisease{Dist}-{suffix}.pdf')
    os.rename(f'figures/umapFig2-SNPsDisease{Dist}-{suffix}.pdf', f'figures/Fig2-SNPsDisease{Dist}-{suffix}.pdf')

    
    
for Phenotype in human.Phenotype.unique():
    plot_umap(Phenotype, Phenotype)

plot_umap('General Phenotype', 'General')

# Compute dot plot data
gene_sets = human.groupby('Phenotype')['Gene'].apply(list).to_dict()
gene_set_matrix = []
gene_set_names = []

# Compute average expression for each gene set
for name, genes in gene_sets.items():
    valid_genes = [g for g in genes if g in adata.var_names]
    if valid_genes:
        expr = adata[:, valid_genes].X
        avg_expr = np.asarray(expr.mean(axis=1)).flatten() if hasattr(expr, 'mean') else expr.mean(axis=1)
        gene_set_matrix.append(avg_expr)
        gene_set_names.append(name)


scp.tl.dendrogram(adata, groupby='annotation')
        
# Create a new AnnData object for dot plot
X = np.vstack(gene_set_matrix).T  # Transpose to (cells, gene_sets)
pseudo_adata = AnnData(
    X=csr_matrix(X),
    obs=adata.obs.copy(),
    var=pd.DataFrame(index=gene_set_names)
)

# Remove phenotype columns from the pseudo_adata to avoid duplication
pseudo_adata.obs = pseudo_adata.obs.drop(columns=[Phenotype for Phenotype in human.Phenotype.unique() if Phenotype in pseudo_adata.obs.columns])

pseudo_adata.uns = adata.uns.copy()

# Plot dot plot
scp.pl.dotplot(
    pseudo_adata,
    var_names=sorted(gene_set_names, key=str),
    groupby='annotation',
    standard_scale='var',
    dendrogram=True,
    save=f'Fig2c-SNPsDisease{Dist}-dotplot.pdf'
)

os.rename(f'figures/dotplot_Fig2c-SNPsDisease{Dist}-dotplot.pdf', f'figures/Fig2c-SNPsDisease{Dist}-dotplot.pdf')


# In[21]:


#NormalVariation old
import re
import pandas as pd
import numpy as np
import os
import scanpy as scp
from anndata import AnnData
from scipy.sparse import csr_matrix

# File paths
file = 'NormalVariation.csv'
Dist = 5000
adata = scp.read_h5ad('DATA/all_samples.h5ad')
human = pd.read_csv(file)

# Check if the 'Segment' column exists
if 'Segment' in human.columns:
    # Remove the 'segment_' prefix and split if necessary
    human['Phenotype'] = human['Segment'].str.replace('segment_', '', regex=False)
    
    # Flatten the list if it's a list of values (explode the DataFrame if needed)
    human = human.explode('Phenotype')
    
    # Convert the 'Phenotype' column to integers after exploding
    human['Phenotype'] = human['Phenotype'].astype(str)
else:
    print("⚠️ Column 'Segment' not found in the DataFrame.")

    
# Step 1: Normalize Phenotype string: replace _ and , with comma, then split
human['Phenotype'] = human['Phenotype'].apply(lambda x: re.split(r'[_\,]+', str(x)))

# Step 2: Explode into separate rows
human = human.explode('Phenotype')

# Step 3: Clean up
human['Phenotype'] = human['Phenotype'].str.strip()  # or .astype(int) if needed  
    

# Pre-filter human data to reduce redundancy
human = human[human['Distance'] < Dist]
Genes = list(set(human.Gene.unique()) & set(adata.var_names))

# Prepare the gene expression matrix
A = pd.DataFrame(adata[:, Genes].X.A, columns=Genes, index=adata.obs.index)

# Initialize new columns for phenotypes
for Phenotype in human.Phenotype.unique():
    
    adata.obs[Phenotype] = 0
    adata.obs[Phenotype + ' Norm'] = 0

# Process each phenotype
for Phenotype in human.Phenotype.unique():
    # Filter genes associated with the current phenotype
    human_Phenotype = human[human.Phenotype == Phenotype]
    human_Phenotype_Genes = list(set(human_Phenotype['Gene'].unique()) & set(adata.var_names))
    
    # Add expression to the phenotype columns
    for gene in human_Phenotype_Genes:
        exp = A[gene]
        adata.obs[Phenotype] += exp
        adata.obs[Phenotype + ' Norm'] += exp / A[gene].sum()

# Compute general phenotype expressions
adata.obs['General Phenotype'] = adata.obs[[Phenotype for Phenotype in human.Phenotype.unique()]].sum(axis=1)
adata.obs['General Phenotype Norm'] = adata.obs[[Phenotype + ' Norm' for Phenotype in human.Phenotype.unique()]].sum(axis=1)

# Plot UMAPs for each phenotype and general phenotype
def plot_umap(phenotype, suffix):
    scp.pl.umap(adata, color=[phenotype], use_raw=False, sort_order=True, size=30,
                legend_loc='on data', save=f'Fig2-SNPsNormalVariation{Dist}-{suffix}.pdf')
    os.rename(f'figures/umapFig2-SNPsNormalVariation{Dist}-{suffix}.pdf', f'figures/Fig2-SNPsNormalVariation{Dist}-{suffix}.pdf')

    
    
for Phenotype in human.Phenotype.unique():
    plot_umap(Phenotype, Phenotype)

plot_umap('General Phenotype', 'General')

# Compute dot plot data
gene_sets = human.groupby('Phenotype')['Gene'].apply(list).to_dict()
gene_set_matrix = []
gene_set_names = []

# Compute average expression for each gene set
for name, genes in gene_sets.items():
    valid_genes = [g for g in genes if g in adata.var_names]
    if valid_genes:
        expr = adata[:, valid_genes].X
        avg_expr = np.asarray(expr.mean(axis=1)).flatten() if hasattr(expr, 'mean') else expr.mean(axis=1)
        gene_set_matrix.append(avg_expr)
        gene_set_names.append(name)


scp.tl.dendrogram(adata, groupby='annotation')
        
# Create a new AnnData object for dot plot
X = np.vstack(gene_set_matrix).T  # Transpose to (cells, gene_sets)
pseudo_adata = AnnData(
    X=csr_matrix(X),
    obs=adata.obs.copy(),
    var=pd.DataFrame(index=gene_set_names)
)

# Remove phenotype columns from the pseudo_adata to avoid duplication
pseudo_adata.obs = pseudo_adata.obs.drop(columns=[Phenotype for Phenotype in human.Phenotype.unique() if Phenotype in pseudo_adata.obs.columns])

pseudo_adata.uns = adata.uns.copy()

# Plot dot plot
scp.pl.dotplot(
    pseudo_adata,
    var_names=sorted(gene_set_names, key=str),
    groupby='annotation',
    standard_scale='var',
    dendrogram=True,
    save=f'Fig2c-SNPsNormalVariation{Dist}-dotplot.pdf'
)

os.rename(f'figures/dotplot_Fig2c-SNPsNormalVariation{Dist}-dotplot.pdf', f'figures/Fig2c-SNPsNormalVariation{Dist}-dotplot.pdf')


# In[ ]:


#scp.pl.umap(adata, color=['General Phenotype', 'General Phenotype' + ' Norm'],use_raw=False,sort_order=True,size=30,
#            legend_loc='on data'
#           )
#
#
#human = pd.read_csv('FaceGenes.csv')
#for Phenotype in human.Phenotype.unique():
#    scp.pl.umap(adata, color=[Phenotype, Phenotype + ' Norm'],use_raw=False,sort_order=True,size=30,
#                legend_loc='on data'#, save = "Face-All-highly_variable-NO-harmony-" + Phenotype + ".png"
#               )

