#!/usr/bin/env python
# coding: utf-8

# In[1]:


import cellrank
import scvelo as scv
import os
os.chdir("/home/yakov/face_data/scFates")
path = os.getcwd()
print(path)
import platform
print(platform.python_version())


# In[2]:


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


# In[3]:


plt.rcParams['figure.figsize'] = [12, 8]
plt.rcParams['figure.dpi'] = 100 # 200 e.g. is really fine, but slower
scp.set_figure_params(scanpy=True, fontsize=14 )


# In[4]:


import os, sys
os.environ['R_HOME'] = sys.exec_prefix+"/lib/R/"
import scFates as scf
import warnings
warnings.filterwarnings("ignore")
import os, sys
# to avoid any possible jupyter crashes due to rpy2 not finding the R install on conda
os.environ['R_HOME'] = sys.exec_prefix+"/lib/R/"
from anndata import AnnData
import numpy as np
import pandas as pd
import scanpy as sc
import scFates as scf
import palantir
import matplotlib.pyplot as plt
sc.settings.verbosity = 3
sc.settings.logfile = sys.stdout

## fix palantir breaking down some plots
import seaborn
seaborn.reset_orig()
get_ipython().run_line_magic('matplotlib', 'inline')
sc.set_figure_params()
scf.set_figure_pubready()


# In[63]:


adata = scp.read_h5ad('jaw_Sample1.h5ad')


# In[64]:


norm_df = sc.pp.normalize_per_cell(adata,copy=True)
norm_df = palantir.preprocess.log_transform(norm_df)


# In[65]:


sc.pp.highly_variable_genes(adata, n_top_genes=2500, flavor='cell_ranger')
sc.pp.pca(adata)


# In[66]:


pca_projections = pd.DataFrame(adata.obsm["X_pca"],index=adata.obs_names)
dm_res = palantir.utils.run_diffusion_maps(pca_projections)
ms_data = palantir.utils.determine_multiscale_space(dm_res,n_eigs=4)


# In[67]:


# generate neighbor draph in multiscale diffusion space
adata.obsm["X_palantir"]=ms_data.values
sc.pp.neighbors(adata,n_neighbors=30,use_rep="X_palantir")


# In[68]:


# draw ForceAtlas2 embedding using 2 first PCs as initial positions
adata.obsm["X_pca2d"]=adata.obsm["X_pca"][:,:2]
sc.tl.draw_graph(adata,init_pos='X_pca2d')


# In[69]:


sc.pl.draw_graph(adata,color_map="RdBu_r")


# In[123]:


adata.write('jaw.h5ad')


# In[5]:


adata = scp.read_h5ad('jaw.h5ad')


# In[6]:


## Load CSV file with two new columns (must have matching index)
#csv_data = pd.read_csv('all_samples_55clusters.csv', index_col=0)  # Make sure to load index

## Select the two columns you want to add
#columns_to_add = csv_data[['general_annotation', 'specific_annotation']]

## Join based on index
#adata.obs = adata.obs.join(columns_to_add)


# In[7]:


sc.settings.set_figure_params(dpi=80, frameon=True, figsize=(12, 12))
sc.pl.umap(adata, color='general_annotation', size = 100 )
#sc.pl.draw_graph(adata,color_map="RdBu_r", color='general_annotation', legend_loc="on data", legend_fontsize='xx-small', size = 10)


# In[8]:


sc.pl.umap(adata, color='general_annotation', groups=['progenitor'], size = 50)


# In[9]:


sc.pl.draw_graph(adata,color_map="RdBu_r", color='general_annotation', size = 10)
#sc.pl.draw_graph(adata,color_map="RdBu_r", color='Age', size = 10)


# In[10]:


#n=400
#l=1
#s=.01

sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(12, 12))

#scf.tl.tree(adata,method="ppt",Nodes=n,use_rep="palantir",plot=True, basis="draw_graph_fa",
#device="cpu",seed=1,ppt_lambda=l,ppt_sigma=s,ppt_nsteps=100,node_fontsize=1,show_node_labels=False)


n=500
l=10
s=0.4



scf.tl.tree(adata,method="ppt",Nodes=n,use_rep="umap",plot=True, basis="umap", 
device="cpu",seed=1,ppt_lambda=l,ppt_sigma=s,ppt_nsteps=300,node_fontsize=1,show_node_labels=False)


# In[11]:


scf.tl.cleanup(adata, minbranchlength = 3)


# In[12]:


scf.pl.graph(adata, size_nodes = 0, tips = True,forks=False, basis = 'umap')


# In[16]:


scf.pl.graph(adata,tips = True, forks=False, basis = 'umap')
scf.tl.cleanup(adata,minbranchlength = 0, leaves=[18,31,102,127,145,157,178,192,200,202,210,227,318,400,417,457])
scf.pl.graph(adata,tips = True, forks=False, basis = 'umap')


# In[13]:


#sc.pl.embedding(adata, basis='X_draw_graph_fa', color='TNMD')
#sc.pl.embedding(adata, basis='X_draw_graph_fa', color=['RUNX2'], size=10)
#sc.pl.embedding(adata, basis='X_draw_graph_fa', color=['SP7'], size=10)
#sc.pl.embedding(adata, basis='X_draw_graph_fa', color=['ALPL'], size=10)
#sc.pl.embedding(adata, basis='X_draw_graph_fa', color=['MSX1'], size=10)


# In[17]:


scf.tl.root(adata,150)
scf.tl.pseudotime(adata,n_jobs=70,n_map=100,seed=45)
scf.pl.trajectory(adata, basis = 'umap')


# In[19]:


sc.pl.umap(adata, color=["seg", "milestones"],legend_loc = 'on data')


# In[43]:


scf.tl.rename_milestones(adata,["Perivascular","Dermal","130","143","Root","154","168","204","Submandibular_1","219","225","232","Meningeal","Bone","239","255","256","Stroma","Perimysium","Tendon","356","Submandibular_2","Cartilage","45","6","68","85","88"])
sc.pl.umap(adata, color=["seg", "milestones"],legend_loc = 'on data')
# we change the color of the root milestone for better visualisations
adata.uns["milestones_colors"][3]="#17bece"


# In[44]:


scf.pl.milestones(adata,annotate=True, basis = 'umap')


# In[45]:


sc.set_figure_params()
fig, axs=plt.subplots(2,2,figsize=(8,8))
axs=axs.ravel()
scf.pl.graph(adata,basis="umap",show=False,ax=axs[0],legend_fontsize
='xx-smal')
scf.pl.trajectory(adata,basis="umap",show=False,ax=axs[1])
sc.pl.umap(adata,color=["seg"],legend_loc="on data",show=False,ax=axs[2],legend_fontoutline=True)
scf.pl.milestones(adata,ax=axs[3],show=False,annotate=True, basis = 'umap')
plt.savefig("figures/A-jaw-umap.pdf",dpi=300)


# In[46]:


scf.tl.dendrogram(adata)


# In[47]:


scf.pl.dendrogram(adata,color="seg")


# In[48]:


sc.set_figure_params(figsize=(15,15),frameon=False,dpi_save=300)
scf.pl.dendrogram(adata,color="t",show_info=False,save="B1-jaw-umap.pdf",cmap="viridis")
scf.pl.dendrogram(adata,color="milestones",legend_loc="on data",color_milestones=True,legend_fontoutline=True,save="B2-jaw-umap.pdf")


# In[49]:


adata.write('jaw_Sample1b.h5ad')


# In[50]:


adata = scp.read_h5ad('jaw_Sample1b.h5ad')


# In[51]:


scf.tl.test_association(adata,n_jobs=60)


# In[52]:


sc.set_figure_params()
scf.pl.test_association(adata)
plt.savefig("figures/C-jaw-umap.pdf",dpi=300)


# In[53]:


adata.write('jaw_Sample1c.h5ad')


# In[54]:


adata = scp.read_h5ad('jaw_Sample1c.h5ad')


# In[55]:


scf.tl.fit(adata,n_jobs=80)


# In[56]:


adata.write('jaw_Sample1d.h5ad')


# In[5]:


adata = scp.read_h5ad('jaw_Sample1d.h5ad')


# In[18]:


for k in adata.uns.keys():
    print(k)


# In[20]:


# Get the result table
assoc = adata.uns["stat_assoc_list"]

print(assoc)


# In[9]:


import os
import scanpy as sc
import scFates as scf
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap

# --- Define genes of interest ---
genes_of_interest = [
    "ALX1", "BARX1", "DLX5",
    "EBF3", "EGR1", "EMX2", "FOXC1", "FOXP2", "GATA2", "ID2", "ID3",
    "MEOX2", "MKX", "MOXD1", 
    "OSR2", "PAX3", "PAX9", 
    "PRRX1", "SIX1", "TBX3", "TFAP2B",
    "PAX7", 
    "SHOX2", 
    "NR2F1", "TRPS1", "HDAC9", "DACH1", "VGLL3", "MAB21L2", "FOXD1", "FOXF2", "MSX1", "OSR1", "TWIST2", "ZNF804A"
]

# --- Filter for valid genes ---
valid_genes = [g for g in genes_of_interest if g in adata.var_names]
missing_genes = [g for g in genes_of_interest if g not in adata.var_names]

print(f"✅ Found {len(valid_genes)} valid genes in adata.")
if missing_genes:
    print(f"⚠️ {len(missing_genes)} missing: {missing_genes[:5]}{'...' if len(missing_genes) > 5 else ''}")

# --- Set figure settings ---
sc.set_figure_params(figsize=(4, 4), frameon=False, dpi=100)

# --- Define colormap ---
custom_cmap = LinearSegmentedColormap.from_list(
    "lightgrey_to_darkred",
    ["#f0f0f0", "#ffaaaa", "#8b0000"]
)

# --- Create multi-plot pages ---
output_pdf = "PSEUDOTIME_PART1_jaw.pdf"
plots_per_page = 35
rows, cols = 5, 7

with PdfPages(output_pdf) as pdf:
    fig, axes = plt.subplots(rows, cols, figsize=(21, 15))  # Adjust size to fit 5x7 grid
    axes = axes.flatten()
    plot_count = 0

    for idx, gene in enumerate(valid_genes):
        ax = axes[plot_count % plots_per_page]
        try:
            print(f"🔹 Plotting {gene}...")
            scf.pl.dendrogram(
                adata,
                color=gene,
                legend_loc="on data",
                color_milestones=False,
                legend_fontoutline=True,
                size=15,  # Smaller point size for dense layout
                cmap=custom_cmap,
                show=False,
                ax=ax
            )
            ax.set_title(gene, fontsize=12, fontweight='bold')  # Larger, bold gene name
        except Exception as e:
            print(f"❌ Failed to plot {gene}: {e}")
            ax.axis("off")

        plot_count += 1

        # Save and reset every 35 plots
        if plot_count % plots_per_page == 0 or gene == valid_genes[-1]:
            # Hide unused subplots on last page
            if plot_count % plots_per_page != 0:
                for j in range(plot_count % plots_per_page, plots_per_page):
                    axes[j].axis("off")
            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

            # Prepare next page (if not last gene)
            if gene != valid_genes[-1]:
                fig, axes = plt.subplots(rows, cols, figsize=(21, 15))
                axes = axes.flatten()
                plot_count = 0

print(f"\n✅ PDF saved: {output_pdf}")


# In[6]:


import os
import scanpy as sc
import scFates as scf
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap

# --- Define genes of interest ---
genes_of_interest = [
    "ALX1", "BARX1", "C1orf61", "DLX2", "DLX3", "DLX5", "DLX6-AS1", "SATB2-AS1",
    "EBF3", "EGR1", "EMX2", "FOXC1", "FOXC2", "FOXD1", "FOXE3", "FOXF2", "FOXG1",
    "FOXL2", "FOXP2", "GATA2", "GATA3", "HAND1", "HIC1", "HMX1", "ID2", "ID3",
    "IRX3", "IRX5", "MAFB", "MEOX1", "MEOX2", "MKX", "MOXD1", "MSX1", "NR5A2",
    "OSR1", "OSR2", "PAX1", "PAX3", "PAX9", "PHOX2B", "PITX1", "PITX2", "POU3F3",
    "PROX1", "PRRX1", "SCX", "SHOX", "SIX1", "SIX2", "SOX2", "TBX3", "TFAP2B",
    "TWIST2", "ZIC1", "ZIC2", "ZNF503", "ZNF804A", "PAX6", "PAX7", "HAND2", "IRX1",
    "MSX2", "NEUROD1", "SHOX2", "INSM2", "MYOG", "NKX3-2", "POU3F4", "TLX3",
    "NR2F1", "TRPS1", "HDAC9", "H3F3B", "PRDM6", "SIM2", "DACH1", "VGLL3", "MAB21L2"
]

# --- Filter for valid genes ---
valid_genes = [g for g in genes_of_interest if g in adata.var_names]
missing_genes = [g for g in genes_of_interest if g not in adata.var_names]

print(f"✅ Found {len(valid_genes)} valid genes in adata.")
if missing_genes:
    print(f"⚠️ {len(missing_genes)} missing: {missing_genes[:5]}{'...' if len(missing_genes) > 5 else ''}")

# --- Set figure settings ---
sc.set_figure_params(figsize=(15, 15), frameon=False, dpi=100)

# --- Define colormap ---
custom_cmap = LinearSegmentedColormap.from_list(
    "lightgrey_to_darkred",
    ["#f0f0f0", "#ffaaaa", "#8b0000"]
)

# --- Create multipage PDF ---
output_pdf = "genes_of_interest_dendrograms2_jaw.pdf"
with PdfPages(output_pdf) as pdf:
    for gene in valid_genes:
        try:
            print(f"🔹 Plotting {gene}...")
            scf.pl.dendrogram(
                adata,
                color=gene,
                legend_loc="on data",
                color_milestones=False,
                legend_fontoutline=True,
                size=150,
                cmap=custom_cmap,
                show=False
            )
            fig = plt.gcf()
            pdf.savefig(fig)      # ✅ Save BEFORE show
            plt.show()            # 👀 Display on screen
            plt.close(fig)        # Close to avoid overlap
        except Exception as e:
            print(f"❌ Failed to plot {gene}: {e}")

print(f"\n✅ PDF saved: {output_pdf}")


# In[6]:


import pandas as pd
import scanpy as sc
import scFates as scf

# Load gene list from CSV
gene_df = pd.read_csv("master_spatial_n2.csv", index_col=0)
genes = gene_df.index.tolist()

# List of genes of interest
genes = [
    "ALX1", "DLX5", "MSX1"
]

# Filter genes to those present in adata.var_names
valid_genes = [g for g in genes if g in adata.var_names]
print(f"Loaded {len(genes)} genes, using {len(valid_genes)} valid ones found in adata.")

# Set figure parameters
sc.set_figure_params(figsize=(15, 15), frameon=False, dpi_save=300)

# Loop and plot
for gene in valid_genes:
    print(f"Plotting {gene}...")
    scf.pl.dendrogram(
        adata,
        color=gene,
        legend_loc="on data",
        color_milestones=False,
        legend_fontoutline=True,
        save=f"-B2-jaw-umap-{gene}.pdf"
    )

print("✅ Done.")


# In[39]:


import matplotlib.pyplot as plt
import numpy as np
import os

# Output folder
output_dir = "dendro_genes_overlay"
os.makedirs(output_dir, exist_ok=True)

# Genes of interest
genes_of_interest = ["ALX1", "DLX5", "MSX1"]

# Filter valid genes
valid_genes = [g for g in genes_of_interest if g in adata.var_names]

# Extract dendrogram coordinates
dendro_coords = adata.obsm['X_dendro']

for gene in valid_genes:
    try:
        print(f"🔹 Plotting {gene}...")

        # --- Create dendrogram edges first ---
        fig, ax = plt.subplots(figsize=(6, 5))
        scf.pl.dendrogram(
            adata,
            color=None,        # just draw lines
            show=False,
            ax=ax
        )

        # --- Overlay points colored by gene expression ---
        expr = adata[:, gene].X.toarray().flatten() if hasattr(adata[:, gene].X, "toarray") else adata[:, gene].X.flatten()

        # Use scFates function to get the same colormap as scFates uses
        from matplotlib.colors import Normalize
        from matplotlib.cm import ScalarMappable

        norm = Normalize(vmin=expr.min(), vmax=expr.max())
        sm = ScalarMappable(norm=norm, cmap='Reds')  # match scFates default for continuous genes
        colors = sm.to_rgba(expr)

        ax.scatter(
            dendro_coords[:, 0],
            dendro_coords[:, 1],
            c=colors,
            s=0.01,       # very small dots
            alpha=0.2     # semi-transparent
        )

        ax.set_title(f"{gene}", fontsize=12)
        ax.axis('off')

        # Add colorbar
        cbar = fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label('Expression', fontsize=10)
        cbar.ax.tick_params(labelsize=10)

        # Save as high-resolution PNG
        plt.savefig(os.path.join(output_dir, f"{gene}_overlay_jaw.png"), dpi=300, bbox_inches='tight')
        plt.close(fig)

    except Exception as e:
        print(f"❌ Failed for {gene}: {e}")

print(f"🎉 All gene plots saved as high-resolution PNGs in '{output_dir}'")


# In[40]:


import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

# Output folder
output_dir = "dendro_genes_overlay"
os.makedirs(output_dir, exist_ok=True)

# Genes of interest
genes_of_interest = ["ALX1", "DLX5", "MSX1"]

# Filter valid genes
valid_genes = [g for g in genes_of_interest if g in adata.var_names]

# Extract dendrogram coordinates
dendro_coords = adata.obsm['X_dendro']

for gene in valid_genes:
    try:
        print(f"🔹 Plotting {gene}...")

        # --- Create dendrogram edges first ---
        fig, ax = plt.subplots(figsize=(6, 5))
        scf.pl.dendrogram(
            adata,
            color=None,          # just draw lines
            show=False,
            ax=ax,
            linewidth=0.5,       # thinner lines
            edgecolor='gray'     # lighter lines
        )

        # --- Overlay points colored by gene expression ---
        expr = adata[:, gene].X.toarray().flatten() if hasattr(adata[:, gene].X, "toarray") else adata[:, gene].X.flatten()

        norm = Normalize(vmin=expr.min(), vmax=expr.max())
        sm = ScalarMappable(norm=norm, cmap='Reds')
        colors = sm.to_rgba(expr)

        ax.scatter(
            dendro_coords[:, 0],
            dendro_coords[:, 1],
            c=colors,
            s=5,         # slightly larger dots
            alpha=0.4    # slightly more opaque
        )

        ax.set_title(f"{gene}", fontsize=12)
        ax.axis('off')

        # Add colorbar
        cbar = fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label('Expression', fontsize=10)
        cbar.ax.tick_params(labelsize=10)

        # Save as high-resolution PNG
        plt.savefig(os.path.join(output_dir, f"{gene}_overlay_jaw.png"), dpi=300, bbox_inches='tight')
        plt.close(fig)

    except Exception as e:
        print(f"❌ Failed for {gene}: {e}")

print(f"🎉 All gene plots saved as high-resolution PNGs in '{output_dir}'")


# In[ ]:





# In[42]:


adata.var[adata.var["signi"]].to_csv('jaw significant features', index=True)


# In[58]:


scf.tl.test_fork(adata,root_milestone="Root",milestones = ["109","Root","119","157","179","36","37","40","70","85"],n_jobs=40,rescale=True)


# In[ ]:





# In[ ]:





# In[ ]:





# In[156]:


adata


# In[ ]:




