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


# In[4]:


adata = scp.read_h5ad('nose_Sample1.h5ad')


# In[6]:


norm_df = sc.pp.normalize_per_cell(adata,copy=True)
norm_df = palantir.preprocess.log_transform(norm_df)


# In[7]:


sc.pp.highly_variable_genes(adata, n_top_genes=1500, flavor='cell_ranger')
sc.pp.pca(adata)


# In[8]:


pca_projections = pd.DataFrame(adata.obsm["X_pca"],index=adata.obs_names)
dm_res = palantir.utils.run_diffusion_maps(pca_projections)
ms_data = palantir.utils.determine_multiscale_space(dm_res,n_eigs=4)


# In[9]:


# generate neighbor draph in multiscale diffusion space
adata.obsm["X_palantir"]=ms_data.values
sc.pp.neighbors(adata,n_neighbors=30,use_rep="X_palantir")


# In[10]:


# draw ForceAtlas2 embedding using 2 first PCs as initial positions
adata.obsm["X_pca2d"]=adata.obsm["X_pca"][:,:2]
sc.tl.draw_graph(adata,init_pos='X_pca2d')


# In[43]:


sc.pl.draw_graph(adata,color_map="RdBu_r", color='annotation')


# In[44]:


sc.pl.draw_graph(adata,color_map="RdBu_r", color='Age')


# In[9]:


adata.write('nose_Sample1a.h5ad')


# In[38]:


adata = scp.read_h5ad('nose_Sample1a.h5ad')
#adata = scp.read_h5ad('nose_Sample.h5ad')


# In[39]:


sc.settings.set_figure_params(dpi=80, frameon=True, figsize=(12, 12))


# In[40]:


## Load CSV file with two new columns (must have matching index)
#csv_data = pd.read_csv('all_samples_55clusters.csv', index_col=0)  # Make sure to load index

## Select the two columns you want to add
#columns_to_add = csv_data[['general_annotation', 'specific_annotation']]

## Join based on index
#adata.obs = adata.obs.join(columns_to_add)


# In[41]:


sc.settings.set_figure_params(dpi=80, frameon=True, figsize=(12, 12))
sc.pl.draw_graph(adata,color_map="RdBu_r", color='general_annotation', legend_loc="on data", legend_fontsize
='xx-small', size = 10)


# In[42]:


scp.pp.neighbors(adata, n_neighbors = 40, use_rep='X_pca_harmony')
scp.tl.leiden(adata, resolution=0.4)
scp.tl.umap(adata, n_components=2)


# In[43]:


#200 principal points, sigma = 0.001, lambda = 10
#n=400
#l=1
#s=.01
#n=400
#l=10
#s=1
n=200
l=100
s=.5


#scf.tl.tree(adata,method="ppt",Nodes=n,use_rep="palantir",plot=True, basis="draw_graph_fa", 
#device="cpu",seed=1,ppt_lambda=l,ppt_sigma=s,ppt_nsteps=100,node_fontsize=1,show_node_labels=False)
scf.tl.tree(adata,method="ppt",Nodes=n,use_rep="umap",plot=True, basis="umap", 
device="cpu",seed=1,ppt_lambda=l,ppt_sigma=s,ppt_nsteps=200,node_fontsize=1,show_node_labels=False)


# In[44]:


scf.pl.graph(adata, size_nodes = 0, tips = True,forks=False, basis = 'umap')


# In[45]:


scf.pl.graph(adata,tips = True, forks=False, basis = 'umap')
scf.tl.cleanup(adata,minbranchlength = 0, leaves=[23,110,37,110,122,174,160,69,82,91,90,164,49,109])
scf.tl.cleanup(adata,minbranchlength = 0, leaves=[122])
scf.pl.graph(adata,tips = True, forks=False, basis = 'umap')


# In[46]:


#sc.pl.embedding(adata, basis='X_draw_graph_fa', color='TNMD')
#sc.pl.embedding(adata, basis='X_draw_graph_fa', color=['RUNX2'], size=10)
#sc.pl.embedding(adata, basis='X_draw_graph_fa', color=['SP7'], size=10)
#sc.pl.embedding(adata, basis='X_draw_graph_fa', color=['ALPL'], size=10)
#sc.pl.embedding(adata, basis='X_draw_graph_fa', color=['MSX1'], size=10)
#sc.pl.embedding(adata, basis='X_dendro', color=['FGF8','SHH','HAND2'], size=10)


# In[47]:


#scf.tl.root(adata,13)
#scf.tl.pseudotime(adata,n_jobs=60,n_map=100,seed=42)


# In[48]:


#adata.obs.seg = adata.obs.seg.astype("category")


# In[49]:


#scf.pl.graph(adata,tips = True, forks=True, color_cells="milestones")


# In[50]:


#scf.tl.subset_tree(adata, root_milestone="13", milestones = ["55"], mode='substract' )


# In[51]:


scf.tl.merge_n_simplify(adata)
scf.pl.graph(adata,tips = True, forks=False, basis = 'umap')


# In[72]:


scf.tl.root(adata,123)
scf.tl.pseudotime(adata,n_jobs=60,n_map=100,seed=42)
scf.pl.trajectory(adata, basis = 'umap')


# In[73]:


sc.pl.umap(adata, color=["seg", "milestones"],legend_loc = 'on data')


# In[74]:


scf.tl.rename_milestones(adata,['Bone','Stromal','Cartilage','114','119','Root','127','dermal','139','3','Perichondrial','45','46','Dental-like','Palatal-like','Teno','Perich','77','88','Perivascular'])
sc.pl.umap(adata, color=["seg", "milestones"],legend_loc = 'on data')
# we change the color of the root milestone for better visualisations
#adata.uns["milestones_colors"][3]="#17bece"


# In[75]:


scf.pl.milestones(adata,annotate=True, basis = 'umap')


# In[76]:


sc.set_figure_params()
fig, axs=plt.subplots(2,2,figsize=(15,15))
axs=axs.ravel()
scf.pl.graph(adata,basis="umap",show=False,ax=axs[0])
scf.pl.trajectory(adata,basis="umap",show=False,ax=axs[1])
sc.pl.umap(adata,color=["seg"],legend_loc="on data",show=False,ax=axs[2],legend_fontoutline=True)
scf.pl.milestones(adata,ax=axs[3],show=False,annotate=True, basis = 'umap')
plt.savefig("figures/A-nose-umap.pdf",dpi=300)


# In[77]:


scf.tl.dendrogram(adata)


# In[78]:


scf.pl.dendrogram(adata,color="seg")


# In[80]:


sc.set_figure_params(figsize=(15,15),frameon=False,dpi_save=300)
scf.pl.dendrogram(adata,color="t",show_info=False,save="B1-nose-umap.pdf",cmap="viridis")
scf.pl.dendrogram(adata,color="milestones",legend_loc="on data",color_milestones=True,legend_fontoutline=True,save="B2-nose-umap.pdf")


# In[81]:


adata.write('nose_Sample1b_140825.h5ad')


# In[82]:


adata = scp.read_h5ad('nose_Sample1b_140825.h5ad')


# In[83]:


scf.tl.test_association(adata,n_jobs=60)


# In[84]:


scf.pl.test_association(adata)
plt.savefig("figures/C-nose-umap.pdf",dpi=300)


# In[85]:


adata.write('nose_Sample1c_140825.h5ad')


# In[86]:


adata = scp.read_h5ad('nose_Sample1c_140825.h5ad')


# In[87]:


scf.tl.fit(adata,n_jobs=60)


# In[88]:


adata.write('nose_Sample1d_140825.h5ad')


# In[5]:


adata = scp.read_h5ad('nose_Sample1d_140825.h5ad')


# In[13]:


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

# --- High DPI settings ---
sc.set_figure_params(figsize=(4, 4), frameon=False, dpi=300)

# --- Define high-quality colormap ---
custom_cmap = LinearSegmentedColormap.from_list(
    "lightgrey_to_darkred",
    ["#f0f0f0", "#ffaaaa", "#8b0000"]
)

# --- Create high-quality multi-plot PDF ---
output_pdf = "PSEUDOTIME_PART1_nose.pdf"
plots_per_page = 35
rows, cols = 5, 7

with PdfPages(output_pdf) as pdf:
    fig, axes = plt.subplots(rows, cols, figsize=(21, 15), constrained_layout=True)
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
                size=15,
                cmap=custom_cmap,
                show=False,
                ax=ax
            )
            ax.set_title(gene, fontsize=12, fontweight='bold')
            ax.set_rasterized(True)  # ✅ Makes dense plots more efficient in vector PDFs
        except Exception as e:
            print(f"❌ Failed to plot {gene}: {e}")
            ax.axis("off")

        plot_count += 1

        # Save and reset every 35 plots
        if plot_count % plots_per_page == 0 or gene == valid_genes[-1]:
            # Hide unused subplots
            if plot_count % plots_per_page != 0:
                for j in range(plot_count % plots_per_page, plots_per_page):
                    axes[j].axis("off")
            pdf.savefig(fig, dpi=300)  # ✅ High-resolution export
            plt.close(fig)

            # Prepare new page if needed
            if gene != valid_genes[-1]:
                fig, axes = plt.subplots(rows, cols, figsize=(21, 15), constrained_layout=True)
                axes = axes.flatten()
                plot_count = 0

print(f"\n✅ High-resolution PDF saved: {output_pdf}")


# In[10]:


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
sc.set_figure_params(figsize=(4, 4), frameon=False, dpi=300)

# --- Define colormap ---
custom_cmap = LinearSegmentedColormap.from_list(
    "lightgrey_to_darkred",
    ["#f0f0f0", "#ffaaaa", "#8b0000"]
)

# --- Create multi-plot pages ---
output_pdf = "PSEUDOTIME_PART1_nose.pdf"
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


# In[90]:


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
output_pdf = "genes_of_interest_dendrograms2.pdf"
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


# In[203]:


#adata.var[adata.var["signi"]].to_csv('nose significant features', index=True)


# In[214]:


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
        plt.savefig(os.path.join(output_dir, f"{gene}_overlay_nose.png"), dpi=300, bbox_inches='tight')
        plt.close(fig)

    except Exception as e:
        print(f"❌ Failed for {gene}: {e}")

print(f"🎉 All gene plots saved as high-resolution PNGs in '{output_dir}'")


# In[18]:


scf.tl.test_fork(adata,root_milestone="Root",milestones = ['Bone','Stromal','Cartilage','Perichondrial','Perivascular'],n_jobs=60,rescale=True)
#scf.tl.test_fork(adata,root_milestone="Root",milestones = ["m170","m90"],n_jobs=60,rescale=True)


# In[19]:


scf.tl.branch_specific(adata,root_milestone="Root",milestones= ['Bone','Stromal','Cartilage','Perichondrial','Perivascular'],effect=.5)
#scf.tl.branch_specific(adata,root_milestone="Root",milestones= ["m170","m90"],effect=0.5)


# In[22]:


sc.set_figure_params(figsize=(15,15),frameon=False,dpi_save=300)
g1=scf.pl.trends(adata,
                 root_milestone="Root",
                 milestones= ['Bone','Stromal','Cartilage','Perichondrial','Perivascular'],
                 branch="Bone", n_features = 20, fontsize = 5 ,basis ='umap',
                 plot_emb=True,ordering="max",return_genes=True, save = "Bone.pdf" )


# In[23]:


sc.set_figure_params(figsize=(15,15),frameon=False,dpi_save=300)
g1=scf.pl.trends(adata,
                 root_milestone="Root",
                 milestones= ['Bone','Stromal','Cartilage','Perichondrial','Perivascular'],
                 branch="Stromal", n_features = 20, fontsize = 5 ,basis ='umap',
                 plot_emb=True,ordering="max",return_genes=True, save = "Stromal.pdf" )


# In[24]:


sc.set_figure_params(figsize=(15,15),frameon=False,dpi_save=300)
g1=scf.pl.trends(adata,
                 root_milestone="Root",
                 milestones= ['Bone','Stromal','Cartilage','Perichondrial','Perivascular'],
                 branch="Cartilage", n_features = 20, fontsize = 5 ,basis ='umap',
                 plot_emb=True,ordering="max",return_genes=True, save = "Cartilage.pdf" )


# In[25]:


sc.set_figure_params(figsize=(15,15),frameon=False,dpi_save=300)
g1=scf.pl.trends(adata,
                 root_milestone="Root",
                 milestones= ['Bone','Stromal','Cartilage','Perichondrial','Perivascular'],
                 branch="Perichondrial", n_features = 20, fontsize = 5 ,basis ='umap',
                 plot_emb=True,ordering="max",return_genes=True, save = "Perichondrial.pdf" )


# In[26]:


sc.set_figure_params(figsize=(15,15),frameon=False,dpi_save=300)
g1=scf.pl.trends(adata,
                 root_milestone="Root",
                 milestones= ['Bone','Stromal','Cartilage','Perichondrial','Perivascular'],
                 branch="Perivascular", n_features = 20, fontsize = 5 ,basis ='umap',
                 plot_emb=True,ordering="max",return_genes=True, save = "Perivascular.pdf" )


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[148]:


g1=scf.pl.trends(adata,
                 root_milestone="Root",
                 milestones=["m170", "m210"],
                 branch="m170",
                 plot_emb=True,ordering="max",return_genes=True)


# In[136]:


g1=scf.pl.trends(adata,
                 root_milestone="Root",
                 milestones=["m170", "m210"],
                 branch="m210",
                 plot_emb=True,ordering="max",return_genes=True)


# In[ ]:





# In[ ]:





# In[40]:


milestones = list(adata.obs["seg"].cat.categories)


# In[41]:


import matplotlib.pyplot as plt

milestones = list(adata.obs["seg"].cat.categories)

cmap = plt.get_cmap("tab20")
colors_array = ['#%02x%02x%02x' % tuple(int(255*x) for x in cmap(i)[:3]) for i in range(len(milestones))]

seg_colors = dict(zip(milestones, colors_array))

adata.uns['seg_colors'] = seg_colors

print(seg_colors)


# In[ ]:





# In[ ]:





# In[10]:


adata


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[32]:


scf.__version__


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[13]:


milestones_of_interest = ["bif-103", "bif-107", "bif-145", "bif-205", "bif-230", "bif-268", "bif-282", "Root", "Bif-6", "Bif-62"]

# Subset the AnnData object to cells that are in any of the selected milestones
milestone_cells = adata[adata.obs["milestones"].isin(milestones_of_interest)]


# In[14]:


# Find duplicated index values (cell names)
duplicates = pd.Series(milestone_cells.obs_names).duplicated(keep=False)

# Show only the duplicated ones
duplicated_cells = milestone_cells.obs_names[duplicates]

# Print the duplicated cell names
print(duplicated_cells)


# In[15]:


milestone_cells


# In[ ]:





# In[ ]:





# In[ ]:




