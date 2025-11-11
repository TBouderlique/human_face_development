# Argument 1: pathway
# Argument 2: sample id ex: 0, 1, 2
# Use commot env. 

# Imports
import pickle
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
import anndata as ad
import tqdm as tqdm
import json
import numpy as np
import os
import statistics
from scipy import stats
from collections import Counter
from scipy.sparse import csr_matrix
import commot as ct
import sys

DATA_PATH='/home/felix/projects/facial/felix/data/reprocessed_data/'
PLOT_PATH="/home/felix/projects/facial/felix/plots"


db="secreted_signaling"

samples = ["A02993E1",
           "A02989E1",
           "A02994D6",
           "A02994E6",
           "C01939A4",
           "C01939A5",
           "C01939A6",
           "C01939B2", 
          ]

sample=samples[int(sys.argv[2])]
pathway=sys.argv[1]

adata = sc.read_h5ad(f"{DATA_PATH}/annotated_stereo/{sample}_with_uncertain_50.h5ad")

color_dict = {
    'chondrocytes': '#FF0000', 
    'dermal': '#00FF00',  
    'endothelia': '#0000FF', 
    'epithelia': '#000000', 
    'erythrocytes': '#9c3a2c',  
    'glia': '#f6fa00', 
    'immune': '#2f6e30',  
    'interneurons': '#00ffee',  
    'melanocytes': '#ff00d0', 
    'mesenchyme A(near-dermal)': '#4b6e00',  
    'mesenchyme B (middle)': '#015361', 
    'mesenchyme C (ebf)': '#d12888', 
    'mesenchyme D (young)': '#f2a7a7', 
    'mesenchyme E (fox)': '#b4a3ff', 
    'mesenchyme F (middle)': '#adffca',
    'mesengyme G (wierd)': '#ebdada', 
    'myoblasts': '#8b3ef0', 
    'neruons': '#dbfc8d', 
    'neural progenitors': '#ff8400',  
    'osteogenic': '#007862', 
    'perivascular': '#85661d',  
    'photoreceptors': '#5c0000',  
    'Uncertain': '#968d8d', 
}


felix_data_path="/home/felix/projects/facial/felix/data/reprocessed_data/"
adata_dis500 = sc.read_h5ad(f"{felix_data_path}{sample}_adata_ccc_dis500_secreted_signaling.h5ad") # Switch here


ct.tl.communication_direction(adata_dis500, database_name='cellchat', pathway_name=pathway, k=5)

# Import directly commot
from typing import Optional, Union
import ot
import sys
import anndata
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D
import plotly
import seaborn as sns
from scipy import sparse
from scipy.spatial import distance_matrix
from scipy.stats import spearmanr, pearsonr
from sklearn.preprocessing import normalize
from sklearn.decomposition import PCA

from commot._utils import plot_cluster_signaling_chord
from commot._utils import plot_cluster_signaling_network
from commot._utils import plot_cluster_signaling_network_multipair
from commot._utils import plot_cell_signaling
from commot._utils import plot_cell_signaling_compare
from commot._utils import get_cmap_qualitative
from commot._utils import leiden_clustering
def plot_cell_communication(
    adata: anndata.AnnData,
    database_name: str = None,
    pathway_name: str = None,
    lr_pair = None,
    keys = None,
    plot_method: str = "cell",
    background: str = "summary",
    background_legend: bool=False,
    clustering: str = None,
    summary: str = "sender",
    cmap: str = "coolwarm",
    cluster_cmap: dict = None,
    pos_idx: np.ndarray = np.array([0,1],int),
    ndsize: float = 1,
    scale: float = 1.0,
    normalize_v: bool = False,
    normalize_v_quantile: float = 0.95,
    arrow_color: str = "#333333",
    grid_density: float = 1.0,
    grid_knn: int = None,
    grid_scale: float = 1.0,
    grid_thresh: float = 1.0,
    grid_width: float = 0.005,
    stream_density: float = 1.0,
    stream_linewidth: float = 1,
    stream_cutoff_perc: float = 5,
    filename: str = None,
    ax: Optional[mpl.axes.Axes] = None
):
    """
    Plot spatial directions of cell-cell communication.
    
    .. image:: cell_communication.png
        :width: 500pt

    The cell-cell communication should have been computed by the function :func:`commot.tl.spatial_communication`.
    The cell-cell communication direction should have been computed by the function :func:`commot.tl.communication_direction`.

    Parameters
    ----------
    adata
        The data matrix of shape ``n_obs`` × ``n_var``.
        Rows correspond to cells or positions and columns to genes.
    database_name
        Name of the ligand-receptor interaction database. 
    pathway_name
        Name of the signaling pathway to be plotted. Will be used only when ``lr_pair`` and ``keys`` are None.
        If none of ``pathway_name``, ``lr_pair``, and ``keys`` are given, the total signaling through all pairs will be plotted.
    lr_pair
        A tuple of ligand name and receptor name. If given, ``pathway_name`` will be ignored. Will be ignored if ``keys`` is given.
    keys
        A list of keys for example 'ligA-recA' for a LR pair or 'pathwayX' for a signaling pathway 'pathwayX'. 
        If given, pathway_name and lr_pair will be ignored.
        If more than one is given, the average will be plotted.
    plot_method
        'cell' plot vectors on individual cells. 
        'grid' plot interpolated vectors on regular grids.
        'stream' streamline plot.
    background
        'summary': scatter plot with color representing total sent or received signal.
        'image': the image in Visium data.
        'cluster': scatter plot with color representing cell clusters.
    background_legend
        Whether to include the background legend when background is set to `summary` or `cluster`.
    clustering
        The key for clustering result. Needed if background is set to `cluster`.
        For example, if ``clustering=='leiden'``, the clustering result should be available in ``.obs['leiden']``.
    summary
        If background is set to 'summary', the numerical value to plot for background.
        'sender': node color represents sender weight.
        'receiver': node color represents receiver weight.
    cmap
        matplotlib colormap name for node summary if numerical (background set to 'summary'), e.g., 'coolwarm'.
        plotly colormap name for node color if summary is (background set to 'cluster'). e.g., 'Alphabet'.
    cluster_cmap
        A dictionary that maps cluster names to colors when setting background to 'cluster'. If given, ``cmap`` will be ignored.
    pos_idx
        The coordinates to use for plotting (2D plot).
    ndsize
        The node size of the spots.
    scale
        The scale parameter passed to the matplotlib quiver function :func:`matplotlib.pyplot.quiver` for vector field plots.
        The smaller the value, the longer the arrows.
    normalize_v
        Whether the normalize the vector field to uniform lengths to highlight the directions without showing magnitudes.
    normalize_v_quantile
        The vector length quantile to use to normalize the vector field.
    arrow_color
        The color of the arrows.
    grid_density
        The density of grid if ``plot_method=='grid'``.
    grid_knn
        If ``plot_method=='grid'``, the number of nearest neighbors to interpolate the signaling directions from spots to grid points.
    grid_scale
        The scale parameter (relative to grid size) for the kernel function of mapping directions of spots to grid points.
    grid_thresh
        The threshold of interpolation weights determining whether to include a grid point. A smaller value gives a tighter cover of the tissue by the grid points.
    grid_width
        The value passed to the ``width`` parameter of the function :func:`matplotlib.pyplot.quiver` when ``plot_method=='grid'``.
    stream_density
        The density of stream lines passed to the ``density`` parameter of the function :func:`matplotlib.pyplot.streamplot` when ``plot_method=='stream'``.
    stream_linewidth
        The width of stream lines passed to the ``linewidth`` parameter of the function :func:`matplotlib.pyplot.streamplot` when ``plot_method=='stream'``.
    stream_cutoff_perc
        The quantile cutoff to ignore the weak vectors. Default to 5 that the vectors shorter than the 5% quantile will not be plotted.
    filename
        If given, save to the filename. For example 'ccc_direction.pdf'.
    ax
        An existing matplotlib ax (`matplotlib.axis.Axis`).

    Returns
    -------
    ax : matplotlib.axis.Axis
        The matplotlib ax object of the plot.

    """

    if not keys is None:
        ncell = adata.shape[0]
        V = np.zeros([ncell, 2], float)
        signal_sum = np.zeros([ncell], float)
        for key in keys:
            if summary == 'sender':
                V = V + adata.obsm['commot_sender_vf-'+database_name+'-'+key][:,pos_idx]
                signal_sum = signal_sum + adata.obsm['commot-'+database_name+"-sum-sender"]['s-'+key]
            elif summary == 'receiver':
                V = V + adata.obsm['commot_receiver_vf-'+database_name+'-'+key][:,pos_idx]
                signal_sum = signal_sum + adata.obsm['commot-'+database_name+"-sum-receiver"]['r-'+key]
        V = V / float( len( keys ) )
        signal_sum = signal_sum / float( len( keys ) )
    elif keys is None:
        if not lr_pair is None:
            vf_name = database_name+'-'+lr_pair[0]+'-'+lr_pair[1]
            sum_name = lr_pair[0]+'-'+lr_pair[1]
        elif not pathway_name is None:
            vf_name = database_name+'-'+pathway_name
            sum_name = pathway_name
        else:
            vf_name = database_name+'-total-total'
            sum_name = 'total-total'
        if summary == 'sender':
            V = adata.obsm['commot_sender_vf-'+vf_name][:,pos_idx]
            signal_sum = adata.obsm['commot-'+database_name+"-sum-sender"]['s-'+sum_name]
        elif summary == 'receiver':
            V = adata.obsm['commot_receiver_vf-'+vf_name][:,pos_idx]
            signal_sum = adata.obsm['commot-'+database_name+"-sum-receiver"]['r-'+sum_name]

    if background=='cluster' and not cmap in ['Plotly','Light24','Dark24','Alphabet']:
        cmap='Alphabet'
    if ax is None:
        fig, ax = plt.subplots()
    if normalize_v:
        V = V / np.quantile(np.linalg.norm(V, axis=1), normalize_v_quantile)
    plot_cell_signaling(
        adata.obsm["spatial"][:,pos_idx],
        V,
        signal_sum,
        cmap = cmap,
        cluster_cmap = cluster_cmap,
        plot_method = plot_method,
        background = background,
        clustering = clustering,
        background_legend = background_legend,
        adata = adata,
        summary = summary,
        scale = scale,
        ndsize = ndsize,
        filename = filename,
        arrow_color = arrow_color,
        grid_density = grid_density,
        grid_knn = grid_knn,
        grid_scale = grid_scale,
        grid_thresh = grid_thresh,
        grid_width = grid_width,
        stream_density = stream_density,
        stream_linewidth = stream_linewidth,
        stream_cutoff_perc = stream_cutoff_perc,
        ax = ax,
        # fig = fig,
    )
    return ax

# Transfer labels of adata to adata_dis500
adata_dis500.obs["labels"] = adata.obs.labels

adata_dis500.obs["labels"].unique().tolist()

# Just set mesenchymes to one mesenchyme
map_dict = {}
for cell in tqdm.tqdm(adata_dis500.obs_names):
    if "mesen" in adata_dis500[cell,].obs.labels.item():
        map_dict[cell] = "mesenchyme"
    else:
        map_dict[cell] = adata_dis500[cell,].obs.labels.item()
        
        
adata_dis500.obs["simplified_labels"] = adata_dis500.obs_names.map(map_dict)

# Select certain labels
colors = ['red', 'blue', 'green', 'orange']
selected_labels = {
    "neruons":"red",
    "osteogenic":"brown",
    "chondrocytes": "blue",     
    "neural progenitors":"green" ,
    "interneurons":"purple",
    "mesenchyme":"orange",
    "other cells":"lightgray"
}


map_dict = {}
for cell in tqdm.tqdm(adata_dis500.obs_names):
    if adata_dis500[cell,].obs.simplified_labels.item() in selected_labels:
        map_dict[cell] = adata_dis500[cell,].obs.simplified_labels.item()
       
    else:
        map_dict[cell] = "other cells"

adata_dis500.obs["selected_labels"] = adata_dis500.obs_names.map(map_dict)

plt.figure(figsize=(30, 8))


ax = plot_cell_communication(adata_dis500,
                                      database_name='cellchat',
                                      pathway_name=pathway,
                                      plot_method='grid', # cells, stream, grid
                                      background_legend=True,
                                      scale=0.00003,
                                      ndsize=3, # used to be 8
                                      grid_density=0.6, # Set number of points
                                      summary='sender', # Or receiver
                                      background='cluster',
                                      clustering='selected_labels',
                                      cluster_cmap=selected_labels,
                                      normalize_v = True,
                                      normalize_v_quantile=0.995,
                                      arrow_color="black",
                                     )
plt.title(f"{pathway} signaling in {sample}")
plt.savefig(f"{pathway}_signaling_in_{sample}.svg")