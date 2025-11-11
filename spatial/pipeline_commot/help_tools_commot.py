import anndata
import commot as ct
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

os.environ['R_HOME'] = '/home/felix/miniconda3/envs/commot/lib/R'
os.environ['R_USER'] = '/home/felix/miniconda3/envs/commot/lib/python3.7/site-packages/rpy2'


def communication_deg_detection_modified(
    adata: anndata.AnnData,
    n_var_genes: int = None,
    var_genes = None,
    database_name: str = None,
    pathway_name: str = None,
    summary: str = 'receiver',
    lr_pair: tuple = ('total','total'),
    nknots: int = 6,
    n_deg_genes: int = None,
    n_points: int = 50,
    deg_pvalue_cutoff: float = 0.05,
    batch=None,
    pathway=None
):
    """
    Identify signaling dependent genes
    Slightly modified version from the original which solved some crashes.

    This function depends on tradeSeq [Van_den_Berge2020]_. Currently, tradeSeq version 1.0.1 with R version 3.6.3 has been tested to work.
    For the R-python interface, rpy2==3.4.2 and anndata2ri==1.0.6 have been tested to work.

    Here, the total received or sent signal for the spots are considered as a "gene expression" where tradeSeq is used to find the correlated genes.

    Parameters
    ----------
    adata
        The data matrix of shape ``n_obs`` × ``n_var``.
        Rows correspond to cells or positions and columns to genes.
        The count data should be available through adata.layers['count'].
        For example, when examining the received signal through the ligand-receptor pair "ligA" and "RecA" infered with the LR database "databaseX", 
        the signaling inference result should be available in 
        ``adata.obsm['commot-databaseX-sum-receiver']['r-ligA-recA']``
    n_var_genes
        The number of most variable genes to test.
    var_genes
        The genes to test. n_var_genes will be ignored if given.
    n_deg_genes
        The number of top deg genes to evaluate yhat.
    pathway_name
        Name of the signaling pathway (choose from the third column of ``.uns['commot-databaseX-info']['df_ligrec']``).
        If ``pathway_name`` is specified, ``lr_pair`` will be ignored.
    summary
        'sender' or 'receiver'
    lr_pair
        A tuple of the ligand-receptor pair.
        If ``pathway_name`` is specified, ``lr_pair`` will be ignored.
    nknots
        Number of knots in spline when constructing GAM.
    n_points
        Number of points on which to evaluate the fitted GAM 
        for downstream clustering and visualization.
    deg_pvalue_cutoff
        The p-value cutoff of genes for obtaining the fitted gene expression patterns.

    Returns
    -------
    df_deg: pd.DataFrame
        A data frame of deg analysis results, including Wald statistics, degree of freedom, and p-value.
    df_yhat: pd.DataFrame
        A data frame of smoothed gene expression values.
    
    References
    ----------

    .. [Van_den_Berge2020] Van den Berge, K., Roux de Bézieux, H., Street, K., Saelens, W., Cannoodt, R., Saeys, Y., ... & Clement, L. (2020). 
        Trajectory-based differential expression analysis for single-cell sequencing data. Nature communications, 11(1), 1-13.

    """
    # setup R environment
    # !!! anndata2ri works only with 3.6.3 on the tested machine
    import rpy2
    import anndata2ri
    import rpy2.robjects as ro
    from rpy2.robjects.conversion import localconverter
    import rpy2.rinterface_lib.callbacks
    import logging
    rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)

    ro.r('library(tradeSeq)')
    ro.r('library(clusterExperiment)')
    anndata2ri.activate()
    ro.numpy2ri.activate()
    ro.pandas2ri.activate()

    # prepare input adata for R
    adata_deg = anndata.AnnData(
        X = adata.layers['counts'],
        var = pd.DataFrame(index=list(adata.var_names)),
        obs = pd.DataFrame(index=list(adata.obs_names)))
    adata_deg_var = adata_deg.copy()
    sc.pp.filter_genes(adata_deg_var, min_cells=3)
    sc.pp.filter_genes(adata_deg, min_cells=3)
    sc.pp.normalize_total(adata_deg_var, target_sum=1e4)
    sc.pp.log1p(adata_deg_var)
    if n_var_genes is None:
        sc.pp.highly_variable_genes(adata_deg_var, min_mean=0.0125, max_mean=3, min_disp=0.5)
    elif not n_var_genes is None:
        sc.pp.highly_variable_genes(adata_deg_var, n_top_genes=n_var_genes)
    if var_genes is None:
        adata_deg = adata_deg[:, adata_deg_var.var.highly_variable]
    else:
        adata_deg = adata_deg[:, var_genes]
    del adata_deg_var

    summary_name = 'commot-'+database_name+'-sum-'+summary
    if summary == 'sender':
        summary_abrv = 's'
    else:
        summary_abrv = 'r'
    if not pathway_name is None:
        comm_sum = adata.obsm[summary_name][summary_abrv+'-'+pathway_name].values.reshape(-1,1)
    elif pathway_name is None:
        comm_sum = adata.obsm[summary_name][summary_abrv+'-'+lr_pair[0]+'-'+lr_pair[1]].values.reshape(-1,1)
    cell_weight = np.ones_like(comm_sum).reshape(-1,1)

    # send adata to R
    adata_r = anndata2ri.py2rpy(adata_deg)
    ro.r.assign("adata", adata_r)
    ro.r("X <- as.matrix( assay( adata, 'X') )")
    ro.r.assign("pseudoTime", comm_sum)
    ro.r.assign("cellWeight", cell_weight)
        
    # perform analysis (tradeSeq-1.0.1 in R-3.6.3)
    string_fitGAM = 'sce <- fitGAM(counts=X, pseudotime=pseudoTime[,1], cellWeights=cellWeight[,1], nknots=%d, verbose=TRUE)' % nknots
    ro.r(string_fitGAM)
    ro.r('assoRes <- data.frame( associationTest(sce, global=FALSE, lineage=TRUE) )')
    ro.r('assoRes[is.nan(assoRes[,"waldStat_1"]),"waldStat_1"] <- 0.0')
    ro.r('assoRes[is.nan(assoRes[,"df_1"]),"df_1"] <- 0.0')
    ro.r('assoRes[is.nan(assoRes[,"pvalue_1"]),"pvalue_1"] <- 1.0')
    
    
    with localconverter(ro.pandas2ri.converter):
        df_assoRes = ro.r['assoRes']
    ro.r('assoRes = assoRes[assoRes[,"pvalue_1"] <= %f,]' % deg_pvalue_cutoff) # This creates NA entries.
    
    ro.r('assoRes <- na.omit(assoRes)') # Removes NAN entries
    
    ro.r('oAsso <- order(assoRes[,"waldStat_1"], decreasing=TRUE)')
    if n_deg_genes is None:
        n_deg_genes = df_assoRes.shape[0]
        
    # The code below is modified from the original. 
    with localconverter(ro.pandas2ri.converter):
        oAsso = ro.r['oAsso']
    
    if oAsso.shape[0] <= 2:
        print(f"Actual number of n_deg is low")
        with open(f"ALERTS/{batch}_{pathway}.txt", "w") as text_file:
            text_file.write(f"Sample: {batch} with pathway {pathway} had {oAsso.shape[0]} deg genes")
        
        df_deg = df_assoRes.rename(columns={'waldStat_1':'waldStat', 'df_1':'df', 'pvalue_1':'pvalue'})
        idx = np.argsort(-df_deg['waldStat'].values)
        df_deg = df_deg.iloc[idx]
        df_deg = df_assoRes.rename(columns={'waldStat_1':'waldStat', 'df_1':'df', 'pvalue_1':'pvalue'})
        
        return df_deg, None
        
        
    
    if oAsso.shape[0] < n_points:
        n_points = oAsso.shape[0]
        print(f"Set n_points to {n_points}")


        # Reducing nredcuedDims version
        string_cluster = 'clusPat <- clusterExpressionPatterns(sce, nPoints = %d,' % n_points\
        + 'verbose=TRUE, genes = rownames(assoRes)[oAsso][1:min(%d,length(oAsso))],' % n_deg_genes \
        + ' nReducedDims=%d,' % int(n_points/2) \
        + ' minSizes = 2,' \
        + ' k0s=4:5, alphas=c(0.1))'


    else:
        string_cluster = 'clusPat <- clusterExpressionPatterns(sce, nPoints = %d,' % n_points\
        + 'verbose=TRUE, genes = rownames(assoRes)[oAsso][1:min(%d,length(oAsso))],' % n_deg_genes \
        + ' k0s=4:5, alphas=c(0.1))'
    
    
    ro.r(string_cluster)
    ro.r('yhatScaled <- data.frame(clusPat$yhatScaled)')
    with localconverter(ro.pandas2ri.converter):
        yhat_scaled = ro.r['yhatScaled']

    df_deg = df_assoRes.rename(columns={'waldStat_1':'waldStat', 'df_1':'df', 'pvalue_1':'pvalue'})
    idx = np.argsort(-df_deg['waldStat'].values)
    df_deg = df_deg.iloc[idx]
    df_yhat = yhat_scaled

    anndata2ri.deactivate()
    ro.numpy2ri.deactivate()
    ro.pandas2ri.deactivate()

    return df_deg, df_yhat