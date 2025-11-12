#!/usr/bin/env python
# coding: utf-8

# In[1]:


import cellrank
import scvelo as scv
import os
os.chdir("/home/yakov/face_data")
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


adata_ear_11p5 = scp.read_h5ad('DATA/ear_11p5.h5ad')
adata_ear_6_5 = scp.read_h5ad('DATA/ear_6_5.h5ad')
adata_ear_9 = scp.read_h5ad('DATA/ear_9.h5ad')
adata_eyes_10_5 = scp.read_h5ad('DATA/eyes_10_5.h5ad')
adata_eye_11 = scp.read_h5ad('DATA/eye_11.h5ad')
adata_eye_11p5 = scp.read_h5ad('DATA/eye_11p5.h5ad')
adata_eye_6_5 = scp.read_h5ad('DATA/eye_6_5.h5ad')
adata_eye_7 = scp.read_h5ad('DATA/eye_7.h5ad')
adata_eye_9 = scp.read_h5ad('DATA/eye_9.h5ad')
adata_F1_6_5 = scp.read_h5ad('DATA/F1_6_5.h5ad')
adata_F2_6_5 = scp.read_h5ad('DATA/F2_6_5.h5ad')
adata_half_11_5 = scp.read_h5ad('DATA/half_11_5.h5ad')
adata_half_9_5 = scp.read_h5ad('DATA/half_9_5.h5ad')
adata_jaw_10_5 = scp.read_h5ad('DATA/jaw_10_5.h5ad')
adata_jaw_11 = scp.read_h5ad('DATA/jaw_11.h5ad')
adata_jaw_11p5 = scp.read_h5ad('DATA/jaw_11p5.h5ad')
adata_jaw_6_5 = scp.read_h5ad('DATA/jaw_6_5.h5ad')
adata_jaw_7 = scp.read_h5ad('DATA/jaw_7.h5ad')
adata_jaw_9 = scp.read_h5ad('DATA/jaw_9.h5ad')
adata_nose_10_5 = scp.read_h5ad('DATA/nose_10_5.h5ad')
adata_nose_11 = scp.read_h5ad('DATA/nose_11.h5ad')
adata_nose_11p5 = scp.read_h5ad('DATA/nose_11p5.h5ad')
adata_nose_6_5 = scp.read_h5ad('DATA/nose_6_5.h5ad')
adata_nose_7 = scp.read_h5ad('DATA/nose_7.h5ad')
adata_nose_9 = scp.read_h5ad('DATA/nose_9.h5ad')
adata_rest_11 = scp.read_h5ad('DATA/rest_11.h5ad')
adata_rest_7 = scp.read_h5ad('DATA/rest_7.h5ad')
adata_suture_10_5 = scp.read_h5ad('DATA/suture_10_5.h5ad')
adata_twin1_7_5 = scp.read_h5ad('DATA/twin1_7_5.h5ad')
adata_twin2_7_5 = scp.read_h5ad('DATA/twin2_7_5.h5ad')
adata_sn_eye_10 = scp.read_h5ad('DATA/single-nucleus/sn_eye_10.h5ad')
adata_sn_eye_11 = scp.read_h5ad('DATA/single-nucleus/sn_eye_11.h5ad')
adata_sn_eye_6 = scp.read_h5ad('DATA/single-nucleus/sn_eye_6.h5ad')
adata_sn_jaw_10 = scp.read_h5ad('DATA/single-nucleus/sn_jaw_10.h5ad')
adata_sn_jaw_11 = scp.read_h5ad('DATA/single-nucleus/sn_jaw_11.h5ad')
adata_sn_jaw_6 = scp.read_h5ad('DATA/single-nucleus/sn_jaw_6.h5ad')
adata_sn_nose_10 = scp.read_h5ad('DATA/single-nucleus/sn_nose_10.h5ad')
adata_sn_nose_11 = scp.read_h5ad('DATA/single-nucleus/sn_nose_11.h5ad')
adata_sn_nose_6 = scp.read_h5ad('DATA/single-nucleus/sn_nose_6.h5ad')
adata_sn_rest_11 = scp.read_h5ad('DATA/single-nucleus/sn_rest_11.h5ad')
adata_sn_anterior_6_5 = scp.read_h5ad('DATA/sn_anterior_6_5.h5ad')
adata_sn_posterior_6_5 = scp.read_h5ad('DATA/sn_posterior_6_5.h5ad')


# In[5]:


adata_ear_11p5.obs["Sample"] = "ear_11p5"
adata_ear_6_5.obs["Sample"] = "ear_6_5"
adata_ear_9.obs["Sample"] = "ear_9"
adata_eyes_10_5.obs["Sample"] = "eyes_10_5"
adata_eye_11.obs["Sample"] = "eye_11"
adata_eye_11p5.obs["Sample"] = "eye_11p5"
adata_eye_6_5.obs["Sample"] = "eye_6_5"
adata_eye_7.obs["Sample"] = "eye_7"
adata_eye_9.obs["Sample"] = "eye_9"
adata_F1_6_5.obs["Sample"] = "F1_6_5"
adata_F2_6_5.obs["Sample"] = "F2_6_5"
adata_half_11_5.obs["Sample"] = "half_11_5"
adata_half_9_5.obs["Sample"] = "half_9_5"
adata_jaw_10_5.obs["Sample"] = "jaw_10_5"
adata_jaw_11.obs["Sample"] = "jaw_11"
adata_jaw_11p5.obs["Sample"] = "jaw_11p5"
adata_jaw_6_5.obs["Sample"] = "jaw_6_5"
adata_jaw_7.obs["Sample"] = "jaw_7"
adata_jaw_9.obs["Sample"] = "jaw_9"
adata_nose_10_5.obs["Sample"] = "nose_10_5"
adata_nose_11.obs["Sample"] = "nose_11"
adata_nose_11p5.obs["Sample"] = "nose_11p5"
adata_nose_6_5.obs["Sample"] = "nose_6_5"
adata_nose_7.obs["Sample"] = "nose_7"
adata_nose_9.obs["Sample"] = "nose_9"
adata_rest_11.obs["Sample"] = "rest_11"
adata_rest_7.obs["Sample"] = "rest_7"
adata_suture_10_5.obs["Sample"] = "suture_10_5"
adata_twin1_7_5.obs["Sample"] = "twin1_7_5"
adata_twin2_7_5.obs["Sample"] = "twin2_7_5"


adata_sn_eye_10.obs["Sample"] = "eye_10"
adata_sn_eye_11.obs["Sample"] = "eye_11"
adata_sn_eye_6.obs["Sample"] = "eye_6"
adata_sn_jaw_10.obs["Sample"] = "jaw_10"
adata_sn_jaw_11.obs["Sample"] = "jaw_11"
adata_sn_jaw_6.obs["Sample"] = "jaw_6"
adata_sn_nose_10.obs["Sample"] = "nose_10"
adata_sn_nose_11.obs["Sample"] = "nose_11"
adata_sn_nose_6.obs["Sample"] = "nose_6"
adata_sn_rest_11.obs["Sample"] = "rest_11"
adata_sn_anterior_6_5.obs["Sample"] = "anterior_6_5"
adata_sn_posterior_6_5.obs["Sample"] = "posterior_6_5"

adata_ear_11p5.obs["Sample_type"] = "new"
adata_ear_6_5.obs["Sample_type"] = "old"
adata_ear_9.obs["Sample_type"] = "mid"
adata_eyes_10_5.obs["Sample_type"] = "old"
adata_eye_11.obs["Sample_type"] = "new"
adata_eye_11p5.obs["Sample_type"] = "mid"
adata_eye_6_5.obs["Sample_type"] = "old"
adata_eye_7.obs["Sample_type"] = "new"
adata_eye_9.obs["Sample_type"] = "mid"
adata_F1_6_5.obs["Sample_type"] = "new"
adata_F2_6_5.obs["Sample_type"] = "mid"
adata_half_11_5.obs["Sample_type"] = "old"
adata_half_9_5.obs["Sample_type"] = "old"
adata_twin1_7_5.obs["Sample_type"] = "old"
adata_twin2_7_5.obs["Sample_type"] = "mid"
adata_rest_11.obs["Sample_type"] = "new"
adata_rest_7.obs["Sample_type"] = "new"
adata_jaw_10_5.obs["Sample_type"] = "old"
adata_jaw_11.obs["Sample_type"] = "new"
adata_jaw_11p5.obs["Sample_type"] = "mid"
adata_jaw_6_5.obs["Sample_type"] = "old"
adata_jaw_7.obs["Sample_type"] = "mid"
adata_jaw_9.obs["Sample_type"] = "new"
adata_nose_10_5.obs["Sample_type"] = "old"
adata_nose_11.obs["Sample_type"] = "new"
adata_nose_11p5.obs["Sample_type"] = "mid"
adata_nose_6_5.obs["Sample_type"] = "old"
adata_nose_7.obs["Sample_type"] = "new"
adata_nose_9.obs["Sample_type"] = "mid"
adata_suture_10_5.obs["Sample_type"] = "old"

adata_sn_eye_10.obs["Sample_type"] = "sn"
adata_sn_eye_11.obs["Sample_type"] = "sn"
adata_sn_eye_6.obs["Sample_type"] = "sn"
adata_sn_jaw_10.obs["Sample_type"] = "sn"
adata_sn_jaw_11.obs["Sample_type"] = "sn"
adata_sn_jaw_6.obs["Sample_type"] = "sn"
adata_sn_nose_10.obs["Sample_type"] = "sn"
adata_sn_nose_11.obs["Sample_type"] = "sn"
adata_sn_nose_6.obs["Sample_type"] = "sn"
adata_sn_rest_11.obs["Sample_type"] = "sn"
adata_sn_anterior_6_5.obs["Sample_type"] = "sn"
adata_sn_posterior_6_5.obs["Sample_type"] = "sn"

'''
adata_ear_11p5.obs["Sample_type"] = "new"
adata_ear_6_5.obs["Sample_type"] = "old"
adata_ear_9.obs["Sample_type"] = "new"
adata_eyes_10_5.obs["Sample_type"] = "old"
adata_eye_11.obs["Sample_type"] = "new"
adata_eye_11p5.obs["Sample_type"] = "new"
adata_eye_6_5.obs["Sample_type"] = "old"
adata_eye_7.obs["Sample_type"] = "new"
adata_eye_9.obs["Sample_type"] = "new"
adata_F1_6_5.obs["Sample_type"] = "new"
adata_F2_6_5.obs["Sample_type"] = "new"
adata_half_11_5.obs["Sample_type"] = "old"
adata_half_9_5.obs["Sample_type"] = "old"
adata_jaw_10_5.obs["Sample_type"] = "old"
adata_jaw_11.obs["Sample_type"] = "new"
adata_jaw_11p5.obs["Sample_type"] = "new"
adata_jaw_6_5.obs["Sample_type"] = "old"
adata_jaw_7.obs["Sample_type"] = "new"
adata_jaw_9.obs["Sample_type"] = "new"
adata_nose_10_5.obs["Sample_type"] = "old"
adata_nose_11.obs["Sample_type"] = "new"
adata_nose_11p5.obs["Sample_type"] = "new"
adata_nose_6_5.obs["Sample_type"] = "old"
adata_nose_7.obs["Sample_type"] = "new"
adata_nose_9.obs["Sample_type"] = "new"
adata_rest_11.obs["Sample_type"] = "new"
adata_rest_7.obs["Sample_type"] = "new"
adata_suture_10_5.obs["Sample_type"] = "old"
adata_twin1_7_5.obs["Sample_type"] = "old"
adata_twin2_7_5.obs["Sample_type"] = "old"

adata_sn_eye_10.obs["Sample_type"] = "sn"
adata_sn_eye_11.obs["Sample_type"] = "sn"
adata_sn_eye_6.obs["Sample_type"] = "sn"
adata_sn_jaw_10.obs["Sample_type"] = "sn"
adata_sn_jaw_11.obs["Sample_type"] = "sn"
adata_sn_jaw_6.obs["Sample_type"] = "sn"
adata_sn_nose_10.obs["Sample_type"] = "sn"
adata_sn_nose_11.obs["Sample_type"] = "sn"
adata_sn_nose_6.obs["Sample_type"] = "sn"
adata_sn_rest_11.obs["Sample_type"] = "sn"
adata_sn_anterior_6_5.obs["Sample_type"] = "sn"
adata_sn_posterior_6_5.obs["Sample_type"] = "sn"
'''

adata_ear_11p5.obs["Age"] = "11"
adata_ear_6_5.obs["Age"] = "6.5"
adata_ear_9.obs["Age"] = "9"
adata_eyes_10_5.obs["Age"] = "10.5"
adata_eye_11.obs["Age"] = "11"
adata_eye_11p5.obs["Age"] = "11"
adata_eye_6_5.obs["Age"] = "6.5"
adata_eye_7.obs["Age"] = "7"
adata_eye_9.obs["Age"] = "9"
adata_F1_6_5.obs["Age"] = "6.5"
adata_F2_6_5.obs["Age"] = "6.5"
adata_half_11_5.obs["Age"] = "11.5"
adata_half_9_5.obs["Age"] = "9.5"
adata_jaw_10_5.obs["Age"] = "10.5"
adata_jaw_11.obs["Age"] = "11"
adata_jaw_11p5.obs["Age"] = "11"
adata_jaw_6_5.obs["Age"] = "6.5"
adata_jaw_7.obs["Age"] = "7"
adata_jaw_9.obs["Age"] = "9"
adata_nose_10_5.obs["Age"] = "10.5"
adata_nose_11.obs["Age"] = "11"
adata_nose_11p5.obs["Age"] = "11"
adata_nose_6_5.obs["Age"] = "6.5"
adata_nose_7.obs["Age"] = "7"
adata_nose_9.obs["Age"] = "9"
adata_rest_11.obs["Age"] = "11"
adata_rest_7.obs["Age"] = "7"
adata_suture_10_5.obs["Age"] = "10.5"
adata_twin1_7_5.obs["Age"] = "7.5"
adata_twin2_7_5.obs["Age"] = "7.5"

adata_sn_eye_10.obs["Age"] = "10"
adata_sn_eye_11.obs["Age"] = "11"
adata_sn_eye_6.obs["Age"] = "6"
adata_sn_jaw_10.obs["Age"] = "10"
adata_sn_jaw_11.obs["Age"] = "11"
adata_sn_jaw_6.obs["Age"] = "6"
adata_sn_nose_10.obs["Age"] = "10"
adata_sn_nose_11.obs["Age"] = "11"
adata_sn_nose_6.obs["Age"] = "6"
adata_sn_rest_11.obs["Age"] = "11"
adata_sn_anterior_6_5.obs["Age"] = "6.5"
adata_sn_posterior_6_5.obs["Age"] = "6.5"

adata_ear_11p5.obs["Origin"] = "ear"
adata_ear_6_5.obs["Origin"] = "ear"
adata_ear_9.obs["Origin"] = "ear"
adata_eyes_10_5.obs["Origin"] = "eye"
adata_eye_11.obs["Origin"] = "eye"
adata_eye_11p5.obs["Origin"] = "eye"
adata_eye_6_5.obs["Origin"] = "eye"
adata_eye_7.obs["Origin"] = "eye"
adata_eye_9.obs["Origin"] = "eye"
adata_F1_6_5.obs["Origin"] = "all"
adata_F2_6_5.obs["Origin"] = "all"
adata_half_11_5.obs["Origin"] = "all"
adata_half_9_5.obs["Origin"] = "all"
adata_jaw_10_5.obs["Origin"] = "jaw"
adata_jaw_11.obs["Origin"] = "jaw"
adata_jaw_11p5.obs["Origin"] = "jaw"
adata_jaw_6_5.obs["Origin"] = "jaw"
adata_jaw_7.obs["Origin"] = "jaw"
adata_jaw_9.obs["Origin"] = "jaw"
adata_nose_10_5.obs["Origin"] = "nose"
adata_nose_11.obs["Origin"] = "nose"
adata_nose_11p5.obs["Origin"] = "nose"
adata_nose_6_5.obs["Origin"] = "nose"
adata_nose_7.obs["Origin"] = "nose"
adata_nose_9.obs["Origin"] = "nose"
adata_rest_11.obs["Origin"] = "rest"
adata_rest_7.obs["Origin"] = "rest"
adata_suture_10_5.obs["Origin"] = "suture"
adata_twin1_7_5.obs["Origin"] = "all"
adata_twin2_7_5.obs["Origin"] = "all"

adata_sn_eye_10.obs["Origin"] = "eye"
adata_sn_eye_11.obs["Origin"] = "eye"
adata_sn_eye_6.obs["Origin"] = "eye"
adata_sn_jaw_10.obs["Origin"] = "jaw"
adata_sn_jaw_11.obs["Origin"] = "jaw"
adata_sn_jaw_6.obs["Origin"] = "jaw"
adata_sn_nose_10.obs["Origin"] = "nose"
adata_sn_nose_11.obs["Origin"] = "nose"
adata_sn_nose_6.obs["Origin"] = "nose"
adata_sn_rest_11.obs["Origin"] = "rest"
adata_sn_anterior_6_5.obs["Origin"] = "anterior"
adata_sn_posterior_6_5.obs["Origin"] = "posterior"


# In[6]:


import pandas as pd

# List all your AnnData objects here — make sure they're loaded in your environment
adatas = [
    adata_ear_11p5, adata_ear_6_5, adata_ear_9, adata_eyes_10_5,
    adata_eye_11, adata_eye_11p5, adata_eye_6_5, adata_eye_7, adata_eye_9,
    adata_F1_6_5, adata_F2_6_5, adata_half_11_5, adata_half_9_5,
    adata_jaw_10_5, adata_jaw_11, adata_jaw_11p5, adata_jaw_6_5, adata_jaw_7,
    adata_jaw_9, adata_nose_10_5, adata_nose_11, adata_nose_11p5, adata_nose_6_5,
    adata_nose_7, adata_nose_9, adata_rest_11, adata_rest_7, adata_suture_10_5,
    adata_twin1_7_5, adata_twin2_7_5,
    adata_sn_eye_10, adata_sn_eye_11, adata_sn_eye_6,
    adata_sn_jaw_10, adata_sn_jaw_11, adata_sn_jaw_6,
    adata_sn_nose_10, adata_sn_nose_11, adata_sn_nose_6,
    adata_sn_rest_11, adata_sn_anterior_6_5, adata_sn_posterior_6_5
]

results = []
for ad in adatas:
    sample = ad.obs["Sample"].iloc[0]
    sample_type = ad.obs["Sample_type"].iloc[0]
    n_genes = ad.n_vars
    results.append({
        "Sample": sample,
        "Sample_type": sample_type,
        "Genes": n_genes
    })

df_gene_counts = pd.DataFrame(results)

# Sort by Sample_type then Sample
df_gene_counts_sorted = df_gene_counts.sort_values(by=["Sample_type", "Sample"]).reset_index(drop=True)

# Display the table
print(df_gene_counts_sorted)


# In[12]:


adata_sn = adata_sn_eye_10.concatenate(adata_sn_eye_11, adata_sn_eye_6,
adata_sn_jaw_10, adata_sn_jaw_11, adata_sn_jaw_6,
adata_sn_nose_10, adata_sn_nose_11, adata_sn_nose_6,
adata_sn_rest_11, adata_sn_anterior_6_5, adata_sn_posterior_6_5, join='outer')


# In[6]:


#adata = adata_ear_11p5.concatenate(adata_ear_6_5,adata_ear_9,adata_eyes_10_5,adata_eye_11,adata_eye_11p5,adata_eye_6_5,adata_eye_7,adata_eye_9,adata_F1_6_5,adata_F2_6_5,adata_half_11_5,adata_half_9_5,adata_jaw_10_5,adata_jaw_11,adata_jaw_11p5,adata_jaw_6_5,adata_jaw_7,adata_jaw_9,adata_nose_10_5,adata_nose_11,adata_nose_11p5,adata_nose_6_5,adata_nose_7,adata_nose_9,adata_rest_11,adata_rest_7,adata_suture_10_5,adata_twin1_7_5,adata_twin2_7_5,adata_sn_eye_10, adata_sn_eye_11, adata_sn_eye_6, adata_sn_jaw_10, adata_sn_jaw_11, adata_sn_jaw_6, adata_sn_nose_10, adata_sn_nose_11, adata_sn_nose_6, adata_sn_rest_11, adata_sn_anterior_6_5, adata_sn_posterior_6_5, join='outer')
adata_old_anterior_posterior = adata_ear_6_5.concatenate(adata_eyes_10_5,adata_eye_6_5,adata_half_11_5,adata_half_9_5,adata_jaw_10_5,adata_jaw_6_5,adata_nose_10_5,adata_nose_6_5,adata_suture_10_5,adata_twin1_7_5,adata_twin2_7_5, adata_sn_anterior_6_5, adata_sn_posterior_6_5, join='outer')
adata_new = adata_ear_11p5.concatenate(adata_ear_9,adata_eye_11,adata_eye_11p5,adata_eye_7,adata_eye_9,adata_F1_6_5,adata_F2_6_5,adata_jaw_11,adata_jaw_11p5,adata_jaw_7,adata_jaw_9,adata_nose_11,adata_nose_11p5,adata_nose_7,adata_nose_9, adata_rest_11,adata_rest_7, adata_sn_eye_10, adata_sn_eye_11, adata_sn_eye_6, adata_sn_jaw_10, adata_sn_jaw_11, adata_sn_jaw_6, adata_sn_nose_10, adata_sn_nose_11, adata_sn_nose_6, adata_sn_rest_11, join='outer')
adata = adata_old_anterior_posterior.concatenate(adata_new, join='inner')


# In[7]:


adata_old_anterior_posterior


# In[8]:


adata_new


# In[9]:


adata


# In[13]:


adata_sn


# In[14]:


import pandas as pd

# Export genes to CSV
genes = pd.DataFrame(adata_sn.var_names, columns=["Gene"])
genes.to_csv("adata_sn_genes.csv", index=False)

print("✅ Saved gene list to adata_sn_genes.csv")


# In[7]:


adata.write('DATA/all_samples-raw-count-matrix.h5ad')


# In[8]:


adata.X.A


# In[7]:


adata.var_names_make_unique()
adata.obs_names_make_unique
    
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
scp.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

fig, axs = plt.subplots(ncols=3, nrows=1, figsize=(22, 4)) # Создаём сетку, чтобы отрисовать три графика

# Отрисовываем распределение числа генов по клеткам
sns.histplot(adata.obs["n_genes_by_counts"], kde=False, ax=axs[1], stat="density")
axs[1].set_xlabel("Genes per cell")
axs[1].set_ylabel("Cells number")
axs[1].set_title("n_genes_by_counts")

    # Отрисовываем распределение общего числа UMI
sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0], stat="density")
axs[0].set_xlabel("UMI per cell, total_counts")
axs[0].set_ylabel("Cells number")
axs[0].set_title("UMI - total_counts")


    # Отрисовываем распределение митохондриальной экспрессии
sns.histplot(adata.obs["pct_counts_mt"], kde=False, ax=axs[2], stat="density")
axs[2].set_xlabel("pct_counts_mt")
axs[2].set_ylabel("Cells number")
axs[2].set_title("pct_counts_mt")

# annotate the group of mitochondrial genes as "mt"
#adata.var["mt"] = adata.var_names.str.startswith("MT-")
#scp.pp.calculate_qc_metrics(
#    adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
#)
#scp.pl.violin(
#    adata,
#    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
#    jitter=0.4,
#    multi_panel=True,
#)

#scp.pl.scatter(adata, x="total_counts", y="pct_counts_mt")
#scp.pl.scatter(adata, x="total_counts", y="n_genes_by_counts")
adata = adata[adata.obs.n_genes_by_counts < 7500, :].copy()
adata = adata[adata.obs.n_genes_by_counts > 250, :].copy()
adata = adata[adata.obs.total_counts_mt < 750, :].copy()
adata = adata[adata.obs.total_counts < 25000, :].copy()
adata = adata[adata.obs.pct_counts_mt < 4.5, :].copy()
adata = adata[adata.obs.doublet_score < 0.15, :].copy()

#scp.pl.scatter(adata, x="total_counts", y="pct_counts_mt")
#scp.pl.scatter(adata, x="total_counts", y="n_genes_by_counts")

#adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
#scp.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

#fig, axs = plt.subplots(ncols=3, nrows=1, figsize=(22, 4)) # Создаём сетку, чтобы отрисовать три графика

# Отрисовываем распределение числа генов по клеткам
#sns.histplot(adata.obs["n_genes_by_counts"], kde=False, ax=axs[1], stat="density")
#axs[1].set_xlabel("Genes per cell")
#axs[1].set_ylabel("Cells number")
#axs[1].set_title("n_genes_by_counts")

    # Отрисовываем распределение общего числа UMI
#sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0], stat="density")
#axs[0].set_xlabel("UMI per cell, total_counts")
#axs[0].set_ylabel("Cells number")
#axs[0].set_title("UMI - total_counts")


    # Отрисовываем распределение митохондриальной экспрессии
#sns.histplot(adata.obs["pct_counts_mt"], kde=False, ax=axs[2], stat="density")
#axs[2].set_xlabel("pct_counts_mt")
#axs[2].set_ylabel("Cells number")
#axs[2].set_title("pct_counts_mt")


#adata.write('DATA/' + file +'.h5ad')

for col in adata.var:
    del adata.var[col]
    
all_genes = [x for x in adata.var_names]
for gene in ['XIST', 'AD16','AIC','APOO','ARMCX6','BEX1','BEX2','BEX4','CCDC120','CCDC22','CD99L2','CDR1-AS','CHRDL1','CMTX2','CMTX3','CT45A5','CT55','CXorf36','CXorf57','CXorf40A','CXorf49','CXorf66','CXorf67','DACH2','EFHC2','ERCC6L','F8A1','FAM104B','FAM120C','FAM122B','FAM122C','FAM127A','FAM50A','FATE1','FMR1-AS1','FRMPD3','FRMPD4','FUNDC1','FUNDC2','GAGE12F','GAGE2A','GATA1','GNL3L','GPRASP2','GRIPAP1','GRDX','HDHD1A','HS6ST2','ITM2A','LAS1L','LINC01420','LOC101059915','MAGEA2','MAGEA5','MAGEA8','MAGED4B','MAGT1','MAGED4','MAP3K15','MBNL3','MBTPS2','MCT-1','MIR106A','MIR222','MIR361','MIR503','MIR6087','MIR660','MIRLET7F2','MORF4L2','MOSPD1','MOSPD2','NAP1L3','NKRF','NRK','OTUD5','PASD1','PAGE1','PAGE2B','PBDC1','PCYT1B','PIN4','PLAC1','PLP2','RPA4','RPS6KA6','RRAGB','RTL3','SFRS17A','SLC38A5','SLITRK2','SMARCA1','SMS','SPANXN1','SPANXN5','SPG16','SSR4','TAF7L','TCEAL1','TCEAL4','TENT5D','TEX11','THOC2','TMEM29','TMEM47','TMLHE','TRAPPC2P1','TREX2','TRO','TSPYL2','TTC3P1','USP51','VSIG1','YIPF6','ZC3H12B','ZCCHC18','ZFP92','ZMYM3','ZNF157','ZNF182','ZNF275','ZNF674', 'SRY','ZFY','RPS4Y1','AMELY','TBL1Y','PCDH11Y','TGIF2LY','TSPY1','AZFa','USP9Y','DDX3Y','UTY','TB4Y','AZFb','CYorf15','RPS4Y2','EIF1AY','KDM5D','XKRY','HSFY1','PRY','RBMY1A1','AZFc','DAZ1','CDY1','VCY1']:
    if gene in all_genes:    
        all_genes.remove(gene)
adata = adata[:, all_genes].copy()

cell_cycle_genes = pd.read_csv('../cell_cycle_genes/all_human_genes.csv')
s_genes = cell_cycle_genes.all_genes[:43]
g2m_genes = cell_cycle_genes.all_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes.all_genes if x in adata.var_names]
scp.pp.filter_cells(adata, min_genes=200)
scp.pp.filter_genes(adata, min_cells=3)
scp.pp.normalize_total(adata, target_sum=1e4)
scp.pp.log1p(adata)
scp.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
scp.pl.highly_variable_genes(adata)
adata.raw = adata.copy()
scp.pp.scale(adata)
scp.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
#adata_cc_genes = adata[:, cell_cycle_genes]
#scp.tl.pca(adata_cc_genes)
#scp.pl.pca_scatter(adata_cc_genes, color='phase')
scp.pp.regress_out(adata, ['S_score', 'G2M_score'])
#scp.pp.scale(adata)
#adata_cc_genes = adata[:, cell_cycle_genes]
#scp.tl.pca(adata_cc_genes)
#scp.pl.pca_scatter(adata_cc_genes, color='phase')
adata.write('DATA/all_samples-no_cell_cycle_260325.h5ad')


# In[8]:


adata = scp.read_h5ad('DATA/all_samples-no_cell_cycle_260325.h5ad')
#adata_old = scp.read_h5ad('DATA/all_old-no_cell_cycle.h5ad')
#adata_new = scp.read_h5ad('DATA/all_new-no_cell_cycle.h5ad')
#adata = adata_old.concatenate(adata_new, join='outer')
#scp.pp.scale(adata, max_value=10)
#adata.raw=adata
scp.tl.pca(adata,n_comps = 30, svd_solver='arpack',use_highly_variable=True)
if adata.raw is not None:
    adata = adata.raw.to_adata()
    del adata.raw
sce.pp.harmony_integrate(adata, 'Sample_type')
scp.pp.neighbors(adata, n_neighbors = 15, n_pcs = 25, use_rep = 'X_pca_harmony' )
scp.tl.leiden(adata, resolution = 0.4)
scp.tl.umap(adata,n_components = 2)
scp.pl.umap(adata, color=['leiden'], save = '-all-harmony4.png')
#for obsp_key in list(adata.obsp.keys()):
#    del adata.obsp[obsp_key]
#for varm_key in list(adata.varm.keys()):
#    del adata.varm[varm_key]
#for uns_key in list(adata.uns.keys()):
#    del adata.uns[uns_key]
#if adata.raw is not None:
#    adata = adata.raw.to_adata()
#    del adata.raw
adata.write('DATA/all_samples-harmony4.h5ad')


# In[11]:


scp.pp.neighbors(adata, n_neighbors = 30, n_pcs = 25, use_rep = 'X_pca_harmony' )
scp.tl.leiden(adata, resolution = 0.4)
scp.tl.umap(adata,n_components = 2)
scp.pl.umap(adata, color=['leiden'], save = '-all-harmony5.png')
#for obsp_key in list(adata.obsp.keys()):
#    del adata.obsp[obsp_key]
#for varm_key in list(adata.varm.keys()):
#    del adata.varm[varm_key]
#for uns_key in list(adata.uns.keys()):
#    del adata.uns[uns_key]
#if adata.raw is not None:
#    adata = adata.raw.to_adata()
#    del adata.raw
adata.write('DATA/all_samples-harmony5.h5ad')


# In[9]:


adata.X


# In[12]:


scp.pl.umap(adata, color=['leiden','HBB','RUNX2','SOX9','SP7','DMP1', 'DSPP','BGLAP'], legend_loc="on data")


# In[ ]:




