# Annotate the stereoseq data by label transfer.

We use celltypist to transfer the data from the full RNA-seq dataset annotations. Use the celltypist.yml conda environemnt for this. 


### 1. First train different models using the first file and change the input as described in the python file.


```
max_iter_1=int(sys.argv[1]) # Original 5
max_iter_2=int(sys.argv[2]) # Original 100
n_genes=int(sys.argv[3]) # Original 100
save_suffix=sys.argv[4]
```

### 2. Compare the different models using the second file and look at the resulting txt file.

### 3. Set the correct models in the third file and transfer annotations.

### 4. Set uncertain annotations using the notebook in the end. 



