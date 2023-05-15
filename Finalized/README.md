##### Please follow this step-by-step guide to preprocess and create objects for different cell types. 

### PREPROCESSING
*Before starting, please upload the .h5ad object into your Google Drive for easier access*

1. Go to **Preprocessing_step1.ipynb**. Read in the desired .h5ad AnnData in this line: `adata = sc.read_h5ad("drive/MyDrive/scRNA ML classifier/data_objects_May_2022/PD_mg.h5ad") `
2. Run **Data Selection** Module. This module drops out cells of LBD.
3. Run **making new AnnData object**. This module creates a processed data.
4. Click the folder label on the leftmost panel --> click "write" --> download the new anndata object you just created, then upload it into your Google Drive. We will use this object for the second preprocessing.
5. Once uploaded, open **Preprocessing_step2.ipynb**. Run all cells.
      - Since we discussed that filtering is not necessary, only `HVGs`, `HVGs_2500` and `filter_by_HVGs` are the necessary functions to run when creating the new objects
      - Depending on the object you want to create, you can customize the functions to run. (for example, method 3 doesn't require PCAs so you don't need to run the `PCA` function. 
6. Once the three objects for three methods are created, download them, then upload to Google Drive. 


### Running Models
Due to the looping nature, you need to modify the path pass in correct .h5ad file for each preprocessing method. 
```
from google.colab import drive
drive.mount('/content/drive')

adata_m1 = sc.read_h5ad("path")
adata_m2 = sc.read_h5ad("path")
adata_m3 = sc.read_h5ad("path")

data_list = [adata_m1, adata_m2, adata_m3]
```
**LOgistic Regression**
- Run all cells, most importantly the mega-loop. Then paste the result dataframe into an excel sheet. 
- Run Create plots to get the confusion matrix. 
**RAndom Forest Classifier**
- Run all cells, most importantly the mega-loop. Then paste the result dataframe into an excel sheet. 
- Run Create plots to get the confusion matrix. 
**DEep Neural Networks**

**SUpport Vector Machine?** 
- Run all cells, most importantly the mega-loop. Then paste the result dataframe into an excel sheet. 


**Missing**
RAndom Forest for both DAN and Astro
heat map plotting script for Bootstrapping 
Bootstrapping for method 1 and method 3
