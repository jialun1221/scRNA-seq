import numpy as np
import pandas as pd
import scanpy as sc
import sys
from scipy.sparse import csr_matrix
import hdf5plugin
import os

def make_anndata(adata, filename):
    #adata = the data object you want to make the new one from
    #file name = a string with file name ended in .h5ad
    file_obj = filename

    new = sc.AnnData(X=adata.X,
                     obs=adata.obs.copy(),
                     var=adata.var.copy(),
                     uns=adata.uns.copy(),
                     obsm=adata.obsm.copy(),
                     varm=adata.varm.copy(),
                     layers=adata.layers.copy(),
                     raw=adata.raw.copy(),
                     dtype="float32",
                     shape=None,
                     # filename = NULL,
                     # filemode = NULL,
                     obsp=adata.obsp.copy(),
                     varp=adata.varp
                     )

    # A random line that I found necessary for the object to work.
    new.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(
        columns={'_index': 'features'})
    new.write(file_obj)
