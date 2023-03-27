import numpy as np
import pandas as pd
import scanpy as sc
import sys
from scipy.sparse import csr_matrix
import hdf5plugin
import os
from numpy import inf

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

def make_new_adata():
  #open a file to store the output AnnData
  file_obj = 'DAN.h5ad'
  adata = sc.read('/Users/lily/Compute Canada/DAN.h5ad')
  print(adata)

  adata.obs = adata.obs.reset_index()  # Set index for the labels
 #k = adata.obs  # create a variable for further uses (a DataFrame)
  lst = adata.obs.index[adata.obs['disease__ontology_label'] == 'Lewy body dementia'].tolist()
  arr = adata.X.toarray()
  arr = np.delete(arr, obj=lst,
                axis=0)  # delete rows that contain Lewy Body Dementia according to the previously generated index stored in y

  # drop command for adata.obs
  adata.obs.drop(adata.obs.index[adata.obs['disease__ontology_label'] == 'Lewy body dementia'], inplace=True)
  #Command for making a new AnnData object. For each parameter, need to make a deep copy of the original object.
  new = sc.AnnData(X = arr,
    obs = adata.obs.copy(),
    var = adata.var.copy(),
    uns = adata.uns.copy(),
    obsm = adata.obsm.copy(),
    varm = adata.varm.copy(),
    layers = adata.layers.copy(),
    raw = adata.raw.copy(),
    dtype = "float32",
    shape = None,
    #filename = NULL,
    #filemode = NULL,
    obsp = adata.obsp.copy(),
    varp = adata.varp
  )

  # A random line that I found necessary for the object to work.
  new.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
  new.write(file_obj)