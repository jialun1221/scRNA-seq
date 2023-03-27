import numpy as np
import pandas as pd
import scanpy as sc
import sys
import seaborn as sns

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

no_PC = "no_PC_HVGs.h5ad"
PC_in_obsm = 'PC_all_genes.h5ad'
PC_of_HVGs = "PC_HVGs.h5ad"

adata = sc.read('/Users/lily/PycharmProjects/scRNA/new_anndata.h5ad')

def filtering(adata):
    adata.var_names_make_unique()
    sc.pp.filter_cells(adata,min_genes=200)  # Minimum number of genes expressed required for a cell to pass filtering.
    sc.pp.filter_genes(adata,min_cells=3)  # Minimum number of cells expressed required for a gene to pass filtering.
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    sc.pl.highest_expr_genes(adata, n_top=5, )
    return adata

def qc_metrics_plot(adata):
    adata.var["mito"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], inplace=True)
    sns.jointplot(
        data=adata.obs,
        x="log1p_total_counts",
        y="log1p_n_genes_by_counts",
        kind="hex",
    )
    sns.histplot(adata.obs["pct_counts_mito"])

def violin_plot(adata):
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                 jitter=0.4, multi_panel=True)

def scatter_plot(adata):
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')

def mito_genes_filter(adata, num_genes, pct_mt):
    '''filters out the unnecessary genes, keep at 2500 HVGs when wanted'''
    #normally, num_genes = 2500, pct_mt = 5
    adata = adata[adata.obs.n_genes_by_counts < num_genes, :]  # wants all the rows with number of counts <2500
    adata = adata[adata.obs.pct_counts_mt < pct_mt, :]  # also has to have mt<5
    return adata

def normalize(adata):
    sc.pp.normalize_total(adata, target_sum=1e4)
    return adata

def logarithmize(adata):
    sc.pp.log1p(adata)
    return adata

def HVGs(adata):
    #min/max params can be passed into the function
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata.raw = adata
    return adata

def HVGs_2500(adata):
    sc.pp.highly_variable_genes(adata, n_top_genes=2500, min_mean=0.0125, max_mean=3, min_disp=0.5)  # B
    adata.raw = adata
    return adata

def plotHVGs(adata):
    sc.pl.highly_variable_genes(adata)

def filter_by_HVGs(adata):
    adata = adata[:, adata.var.highly_variable]  # 2500 highly variable genes
    print("dimension check:" + str(adata))
    return adata

def regress_out(adata):
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    return adata

def scale(adata):
    #max_value defulat = 10
    sc.pp.scale(adata, max_value=10)
    return adata

def PCA(adata, n_comps, filename):
    #set the number of PCs you want in n_comps
    sc.tl.pca(adata, n_comps, svd_solver='arpack')
    pca = pd.DataFrame(adata.obsm['X_pca']) #assign pca data into X_pca
    #adata = adata.uns['pca']
    print('here')
    print(pca)
    adata.obsm['X_pca'] = sc.AnnData(pca)
    print("pca check: " + str(adata.obs))
    adata.write(filename)
    return adata