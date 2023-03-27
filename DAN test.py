import scanpy as sc
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

# adata = sc.read('/Users/lily/Compute Canada/DAN.h5ad')
# print(adata)
#
# b  = sc.read('/Users/lily/PycharmProjects/scRNA/DAN.h5ad')
# print(b)

print("DAN related")
DAN_PC_all_genes = sc.read('/Users/lily/PycharmProjects/scRNA/DAN_PC_all_genes.h5ad')
print(DAN_PC_all_genes.shape)
# DAN_no_PC_HVGs = sc.read('/Users/lily/PycharmProjects/scRNA/DAN_no_PC_HVGs.h5ad')
# print(DAN_no_PC_HVGs)
# DAN_PC_HVGs = sc.read('/Users/lily/PycharmProjects/scRNA/DAN_PC_HVGs.h5ad')
# print(DAN_PC_HVGs)
# print(DAN_PC_all_genes.obsm['X_pca'].X.shape)
#
# if 'X_pca' in DAN_PC_all_genes.obsm:
#
#     X_train = DAN_PC_all_genes.obsm['X_pca'].X
#     print('reached here', X_train, '\n','nex')
# else:
#     print('reached here istead!')
#     X_train = DAN_PC_all_genes.X

m2 =  sc.read('/Users/lily/PycharmProjects/scRNA/PC_HVGs.h5ad')
m3 = sc.read('/Users/lily/PycharmProjects/scRNA/no_PC_HVGs.h5ad')
print(m2)
print(m3)