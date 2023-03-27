from preprocessing import make_new_adata
from preprocessing2 import *

def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press ⌘F8 to toggle the breakpoint.

def make_no_PC_HVGs(new):
    print("before gen filter" + str(new.obs.shape))
    print("before gen filter" + str(new.X.shape))
    new = filtering(new)
    print("after gen filter" + str(new.obs.shape))
    print("after gen filter" + str(new.X.shape))
    # qc_metrics_plot(new)
    # violin_plot(new)
    # scatter_plot(new)
    new = mito_genes_filter(new, 2500, 5)
    print("after mito filter" + str(new.obs.shape))
    print("after mito filter" + str(new.X.shape))

    new = normalize(new)
    new = logarithmize(new)
    new = HVGs_2500(new)
    new = filter_by_HVGs(new)
    new = regress_out(new)
    new = scale(new)
    # Obj created for pure 2500 HVGs
    obj = sc.AnnData(X=new.X.copy(),
                     obs=new.obs.copy(),
                     var=new.var.copy(),
                     uns=new.uns.copy(),
                     obsm=new.obsm.copy(),
                     varm=new.varm.copy(),
                     layers=new.layers.copy(),
                     raw=new.raw.copy(),
                     dtype="float32",
                     shape=None,
                     # filename = NULL,
                     # filemode = NULL,
                     obsp=new.obsp.copy(),
                     varp=new.varp
                     )
    # A random line that I found necessary for the object to work.
    obj.__dict__['_raw'].__dict__['_var'] = new.__dict__['_raw'].__dict__['_var'].rename(
        columns={'_index': 'features'})
    obj.write(no_PC)

def make_PC_all_genes(new):
    print("before gen filter" + str(new.obs.shape))
    print("before gen filter" + str(new.X.shape))
    new = filtering(new)
    print("after gen filter" + str(new.obs.shape))
    print("after gen filter" + str(new.X.shape))
    # qc_metrics_plot(new)
    # violin_plot(new)
    # scatter_plot(new)
    new = mito_genes_filter(new, 2500, 5)
    print("after mito filter" + str(new.obs.shape))
    print("after mito filter" + str(new.X.shape))

    new = normalize(new)
    new = logarithmize(new)
    new = HVGs(new)
    #plotHVGs(new)
    new = regress_out(new)
    new = scale(new)
    new = PCA(new, 50,PC_all_genes)
    obj = sc.AnnData(X=new.X.copy(),
                     obs=new.obs.copy(),
                     var=new.var.copy(),
                     uns=new.uns.copy(),
                     obsm=new.obsm.copy(),
                     varm=new.varm.copy(),
                     layers=new.layers.copy(),
                     #raw=new.raw.copy(),
                     dtype="float32",
                     shape=None,
                     # filename = NULL,
                     # filemode = NULL,
                     obsp=new.obsp.copy(),
                     varp=new.varp
                     )
    # A random line that I found necessary for the object to work.
    #obj.__dict__['_raw'].__dict__['_var'] = new.__dict__['_raw'].__dict__['_var'].rename(
        #columns={'_index': 'features'})
    obj.write(PC_all_genes)

def make_PC_HVGs(new):
    print("before gen filter" + str(new.obs.shape))
    print("before gen filter" + str(new.X.shape))
    new = filtering(new)
    print("after gen filter" + str(new.obs.shape))
    print("after gen filter" + str(new.X.shape))
    # qc_metrics_plot(new)
    # violin_plot(new)
    # scatter_plot(new)
    new = mito_genes_filter(new, 2500, 5)
    # print("after mito filter" + str(new.obs.shape))
    # print("after mito filter" + str(new.X.shape))

    new = normalize(new)
    new = logarithmize(new)
    new = HVGs_2500(new)
    new = filter_by_HVGs(new)
    new = regress_out(new)
    new = scale(new)
    new = PCA(new, 50, PC_HVGs)
    obj = sc.AnnData(X=new.X.copy(),
                     obs=new.obs.copy(),
                     var=new.var.copy(),
                     uns=new.uns.copy(),
                     obsm=new.obsm.copy(),
                     varm=new.varm.copy(),
                     layers=new.layers.copy(),
                     #raw=new.raw.copy(),
                     dtype="float32",
                     shape=None,
                     # filename = NULL,
                     # filemode = NULL,
                     obsp=new.obsp.copy(),
                     varp=new.varp
                     )
    # A random line that I found necessary for the object to work.
    #obj.__dict__['_raw'].__dict__['_var'] = new.__dict__['_raw'].__dict__['_var'].rename(
        #columns={'_index': 'features'})
    obj.write(PC_HVGs)
    return new

if __name__ == '__main__':
    no_PC_HVGs = "no_PC_HVGs.h5ad"
    # PC_all_genes = 'PC_all_genes.h5ad'
    # PC_HVGs = "PC_HVGs.h5ad" #only this working fine
    #make_new_adata()
    new = sc.read('/Users/lily/PycharmProjects/scRNA/new_anndata.h5ad')
    make_no_PC_HVGs(new)
    # make_PC_all_genes(new)
    # make_PC_HVGs(new)
    # print("this is PC_HVGs")
    # print(PC_HVGs)
    '''
    March 19, 2023
    DAN data object
    '''
    #file_obj = 'DAN.h5ad'
    # DAN = sc.read('/Users/lily/PycharmProjects/scRNA/DAN.h5ad')
    # print(DAN)
    # #make_new_adata()
    #
    # no_PC_HVGs = "DAN_no_PC_HVGs.h5ad"
    # PC_all_genes = 'DAN_PC_all_genes.h5ad'
    # PC_HVGs = "DAN_PC_HVGs.h5ad" #only this working fine

    #make_no_PC_HVGs(DAN)
    # make_PC_all_genes(DAN)
    # make_PC_HVGs(DAN)

