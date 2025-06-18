import numpy as np
import pandas as pd
import umap
import argparse
import anndata
import scipy.sparse as sp
import sys
sys.path.append("..")
from src.pystime.read import init_cell_anndata_from_csv
from src.pystime.stime import STIME

def main(adata,species):
    

    stime = STIME(adata=adata,species=species, search_method='neighbor', n_neighs=6)
    idata = stime.idata
    pdata = stime.pdata

    # save data
    mat = idata.X.todense()
    np.savetxt("../tmp/N_M_mat.csv", mat)
    idata.obs.to_csv("../tmp/interface_meta.csv",index=False)
    idata.var.to_csv("../tmp/interface_LR_meta.csv",index=False)
    csr_matrix = sp.csr_matrix(idata.X)
    idata.X = csr_matrix
    idata.X = np.nan_to_num(idata.X)
    idata.layers['MCMF-LRTf'] = sp.csr_matrix(idata.layers['MCMF-LRTf'])
    idata.write("../tmp/idata.h5ad")

if __name__ == "__main__":
    row_col_order ='cell_gene'
    parser = argparse.ArgumentParser()

    parser.add_argument("--mode", type=str, required=True, help="running mode")
    parser.add_argument("--adata", type=str, required=False, help="directory of h5ad file")
    parser.add_argument("--count", type=str, required=False, help="cell gene expression count matrix absolute directory")
    parser.add_argument("--cellmeta", type=str, required=False, help="cellmeta csv file absolute directory")
    parser.add_argument("--genemeta", type=str, required=False,default=None, help="genemeta csv file absolute directory")
    parser.add_argument("--x_col", type=str, required=False,default="x", help="x coordinate of cell")
    parser.add_argument("--y_col", type=str, required=False,default="y", help="y coordinate of cell")
    parser.add_argument("--cell_id_col", type=str, required=False,default="barcode", help="name of cell id column")
    parser.add_argument("--gene_id_col", type=str, required=False,default="gene_id", help="name of gene id column")
    parser.add_argument("--gene_symbol_col", type=str, required=False,default="gene_name", help="name of gene name column")
    parser.add_argument("--species",type=str, required=False,default="Human", help="Human or Mouse")

    args = parser.parse_args()

    mode = args.mode
    if mode != 'csv' and mode != 'adata':
        print("ERROR! Wrong mode.")
        exit(0)
    x_col = args.x_col
    y_col = args.y_col
    cell_id_col = args.cell_id_col
    gene_id_col = args.gene_id_col
    gene_symbol_col = args.gene_symbol_col
    if mode =="csv":
        count_fn = args.count
        cell_meta_fn = args.cellmeta
        if args.genemeta:
            gene_meta_fn = args.genemeta
        else:
            gene_meta_fn = None
        adata = init_cell_anndata_from_csv(count_fn, cell_meta_fn, gene_meta_fn=gene_meta_fn, row_col_order=row_col_order,
                               cell_id_col=cell_id_col, x_col=x_col, y_col=y_col,
                               gene_id_col=gene_id_col, gene_symbol_col=gene_symbol_col,
                               delimiter=','
                              )
    elif mode=="adata":
        # create csv files
        if args.adata:
            adata = anndata.read_h5ad(args.adata)
        else:
            print("ERROR! Need adata dir")
            exit(0)
    main(adata,species=args.species)
    
