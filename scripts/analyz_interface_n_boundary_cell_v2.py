import pdb
import anndata
import pandas as pd
from scipy.sparse import csr_matrix
import os
import sys

sys.path.append("..")
from src.analyz_interface import Analysis_Interface
from src.analyz_boundary_cell import Analysis_Boundary_Cell


def analyz_infc(idata_file, clus_result_df, lr_pathway_db_file, save_root_dir, plot_volcano=False):
    idata = anndata.read_h5ad(idata_file)
    csr_sparse_matrix = csr_matrix(idata.X)
    idata.X = csr_sparse_matrix
    analyz = Analysis_Interface(idata_file=idata, result_file_dir=clus_result_df,
                                lr_pathway_db_file_dir=lr_pathway_db_file, plot_volcano=plot_volcano)

    # Create directories based on the root directory
    lr_enrich_dir = os.path.join(save_root_dir, "lr")
    os.makedirs(lr_enrich_dir, exist_ok=True)
    analyz.proc_lr(save_lr_enrich=True, lr_enrich_dir=lr_enrich_dir)

    pathway_enrich_dir = os.path.join(save_root_dir, "pathway")
    os.makedirs(pathway_enrich_dir, exist_ok=True)
    analyz.proc_pathway(save_pathway_enrich=True, pathway_enrich_dir=pathway_enrich_dir)


def analyz_bc(idata, adata, clus_result_df, lr_pathway_db_file, celltype_key="cell_type", x_row="x", y_col="y", save_df=True, save_df_dir=None, plot_volcano=False, save_pathway_enrich=None, plot_distribution=False, save_distribution_dir=None, save_distribution_type_dir=None, save_gene_exp=False):
    analyz = Analysis_Boundary_Cell(idata_file=idata, adata_file=adata, result_file_dir=clus_result_df, lr_pathway_db_file_dir=lr_pathway_db_file, celltype_key=celltype_key, x_row=x_row, y_col=y_col, save_df=save_df, save_df_dir=save_df_dir, plot_volcano=plot_volcano, save_pathway_enrich=save_pathway_enrich)
    
    if plot_distribution:
        bcf = analyz.bc_finder
        bcf.plot_boundary_cell_cluster(save_figs=True, fig_dir=save_distribution_dir)
        bcf.plot_boundary_cell_celltype(save_figs=True, fig_dir=save_distribution_type_dir)
    
    if save_gene_exp:
        analyz.plot_gene_heatmap(mode="avg")


def main_proc(idata, adata, clus_result_df, save_root_dir, lr_pathway_db_file, celltype_key="ground_truth", x_row="x", y_col="y", save_df=True, save_df_dir=None, plot_volcano=True, plot_distribution=True, save_gene_exp=False):
    analyz_infc(idata_file=idata, clus_result_df=clus_result_df, lr_pathway_db_file=lr_pathway_db_file, save_root_dir=save_root_dir, plot_volcano=plot_volcano)

    # Create additional directories based on the root directory
    save_bc_df_dir = os.path.join(save_root_dir, "bc_df")
    save_pathway_dir_bc = os.path.join(save_root_dir, "pathway_bc")
    save_distribution_dir = os.path.join(save_root_dir, "distribution_cluster_id")
    save_distribution_type_dir = os.path.join(save_root_dir, "distribution_celltype")

    analyz_bc(idata, adata, clus_result_df, lr_pathway_db_file,
              celltype_key=celltype_key, x_row=x_row, y_col=y_col, save_df=save_df, save_df_dir=save_bc_df_dir, plot_volcano=plot_volcano, save_pathway_enrich=save_pathway_dir_bc, plot_distribution=plot_distribution, save_distribution_dir=save_distribution_dir, save_distribution_type_dir=save_distribution_type_dir, save_gene_exp=save_gene_exp)


lr_pathway_db_file_r = "../data/lr_pathway_meta.csv"

if len(sys.argv) < 6:
    print("Usage: python analyz_interface_n_boundary_cell_v2.py <idata name> <adata name> <cluster result file> <save root dir> <cell type key(optional)> <save gene exp map flag(optional)>")
    sys.exit(1)

celltype_key = "ground_truth"
if len(sys.argv) >= 6:
    celltype_key = sys.argv[5]

save_gene_map_f = False
if len(sys.argv) >= 7:
    save_gene_map_f = sys.argv[6].lower() == "true"

main_proc(idata=sys.argv[1], adata=sys.argv[2], clus_result_df=sys.argv[3], save_root_dir=sys.argv[4], lr_pathway_db_file=lr_pathway_db_file_r, celltype_key=celltype_key, save_gene_exp=save_gene_map_f)

print("Over..")