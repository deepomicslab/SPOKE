"""
python "/home/grads/qiusliang2/sctipts6/analyz_interface_n_boundary_cell_t.py" "/public/qiusliang2/benchmark_st_br_cancer/pc7/proc_interface/H1_removed_20240611155949/idata.h5ad" "/public/qiusliang2/benchmark_st_br_cancer/pc7/H1_removed_20240611155949_recons.h5ad" "/public/qiusliang2/benchmark_st_br_cancer/pc7/proc_interface/H1_removed_20240611155949/n_neigh=180_eta=4_opK=6_clubs=15.csv" "/home/grads/qiusliang2/result8/benchmark_brc/pc7/infc/lr/" "/home/grads/qiusliang2/result8/benchmark_brc/pc7/infc/pathway/" "/home/grads/qiusliang2/result8/benchmark_brc/pc7/bc/" "/home/grads/qiusliang2/result8/benchmark_brc/pc7/bc/pathway/" "/home/grads/qiusliang2/result8/benchmark_brc/pc7/bc/distribution/cluster_id/" "/home/grads/qiusliang2/result8/benchmark_brc/pc7/bc/distribution/celltype/" ground_truth True
"""
import pdb

import anndata
import pandas as pd
from scipy.sparse import csr_matrix
import os
import sys
sys.path.append("..")
from src.analyz_interface import Analysis_Interface
from src.analyz_boundary_cell import Analysis_Boundary_Cell


def analyz_infc(idata_file, clus_result_df, lr_pathway_db_file, save_lr_dir, save_pathway_dir, plot_volcano=False):
    idata = anndata.read_h5ad(idata_file)
    csr_sparse_matrix = csr_matrix(idata.X)
    idata.X = csr_sparse_matrix
    analyz = Analysis_Interface(idata_file=idata, result_file_dir=clus_result_df,
                                lr_pathway_db_file_dir=lr_pathway_db_file, plot_volcano=plot_volcano)
    # 1 for lr pairs
    lr_enrich_dir = save_lr_dir + "/"
    if not lr_enrich_dir is None and not os.path.exists(lr_enrich_dir):
        os.makedirs(lr_enrich_dir, exist_ok=True)
    analyz.proc_lr(save_lr_enrich=True, lr_enrich_dir=lr_enrich_dir)
    # 2 for pathways
    pathway_enrich_dir = save_pathway_dir + "/"
    if not pathway_enrich_dir is None and not os.path.exists(pathway_enrich_dir):
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


def main_proc(idata, adata, clus_result_df, save_lr_dir, save_pathway_dir_infc, lr_pathway_db_file, celltype_key="ground_truth", x_row="x", y_col="y", save_df=True,
              save_df_dir=None, plot_volcano=True, save_pathway_dir_bc=None,plot_distribution=True, save_distribution_dir=None, save_distribution_type_dir=None, save_gene_exp=False):
    analyz_infc(idata_file=idata, clus_result_df=clus_result_df, lr_pathway_db_file=lr_pathway_db_file, save_lr_dir=save_lr_dir, save_pathway_dir=save_pathway_dir_infc, plot_volcano=plot_volcano)
    analyz_bc(idata, adata, clus_result_df, lr_pathway_db_file,
              celltype_key=celltype_key, x_row=x_row, y_col=y_col, save_df=save_df, save_df_dir=save_df_dir, plot_volcano=plot_volcano, save_pathway_enrich=save_pathway_dir_bc, plot_distribution=plot_distribution, save_distribution_dir=save_distribution_dir, save_distribution_type_dir=save_distribution_type_dir, save_gene_exp=save_gene_exp)
    return


lr_pathway_db_file_r = "../data/lr_pathway_meta.csv"

if len(sys.argv) < 11:
    print("Usage: python analyz_interface_n_boundary_cell.py <idata name> <adata name> <cluster result file> <save lr dir> <save_pathway_dir_infc> <save_bc_df_dir> <save_pathway_dir_bc> <save bc distribution (cluster id combination)> <save bc distribution (cluster id combination with cell type)> <cell type key(optional)> <save gene exp map flag(optional)>")
    sys.exit(1)

celltype_key="ground_truth"
if len(sys.argv) >= 11:
    celltype_key = sys.argv[10]
save_gene_map_f = False
if len(sys.argv) >= 12:
    save_gene_map_f = sys.argv[11]
    save_gene_map_f = save_gene_map_f.lower() == "true"

main_proc(idata=sys.argv[1], adata=sys.argv[2], clus_result_df=sys.argv[3], save_lr_dir=sys.argv[4], save_pathway_dir_infc=sys.argv[5], save_df_dir=sys.argv[6], save_pathway_dir_bc=sys.argv[7], save_distribution_dir=sys.argv[8], save_distribution_type_dir=sys.argv[9],celltype_key=celltype_key, save_gene_exp=save_gene_map_f,lr_pathway_db_file=lr_pathway_db_file_r)

print("Over..")

