"""
python "/home/grads/qiusliang2/sctipts6/analyz_color_motif.py" "/public/qiusliang2/benchmark_st_data/pc1/type_def_GSM3036911_spatial_transcriptomics_20240607225416_re01.h5ad" "/public/qiusliang2/benchmark_st_data/pc1/idata.h5ad" Region "/home/grads/qiusliang2/sctipts6/colormotifs/fanmod_detect/" "/home/grads/qiusliang2/sctipts6/colormotifs/fanmod_detect/" False
"""
import sys
sys.path.append("..")
from src.reformat_adj import color_motif


def main_proc(adata_file, idata_file, celltype_key, save_adj_dir, output_dir, n_neighbor_v, delete_adj='True', pick_cells=True):
    analyz = color_motif(adata_file=adata_file, celltype_key=celltype_key, save_adj_dir=save_adj_dir, delete_adj=delete_adj, output_dir=output_dir,idata_file=idata_file, n_neighbor_v=n_neighbor_v)
    analyz.plot_motif_counts_distribution()  # same as all cells in motifs
    analyz.analyz_only_cells_in_motifs()
    analyz.analyz_all_cells_in_motifs()
    analyz.analyz_only_cells_in_motifs_cell_types()



if len(sys.argv) < 5:
    print("Usage: python analyz_color_motif.py <adata name> <idata name> <cell type key> <save dump files dir> <delete adj file or not>")
    sys.exit(1)

del_adj = 'False'
if len(sys.argv) >= 6:
    del_adj = sys.argv[5]

n_neighbor = None
if len(sys.argv) >= 7:
    n_neighbor = int(sys.argv[6])

main_proc(adata_file=sys.argv[1], idata_file=sys.argv[2], celltype_key=sys.argv[3], save_adj_dir=sys.argv[4], delete_adj=del_adj,n_neighbor_v=n_neighbor)
print("Over..")
