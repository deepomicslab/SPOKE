"""
2024/06/25
---
to get stague embedding:
conda activate stague5
cd STAGUE
python main.py --adata_file ./data/V1/raw/adata.h5ad --output_dir ./result/V1 --n_clusters 7 --adj_weight 0.3 --k 25 --clu_model kmeans
---
import anndata
import umap
import numpy as np
adata = anndata.read_h5ad("/public/qiusliang2/benchmark_st_data/pc1/stague_cell_gene/adata_processed.h5ad")
umap_embed = umap.UMAP(random_state=0).fit_transform(adata.X.todense().A)
np.savetxt("/public/qiusliang2/benchmark_st_data/pc1/stague_cell_gene/umap_aff_m_cell_gene.txt",umap_embed)
---
about co-occurrence:
import squidpy as sq
sq.datasets.imc()
sq.gr.co_occurrence(adata, cluster_key="optimal_cluster")
sq.pl.spatial_scatter(adata, color="optimal_cluster", size=20, shape=None)
sq.pl.co_occurrence(adata, cluster_key="optimal_cluster",clusters=0,save=)
"""
import numpy as np
import pandas as pd
import pdb
import os
import anndata
from pandas.core.frame import DataFrame
import sys
import time
import subprocess
import matplotlib.pyplot as plt
import re
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
import random
from scipy.sparse import csr_matrix
from collections import OrderedDict
import scanpy as sc
import gseapy as gp
import copy
import math
from itertools import combinations
import sys
sys.path.append("..")
from src.TEC import TEC


def remove_duplicates(text):
    words = text.split(',')
    words = list(OrderedDict.fromkeys(words))
    return ','.join(words)


class color_motif:
    """
    1. construct input adj matrix with cell type
    2. detect color motifs
    3. construct color_weight_graph
    4. cluster
    5. analyz
    (delete adj file or not)
    """
    def __init__(self, adata_file, idata_file, celltype_key, save_adj_dir, save_colored_df_pre=None, weight_mat_dir_path=None, delete_adj="false", fanmod_dir="../Network-Motif/fanmod/", output_dir="color_motifs", x_row="x", y_col="y", n_neighbor_v=None, eta_v=None, spatial_weight=None, save_color_graph_result_path=None, p_value_thres=0.05, shuffle_t=1000, save_enrich_df_path=None, lr_pathway_db_file="../data/lr_pathway_meta.csv", save_motif_count_dir=None, pvalue_key='pvalue', save_plots_dir=None, dot_size=8, pvals_enrich='pvals', save_enrich_dir=None, plot_volcano=True):
        if isinstance(adata_file, str):
            self.adata = anndata.read_h5ad(adata_file)
        else:
            self.adata = adata_file
        if isinstance(idata_file, str):
            self.idata = anndata.read_h5ad(idata_file)
        else:
            self.idata = idata_file
        self.cell_type = celltype_key
        self.save_adj_dir_root = save_adj_dir
        self.save_adj_dir = self.save_adj_dir_root + 'output_file.txt'
        self.fanmod_dir = fanmod_dir
        self.color_graph = None
        self.output_dir = output_dir + "/result"
        if not output_dir is None and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        self.x_row = x_row
        self.y_col = y_col
        self.del_adj = delete_adj.lower() == "true"
        self.save_colored_df = not self.del_adj
        self.save_color_graph = not self.del_adj
        self.save_colored_df_pre = self.save_adj_dir_root if save_colored_df_pre is None else save_colored_df_pre
        self.weight_mat_dir_path = self.save_adj_dir_root if weight_mat_dir_path is None else weight_mat_dir_path
        self.save_color_graph_result_path = self.save_adj_dir_root if save_color_graph_result_path is None else save_color_graph_result_path
        self.save_motif_count_dir = self.save_adj_dir_root if save_motif_count_dir is None else save_motif_count_dir
        self.save_plots_dir = self.save_adj_dir_root if save_plots_dir is None else save_plots_dir
        self.save_enrich_dir = self.save_adj_dir_root if save_enrich_dir is None else save_enrich_dir
        self.weight_mat_dir_path += "/color_graph.csv"
        num_celltype = len(self.adata.obs[self.cell_type].unique())
        t_round = (len(self.adata.obs)//(num_celltype+1))/10
        t_round = math.floor(t_round / 10) if math.floor(t_round / 10) >=1  else int(t_round)
        t_round = t_round * 10
        self.n_neighbor_v = t_round if n_neighbor_v is None else n_neighbor_v
        self.eta = 4 if eta_v is None else eta_v
        self.spatial_weight = 0.1 if spatial_weight is None else spatial_weight
        self.spatial_embedding = self.adata.obsm['X_stague_umap']
        self.save_color_graph_result = not self.del_adj
        self.cells_in_motifs = None
        self.mapping_cells = None
        self.p_value_thres = p_value_thres
        self.shuffle = shuffle_t
        self.save_enrich_df_path = self.save_adj_dir_root if save_enrich_df_path is None else save_enrich_df_path
        self.lr_pathway_db = lr_pathway_db_file
        self.pvalue_key = pvalue_key    # for analysis on motifs
        self.vertices_header = None
        self.pval_key_enrich = pvals_enrich   # for cluster on color map
        self.gene_sets = ['KEGG_2016', 'KEGG_2013']
        self.pval_key_cm = 'pvals'
        self.lr_db_cm = None

        self.dot_size = dot_size
        self.min_x = int(self.adata.obs['x'].min())
        self.max_x = int(self.adata.obs['x'].max() + 0.5)
        self.min_y = int(self.adata.obs['y'].min())
        self.max_y = int(self.adata.obs['y'].max() + 0.5)

        # plot paras for color motifs enrich
        self.plot_volcano = plot_volcano
        self.x_threshold_left = -1
        self.x_threshold_right = 1
        self.y_threshold = 1
        self.xmin = -6
        self.xmax = 6 + 2
        self.ymin = -1
        self.ymax = 3 + 2

        self.lr_db_infc = self.constrcut_infc_lr_meta_info()

        # main proc
        self.construct_adj()
        self.detect_cm()  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        self.out_df = self.construct_color_graph()
        self.cluster_model, self.color_graph_result = self.run_clustering_color_graph(pick_cells=False)
        self.cluster_model_picked, self.color_graph_result_picked = self.run_clustering_color_graph(
            pick_cells=True)
        # self.enrich_df, self.enrich_df_filtered = self.analyz_pathway_enrich(computation="man", all_path=True, shuffle_t=self.shuffle,
        #                       save_df=True, save_dir=self.save_enrich_df_path, p_thres=self.p_value_thres)
        self.proc_n_annotate(pick_cells=True)

        # after proc
        if self.del_adj:
            if os.path.exists(self.save_adj_dir):
                os.remove(self.save_adj_dir)

    def constrcut_infc_lr_meta_info(self):
        """
            pathways can be duplicated for each cell
        """
        mat = self.idata.X.todense().A
        lr_meta = self.idata.var
        lr_pathway_db = pd.read_csv(self.lr_pathway_db, index_col=False)
        interface_cell = self.idata.obs
        lr_pathway_db2 = lr_pathway_db.drop(["dest", "dest_i", "src_dest_i", "i"],
                                            axis=1)
        # lr_pathway_db3 = lr_pathway_db2.drop_duplicates()  # to get duplicates representing counts
        lr_pathway_db4 = lr_pathway_db2.groupby(["src"]).agg(
            {'type': lambda x: ','.join(x.unique()),
             'pathway': lambda x: ','.join(x.unique()),
             'source': lambda x: ','.join(x.unique())}).reset_index()

        def remove_duplicates(text):
            words = text.split(',')
            words = list(OrderedDict.fromkeys(words))
            return ','.join(words)

        lr_pathway_db4['type'] = lr_pathway_db4['type'].apply(remove_duplicates)
        # lr_pathway_db4['pathway'] = lr_pathway_db4['pathway'].apply(remove_duplicates)
        lr_pathway_db4['source'] = lr_pathway_db4['source'].apply(remove_duplicates)

        # change to lower case simutaneously
        lr_meta['receptor'] = lr_meta['receptor'].str.lower()
        lr_pathway_db4['src'] = lr_pathway_db4['src'].str.lower()
        lr_meta_pathway = pd.merge(lr_meta, lr_pathway_db4, how='left',
                                   left_on='receptor',
                                   right_on='src')
        lr_meta_pathway = lr_meta_pathway.drop("i", axis=1)
        # lr_meta_pathway.to_csv(infc_LR_pathway_file_name, index=False)

        lr_pos = mat != 0
        # nonzero_indices = np.nonzero(mat)
        # pathways = lr_meta_pathway.loc[nonzero_indices][
        #     ['ligand', 'receptor', 'ligand_i', 'receptor_i', 'pathway', 'source', 'type']]
        lr_pathways_df = None
        for row in lr_pos:
            lr_pairs = lr_meta_pathway[row]
            lr_pathways = lr_pairs[["type", "pathway", "source"]].apply(
                lambda x: ','.join(x.dropna()), axis=0)
            lr_pathways['type'] = remove_duplicates(lr_pathways['type'])
            # lr_pathways['pathway'] = remove_duplicates(lr_pathways['pathway'])
            lr_pathways['source'] = remove_duplicates(lr_pathways['source'])
            lr_pathways_df = pd.concat([lr_pathways_df, lr_pathways], axis=1)
        lr_pathways_df_t_out = lr_pathways_df.transpose()
        lr_pathways_df_t_out.index = interface_cell.index
        interface_info_merged = pd.concat(
            [interface_cell, lr_pathways_df_t_out], axis=1)
        interface_info_merged = interface_info_merged.rename(
            columns={"i": "interface_i"})
        return interface_info_merged

    def constrcut_cm_lr_meta_info(self, adata):
        """
            pathways can be duplicated for each cell
        """
        mat = adata.X
        lr_meta = adata.var
        lr_pathway_db = pd.read_csv(self.lr_pathway_db, index_col=False)
        interface_cell = adata.obs
        lr_pathway_db2 = lr_pathway_db.drop(["dest", "dest_i", "src_dest_i", "i"],
                                            axis=1)
        # lr_pathway_db3 = lr_pathway_db2.drop_duplicates()  # to get duplicates representing counts
        lr_pathway_db4 = lr_pathway_db2.groupby(["src"]).agg(
            {'type': lambda x: ','.join(x.unique()),
             'pathway': lambda x: ','.join(x.unique()),
             'source': lambda x: ','.join(x.unique())}).reset_index()

        def remove_duplicates(text):
            words = text.split(',')
            words = list(OrderedDict.fromkeys(words))
            return ','.join(words)

        lr_pathway_db4['type'] = lr_pathway_db4['type'].apply(remove_duplicates)
        # lr_pathway_db4['pathway'] = lr_pathway_db4['pathway'].apply(remove_duplicates)
        lr_pathway_db4['source'] = lr_pathway_db4['source'].apply(remove_duplicates)

        # change to lower case simutaneously
        lr_meta['receptor'] = lr_meta['receptor'].str.lower()
        lr_pathway_db4['src'] = lr_pathway_db4['src'].str.lower()
        lr_meta_pathway = pd.merge(lr_meta, lr_pathway_db4, how='left',
                                   left_on='receptor',
                                   right_on='src')
        lr_meta_pathway = lr_meta_pathway.drop("i", axis=1)
        # lr_meta_pathway.to_csv(infc_LR_pathway_file_name, index=False)

        lr_pos = mat != 0
        # nonzero_indices = np.nonzero(mat)
        # pathways = lr_meta_pathway.loc[nonzero_indices][
        #     ['ligand', 'receptor', 'ligand_i', 'receptor_i', 'pathway', 'source', 'type']]
        lr_pathways_df = None
        for row in lr_pos:
            lr_pairs = lr_meta_pathway[row]
            lr_pathways = lr_pairs[["type", "pathway", "source"]].apply(
                lambda x: ','.join(x.dropna()), axis=0)
            lr_pathways['type'] = remove_duplicates(lr_pathways['type'])
            # lr_pathways['pathway'] = remove_duplicates(lr_pathways['pathway'])
            lr_pathways['source'] = remove_duplicates(lr_pathways['source'])
            lr_pathways_df = pd.concat([lr_pathways_df, lr_pathways], axis=1)
        lr_pathways_df_t_out = lr_pathways_df.transpose()
        lr_pathways_df_t_out.index = interface_cell.index
        cm_info_merged = pd.concat(
            [interface_cell, lr_pathways_df_t_out], axis=1)
        cm_info_merged = cm_info_merged.rename(
            columns={"i": "cm_i"})
        return cm_info_merged

    def construct_cm_lr_db(self):
        # generate pathway info
        lr_pathway_db = pd.read_csv(self.lr_pathway_db,index_col=False)
        lr_pathway_db2 = lr_pathway_db.drop(["dest", "dest_i", "src_dest_i", "i"],
                                            axis=1)
        lr_pathway_db3 = lr_pathway_db2.drop_duplicates()
        lr_pathway_db4 = lr_pathway_db3.groupby(["src"]).agg(
            {'type': lambda x: ','.join(x.unique()),
             'pathway': lambda x: ','.join(x.unique()),
             'source': lambda x: ','.join(x.unique())}).reset_index()
        lr_pathway_db4['type'] = lr_pathway_db4['type'].apply(remove_duplicates)
        lr_pathway_db4['pathway'] = lr_pathway_db4['pathway'].apply(remove_duplicates)
        lr_pathway_db4['source'] = lr_pathway_db4['source'].apply(remove_duplicates)
        var_lower = self.idata.var.copy()
        var_lower['receptor'] = var_lower['receptor'].str.lower()
        lr_pathway_db4_lower = lr_pathway_db4.copy()
        lr_pathway_db4_lower['src'] = lr_pathway_db4_lower['src'].str.lower()
        lr_meta_pathway = pd.merge(var_lower, lr_pathway_db4_lower, how='left',
                                   left_on='receptor', right_on='src')
        lr_meta_pathway['pathway'].fillna('nan').tolist()
        return lr_meta_pathway

    def construct_adj(self):
        mat = self.adata.obsp['spatial_connectivities']
        if isinstance(mat, csr_matrix):
            mat = mat.toarray()
        lines = np.where(mat == 1)
        df = pd.DataFrame(lines).T
        dff = copy.deepcopy(self.adata.obs)
        dff.reset_index(drop=True, inplace=True)
        df[2] = df[0].map(dff[self.cell_type])
        df[3] = df[1].map(dff[self.cell_type])
        mapping = {val: idx for idx, val in enumerate(df[2].unique())}
        self.mapping_cells = mapping
        df[2] = df[2].map(mapping)
        df[3] = df[3].map(mapping)
        df.to_csv(self.save_adj_dir, sep='\t', index=False)

    def detect_cm(self):
        """
        color motifs (cm)
        """
        os.chdir(self.fanmod_dir)
        command = "./fanmod_command_line_linux"
        args = ["3", "100000", "1", self.save_adj_dir, "0", "1", "0", "2", "1", "0", "1",
                "1000", "3", "3", self.output_dir+".csv", "0", "1"]

        result = subprocess.run([command] + args, stdout=subprocess.PIPE)

        print(result.stdout.decode('utf-8'))

    def construct_color_graph(self):
        """
        process dump file and csv file to construct the color graph
        """
        print("Constructing color graph...")
        result_df, column_names = self.extract_table_from_text_file(self.output_dir+".csv")
        result_df = self.proc_df(df=result_df, header=column_names)
        result_df, vertices_header_name = self.map_n_enrich(df=result_df, dump=self.output_dir+".csv.dump")
        self.vertices_header = vertices_header_name
        result_df, out_df = self.generate_color_df(df=result_df, save_df=self.save_colored_df,
                                                 save_dir=self.save_colored_df_pre,
                                                 save_vertices_merged=True)
        # plot_motif_distribution_colors(color_df=result_df, save_figs=True, fig_dir=save_figs_dir) # show the result of motif detection
        weighted_mat = self.convert_2_cell_weight_graph_v2(df=out_df,
                                                   save_mat=self.save_color_graph,
                                                   weight_mat_dir=self.weight_mat_dir_path)
        self.color_graph = weighted_mat
        self.cells_in_motifs = self.cell_rows_related(df=out_df)
        return out_df

    def analyz_cm_enrich(self, df, save_dir, mode="avg"):
        """
        2024/08/18: analyze the enrichment of color motifs using the infc-lr exp,
        construct adata for color motifs * lr data,
        and analyze
        """
        X_infc = self.idata.X.toarray()
        X_infc= np.where(X_infc < 0, 0, X_infc)
        col_with_color = df.columns[df.columns.str.contains('color')].tolist()
        vertices_col = [s.replace('_color', '') for s in col_with_color]
        if 'color_combination' in vertices_col:
            vertices_col.remove('color_combination')
        infc_meta = self.idata.obs
        X = np.zeros(shape=(len(df),X_infc.shape[1]))
        infc_meta['sender'] = infc_meta['sender'].astype(str)
        infc_meta['receiver'] = infc_meta['receiver'].astype(str)
        df = df.reset_index(drop=True)
        for index, cells in df[vertices_col].iterrows():
            c_n = len(cells)  # number of cells
            cells = cells.astype(str)
            pth_info = pd.DataFrame()
            for i in range(c_n - 1):
                c1 = cells[i]
                for j in range(i + 1, c_n):
                    c2 = cells[j]
                    # note: double directions: c1->c2, c2->c1
                    pth = infc_meta[(infc_meta['sender'] == c1) & (
                            infc_meta['receiver'] == c2) |
                                       (infc_meta['sender'] == c2) & (
                                               infc_meta['receiver'] == c1)]
                    if not pth.empty:
                        pth_info = pd.concat([pth_info, pth], axis=0)
            if not pth_info.empty:
                # get all the row indices of X
                rows = pth_info['i']
                X_t = X_infc[rows,:]
                if mode == "avg":
                    X[index, :] += np.mean(X_t, axis=0)
                else:
                    X[index,:] += np.sum(X_t, axis=0)
        adata = anndata.AnnData(X=X, var=self.idata.var, obs=df)
        # analyz lr enrich
        self.analyz_cm_lr(adata=adata, save_dir=save_dir)
        # analyz pathway enrich
        self.analyz_cm_pathway(adata_lr=adata, save_dir=save_dir)
        return

    def analyz_cm_lr(self, adata, save_dir, sta_method='t-test'):
        # process
        adata.obs['color_combination'] = adata.obs['color_combination'].astype('category')
        sc.pp.normalize_total(adata, target_sum=1e4)  # target_sum can be varied
        sc.pp.log1p(adata)
        sc.tl.rank_genes_groups(adata, 'color_combination', method=sta_method,
                                key_added=sta_method)
        self.lr_db_cm = self.construct_cm_lr_db()
        lr_meta_pathway2 = self.lr_db_cm.set_index(adata.var.index)  # set index to synchronize
        adata.var['pathway'] = lr_meta_pathway2['pathway']
        key_ids = [item[0] for item in
                   adata.uns[sta_method]['pvals'].dtype.descr]
        save_dir += "/cm_lrs/"
        if not save_dir is None and not os.path.exists(save_dir):
            os.makedirs(save_dir, exist_ok=True)
        for key_id in key_ids:
            lr_t = lr_meta_pathway2.loc[
                list(adata.uns[sta_method]['names'][key_id])]
            lr_t['score'] = adata.uns[sta_method]['scores'][key_id]
            lr_t['pvals'] = adata.uns[sta_method]['pvals'][key_id]
            lr_t['pvals_adj'] = adata.uns[sta_method]['pvals_adj'][key_id]
            lr_t['logfoldchanges'] = \
            adata.uns[sta_method]['logfoldchanges'][key_id]
            if save_dir:
                lr_t.to_csv(save_dir + '/' + "cms_enrich_" + str(key_id) + ".csv",
                            index=False)

            if self.plot_volcano:
                # plot volcano plot
                result = pd.DataFrame()
                result['x'] = lr_t['logfoldchanges']
                result['y'] = lr_t[self.pval_key_cm]
                result['y'] = -np.log10(result['y'])
                # result['x'] = result['x'].apply(lambda x: np.exp2(x))

                result['group'] = 'black'
                result.loc[(result.x > self.x_threshold_right) & (
                        result.y > self.y_threshold), 'group'] = 'tab:red'  # x=-+x_threshold直接截断
                result.loc[(result.x < self.x_threshold_left) & (
                        result.y > self.y_threshold), 'group'] = 'tab:blue'  # x=-+x_threshold直接截断
                result.loc[result.y < self.y_threshold, 'group'] = 'dimgrey'

                fig = plt.figure(figsize=plt.figaspect(7 / 6))
                ax = fig.add_subplot()
                ax.set(xlim=(self.xmin, self.xmax), ylim=(self.ymin, self.ymax),
                       title='')
                ax.scatter(result['x'], result['y'], s=2, c=result['group'])
                ax.set_ylabel('-Log10(p value)', fontweight='bold')
                ax.set_xlabel('Log2 (fold change)', fontweight='bold')
                ax.spines['right'].set_visible(False)
                ax.spines['top'].set_visible(False)

                # 水平和竖直线
                ax.vlines(self.x_threshold_left, self.ymin, self.ymax, color='dimgrey',
                          linestyle='dashed',
                          linewidth=1)
                ax.vlines(self.x_threshold_right, self.ymin, self.ymax, color='dimgrey',
                          linestyle='dashed',
                          linewidth=1)
                ax.hlines(self.y_threshold, self.xmin, self.xmax, color='dimgrey',
                          linestyle='dashed',
                          linewidth=1)

                ax.set_xticks(range(self.xmin, self.xmax, 2))
                ax.set_yticks(range(self.ymin, self.ymax, 2))
                plt.title(str(key_id))
                plt.savefig(save_dir + '/' + "cms_enrich_lr_" + str(key_id) + ".png",
                            format="png")

    def analyz_cm_pathway(self, adata_lr, save_dir, save_pathway_enrich=True, sta_method='t-test'):
        pathway_meta = self.lr_db_cm
        lr_meta_pathway2 = pathway_meta.set_index(
            adata_lr.var.index)  # set index to synchronize
        adata = copy.deepcopy(adata_lr)
        adata.var['pathway'] = lr_meta_pathway2['pathway']
        X_n = np.zeros((len(adata.obs), len(adata.var['pathway'].unique())))
        pathways = adata.var['pathway'].unique()
        for i, pathway in enumerate(pathways):
            rows = adata.var['pathway'] == pathway
            X_n[:, i] = adata.X[:, adata.var[rows]['i']].sum(axis=1)
        adata1 = anndata.AnnData(X=X_n, var=pd.DataFrame({'pathway': pathways}),
                                 obs=adata.obs)
        adata = adata1
        save_dir += "/cm_pathways/"

        # process
        adata.obs['color_combination'] = adata.obs['color_combination'].astype('category')
        sc.pp.normalize_total(adata, target_sum=1e4)  # target_sum can be varied
        sc.pp.log1p(adata)
        sc.tl.rank_genes_groups(adata, 'color_combination', method=sta_method,
                                key_added=sta_method)
        key_ids = [item[0] for item in adata.uns[sta_method]['pvals'].dtype.descr]
        if not save_dir is None and not os.path.exists(save_dir):
            os.makedirs(save_dir, exist_ok=True)
        for key_id in key_ids:
            pathway_t = adata.var.loc[list(adata.uns[sta_method]['names'][key_id])]
            pathway_t['score'] = adata.uns[sta_method]['scores'][key_id]
            pathway_t['pvals'] = adata.uns[sta_method]['pvals'][key_id]
            pathway_t['pvals_adj'] = adata.uns[sta_method]['pvals_adj'][key_id]
            pathway_t['logfoldchanges'] = adata.uns[sta_method]['logfoldchanges'][
                key_id]
            pathway_t['color_combination'] = key_id
            if save_pathway_enrich:
                pathway_t.to_csv(save_dir + "pathway_" + str(key_id) + ".csv",
                                 index=True)

    def remap_cell_rows(self, z_rows, result, o_rows):
        """
        cells not in motifs are labeled as -1
        """
        re = pd.DataFrame(-1, index=o_rows, columns=result.columns, dtype=int)
        result['index'] = z_rows
        re.loc[z_rows] = result.set_index('index')
        re = re.astype(int)
        return re

    def cell_rows_related(self, df):
        """
        return indices of cells related to motifs
        (remove indices of cells not in color motifs)
        type = arr
        """
        indices = []
        for vertex_index in self.vertices_header:
            indices.extend(list(df[vertex_index]))
        indices = set(indices)
        indices = list(indices)
        indices.sort()
        return indices

    def run_clustering_color_graph(self, pick_cells=False):
        X = self.color_graph
        spatial_embedding = self.spatial_embedding
        n_neighbor_v = self.n_neighbor_v
        if pick_cells:
            org_rows = list(range(X.shape[0]))
            rows = self.cells_in_motifs
            rows_indices = np.ix_(rows, rows)
            X = X[rows_indices]
            spatial_embedding = self.spatial_embedding[rows,:]
            num_celltype = len(self.adata.obs[self.cell_type].unique())
            n_neighbor_t = int(len(rows)/(num_celltype+1))
            n_neighbor_t = round(n_neighbor_t,-1) if round(n_neighbor_t,-1)>=10 else n_neighbor_t
            n_neighbor_v = min(n_neighbor_v, n_neighbor_t)
        bottom_up = TEC(affinity="gaussian_kernel", kernel_gamma=None,
                         sparsification="knn_neighbors",
                         n_neighbors=n_neighbor_v, objective="KL", strategy="bottom_up",
                         eta=self.eta, eta1=1, eta_mode="coefficient", eta2=1,
                         merge_layers=True, plot_cluster_map=True,
                         __verbose__=True, max_k=30)
        bottom_up.fit_v3(X, spatial_embedding, self.spatial_weight)
        clustering_result = {}
        clustering_result["optimal_cluster"] = bottom_up.labels_
        clustering_result['clubs'] = bottom_up.clubs
        df_tmp = DataFrame(clustering_result)
        df_tmp = df_tmp.join(bottom_up.ks_clusters)
        file_name = "/color_graph_result_n_neighbor_"+str(n_neighbor_v)+"_opK=" + str(
                            bottom_up.optimal_k) + "_clubK=" + str(
                            bottom_up.club_k) + ".csv"
        if pick_cells:
            df_tmp = self.remap_cell_rows(z_rows=rows, result=df_tmp, o_rows=org_rows)
            file_name = "/color_graph_result_n_neighbor_"+str(n_neighbor_v)+"_opK=" + str(
                bottom_up.optimal_k) + "_clubK=" + str(
                bottom_up.club_k) + "_cells_in_motifs_picked.csv"
        if self.save_color_graph_result:
            df_tmp.to_csv(self.save_color_graph_result_path + file_name)
        return bottom_up, df_tmp

    def extract_table_from_text_file(self, file_path):
        """
        extract dataframe from the csv file
        """
        with open(file_path, 'r') as file:
            text = file.read()

        result_table = re.search(r'Result overview:(.*?)$', text, re.DOTALL)

        current_row = []
        if result_table:
            table_text = result_table.group(1)
            rows = table_text.strip().split('\n')
            headers1 = rows[0].strip().split(',')
            headers2 = rows[1].strip().split(',')
            header = pd.concat([pd.DataFrame(headers1), pd.DataFrame(headers2)], axis=1)
            header = header.apply(
                lambda x: ' '.join(str(val) for val in x if pd.notnull(val)), axis=1)
            df = header
            for row in rows[3:]:
                if row:
                    current_row.extend(row.split())
                else:
                    line = current_row[0].strip().split(',')
                    adj_new = line[1] + current_row[1] + current_row[2]
                    adj_new = adj_new.replace(",", "")
                    line[1] = adj_new
                    df = pd.concat([df, pd.DataFrame(line)], axis=1)
                    current_row = []
            line = current_row[0].strip().split(',')
            adj_new = line[1] + ''.join(current_row[1:])
            adj_new = adj_new.replace(",", "")
            line[1] = adj_new
            df = pd.concat([df, pd.DataFrame(line)], axis=1)
            df = df.T
            df.columns = df.iloc[0]
            df = df[1:]
            df = df.reset_index()
            df = df.drop("index", axis=1)
            # change types
            df['p-Value'] = df['p-Value'].astype(float)
            return df, header
        return

    def proc_df(self, df, header, pvalue_thr=0.05):
        # sort according to p-value
        df_sorted = df.sort_values(header.iloc[-1])
        df_sorted[header.iloc[-1]] = df_sorted[header.iloc[-1]].astype(float)
        filtered_df = df_sorted[df_sorted[header.iloc[-1]] < pvalue_thr]
        return filtered_df

    def extract_dump(self, dump):
        with open(dump, 'r') as file:
            text = file.read()
        result_table = re.search(
            r'Format: adjacency matrix, <participating vertices>(.*?)$', text,
            re.DOTALL)
        df = pd.DataFrame()
        if result_table:
            table_text = result_table.group(1)
            rows = table_text.strip().split('\n')
            # for row in rows:
            #     row = row.strip().split(',')
            #     df = pd.concat([df, pd.DataFrame([row])])
            # df.reset_index(drop=True, inplace=True)
            # n_col = len(df.columns) - 1
            # col_list = ['Adj_Matrix']
            # col_list.extend(["vertex_" + str(i) for i in range(1, n_col + 1)])
            # df.columns = col_list
            data = [row.strip().split(',') for row in rows]

            n_col = len(data[0]) - 1
            col_list = ['Adj_Matrix'] + [f'vertex_{i}' for i in range(1, n_col + 1)]

            df = pd.DataFrame(data, columns=col_list)
            return df
        return df

    def get_indices_up_triangle(self, vertex_num):
        indices = []
        for i in range(vertex_num):
            for j in range(i + 1, vertex_num):
                indices.append(i * vertex_num + j)
        return indices

    def remove_not_complete(self, adj_matrix, vertex_num):
        """
        only remain motif with all the vertices connected with each other
        ---
        1. transform the list to a matrix
        2. obtain the upper triangle
        3. compute the product, if 0, return True
        ---
        we can also only keep one single type (ID) that correspond to the complete graph
        """
        indices = self.get_indices_up_triangle(vertex_num)
        code_list = list(map(int, adj_matrix))
        arr = np.array(code_list)[indices]
        if np.prod(arr):
            return True
        else:
            return False

    def map_vertices(self, left, right):
        """
        map adj_matrix back to the # of vertices
        """
        merged = left.merge(right, left_on=left.columns[1], right_on='Adj_Matrix',
                            how='left')
        return merged

    def remove_duplicate_rows(self, df, columns, vertex_num):
        """
        remove motif with same colors
        """
        rows_to_remove = []
        for index, row in df.iterrows():
            values_to_compare = row[columns]
            if len(set(values_to_compare)) <= vertex_num - 1:
                rows_to_remove.append(index)
        df = df.drop(rows_to_remove)
        complete_graph_idx = df.apply(
            lambda row, N: self.remove_not_complete(row['Adj-Matrix '], N),
            args=(vertex_num,), axis=1)
        # df.apply(remove_not_complete, args=(vertex_num,), axis=1)
        df = df[complete_graph_idx]
        return df

    def filter_vertices(self, df):
        """
        1. map colors
        2. map xy
        delete motifs with duplicate color(s)
        """
        dff = self.adata.obs
        if 'barcode' in dff.columns:
            dff = dff.drop('barcode', axis=1)
        dff.reset_index(inplace=True)
        colors = dff[[self.cell_type]]
        xy = dff[[self.x_row, self.y_col]]
        c_type = dff[[self.cell_type]]
        vertices_header = df.columns[8:]
        key_w = colors.columns[0]
        vertices_color = []
        vertices_x = []
        vertices_y = []
        for col in vertices_header:
            # map color
            col_name = col + '_color'
            df[col] = df[col].astype(int)
            df[col_name] = df[col].map(colors[key_w])
            vertices_color.append(col_name)
            # map spot
            spot_name_x = col + '_x'
            spot_name_y = col + '_y'
            vertices_x.append(spot_name_x)
            vertices_y.append(spot_name_y)
            df[spot_name_x] = df[col].map(xy[self.x_row])
            df[spot_name_y] = df[col].map(xy[self.y_col])
            # map cell type
            type_col = col + "_type"
            df[type_col] = df[col].map(c_type[c_type.columns[0]])
        # filter motifs with vertices all in the same color

        vertex_num = len(vertices_header)
        df = self.remove_duplicate_rows(df, vertices_color, vertex_num)
        # average spot of vertices
        df['vertices_x_avg'] = df[vertices_x].mean(axis=1)
        df['vertices_y_avg'] = df[vertices_y].mean(axis=1)
        return df, vertices_header

    def map_n_enrich(self, df, dump):
        """
        1. map the motifs to the cells(nodes) in the graph
        2. filter the motifs so that there are at least 2 colors
        3. map cell location info (default: x,y)
        4. analyze the enrichment
        """
        meta_info = self.extract_dump(dump)  # adj_matrix,vertex,vertex,vertex
        dfs = self.map_vertices(left=df, right=meta_info)
        dfs, vertices_color = self.filter_vertices(df=dfs)
        return dfs, vertices_color

    def map_lr_n_pathways(self, df):
        """
        map l-r pairs and pathways enriched in the morif according to the cell-cell combination
        1. read meta file
        2. for each motif, generate cell combination, find l-r pair
        """
        pathway_meta = self.lr_db_infc
        col_with_color = df.columns[df.columns.str.contains('color')].tolist()
        # col_with_color.remove('color_combination')
        vertices_col = [s.replace('_color', '') for s in col_with_color]
        pathway_info_cols = ['type', 'pathway', 'source']
        pth_info_merged = pd.DataFrame()
        pathway_meta['sender'] = pathway_meta['sender'].astype(df[vertices_col[0]].dtype)
        pathway_meta['receiver'] = pathway_meta['receiver'].astype(
            df[vertices_col[0]].dtype)

        pth_info_list = []

        for index, cells in df[vertices_col].iterrows():
            c_n = len(cells)  # number of cells
            pth_info = pd.DataFrame()
            for i in range(c_n - 1):
                c1 = cells[i]
                for j in range(i + 1, c_n):
                    c2 = cells[j]
                    # note: double directions: c1->c2, c2->c1
                    pth = pathway_meta[(pathway_meta['sender'] == c1) & (
                                pathway_meta['receiver'] == c2) |
                                       (pathway_meta['sender'] == c2) & (
                                                   pathway_meta['receiver'] == c1)]
                    if not pth.empty:
                        pth_info_t = pth[pathway_info_cols]
                        pth_info = pd.concat([pth_info, pth_info_t], axis=0)
            # pth_info.drop_duplicates()
            # pth_info.dropna(inplace=True)
            if pth_info.empty or all(pth_info['pathway']==""):
                pth_info = pd.DataFrame(columns=pathway_info_cols)
                # pth_info.loc[0] = {'type': '', 'pathway': '', 'source': ''}
                pth_info.loc[0] = [''] * len(pathway_info_cols)
            else:
                pth_info.dropna(how='all', inplace=True)
                pth_info = pd.DataFrame(pth_info.apply(lambda x: ','.join(x), axis=0)).T

                pth_info.loc[0]['type'] = ','.join(set(pth_info.loc[0]['type'].split(',')))
                pth_info.loc[0]['source'] = ','.join(
                    set(pth_info.loc[0]['source'].split(',')))

            # pth_info_merged = pd.concat([pth_info_merged, pth_info], axis=0)
            pth_info_list.append(pth_info)
        dft = pd.concat(pth_info_list, axis=0).reset_index()
        dft = dft.drop('index', axis=1)
        df = pd.concat([df.reset_index(), dft], axis=1)
        return df

    def map_lr_n_pathways_v2(self, df):
        """
        map l-r pairs and pathways enriched in the morif according to the cell-cell combination
        1. read meta file
        2. for each motif, generate cell combination, find l-r pair
        """
        pathway_meta = self.lr_db_infc
        col_with_color = df.columns[df.columns.str.contains('color')].tolist()
        # col_with_color.remove('color_combination')
        vertices_col = [s.replace('_color', '') for s in col_with_color]
        pathway_info_cols = ['type', 'pathway', 'source']
        pathway_meta['sender'] = pathway_meta['sender'].astype(df[vertices_col[0]].dtype)
        pathway_meta['receiver'] = pathway_meta['receiver'].astype(
            df[vertices_col[0]].dtype)

        pth_info_list = []

        for index, cells in df[vertices_col].iterrows():
            cell_pairs = pd.DataFrame(
                [(c1, c2) for i, c1 in enumerate(cells) for c2 in cells[i + 1:]],
                columns=['cell1', 'cell2'])

            # Merge to find pathways in both directions
            merged1 = cell_pairs.merge(pathway_meta, left_on=['cell1', 'cell2'],
                                       right_on=['sender', 'receiver'], how='left')
            merged2 = cell_pairs.merge(pathway_meta, left_on=['cell1', 'cell2'],
                                       right_on=['receiver', 'sender'], how='left')

            # Combine results from both directions
            merged = pd.concat([merged1, merged2], ignore_index=True)
            pth_info = merged.dropna(subset=pathway_info_cols)

            # Select only the necessary columns
            pth_info = pth_info[pathway_info_cols]
            # pth_info.drop_duplicates()
            # pth_info.dropna(inplace=True)
            if pth_info.empty or all(pth_info['pathway']==""):
                pth_info = pd.DataFrame(columns=pathway_info_cols)
                # pth_info.loc[0] = {'type': '', 'pathway': '', 'source': ''}
                pth_info.loc[0] = [''] * len(pathway_info_cols)
            else:
                pth_info.dropna(how='all', inplace=True)
                pth_info = pd.DataFrame(pth_info.apply(lambda x: ','.join(x), axis=0)).T

                pth_info.loc[0]['type'] = ','.join(set(pth_info.loc[0]['type'].split(',')))
                pth_info.loc[0]['source'] = ','.join(
                    set(pth_info.loc[0]['source'].split(',')))

            # pth_info_merged = pd.concat([pth_info_merged, pth_info], axis=0)
            pth_info_list.append(pth_info)
        dft = pd.concat(pth_info_list, axis=0).reset_index()
        dft = dft.drop('index', axis=1)
        df = pd.concat([df.reset_index(), dft], axis=1)
        return df

    def generate_color_df(self, df, save_df=True, save_dir=None, save_vertices_merged=True):
        """
        add "color_combination" column;
        save indices of vertex not in any motif;
        map pathway function.
        """
        print("Generating color dataframe...")
        col_with_color = df.columns[df.columns.str.contains('color')].tolist()
        # df = self.map_lr_n_pathways_v2(df=df)
        # pathway_info_cols = ['type', 'pathway', 'source']
        # color_df = df[
        #     col_with_color + ['vertices_x_avg',
        #                       'vertices_y_avg'] + pathway_info_cols].copy()
        color_df = df[
            col_with_color + ['vertices_x_avg',
                              'vertices_y_avg'] ].copy()
        color_df['color_combination'] = color_df[col_with_color].apply(
            lambda x: tuple(sorted(x.tolist())), axis=1)
        df['color_combination'] = color_df['color_combination']
        self.analyz_cm_enrich(df=df, mode="sum", save_dir=save_dir)
        if save_df:
            df.to_csv(save_dir + "color_df.csv", index=False)

        if save_vertices_merged:
            vertices = []
            for vertex in col_with_color:
                vertex = vertex.replace("_color", "")
                vertices.extend(list(df[vertex].unique()))
            vertices = list(set(vertices))
            np.savetxt(save_dir + "vertices_in_motifs.txt", vertices)
        return color_df, df

    def generate_color_df_v2(self, df, save_df=True, save_dir=None,
                             save_vertices_merged=True):
        """
        add "color_combination" column;
        save indices of vertex not in any motif;
        not map pathway function.
        """
        col_with_color = df.columns[df.columns.str.contains('color')].tolist()
        color_df = df[
            col_with_color + ['vertices_x_avg',
                              'vertices_y_avg']].copy()
        color_df['color_combination'] = color_df[col_with_color].apply(
            lambda x: tuple(sorted(x.tolist())), axis=1)
        df['color_combination'] = color_df['color_combination']
        if save_df:
            df.to_csv(save_dir + "color_df.csv", index=False)

        if save_vertices_merged:
            vertices = []
            for vertex in col_with_color:
                vertex = vertex.replace("_color", "")
                vertices.extend(list(df[vertex].unique()))
            vertices = list(set(vertices))
            np.savetxt(save_dir + "vertices_in_motifs.txt", vertices)
        return color_df, df

    # def count_cells_with_pathway(df, pathway):
    #     count = 0
    #     for value in df['pathway']:
    #         if pathway in value:
    #             count += 1
    #     return count

    def count_cells_with_pathway_once(self,df, pathway):
        """
        count how many times this pathway merge in df, one motif counts only once
        """
        count = 0
        for index, row in df.iterrows():
            pths = row['pathway']
            pth_list = pths.split(',')
            if pathway in pth_list:
                count += 1
        return count

    def count_cells_with_pathway_all(self, df, pathway):
        """
        count how many times this pathway merge in df, one motif counts all
        """
        count = 0
        for index, row in df.iterrows():
            pths = row['pathway']
            pth_list = pths.split(',')
            count += pth_list.count(pathway)
        return count

    def count_cells_with_pathway_all_v2(self, df, pathway):
        """
        count how many times this pathway merge in df, one motif counts all
        """
        count = df['pathway'].str.split(',').explode().value_counts().get(pathway, 0)
        return count

    def count_cells_with_pathway_pack(self, all_path=True):
        if all_path:
            count_cells_with_pathway = self.count_cells_with_pathway_all_v2
        else:
            count_cells_with_pathway = self.count_cells_with_pathway_once
        return count_cells_with_pathway

    def compute_pathway_n_pvalue_surv_func(self, df, total_df, total_cells, pathways_all,
                                           all_path=True):
        """
        compute pathway enrichment for each color motif
        with scipy.stats.hypergeom.sf version
        for no duplicate pathway in each row
        """
        enrichment_results = pd.DataFrame(columns=['pathway', 'pvalue'])
        pathways = self.get_pathways(df)
        count_cells_with_pathway = self.count_cells_with_pathway_pack(all_path)
        pvalue_z = 1
        for pathway in pathways:
            if pd.isnull(pathway) or pathway == "":
                pvalue = pvalue_z
            else:
                # 计算包含该路径的细胞数（k）
                pathway_cells = count_cells_with_pathway(df, pathway)
                pathway_all_cells = count_cells_with_pathway(total_df, pathway)

                # survival function: smaller p, more possible enrichment
                pvalue = hypergeom.sf(pathway_cells - 1, total_cells,
                                      pathway_all_cells, len(df))

            enrichment_results = pd.concat([enrichment_results, pd.DataFrame(
                {'pathway': [pathway], 'pvalue': [pvalue]})], ignore_index=True)

        # add other pvalues
        pathways_res = list(set(pathways_all) - set(pathways))
        enrichment_results = pd.concat([enrichment_results, pd.DataFrame(
            {'pathway': pathways_res, 'pvalue': [pvalue_z] * len(pathways_res)})],
                                       ignore_index=True)

        # 对p值进行多重校正，例如使用Bonferroni校正
        enrichment_results['adjusted_pvalue'] = \
            multipletests(enrichment_results['pvalue'], method='fdr_bh')[1]
        return enrichment_results

    def get_pvalue(self, total_df, pathway, query, n_motifs, shuffle_t=1000, all_path=True,
                   random_seed=0):
        """
        for each pathway, compute the pvalue of enrichment
        ---
        shuffle_t: the times of random picking
        """
        random.seed(random_seed)
        # p_arr = np.zeros(shuffle_t)
        count_cells_with_pathway = self.count_cells_with_pathway_pack(all_path)
        index_list = list(range(len(total_df)))
        pathway_cells_array = np.zeros(shuffle_t)
        for t in range(shuffle_t):
            pick_motifs = random.sample(index_list, n_motifs)
            pathway_cells = count_cells_with_pathway(total_df.loc[pick_motifs], pathway)
            # if pathway_cells < query:  # smaller pvalue, more possible enrichment
            #     p_arr[t] = 1
            pathway_cells_array[t] = pathway_cells
        p_arr = (pathway_cells_array < query).astype(int)
        pvalue = 1 - p_arr.sum() / shuffle_t
        return pvalue

    def compute_pathway_n_pvalue_man(self,df, total_df, pathways_all, all_path=True,
                                     shuffle_t=1000):
        """
        compute pathway enrichment for each color motif
        with procedure written by ourselves
        """
        pvalue_z = 1  # larger pvalue, more possible enrichment, then pvalue_z=0
        enrichment_results = pd.DataFrame(columns=['pathway', 'pvalue'])
        pathways = self.get_pathways(df)
        count_cells_with_pathway = self.count_cells_with_pathway_pack(all_path)
        for pathway in pathways:
            if pd.isnull(pathway) or pathway == "":
                pvalue = pvalue_z
            else:
                # count cells containing this pathway
                query = count_cells_with_pathway(df, pathway)
                pvalue = self.get_pvalue(total_df=total_df, pathway=pathway, query=query,
                                    n_motifs=len(df), shuffle_t=shuffle_t,
                                    all_path=all_path)

            # 将结果添加到结果DataFrame中
            enrichment_results = pd.concat([enrichment_results, pd.DataFrame(
                {'pathway': [pathway], 'pvalue': [pvalue]})], ignore_index=True)

        # add other pvalues
        pathways_res = list(set(pathways_all) - set(pathways))
        enrichment_results = pd.concat([enrichment_results, pd.DataFrame(
            {'pathway': pathways_res, 'pvalue': [pvalue_z] * len(pathways_res)})],
                                       ignore_index=True)

        # 对p值进行多重校正，例如使用Bonferroni校正
        enrichment_results['adjusted_pvalue'] = \
            multipletests(enrichment_results['pvalue'], method='fdr_bh')[1]
        return enrichment_results

    def get_pathways(self,df):
        """
        return a list of pathways (duplicates are droped)
        """
        if len(df) > 1:
            merged_str = ','.join(df['pathway'])
        else:
            merged_str = df.reset_index().loc[0]['pathway']
        pathways = merged_str.split(',')
        pathways = list(set(pathways))
        return pathways

    def filter_enrich_pathways(self, df, p_thres=0.05):
        """
        only remain the rows with adjust pvalue<=p_thres, p_thres=0.05 by default
        """
        df = df[df[self.pvalue_key] <= p_thres]
        return df

    def analyz_pathway_enrich(self, computation="man", all_path=True, shuffle_t=1000,
                              save_df=True, save_dir=None, p_thres=0.05):
        """
        for each color motif, compute pathway enrichment and pvalue
        ---
        computation: {"lib", "man"}
        """
        df = self.out_df
        color_motifs = df[['color_combination', 'pathway']]
        # color_motifs = color_motifs.groupby('color_combination')['pathway'].agg(','.join).reset_index()
        pathways = self.get_pathways(df)
        total_cells = len(df)
        enrich_df = pd.DataFrame()
        for color_combination, group in color_motifs.groupby('color_combination'):
            if computation == "lib":
                enrich_result = self.compute_pathway_n_pvalue_surv_func(df=group,
                                                                   total_df=df,
                                                                   total_cells=total_cells,
                                                                   pathways_all=pathways,
                                                                   all_path=all_path)
            elif computation == "man":
                enrich_result = self.compute_pathway_n_pvalue_man(df=group, total_df=df,
                                                             all_path=all_path,
                                                             pathways_all=pathways,
                                                             shuffle_t=shuffle_t)
            # combine pvalue results
            enrich_df_tmp = pd.concat([pd.DataFrame(
                {'color_combination': [color_combination] * len(enrich_result)}),
                                       enrich_result], axis=1)
            enrich_df = pd.concat([enrich_df, enrich_df_tmp], axis=0)
        enrich_df = enrich_df.reset_index()
        filtered_enrich_df = self.filter_enrich_pathways(enrich_df, p_thres=p_thres)
        if save_df:
            enrich_df.to_csv(save_dir + "/enriched_df.csv", index=False)
            filtered_enrich_df.to_csv(
                save_dir + "enriched_df_filtered_thes=" + str(p_thres) + ".csv",
                index=False)
        return enrich_df, filtered_enrich_df

    def convert_2_cell_weight_graph(self, df,
                                    save_mat=True, weight_mat_dir=None):
        c_type = self.adata.obs[[self.cell_type]]
        _N_ = len(c_type)
        mat = np.zeros((_N_, _N_), dtype=int)
        df = df[self.vertices_header]
        for point_i in range(len(df.columns) - 1):
            col1 = df.columns[point_i]
            # print("point i:", point_i)
            for point_j in range(point_i + 1, len(df.columns)):
                # print("           point j:", point_j)
                col2 = df.columns[point_j]

                if col1 > col2:
                    col1, col2 = col2, col1

                def update_mat(row):
                    col1_value = row[col1]
                    col2_value = row[col2]
                    mat[col1_value][col2_value] += 1

                df.apply(lambda row: update_mat(row), axis=1)
        mat = mat + mat.T
        if save_mat:
            np.savetxt(weight_mat_dir, mat)
        return mat

    def convert_2_cell_weight_graph_v2(self, df,
                                    save_mat=True, weight_mat_dir=None):
        c_type = self.adata.obs[[self.cell_type]]
        _N_ = len(c_type)
        mat = np.zeros((_N_, _N_), dtype=int)
        df = df[self.vertices_header]

        col_pairs = list(combinations(df.columns, 2))

        for col1, col2 in col_pairs:
            col1_values = df[col1].values
            col2_values = df[col2].values

            # Use np.add.at to efficiently accumulate counts in the matrix
            np.add.at(mat, (col1_values, col2_values), 1)

            # Make the matrix symmetric
        mat = mat + mat.T
        if save_mat:
            np.savetxt(weight_mat_dir, mat)
        return mat

    def plot_motif_distribution(self,df):
        """
        only color type
        """
        color_df = df[
            ['vertex_1_color', 'vertex_2_color', 'vertex_3_color', 'vertices_x_avg',
             'vertices_y_avg']].copy()
        color_df['color_combination'] = color_df[
            ['vertex_1_color', 'vertex_2_color', 'vertex_3_color']].apply(
            lambda x: tuple(set(x)), axis=1)
        fig, ax = plt.subplots()

        # 遍历每个颜色组合
        for color_combination, group in color_df.groupby('color_combination'):
            # 提取x和y坐标
            x = group['vertices_x_avg']
            y = group['vertices_y_avg']

            # 绘制散点图，并使用不同的颜色
            ax.scatter(x, y, label=color_combination, s=self.dot_size)

        # 添加图例
        ax.legend()

        # 显示图形
        # plt.show()
        plt.savefig(
            self.save_plots_dir + "/motifs_distribution_.png",
            format="png")

    def plot_motif_distribution_colors(self,color_df, save_figs=False):
        """
            color type and # of color types
            ---
            save_df = True: save color_df for further use
            save_dir : the directory for saving color_df
            save_figs, fig_dir: to save the figures
        """
        cluster_values = color_df['color_combination'].unique()
        colors = plt.cm.tab20c(np.linspace(0, 1, len(cluster_values)))
        color_iterator = iter(colors)
        for color_combination, group in color_df.groupby('color_combination'):
            fig, ax = plt.subplots()
            x = group['vertices_x_avg']
            y = group['vertices_y_avg']

            # 绘制散点图，并使用不同的颜色
            color = next(color_iterator)
            ax.scatter(x, y, label=color_combination, color=color, s=self.dot_size)

            ax.legend()
            plt.xticks(range(self.min_x, self.max_x, 2))
            plt.yticks(range(self.min_y, self.max_y, 2))
            if save_figs:
                plt.savefig(self.save_plots_dir +"/motifs_colors_"+ str(color_combination) + ".png", format="png")
            # plt.show()
            plt.pause(1)
            plt.clf()

        # 遍历每个颜色组合
        fig, ax = plt.subplots()
        color_iterator = iter(colors)
        for color_combination, group in color_df.groupby('color_combination'):
            x = group['vertices_x_avg']
            y = group['vertices_y_avg']
            color = next(color_iterator)
            ax.scatter(x, y, label=color_combination, color=color, s=self.dot_size)
            plt.xticks(range(self.min_x, self.max_x, 2))
            plt.yticks(range(self.min_y, self.max_y, 2))

        ax.legend()
        if save_figs:
            plt.savefig(self.save_plots_dir +"/motifs_colors_combined.png", format="png")
        # plt.show()

    def plot_motif_distribution_colors_pvalue(self,df):
        """
            color type and # of color types
            ---
            plot the distribution of pvalue for each combination of motif
        """
        col_with_color = df.columns[df.columns.str.contains('color')].tolist()
        color_df = df[
            col_with_color + ['vertices_x_avg',
                              'vertices_y_avg', 'p-Value']].copy()
        color_df['color_combination'] = color_df[col_with_color].apply(
            lambda x: tuple(sorted(x.tolist())), axis=1)

        min_x = int(color_df['vertices_x_avg'].min())
        max_x = int(color_df['vertices_x_avg'].max() + 0.5)
        min_y = int(color_df['vertices_y_avg'].min())
        max_y = int(color_df['vertices_y_avg'].max() + 0.5)
        cluster_values = color_df['color_combination'].unique()
        colors = plt.cm.tab20c(np.linspace(0, 1, len(cluster_values)))
        color_iterator = iter(colors)
        for color_combination, group in color_df.groupby('color_combination'):
            fig, ax = plt.subplots()
            x = group['vertices_x_avg']
            y = group['vertices_y_avg']

            # 绘制散点图，并使用不同的颜色
            color = next(color_iterator)
            plt.scatter(x, y, c=group['p-Value'],
                        cmap='viridis')
            plt.colorbar(label='pvalue')
            plt.title('pvalue Distribution: ' + str(color_combination))
            plt.xlabel('avg_x')
            plt.ylabel('avg_y')
            # plt.gcf().set_size_inches(10, 6)
            plt.tight_layout()
            # plt.show()
            plt.savefig(self.save_plots_dir + "/motif_distribution_pvalue_"+ str(color_combination) +".png",
                        format="png")
            plt.pause(1)
            plt.clf()

        # 遍历每个颜色组合
        fig, ax = plt.subplots()
        color_iterator = iter(colors)
        for color_combination, group in color_df.groupby('color_combination'):
            x = group['vertices_x_avg']
            y = group['vertices_y_avg']
            color = next(color_iterator)
            ax.scatter(x, y, label=color_combination, color=color, s=self.dot_size)
            plt.xticks(range(min_x, max_x, 2))
            plt.yticks(range(min_y, max_y, 2))

        ax.legend()
        # plt.show()
        plt.savefig(self.save_plots_dir + "/motif_distribution_pvalue_combined.png", format="png")

    def plot_motif_distribution_cell_types(self, df):
        """
            color type and # of color types
        """
        color_df = df[
            ['vertex_1_type', 'vertex_2_type', 'vertex_3_type', 'vertices_x_avg',
             'vertices_y_avg']].copy()
        color_df['color_combination'] = color_df[
            ['vertex_1_type', 'vertex_2_type', 'vertex_3_type']].apply(
            lambda x: tuple(sorted(x.tolist())), axis=1)

        for color_combination, group in color_df.groupby('color_combination'):
            fig, ax = plt.subplots()
            x = group['vertices_x_avg']
            y = group['vertices_y_avg']

            # 绘制散点图，并使用不同的颜色
            ax.scatter(x, y, label=color_combination, s=self.dot_size)

            ax.legend()
            plt.xticks(range(self.min_x, self.max_x, 2))
            plt.yticks(range(self.min_y, self.max_y, 2))
            # plt.show()
            plt.savefig(self.save_plots_dir + "/motif_distribution_" + str(color_combination) +".png", format="png")
            plt.pause(1)
            plt.clf()

        fig, ax = plt.subplots()
        # 遍历每个颜色组合
        for color_combination, group in color_df.groupby('color_combination'):
            x = group['vertices_x_avg']
            y = group['vertices_y_avg']

            # 绘制散点图，并使用不同的颜色
            ax.scatter(x, y, label=color_combination, s=self.dot_size)

        # ax.legend()
        # plt.show()
        plt.savefig(self.save_plots_dir + "/motif_distribution_combined.png", format="png")

    def generate_color_label(self, cell_gene_clus, save_file=False):
        result = pd.read_csv(cell_gene_clus, index_col=False)
        df = result['optimal_cluster']
        if save_file:
            df.to_csv(self.save_adj_dir_root + "/color_label.csv", index=False, header=None)
        df = pd.DataFrame(df)
        df.columns = [0]
        return df

    def find_vertices(self, df):
        """
        return a list of vertices column name
        """
        ## for colors
        col_with_color = df.columns[df.columns.str.contains('color')].tolist()
        col_with_color.remove('color_combination')
        vertices_col = [s.replace('_color', '') for s in col_with_color]
        ## for xy
        col_with_x = df.columns[df.columns.str.contains('_x')].tolist()
        col_with_y = df.columns[df.columns.str.contains('_y')].tolist()
        cols = col_with_x + col_with_y
        cols.remove('vertices_x_avg')
        cols.remove('vertices_y_avg')
        spot_header = cols
        return vertices_col, spot_header

    def plot_motif_counts_distribution(self):
        colored_df, color_header, xy_cols, num_clus = self.map_colors()
        self.plot_cell_distribution(df=colored_df, vertices_header=color_header, spot_header=xy_cols, label_num=num_clus, save_figs=True)

    def map_colors(self):
        df = self.out_df
        clus_re = self.color_graph_result
        clus_re2 = self.color_graph_result_picked
        vertices_col, spot_header = self.find_vertices(df)
        vertices_color_cols = []
        for vertex_col in vertices_col:
            col_name = vertex_col + "_fin_color_all"
            df[col_name] = df[vertex_col].map(clus_re['optimal_cluster'])
            col_name2 = vertex_col + "_fin_color_picked"
            df[col_name2] = df[vertex_col].map(clus_re2['optimal_cluster'])
            vertices_color_cols.append(col_name)
        return df, vertices_color_cols, spot_header, clus_re['optimal_cluster'].unique()

    def plot_cell_distribution(self, df, vertices_header, spot_header, label_num,
                               save_figs=True):
        """
            color cell distribution with the color label for each combination of motif
            ---
            loc:
            'best'：自动选择最佳位置，尽量不遮挡数据。
            'upper right'：右上角。
            'upper left'：左上角。
            'lower left'：左下角。
            'lower right'：右下角。
            'right'：右侧中间。
            'center left'：左侧中间。
            'center right'：右侧中间。
            'lower center'：底部中间。
            'upper center'：顶部中间。
            'center'：图表中心。
        """
        x_spot = [string for string in spot_header if '_x' in string]
        y_spot = [string for string in spot_header if '_y' in string]
        colors = plt.cm.tab20c(np.linspace(0, 1, len(label_num)))
        # color_iterator = iter(colors)
        for color_combination, group in df.groupby('color_combination'):
            fig, ax = plt.subplots()
            for x, y, color_label in zip(x_spot, y_spot, vertices_header):
                # 绘制散点图，并使用不同的颜色
                # color = next(color_iterator)
                for label in label_num:
                    xs = group[group[color_label] == label][x]
                    ys = group[group[color_label] == label][y]
                    ax.scatter(xs, ys, color=colors[label], s=self.dot_size)
            # unique_values = np.unique(group[vertices_header])
            #
            # # 创建一个空的图例列表
            # legend_elements = []
            #
            # # 遍历唯一取值，为每个取值创建一个图例元素
            # for value in unique_values:
            #     # 创建一个散点图元素，并设置其颜色为对应取值的颜色
            #     element = plt.Line2D([0], [0], marker='o', color='w', label=value,
            #                          markerfacecolor=scatter.cmap(scatter.norm(value)),
            #                          markersize=10)
            #     # 将图例元素添加到图例列表中
            #     legend_elements.append(element)
            #
            # # 创建图例，并设置其位置和标题
            # plt.legend(handles=legend_elements, loc='upper right')
            # # ax.legend()
            unique_values = np.unique(group[vertices_header])
            legend_elements = []

            for value in unique_values:
                element = plt.Line2D([0], [0], marker='o', color='w', label=value,
                                     markerfacecolor=colors[value],
                                     markersize=10)
                legend_elements.append(element)

            plt.legend(handles=legend_elements, loc='best', title="cluster label")

            plt.title('Color label Distribution: ' + str(color_combination))
            plt.xlabel('x')
            plt.ylabel('y')
            plt.xticks(range(self.min_x, self.max_x, 2))
            plt.yticks(range(self.min_y, self.max_y, 2))
            # plt.gcf().set_size_inches(10, 6)
            plt.tight_layout()
            if save_figs:
                plt.savefig(self.save_motif_count_dir + "/cells_in_motifs_"+str(color_combination) + ".png", format="png")
            # plt.show()
            # plt.pause(1)
            plt.clf()
        plt.close()

        # # 遍历每个颜色组合
        fig, ax = plt.subplots()
        color_iterator = iter(colors)
        for color_combination, group in df.groupby('color_combination'):
            x = group['vertices_x_avg']
            y = group['vertices_y_avg']
            color = next(color_iterator)
            ax.scatter(x, y, label=color_combination, color=color, s=self.dot_size)
            plt.xticks(range(self.min_x, self.max_x, 2))
            plt.yticks(range(self.min_y, self.max_y, 2))

        ax.legend()
        # plt.show()
        if save_figs:
            plt.savefig(self.save_motif_count_dir + "/motif_counts_clustering.png", format="png")

    def plot_only_cells_in_motifs(self, max_color_label, cluster_id, save_figs=True):
        colored_df = self.out_df
        col_with_x = colored_df.columns[colored_df.columns.str.contains('_x')].tolist()
        col_with_y = colored_df.columns[colored_df.columns.str.contains('_y')].tolist()
        col_with_x.remove('vertices_x_avg')
        col_with_y.remove('vertices_y_avg')
        col_with_color = colored_df.columns[
            colored_df.columns.str.contains('fin_color_picked')].tolist()
        cluster_data = colored_df[colored_df['color_combination'] == cluster_id]
        unique_values = np.unique(cluster_data[col_with_color])
        colors = plt.cm.tab20c(np.linspace(0, 1, max_color_label + 1))
        fig, ax = plt.subplots()
        # color_iterator = iter(colors)
        for x, y, color in zip(col_with_x, col_with_y, col_with_color):
            for value in zip(unique_values):
                cluster_data_sub = cluster_data[cluster_data[color] == value]
                ax.scatter(cluster_data_sub[x], cluster_data_sub[y],
                           color=colors[value], s=self.dot_size)

        plt.xticks(range(self.min_x, self.max_x + 1, 2))
        plt.yticks(range(self.min_y, self.max_y + 1, 2))

        legend_elements = []

        for value in unique_values:
            # element = plt.Line2D([0], [0], marker='o', color='w', label=value,
            #                      markerfacecolor=ax.collections[0].cmap(
            #                          ax.collections[0].norm(value)),
            #                      markersize=10)
            element = plt.Line2D([0], [0], marker='o', color='w', label=value,
                                 markerfacecolor=colors[value],
                                 markersize=10)
            legend_elements.append(element)

        plt.legend(handles=legend_elements, loc='best', title='color label')
        plt.title("cells_in_motifs_" + str(cluster_id))
        if save_figs:
            plt.savefig(self.save_plots_dir + "/only_cells_in_motifs_" + str(cluster_id) + ".png", format="png")
        # plt.show()
        plt.pause(1)

    def plot_only_cells_in_motifs_cell_type(self, cluster_id, save_figs=True):
        colored_df = self.out_df
        col_with_x = colored_df.columns[colored_df.columns.str.contains('_x')].tolist()
        col_with_y = colored_df.columns[colored_df.columns.str.contains('_y')].tolist()
        col_with_x.remove('vertices_x_avg')
        col_with_y.remove('vertices_y_avg')
        col_with_type = colored_df.columns[
            colored_df.columns.str.contains('_type')].tolist()
        cluster_data = colored_df[colored_df['color_combination'] == cluster_id]
        unique_values = np.unique(cluster_data[col_with_type])
        fig, ax = plt.subplots()
        types = []
        for col_t in col_with_type:
            types.extend(list(colored_df[col_t].unique()))
        types = list(set(types))
        colors = plt.cm.tab20c(np.linspace(0, 1, len(types)))
        # color_iterator = iter(colors)
        for x, y, type in zip(col_with_x, col_with_y, col_with_type):
            for value in types:
                cluster_data_sub = cluster_data[cluster_data[type] == value]
                ax.scatter(cluster_data_sub[x], cluster_data_sub[y],
                           color=colors[types.index(value)], s=self.dot_size)

        plt.xticks(range(self.min_x, self.max_x + 1, 2))
        plt.yticks(range(self.min_y, self.max_y + 1, 2))

        legend_elements = []

        for value in unique_values:
            # element = plt.Line2D([0], [0], marker='o', color='w', label=value,
            #                      markerfacecolor=ax.collections[0].cmap(
            #                          ax.collections[0].norm(value)),
            #                      markersize=10)
            element = plt.Line2D([0], [0], marker='o', color='w', label=value,
                                 markerfacecolor=colors[types.index(value)],
                                 markersize=10)
            legend_elements.append(element)

        plt.legend(handles=legend_elements, loc='best', title='cell type')
        plt.title("cell type distribution of color labels " + str(cluster_id))
        if save_figs:
            plt.savefig(self.save_plots_dir + "/" + str(cluster_id) + ".png", format="png")
        # plt.show()

    def extract_numbers(self, string):
        """
        extract # of cluster from string, output the list,
        for clusters built motifs
        """
        numbers = re.findall(r'\d+', string)
        return [int(num) for num in numbers]

    def extract_numbers_v2(self, cluster_id_t, keyw="_fin_color_all"):
        """
        extract # of cluster from string, output the list,
        for cell type built motifs
        """
        lines = self.out_df[self.out_df['color_combination']==cluster_id_t]
        headers = lines.columns[lines.columns.str.contains(keyw)].tolist()
        ids = lines[headers].values.ravel()
        ids = list(set(ids))
        return ids

    def extract_vertices_ids(self, cluster_id_t):
        lines = self.out_df[self.out_df['color_combination'] == cluster_id_t]
        vertices_ids = lines[self.vertices_header].values.ravel()
        vertices_ids = list(set(vertices_ids))
        return vertices_ids

    def plot_all_cells_in_motifs(self, max_color_label, cluster_id):
        result = self.color_graph_result
        spot = self.adata.obs[['x', 'y']]
        colors = plt.cm.tab20c(np.linspace(0, 1, max_color_label + 1))

        fig, ax = plt.subplots()
        keyw = "_fin_color_all"
        cluster_ids = self.extract_numbers_v2(cluster_id, keyw=keyw)
        vertices_ids = self.extract_vertices_ids(cluster_id)
        result = result.loc[vertices_ids]
        for cluster_value, color in zip(cluster_ids, colors):
            cluster_data = spot.loc[result[result['optimal_cluster'] == cluster_value].index]
            ax.scatter(cluster_data['x'], cluster_data['y'], color=color,
                       label=cluster_value, s=self.dot_size)

        plt.xticks(range(self.min_x, self.max_x + 1, 2))
        plt.yticks(range(self.min_y, self.max_y + 1, 2))
        ax.legend(loc="best", title="color label")
        # plt.show()
        plt.savefig(self.save_plots_dir + "/all_cells_in_motifs_"+str(cluster_id)+".png", format="png")
        plt.pause(1)

    def analyz_only_cells_in_motifs(self, save_figs=True):
        """
        plot all cells in motifs,
        the cluster result of color graph is labeled with color,
        only cells in motifs are clustered
        """
        result = self.color_graph_result_picked
        colored_df = self.out_df
        max_color_label = result['optimal_cluster'].max()
        cluster_ids = colored_df[
            'color_combination'].unique()  # extract all the cluster_id
        for cluster_id in cluster_ids:
            self.plot_only_cells_in_motifs(max_color_label=max_color_label,
                                      cluster_id=cluster_id, save_figs=save_figs)

    def analyz_all_cells_in_motifs(self):
        """
        plot all cells in motifs,
        the cluster result of color graph is labeled with color,
        all cells are clustered
        """
        result = self.color_graph_result
        colored_df = self.out_df
        max_color_label = result['optimal_cluster'].max()
        cluster_ids = colored_df[
            'color_combination'].unique()  # extract all the cluster_id
        for cluster_id in cluster_ids:
            self.plot_all_cells_in_motifs(
                                     max_color_label=max_color_label,
                                     cluster_id=cluster_id)

    def analyz_only_cells_in_motifs_cell_types(self, save_figs=True):
        """
        plot only cells in motifs, which that means cells not in motifs are ignored there,
        the cell type is labeled with color
        """
        colored_df = self.out_df
        cluster_ids = colored_df[
            'color_combination'].unique()  # extract all the cluster_id
        cluster_ids.sort()
        for cluster_id in cluster_ids:
            self.plot_only_cells_in_motifs_cell_type(cluster_id=cluster_id,
                                                save_figs=save_figs)

    def proc_n_annotate(self, pick_cells=True, plot_volcano=True, save_df=True):
        """
        gene and pathway enrichment of cell clusters of color graph
        only cells in motifs should be considered
        """
        if pick_cells:
            adata = self.adata[self.cells_in_motifs]
            dff = self.color_graph_result_picked.loc[self.cells_in_motifs]['optimal_cluster']
            dff.index = adata.obs.index
            adata.obs = pd.concat([adata.obs,dff],axis=1)
        else:
            adata = self.adata
            dff = self.color_graph_result['optimal_cluster']
            dff.index = adata.obs.index
            adata.obs['optimal_cluster'] = dff
        adata.obs['optimal_cluster'] = adata.obs['optimal_cluster'].astype('category')
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.tl.rank_genes_groups(adata, 'optimal_cluster', method='t-test',
                                key_added="t-test")

        # generate dataframe for pathway enrichment
        key_ids = [item[0] for item in adata.uns['t-test']['pvals'].dtype.descr]

        for key_id in key_ids:
            enrich_df = pd.DataFrame()
            enrich_df['gene_name'] = adata.uns['t-test']['names'][key_id]
            enrich_df['score'] = adata.uns['t-test']['scores'][key_id]
            enrich_df['pvals'] = adata.uns['t-test']['pvals'][key_id]
            enrich_df['pvals_adj'] = adata.uns['t-test']['pvals_adj'][key_id]
            enrich_df['logfoldchanges'] = adata.uns['t-test']['logfoldchanges'][key_id]

            # sort and pick gene set
            # pvals_adj = enrich_df['pvals_adj']
            # names = enrich_df['names']
            # names_sorted_by_pvals = [name for _, name in sorted(zip(pvals_adj, names))]
            gene_list = enrich_df[enrich_df[self.pval_key_enrich] < 0.05]['gene_name'].tolist()
            # if gene_list is empty
            if gene_list == []:
                df_t = pd.DataFrame()
            else:
                enr = gp.enrichr(gene_list=gene_list, gene_sets=self.gene_sets, organism='human')
                df_t = pd.concat(
                    [pd.DataFrame({"cluster id combination": [key_id] * len(enr.res2d)}),
                     enr.res2d],
                    axis=1)
            # save
            if save_df:
                df_t.to_csv(self.save_enrich_dir+'/'+"color_motifs_enrich_" + str(key_id) + ".csv", index=False)
                with open(self.save_enrich_dir+'/'+"color_motifs_enrich_gene_"+ str(key_id) +".txt", 'w') as file:
                    # 将列表中的每个元素写入文件的不同行
                    for item in gene_list:
                        file.write(item + '\n')

            # plot volcano plot
            # plot paras
            if plot_volcano:
                x_threshold_left = -1
                x_threshold_right = 1
                y_threshold = 1
                xmin = -6
                xmax = 6 + 2
                ymin = -1
                ymax = 3 + 2
                result = pd.DataFrame()
                result['x'] = enrich_df['logfoldchanges']
                result['y'] = enrich_df[self.pval_key_enrich]
                result['y'] = -np.log10(result['y'])
                # result['x'] = result['x'].apply(lambda x: np.exp2(x))

                result['group'] = 'black'
                result.loc[(result.x > x_threshold_right) & (
                        result.y > y_threshold), 'group'] = 'tab:red'  # x=-+x_threshold直接截断
                result.loc[(result.x < x_threshold_left) & (
                        result.y > y_threshold), 'group'] = 'tab:blue'  # x=-+x_threshold直接截断
                result.loc[result.y < y_threshold, 'group'] = 'dimgrey'

                fig = plt.figure(figsize=plt.figaspect(7 / 6))
                ax = fig.add_subplot()
                ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax), title='')
                ax.scatter(result['x'], result['y'], s=2, c=result['group'])
                ax.set_ylabel('-Log10(p value)', fontweight='bold')
                ax.set_xlabel('Log2 (fold change)', fontweight='bold')
                ax.spines['right'].set_visible(False)
                ax.spines['top'].set_visible(False)

                # 水平和竖直线
                ax.vlines(x_threshold_left, ymin, ymax, color='dimgrey', linestyle='dashed',
                          linewidth=1)
                ax.vlines(x_threshold_right, ymin, ymax, color='dimgrey',
                          linestyle='dashed',
                          linewidth=1)
                ax.hlines(y_threshold, xmin, xmax, color='dimgrey', linestyle='dashed',
                          linewidth=1)

                ax.set_xticks(range(xmin, xmax, 2))
                ax.set_yticks(range(ymin, ymax, 2))
                plt.title("color motifs enrich: "+str(key_id))
                # plt.show()
                plt.savefig(self.save_enrich_dir+'/'+"color_motifs_enrich_pathway_" + str(key_id) + ".png", format="png")
                plt.close()
