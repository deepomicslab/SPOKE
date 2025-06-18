"""
---
2024/06/12
reformat as pkg version
---
2024/04/17
add: plot volcano
---
2024/04/16
analyze pathway enrichment on boundary cells with scanpy
"""
import pdb
import os
import anndata
from anndata import AnnData
import pandas as pd
import scanpy as sc
import gseapy as gp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import seaborn as sns
import sys
sys.path.append("..")
from src.boundary_cell_finder import Boundary_Cell_Finder


def compute_pathway_enrich(df, group_id='counts'):
    """
    for each cluster, generate the gene_list, and compute pathway enrichment
    """
    enrich_df = pd.DataFrame()
    for g_id, group in df.groupby(group_id):
        gene_list = list(group['gene_name'])
        # sort dataframe according to the "adjusted_pvakue"
        # gene_set: https://maayanlab.cloud/Enrichr/#libraries
        gene_sets = ['KEGG_2016', 'KEGG_2013']
        enr = gp.enrichr(gene_list=gene_list, gene_sets=gene_sets, organism='human')
        df_t = pd.concat([pd.DataFrame({group_id: [g_id] * len(enr.res2d)}),enr.res2d], axis=1)
        enrich_df = pd.concat([enrich_df, df_t], axis=0)
    enrich_df = enrich_df.reset_index()
    return enrich_df


class Analysis_Boundary_Cell:
    """
    analysis class for interface, boundary cell part
    ---
    idata: interface data
    adata: cell data
    """
    def __init__(self, idata_file, adata_file, result_file_dir, lr_pathway_db_file_dir="/public/qiusliang2/sc_br/lr_pathway_meta.csv", celltype_key="cell_type", x_row="x", y_col="y", save_df=True, save_df_dir=None, plot_volcano=False, save_pathway_enrich=None):
        """
        Initialize with finding boundary cells
        """
        if isinstance(idata_file, str):
            self.idata = anndata.read_h5ad(idata_file)
        else:
            self.idata = idata_file
        if isinstance(adata_file, str):
            self.adata = anndata.read_h5ad(adata_file)
        else:
            self.adata = adata_file
        if 'barcode' in self.adata.obs.columns:
            self.adata.obs['barcode1'] = self.adata.obs['barcode']
            self.adata.obs.drop('barcode', axis=1, inplace=True)
        self.adata.obs.reset_index(inplace=True)
        self.save_df = save_df
        self.bc_finder = Boundary_Cell_Finder(result_file_dir, self.adata, idata_file, lr_pathway_db_file_dir=lr_pathway_db_file_dir, celltype_key=celltype_key, x_row=x_row, y_col=y_col, save_df=self.save_df, save_df_dir=save_df_dir)
        self.boundary_cells = self.bc_finder.boundary_cell_df
        self.optional_sort_rules = ['score', 'pvals', 'pvals_adj', 'logfoldchanges',
                               'gene_name']
        self.sort_rule = 'pvals_adj'
        self.pval_key = 'pvals'
        self.gene_sets = ['KEGG_2016', 'KEGG_2013']
        self.plot_volcano = plot_volcano
        self.save_pathway_enrich = save_pathway_enrich
        if not save_pathway_enrich is None and not os.path.exists(save_pathway_enrich):
            os.makedirs(save_pathway_enrich, exist_ok=True)
        self.save_df_dir = save_df_dir
        self.celltype_key = celltype_key
        self.gene_set = []

        # # for debug
        # self.boundary_cells = pd.read_csv(
        #     "/home/grads/qiusliang2/sctipts6/colormotifs_mela2/infc_bc/bc/boundary_cells_df.csv",
        #     index_col=False)
        # lines = [line.strip() for line in open(
        #     "/home/grads/qiusliang2/sctipts6/colormotifs_mela2/infc_bc/bc/output.txt")]
        # self.gene_set = list(set(lines))
        # self.plot_gene_heatmap()

        self.proc_n_annotate()
        self.gene_set = list(set(self.gene_set))
        self.gene_exp_map = None

    def plot_gene_heatmap(self, mode="avg", exp_threshold_por=0.65, count_threshold=3):
        mat = self.adata.X.todense().A
        # get related rows and columns
        bc_ids = self.boundary_cells['cluster_id'].unique()
        col_ids = [list(self.adata.var_names).index(x) for x in self.gene_set]
        row_vectors_list = []
        for bc_id in bc_ids:
            row_ids = self.boundary_cells[self.boundary_cells['cluster_id']==bc_id].index
            mat_c = np.sum(mat[row_ids,:][:,col_ids], axis=0)
            if mode == "avg":
                mat_c = np.mean(mat[row_ids,:][:,col_ids], axis=0)
            row_vectors_list.append(mat_c)
        matrix = np.vstack(row_vectors_list)
        df = pd.DataFrame(matrix.T, columns=bc_ids, index=list(self.gene_set))
        self.gene_exp_map = df
        # PLOT
        # plt.figure(figsize=(10, 20))
        # cols = df.T.columns[(df.T > exp_threshold).sum() >= count_threshold]
        # df_filtered = df.loc[cols]
        exp_threshold = df.max().mean() * exp_threshold_por
        df_filtered = df.loc[df.gt(exp_threshold).any(axis=1), :]
        sns.clustermap(df_filtered, row_cluster=True,col_cluster=True,figsize=(10,20))
        # plt.title('Gene Expression Heatmap')
        # plt.xlabel('Cells')
        # plt.ylabel('Genes')
        plt.savefig(self.save_df_dir+'/'+"bc_enrich_gene_exp"+str(exp_threshold_por)+".png")


    def proc_n_annotate(self):
        self.boundary_cells['barcode'] = self.boundary_cells['cell'].astype(str)
        self.boundary_cells.set_index('barcode', inplace=True)
        self.boundary_cells.index = self.boundary_cells.index.astype(int)
        adata_index1 = self.adata.obs[
            self.adata.obs.index.isin(self.boundary_cells['cell'].tolist())]
        adata = self.adata[adata_index1.index]
        # avoid same col names in adata.obs and self.boundary_cells
        keywords = ['x', 'y', self.celltype_key]
        columns_with_keywords = []
        for keyword in keywords:
            if keyword in adata.obs.columns and keyword in self.boundary_cells.columns:
                columns_with_keywords.append(keyword)
        if len(columns_with_keywords)>0:
            dft = adata.obs.join(self.boundary_cells.drop(columns_with_keywords,axis=1))
        else:
            dft = adata.obs.join(self.boundary_cells)
        adata2 = AnnData(X=adata.X, obs=dft, var=adata.var)
        adata2.uns = adata.uns
        adata2.obsm = adata.obsm
        adata2.varm = adata.varm
        adata2.layers = adata.layers
        adata2.obs['cluster_id'].apply(lambda x: ', '.join(str(i) for i in x))
        adata = adata2
        # fixed: only use boundary cells for further analysis
        # adata.obs = pd.merge(self.boundary_cells, adata.obs, left_index=True, right_index=True, how='left')
        # adata.obs.join(self.boundary_cells)
        # fileter: group of boundary cells with the same cluster id combination contains only one sample
        index_list = list(adata.obs[adata.obs.duplicated(subset='cluster_id', keep=False)].index)
        adata_index = adata.obs[adata.obs.index.isin(index_list)]
        adata = adata[adata_index.index]
        adata.obs['cluster_id'] = adata.obs['cluster_id'].astype('category')
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.tl.rank_genes_groups(adata, 'cluster_id', method='t-test',
                                key_added="t-test")

        # generate dataframe for pathway enrichment
        key_ids = [item[0] for item in adata.uns['t-test']['pvals'].dtype.descr]
        gene_back_set = set()
        for key_id in key_ids:
            gene_back_set.update(set(adata.uns['t-test']['names'][key_id]))
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
            gene_list = enrich_df[enrich_df[self.pval_key] < 0.05]['gene_name'].tolist()
            # if gene_list is empty
            if gene_list == []:
                df_t = pd.DataFrame()
            else:
                self.gene_set.extend(gene_list)
                enr = gp.enrichr(gene_list=gene_list, gene_sets=self.gene_sets, organism='human',background=gene_back_set)
                df_t = pd.concat(
                    [pd.DataFrame({"cluster id combination": [key_id] * len(enr.res2d)}),
                     enr.res2d],
                    axis=1)
            # save
            if self.save_df:
                df_t.to_csv(self.save_df_dir+'/'+"bc_enrich_" + str(key_id) + ".csv", index=False)
                with open(self.save_df_dir+'/'+"bc_enrich_gene_"+ str(key_id) +".txt", 'w') as file:
                    # 将列表中的每个元素写入文件的不同行
                    for item in gene_list:
                        file.write(item + '\n')

            # plot volcano plot
            # plot paras
            if self.plot_volcano:
                x_threshold_left = -1
                x_threshold_right = 1
                y_threshold = 1
                xmin = -6
                xmax = 6 + 2
                ymin = -1
                ymax = 3 + 2
                result = pd.DataFrame()
                result['x'] = enrich_df['logfoldchanges']
                result['y'] = enrich_df['pvals_adj']
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
                plt.title(str(key_id))
                # plt.show()
                plt.savefig(self.save_pathway_enrich+'/'+"bc_enrich_" + str(key_id) + ".png", format="png")
                plt.close()



