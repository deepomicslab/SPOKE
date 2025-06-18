"""
can be imported as package
---
2024/06/10
"""
import pandas as pd
import anndata
import scanpy as sc
from collections import OrderedDict
import numpy as np
import matplotlib.pyplot as plt
import os

import copy
import pdb


def remove_duplicates(text):
    words = text.split(',')
    words = list(OrderedDict.fromkeys(words))
    return ','.join(words)


class Analysis_Interface:
    """
    analysis class for interface
    """
    def __init__(self, idata_file, result_file_dir, lr_pathway_db_file_dir="/public/qiusliang2/sc_br/lr_pathway_meta.csv", sta_method='t-test', plot_volcano=True):
        if isinstance(idata_file, str):
            self.idata = anndata.read_h5ad(idata_file)
        else:
            self.idata = idata_file
        self.pre_proc()
        self.lr_pathway_db_file = pd.read_csv(lr_pathway_db_file_dir, index_col=False)
        self.result = pd.read_csv(result_file_dir, index_col=False)
        self.lr_meta = self.idata.var
        self.sort_rule = 'pvals_adj'
        self.pval_key = 'pvals'
        self.sta_method = sta_method

        # plot paras
        self.plot_volcano = plot_volcano
        self.x_threshold_left = -1
        self.x_threshold_right = 1
        self.y_threshold = 1
        self.xmin = -6
        self.xmax = 6 + 2
        self.ymin = -1
        self.ymax = 3 + 2

        # preparation
        self.lr_meta_pathway = self.prepare_lr_db()
        self.idata.obs['cluster_id'] = self.result.set_index(self.idata.obs.index)['optimal_cluster']
        self.adata_lr = self.compute_for_lr()

    def pre_proc(self):
        mat = self.idata.X.todense().A
        self.idata.X[mat < 0] = 0

    def compute_for_lr(self):
        adata = copy.deepcopy(self.idata)
        # process
        adata.obs['cluster_id'] = adata.obs['cluster_id'].astype('category')
        sc.pp.normalize_total(adata, target_sum=1e4)  # target_sum can be varied
        # sc.pp.log1p(adata)
        sc.tl.rank_genes_groups(adata, 'cluster_id', method=self.sta_method,
                                key_added=self.sta_method)
        return adata

    def prepare_lr_db(self):
        # generate pathway info
        lr_pathway_db = self.lr_pathway_db_file
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
        lr_meta_pathway = pd.merge(self.lr_meta, lr_pathway_db4, how='left',
                                   left_on='receptor', right_on='src')
        lr_meta_pathway['pathway'].fillna('nan').tolist()
        return lr_meta_pathway

    def proc_lr(self, save_lr_enrich=True, lr_enrich_dir=None):
        lr_meta_pathway2 = self.lr_meta_pathway.set_index(
            self.adata_lr.var.index)  # set index to synchronize
        self.adata_lr.var['pathway'] = lr_meta_pathway2['pathway']
        key_ids = [item[0] for item in self.adata_lr.uns[self.sta_method]['pvals'].dtype.descr]
        if not lr_enrich_dir is None and not os.path.exists(lr_enrich_dir):
            os.makedirs(lr_enrich_dir, exist_ok=True)
        for key_id in key_ids:
            lr_t = lr_meta_pathway2.loc[list(self.adata_lr.uns[self.sta_method]['names'][key_id])]
            lr_t['score'] = self.adata_lr.uns[self.sta_method]['scores'][key_id]
            lr_t['pvals'] = self.adata_lr.uns[self.sta_method]['pvals'][key_id]
            lr_t['pvals_adj'] = self.adata_lr.uns[self.sta_method]['pvals_adj'][key_id]
            lr_t['logfoldchanges'] = self.adata_lr.uns[self.sta_method]['logfoldchanges'][key_id]
            if save_lr_enrich:
                lr_t.to_csv(lr_enrich_dir +'/'+"infc_enrich_"+ str(key_id) + ".csv", index=False)

            if self.plot_volcano:
                # plot volcano plot
                result = pd.DataFrame()
                result['x'] = lr_t['logfoldchanges']
                result['y'] = lr_t[self.pval_key]
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
                ax.set(xlim=(self.xmin, self.xmax), ylim=(self.ymin, self.ymax), title='')
                ax.scatter(result['x'], result['y'], s=2, c=result['group'])
                ax.set_ylabel('-Log10(p value)', fontweight='bold')
                ax.set_xlabel('Log2 (fold change)', fontweight='bold')
                ax.spines['right'].set_visible(False)
                ax.spines['top'].set_visible(False)

                # 水平和竖直线
                ax.vlines(self.x_threshold_left, self.ymin, self.ymax, color='dimgrey', linestyle='dashed',
                          linewidth=1)
                ax.vlines(self.x_threshold_right, self.ymin, self.ymax, color='dimgrey',
                          linestyle='dashed',
                          linewidth=1)
                ax.hlines(self.y_threshold, self.xmin, self.xmax, color='dimgrey', linestyle='dashed',
                          linewidth=1)

                ax.set_xticks(range(self.xmin, self.xmax, 2))
                ax.set_yticks(range(self.ymin, self.ymax, 2))
                plt.title(str(key_id))
                plt.savefig(lr_enrich_dir +'/'+"infc_enrich_"+ str(key_id) + ".png", format="png")

    def proc_pathway(self, save_pathway_enrich=True, pathway_enrich_dir=None):
        # for pathway analysis
        # concat cols with same pathway info
        lr_meta_pathway2 = self.lr_meta_pathway.set_index(
            self.adata_lr.var.index)  # set index to synchronize
        adata = copy.deepcopy(self.idata)
        adata.var['pathway'] = lr_meta_pathway2['pathway']
        X_n = np.zeros((len(adata.obs), len(adata.var['pathway'].unique())))
        pathways = adata.var['pathway'].unique()
        for i, pathway in enumerate(pathways):
            rows = adata.var['pathway'] == pathway
            X_n[:, i] = adata.X.todense().A[:, adata.var[rows]['i']].sum(axis=1)
        adata1 = anndata.AnnData(X=X_n, var=pd.DataFrame({'pathway': pathways}),
                                 obs=adata.obs)
        adata = adata1

        # process
        adata.obs['cluster_id'] = adata.obs['cluster_id'].astype('category')
        sc.pp.normalize_total(adata, target_sum=1e4)  # target_sum can be varied
        sc.pp.log1p(adata)
        sc.tl.rank_genes_groups(adata, 'cluster_id', method=self.sta_method,
                                key_added=self.sta_method)
        key_ids = [item[0] for item in adata.uns[self.sta_method]['pvals'].dtype.descr]
        if not pathway_enrich_dir is None and not os.path.exists(pathway_enrich_dir):
            os.makedirs(pathway_enrich_dir,exist_ok=True)
        for key_id in key_ids:
            pathway_t = adata.var.loc[list(adata.uns[self.sta_method]['names'][key_id])]
            pathway_t['score'] = adata.uns[self.sta_method]['scores'][key_id]
            pathway_t['pvals'] = adata.uns[self.sta_method]['pvals'][key_id]
            pathway_t['pvals_adj'] = adata.uns[self.sta_method]['pvals_adj'][key_id]
            pathway_t['logfoldchanges'] = adata.uns[self.sta_method]['logfoldchanges'][key_id]
            pathway_t['cluster_id'] = key_id
            if save_pathway_enrich:
                pathway_t.to_csv(pathway_enrich_dir + "pathway_" + str(key_id) + ".csv",
                             index=True)
