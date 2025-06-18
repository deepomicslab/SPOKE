import pdb
import os
import anndata
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from collections import OrderedDict


def remove_duplicates(text):
    words = text.split(',')
    words = list(OrderedDict.fromkeys(words))
    return ','.join(words)


class Boundary_Cell_Finder:
    """
        to find boundary cells for interface data
    """
    def __init__(self, result_file, adata_file, idata_file, lr_pathway_db_file_dir="/public/qiusliang2/sc_br/lr_pathway_meta.csv", celltype_key="cell_type", x_row="x", y_col="y", save_df=False, save_df_dir=None, plot_unit=5):
        self.x_row = x_row
        self.y_col = y_col
        self.save_df = save_df
        self.save_df_dir = save_df_dir
        if not save_df_dir is None and not os.path.exists(save_df_dir):
            os.makedirs(save_df_dir, exist_ok=True)
        self.result = pd.read_csv(result_file, index_col=False)
        if isinstance(idata_file, str):
            self.idata = anndata.read_h5ad(idata_file)
        else:
            self.idata = idata_file
        self.lr_meta = self.idata.var
        self.lr_pathway_db_file = pd.read_csv(lr_pathway_db_file_dir, index_col=False)
        self.interface_meta = self.prepare_lr_db()

        if isinstance(adata_file, str):
            self.adata = anndata.read_h5ad(adata_file)
        else:
            self.adata = adata_file
        self.cell_type = self.adata.obs[celltype_key]
        self.celltype_key = celltype_key
        self.cell_info = pd.DataFrame(self.adata.obsm['spatial'])
        self.plot_unit = plot_unit
        self.min_x = int(self.cell_info[0].min()) - self.plot_unit
        self.max_x = int(self.cell_info[0].max()) + self.plot_unit
        self.min_y = int(self.cell_info[1].min()) - self.plot_unit
        self.max_y = int(self.cell_info[1].max()) + self.plot_unit
        self.cell_info.columns=[self.x_row, self.y_col]
        self.boundary_cell_df = self.get_boundary_cells_df()
        return

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
        lr_meta_pathway2 = lr_meta_pathway.set_index(
            self.idata.var.index)  # set index to synchronize
        self.idata.var['pathway'] = lr_meta_pathway2['pathway']

        mat = self.idata.X.todense().A
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
            lr_pathways['pathway'] = remove_duplicates(lr_pathways['pathway'])
            lr_pathways['source'] = remove_duplicates(lr_pathways['source'])
            lr_pathways_df = pd.concat([lr_pathways_df, lr_pathways], axis=1)
        lr_pathways_df_t_out = lr_pathways_df.transpose()
        interface_cell = self.idata.obs
        interface_cell.reset_index(inplace=True)
        interface_info_merged = pd.concat(
            [interface_cell, lr_pathways_df_t_out.reset_index(drop=True)], axis=1)
        interface_info_merged = interface_info_merged.rename(
            columns={"i": "interface_i"})
        return interface_info_merged

    def get_boundary_cells_df(self):
        # 2. read interface_meta_file
        df = self.get_boundary_cells()
        # 3. map spots
        df = self.map_boundary_cells(df=df)
        # 4. map cell types
        df = self.map_cell_type(df=df)
        if self.save_df:
            df.to_csv(self.save_df_dir + "boundary_cells_df.csv", index=False)
        return df

    def get_boundary_cells(self):
        """output cells with multiple cluster id """
        self.interface_meta['optimal_cluster'] = self.result['optimal_cluster']
        cells_info = pd.concat([self.interface_meta[['sender', 'optimal_cluster']].rename(
            columns={'sender': 'cell', 'optimal_cluster': 'cluster_id'}),
                                self.interface_meta[['receiver', 'optimal_cluster']].rename(
                                    columns={'receiver': 'cell',
                                             'optimal_cluster': 'cluster_id'})], axis=0)
        df = cells_info.drop_duplicates()
        df.loc[:, 'count'] = df.groupby('cell')['cell'].transform('count')
        # df['count'] = df.groupby('cell')['cell'].transform('count')

        # 过滤掉只在一行中出现的'cell'值
        df = df[df['count'] > 1]

        # 合并相同'cell'值的行，连接'cluster_id'值为tuple
        df['cell'] = df['cell'].astype(int)
        df = df.groupby('cell').agg({'cluster_id': tuple}).reset_index()
        df['cluster_id'] = df['cluster_id'].apply(lambda x: tuple(sorted(x)))  # sort
        df['counts'] = df['cluster_id'].apply(lambda x: len(x))

        # # delete count column
        # df = df.drop('count', axis=1)
        return df

    def map_boundary_cells(self, df):
        """
        map spots to the df
        """
        xy = self.cell_info
        xy = pd.DataFrame(xy).reset_index()
        spot_name_x = "x"
        spot_name_y = "y"
        col = "cell"
        df[col] = df[col].astype(int)
        df[spot_name_x] = df[col].map(xy[self.x_row])
        df[spot_name_y] = df[col].map(xy[self.y_col])
        return df

    def map_cell_type(self, df):
        cell_type_info = self.cell_type
        cell_type_info = pd.DataFrame(cell_type_info).reset_index()
        df['cell'] = df['cell'].astype(int)
        df['cell_type'] = df['cell'].map(cell_type_info[self.celltype_key])
        return df

    def plot_boundary_cells_all(self, save_figs=False, fig_dir=None):
        """
        plot the distribution of boundary cells
        simply plot all the boundary cells
        """
        if not fig_dir is None and not os.path.exists(fig_dir):
            os.makedirs(fig_dir, exist_ok=True)
        df = self.boundary_cell_df
        min_x = self.min_x
        max_x = self.max_x
        min_y = self.min_y
        max_y = self.max_y
        fig, ax = plt.subplots()
        x = df['x']
        y = df['y']
        ax.scatter(x, y)
        plt.xticks(range(min_x, max_x, self.plot_unit))
        plt.yticks(range(min_y, max_y, self.plot_unit))
        if save_figs:
            plt.savefig(fig_dir + "boundary_cells_all.png", format="png")
        plt.show()
        plt.pause(1)
        plt.clf()

    def plot_boundary_cell_bounds(self, save_figs=False, fig_dir=None):
        """
        plot the distribution of boundary cells
        according to the number of cluster id the cell occupies
        """
        if not fig_dir is None and not os.path.exists(fig_dir):
            os.makedirs(fig_dir, exist_ok=True)
        df = self.boundary_cell_df
        min_x = self.min_x
        max_x = self.max_x
        min_y = self.min_y
        max_y = self.max_y
        cluster_values = df['counts'].unique()
        colors = plt.cm.tab20c(np.linspace(0, 1, len(cluster_values) + 1))
        color_iterator = iter(colors)
        for combination, group in df.groupby('counts'):
            fig, ax = plt.subplots()
            x = group['x']
            y = group['y']

            # 绘制散点图，并使用不同的颜色
            color = next(color_iterator)
            ax.scatter(x, y, label=str(combination), color=color)

            ax.legend()
            plt.xticks(range(min_x, max_x, self.plot_unit))
            plt.yticks(range(min_y, max_y, self.plot_unit))
            if save_figs:
                plt.savefig(fig_dir + str(combination) + ".png", format="png")
            plt.show()
            plt.pause(1)
            plt.clf()

        # 遍历每个颜色组合
        fig, ax = plt.subplots()
        color_iterator = iter(colors)
        for color_combination, group in df.groupby('counts'):
            x = group['x']
            y = group['y']
            color = next(color_iterator)
            ax.scatter(x, y, label=color_combination, color=color)
            plt.xticks(range(min_x, max_x, self.plot_unit))
            plt.yticks(range(min_y, max_y, self.plot_unit))

        ax.legend()
        if save_figs:
            plt.savefig(fig_dir + "combine.png", format="png")
        plt.show()
        plt.pause(1)
        plt.clf()

    def plot_boundary_cell_cluster(self, save_figs=False, fig_dir=None, figsize=(10, 6)):
        """
        plot the distribution of boundary cells
        according to the cluster_id combination
        """
        if not fig_dir is None and not os.path.exists(fig_dir):
            os.makedirs(fig_dir, exist_ok=True)
        df = self.boundary_cell_df
        min_x = self.min_x
        max_x = self.max_x
        min_y = self.min_y
        max_y = self.max_y
        cluster_values = df['cluster_id'].unique()
        colors = plt.cm.tab20c(np.linspace(0, 1, len(cluster_values)))
        color_iterator = iter(colors)
        plt.axis('equal')
        for combination, group in df.groupby('cluster_id'):
            fig, ax = plt.subplots()
            x = group['x']
            y = group['y']

            # 绘制散点图，并使用不同的颜色
            color = next(color_iterator)
            ax.scatter(x, y, label=str(combination), color=color)

            ax.legend()
            plt.xticks(range(min_x, max_x, self.plot_unit))
            plt.yticks(range(min_y, max_y, self.plot_unit))
            if save_figs:
                plt.savefig(fig_dir + str(combination) + ".png", format="png")
            plt.clf()

        # 遍历每个颜色组合
        fig, ax = plt.subplots(figsize=figsize)
        plt.axis('equal')
        color_iterator = iter(colors)
        for color_combination, group in df.groupby('cluster_id'):
            x = group['x']
            y = group['y']
            color = next(color_iterator)
            ax.scatter(x, y, label=color_combination, color=color)
        plt.xticks(range(min_x, max_x, self.plot_unit))
        plt.yticks(range(min_y, max_y, self.plot_unit))

        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.tight_layout()
        if save_figs:
            plt.savefig(fig_dir + "combine.png", format="png")
        plt.clf()

        fig, ax = plt.subplots(figsize=figsize)
        plt.axis('equal')
        color_iterator = iter(colors)
        for color_combination, group in df.groupby('cluster_id'):
            x = group['x']
            y = group['y']
            color = next(color_iterator)
            ax.scatter(x, y, label=color_combination, color=color)
        plt.xticks(range(min_x, max_x, self.plot_unit))
        plt.yticks(range(min_y, max_y, self.plot_unit))

        if save_figs:
            plt.savefig(fig_dir + "combine_wo_label.png", format="png")
        # plt.show()
        plt.pause(1)
        plt.clf()

    def plot_boundary_cell_celltype(self, save_figs=False, fig_dir=None):
        """
        plot the distribution of boundary cells
        according to the cluster_id and cell types
        """
        if not fig_dir is None and not os.path.exists(fig_dir):
            os.makedirs(fig_dir, exist_ok=True)
        df = self.boundary_cell_df
        min_x = self.min_x
        max_x = self.max_x
        min_y = self.min_y
        max_y = self.max_y
        cell_types = list(df['cell_type'].unique())
        colors = plt.cm.tab20c(np.linspace(0, 1, len(cell_types)))
        for cluster_id, group in df.groupby('cluster_id'):
            fig, ax = plt.subplots()
            plt.axis('equal')
            for value in cell_types:
                cluster_data_sub = group[group["cell_type"] == value]
                ax.scatter(cluster_data_sub['x'], cluster_data_sub['y'],
                           color=colors[cell_types.index(value)])

            plt.xticks(range(min_x, max_x + 1, self.plot_unit))
            plt.yticks(range(min_y, max_y + 1, self.plot_unit))

            legend_elements = []

            for value in cell_types:
                # element = plt.Line2D([0], [0], marker='o', color='w', label=value,
                #                      markerfacecolor=ax.collections[0].cmap(
                #                          ax.collections[0].norm(value)),
                #                      markersize=10)
                element = plt.Line2D([0], [0], marker='o', color='w', label=value,
                                     markerfacecolor=colors[cell_types.index(value)],
                                     markersize=10)
                legend_elements.append(element)

            plt.legend(handles=legend_elements, loc='best', title='cell type')
            plt.title(
                "cell type distribution of boundary cells with color labels " + str(
                    cluster_id))
            if save_figs:
                plt.savefig(fig_dir + str(cluster_id) + ".png", format="png")
            # plt.show()
            plt.pause(1)
            plt.clf()
