import os
import scanpy as sc
import pandas as pd
import sys
sys.path.append("..")
from src.pystime.spatial_graph import cell_spatial_neighbors
from src.pystime.cell_interface import init_cell_interface_anndata, init_cell_pathway_anndata, LRI_norm


class STIME():
    
    def __init__(self, adata=None, in_dir=None, mtx_fn=None, lr_fn=None, pathway_fn=None,
                 species='Human', library_id='visium', quality='lowres',
                 search_method='neighbor', n_rings=1, n_neighs=6, distance_cutoff=500):


        self.init_adata(adata=adata, in_dir=in_dir, mtx_fn=mtx_fn, library_id=library_id, quality=quality)
        self.init_lrdb(lr_fn=lr_fn, species=species)
        self.init_pathway_df(pathway_fn=pathway_fn, species=species)
        
        cell_spatial_neighbors(adata, search_method=search_method, 
                       n_rings=n_rings, n_neighs=n_neighs,
                       distance_cutoff = distance_cutoff)

        self.idata = init_cell_interface_anndata(adata, self.lr_df)
        print('Successfully construct idata...')
        
        self.pdata = init_cell_pathway_anndata(adata, self.idata, self.pathway_df, species=species)
        print('Successfully construct pdata...')

        LRI_norm(self.idata, adata, self.pdata, method='MCMF-LRTf')
        print('Successfully norm LRI...')        
        
          
    def init_adata(self, adata=None, in_dir=None, mtx_fn=None, library_id='visium', quality='lowres'):
        if adata is None:
            adata = visium(in_dir, counts_file=mtx_fn, library_id=library_id, quality=quality)
            
        adata.obs['i'] = range(adata.shape[0])
        adata.var_names_make_unique()

        #vQC
        sc.pp.calculate_qc_metrics(adata, inplace=True, percent_top=None)
        sc.pp.filter_genes(adata, min_cells=10)
        sc.pp.filter_cells(adata, min_genes=10)
        adata.layers['raw'] = adata.X
        
        sc.pp.normalize_total(adata, inplace=True)
        sc.pp.log1p(adata)
        adata.layers['log1p'] = adata.X
        
        self.adata = adata
        print('Successfully init_adata...')
        
    def init_lrdb(self, lr_fn=None, species='Human'):

        if lr_fn is None:
            lr_df = pd.read_csv(os.path.join(os.path.dirname(__file__), 'data', 'stime-consensus_LRDB.csv'), index_col=0)
        else:
            lr_df = pd.read_csv(lr_fn)

        lr_df = lr_df[lr_df.ligand.isin(self.adata.var_names) & lr_df.receptor.isin(self.adata.var_names)]
        lr_df = lr_df[lr_df.species==species]

        lr_df = lr_df.reset_index()
        del lr_df['index']
        self.lr_df = lr_df
        print('Successfully init_lrdb...')
        
        
    def init_pathway_df(self, pathway_fn=None, species='Human'):

        if pathway_fn is None:
            pathway_df = pd.read_csv(os.path.join(os.path.dirname(__file__), 'data', 'SpaTalk_pathways.csv'), index_col=0)
        else:
            pathway_df = pd.read_csv(pathway_fn)

        self.pathway_df = pathway_df
        print('Successfully init_pathway...')
