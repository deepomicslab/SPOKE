from pulp import LpVariable, lpSum, LpProblem, LpStatus, LpMaximize, LpMinimize, value, PULP_CBC_CMD
from scipy.sparse import lil_matrix, coo_matrix, csr_matrix, triu, tril

import anndata as ad
import pandas as pd
import numpy as np
import time
from multiprocessing import Pool
from functools import partial
import networkx as nx

import pdb

#from stime.back import LRI_norm_by_MCMF_LRT_by_RTpath


def init_cell_interface_anndata(
    adata: ad.AnnData,
    lr_df: pd.DataFrame, 
    adata_layer: str='raw', 
    specie: str='Human', 
    **kwargs
): # Mouse
    
    data_gene2i = {gene:i for i, gene in enumerate(adata.var_names)}
    ligand_is = [data_gene2i[x] for x in lr_df.ligand]
    receptor_is = [data_gene2i[x] for x in lr_df.receptor]
    
    sender_is, receiver_is = adata.obsp["spatial_connectivities"].nonzero()
    sender_ligand_X = adata.layers[adata_layer][sender_is, ][:, ligand_is]
    receiver_receptor_X = adata.layers[adata_layer][receiver_is, ][:, receptor_is]

    idata = ad.AnnData(np.zeros_like(sender_ligand_X.todense()), dtype=sender_ligand_X.dtype)

    idata.obs_names = adata.obs_names[sender_is] + '->' + adata.obs_names[receiver_is]
    idata.obs['i'] = range(idata.obs.shape[0])
    
    idata.layers['sender_ligand_x'] = sender_ligand_X
    idata.layers['receiver_receptor_x'] = receiver_receptor_X
    
    idata.obs['sender'] = adata.obs_names[sender_is]
    idata.obs['receiver'] = adata.obs_names[receiver_is]
    idata.obs['sender_i'] = sender_is
    idata.obs['receiver_i'] = receiver_is
    idata.obs['sender_receiver_i'] = idata.obs['sender_i'].astype(str) + '->' + idata.obs['receiver_i'].astype(str)
    
    idata.obs['receiver_sender_i'] = idata.obs.receiver_i.astype(str) + '->' + idata.obs.sender_i.astype(str)
    '''
    original_index = idata.obs.index
    idata.obs.index = idata.obs.sender_receiver_i
    print(idata.obs.head())
    idata.obs['reverse_i'] = idata[idata.obs.receiver_sender_i, :].obs.i.to_list()
    idata.obs.index = original_index    
    '''
    
    idata.var_names = lr_df.ligand + '->' + lr_df.receptor
    idata.var['i'] = range(idata.var.shape[0])
    idata.var['ligand'] = lr_df.ligand.values
    idata.var['receptor'] = lr_df.receptor.values
    idata.var['ligand_i'] = ligand_is
    idata.var['receptor_i'] = receptor_is
    
    idata.obsm['spatial'] = (adata.obsm['spatial'][idata.obs['sender_i'].values] + adata.obsm['spatial'][idata.obs['receiver_i'].values])/2
    idata.obsm['spatial_distances'] = np.array(adata.obsp['spatial_distances'][idata.obs['sender_i'].values, \
                                idata.obs['receiver_i'].values]).reshape(-1,)
        
    if 'spatial' not in idata.uns.keys():
        idata.uns['spatial'] = adata.uns['spatial']
        
    '''
    cell_ligand_X = adata.X[:, ligand_is]
    cell_receptor_X = adata.X[:, receptor_is]
    lr_coexp = cell_ligand_X.todense() @ cell_receptor_X.todense().T + cell_receptor_X.todense() @ cell_ligand_X.todense().T
    adata.obsp['spatial_lr_coexps'] = csr_matrix(adata.obsp["spatial_connectivities"].multiply(lr_coexp))
    '''
    
    adata.obsp["idata_i"] = csr_matrix((idata.obs.i, (sender_is, receiver_is)), shape=(adata.shape[0], adata.shape[0]))
    
    adata.obsp["triu_spatial_connectivities"] = triu(adata.obsp["spatial_connectivities"], k=1, format="csr")
    adata.obsp["tril_spatial_connectivities"] = tril(adata.obsp["spatial_connectivities"], k=-1, format="csr")

    
    triu_df = pd.DataFrame({i: x for i, x in enumerate(adata.obsp['triu_spatial_connectivities'].nonzero())})
    idata.obs['triu'] = idata.obs['sender_receiver_i'].isin(triu_df[0].astype(str) + '->' + triu_df[1].astype(str))

    tril_df = pd.DataFrame({i: x for i, x in enumerate(adata.obsp['tril_spatial_connectivities'].nonzero())})
    idata.obs['tril'] = idata.obs['sender_receiver_i'].isin(tril_df[0].astype(str) + '->' + tril_df[1].astype(str))
    
    LRI_norm_by_LR_product(idata, inplace=True)
    record_LRI_to_adata(idata, adata, 'LR_product')

    return idata


def init_cell_pathway_anndata(adata, idata, pathway_df, adata_layer='raw', species='Human', depth_limit=None,
                              **kwargs): # Mouse
    
    pdata = init_cell_RTNet(adata, idata, pathway_df, species=species, adata_layer=adata_layer)
    
    if pdata is None:
        print('Pathway data is empty!')
        return None
    
    compute_cell_RTNet_product(adata, idata, pdata, species=species, depth_limit=depth_limit,
                                       adata_layer=adata_layer, idata_layer='LR_product')        
        
    LRI_norm_by_LRTNet_product(idata, adata)        
    record_LRI_to_adata(idata, adata, 'LRTNet_product')
    
    return pdata


def init_cell_RTNet(adata, idata, pathway_df, species='Human', adata_layer='raw'):   
    
    pathway_df = pathway_df[pathway_df.species==species]
    pathway_df = pathway_df[pathway_df.src.isin(adata.var_names) & pathway_df.dest.isin(adata.var_names)]
    pathway_df = pathway_df[~pathway_df.dest.isin(idata.var.receptor)]

    if pathway_df.empty:
        return None
    
    tf_genes = set(pathway_df[(pathway_df.src_tf == 'YES')].src.to_list() + pathway_df[(pathway_df.dest_tf == 'YES')].dest.to_list())
    
    pathway_df = pathway_df.groupby(['src', 'dest']).agg({
        'pathway': lambda x: ','.join(x),
        'source': lambda x: ','.join(x),
        'type': lambda x: ','.join(x),
    }).reset_index()
    data_gene2i = {gene:i for i, gene in enumerate(adata.var_names)}
    pathway_df['src_i'] = pathway_df.src.apply(lambda x: data_gene2i[x])
    pathway_df['dest_i'] = pathway_df.dest.apply(lambda x: data_gene2i[x])
    pathway_df['src_dest_i'] = pathway_df.src_i.astype(str) + '->' + pathway_df.dest_i.astype(str)
    pathway_df['i'] = range(pathway_df.shape[0])
    
    tf_is = np.array([data_gene2i[tf] for tf in tf_genes])
    non_tf_is = np.delete(range(adata.shape[1]), tf_is)
    adata.uns['TF_is'] = tf_is
    adata.uns['non_TF_is'] = non_tf_is    

    A = csr_matrix(([1]*pathway_df.shape[0], (pathway_df.src_i, pathway_df.dest_i)), shape=(adata.shape[1], adata.shape[1]))
    adata.varp['RTNet'] = A
    
    src_x = adata.layers[adata_layer][:, pathway_df.src_i]
    dest_x = adata.layers[adata_layer][:, pathway_df.dest_i]
    
    X = src_x.multiply(dest_x)
    
    pdata = ad.AnnData(X, dtype=X.dtype)
    pdata.layers['SD_product'] = X
    pdata.layers['src_x'] = src_x
    pdata.layers['dest_x'] = dest_x
    pdata.obs = adata.obs
    pdata.var = pathway_df

    return pdata


def compute_cell_RTNet_product(adata, idata, pdata, species='Human', depth_limit=None,
                   adata_layer='raw', idata_layer='LR_product',
                   ): 
    if depth_limit is None:
        if species == 'Human':
            depth_limit = 3
        else:
            depth_limit = 4
        
    forward_cell_RTNet_product_for_T(adata, idata, pdata, depth_limit=depth_limit,
                                     adata_layer=adata_layer, idata_layer=idata_layer)
        
    backward_cell_RTNet_product_for_R(adata, idata, pdata, depth_limit=depth_limit,
                                      adata_layer=adata_layer, idata_layer=idata_layer)
    
    remove_broken_RTNet_edge(pdata, adata, depth_limit, pdata_layer='RTNet_product_sumT')    
    
    pdata_X = csr_matrix(coo_matrix(pdata.shape))        
    pdata.var['depth'] = -1
    for depth in range(1, depth_limit+1):
        X = pdata.layers['RTNet_product_sumT_depth{}'.format(depth)] 
        pdata_X = pdata_X + X
        _, path_is = X.nonzero()
        pdata.var['depth'][pdata.var.i.isin(path_is)] = depth  ## imporve
    pdata.layers['RTNet_product'] = pdata_X    
    
    return pdata


def forward_cell_RTNet_product_for_T(adata, idata, pdata, depth_limit=3,
                   adata_layer='raw', idata_layer='LR_product'):    
    
    interface_is, lr_is = idata.layers[idata_layer].nonzero()

    receiver_is = idata[interface_is, lr_is].obs['receiver_i']
    receptor_is = idata[interface_is, lr_is].var['receptor_i']    
            
    RTNet_product = csr_matrix(coo_matrix(adata.shape))

    prev_signals = csr_matrix(coo_matrix(adata.shape))
    prev_signals[receiver_is, receptor_is] = adata.layers[adata_layer][receiver_is, receptor_is].tolist()[0]

    pdata_X = csr_matrix(coo_matrix(pdata.shape))
        
    depth = 1
    while depth <= depth_limit:

        signals = (prev_signals @ adata.varp['RTNet']).multiply(adata.layers[adata_layer])
        
        prev_cell_is, prev_tf_is = prev_signals[:, adata.uns['TF_is']].nonzero()
        prev_tf_is = adata.uns['TF_is'][prev_tf_is]

        if len(prev_cell_is) > 0 and len(prev_tf_is) > 0:
            tf_signals = csr_matrix(coo_matrix(adata.shape))
            tf_signals[prev_cell_is, prev_tf_is] = prev_signals[prev_cell_is, prev_tf_is].tolist()[0]

            target_signals = (tf_signals @ adata.varp['RTNet']).multiply(adata.layers[adata_layer])
            tmp_cell_is, tmp_target_is = target_signals[prev_cell_is, :][:, adata.uns['non_TF_is']].nonzero()
            current_cell_is = prev_cell_is[tmp_cell_is]
            current_target_is = adata.uns['non_TF_is'][tmp_target_is]

            if len(current_cell_is) > 0 and len(current_target_is) > 0:
                RTNet_product[current_cell_is, current_target_is] += target_signals[current_cell_is, current_target_is].tolist()[0]
                
        X = compute_LRTNet_edge(pdata, adata, prev_signals, signals, adata_layer=adata_layer)
        pdata.layers['RTNet_product_sumT_depth{}'.format(depth)] = X
        pdata_X = pdata_X + X        
            
        prev_signals = signals.copy()
        depth += 1
            
                    
    adata.layers['RTNet_product_sumT'] = RTNet_product

    return pdata_X    

    
def backward_cell_RTNet_product_for_R(adata, idata, pdata, depth_limit=3,
                   adata_layer='raw', idata_layer='LR_product'):    
    
    interface_is, lr_is = idata.layers[idata_layer].nonzero()
    receiver_is = idata[interface_is, lr_is].obs['receiver_i']   
    receptor_is = idata[interface_is, lr_is].var['receptor_i']
    receiver_receptor_from_idata = csr_matrix(coo_matrix(adata.shape))
    receiver_receptor_from_idata[receiver_is, receptor_is] = 1
    
    cell_is, target_is = adata.layers['RTNet_product_sumT'].nonzero()

    
    RTNet_product_sumR = csr_matrix(coo_matrix(adata.shape))

    prev_signals = csr_matrix(coo_matrix(adata.shape))
    prev_signals[cell_is, target_is] = adata.layers[adata_layer][cell_is, target_is].tolist()[0]
    
    pdata_X = csr_matrix(coo_matrix(pdata.shape))

    depth = 1
    while depth <= depth_limit:
        signals = (prev_signals @ adata.varp['RTNet'].T).multiply(adata.layers[adata_layer])
        if depth == 1: # the current gene must be tf
            current_cell_is, current_tf_is = signals[:, adata.uns['TF_is']].nonzero()
            current_tf_is = adata.uns['TF_is'][current_tf_is]
            if len(current_cell_is) > 0 and len(current_tf_is) > 0:
                tf_signals = csr_matrix(coo_matrix(adata.shape))
                tf_signals[current_cell_is, current_tf_is] = signals[current_cell_is, current_tf_is].tolist()[0]
                signals = tf_signals
            else:
                #print('-no tf-')
                break

        receptor_signals = signals.multiply(receiver_receptor_from_idata)
        current_cell_is, current_receptor_is = receptor_signals.nonzero()            

        if len(current_cell_is) > 0 and len(current_receptor_is) > 0:
            RTNet_product_sumR[current_cell_is, current_receptor_is] += signals[current_cell_is, current_receptor_is].tolist()[0]
            
        X = compute_LRTNet_edge(pdata, adata, signals, prev_signals, adata_layer=adata_layer)
        pdata.layers['RTNet_product_sumR_depth{}'.format(depth)] = X            
        pdata_X = pdata_X + X
            
        prev_signals = signals.copy()
        depth += 1
    
    adata.layers['RTNet_product_sumR'] = RTNet_product_sumR

    return pdata_X
    
    
def compute_LRTNet_edge(pdata, adata, prev_signals, signals, 
                       adata_layer='raw'):    
    X = csr_matrix(coo_matrix(pdata.shape))
    
    src_cell_is, src_gene_is = prev_signals.nonzero() 
    dest_cell_is, dest_gene_is = signals.nonzero()
    
    path_is = np.array(pdata.var[pdata.var.src_i.isin(src_gene_is) & pdata.var.dest_i.isin(dest_gene_is)].i.tolist())
    cell_is = np.array(list(set(src_cell_is) & set(dest_cell_is)))

    if cell_is.shape[0] == 0 or path_is.shape[0] == 0:
        return X
            
    path_df = pdata[:, path_is].var[['src_i', 'dest_i', 'i']]

    src_X = prev_signals[cell_is, :][:, path_df.src_i]
    dest_X = adata.layers[adata_layer][cell_is, :][:, path_df.dest_i]

    tmp_cell_is, tmp_path_is = src_X.multiply(dest_X).nonzero()
    cell_is = cell_is[tmp_cell_is]
    path_is = path_is[tmp_path_is]
    
    X[cell_is, path_is] = 1

    return X


def remove_broken_RTNet_edge(pdata, adata, depth_limit, pdata_layer='RTNet_product_sumT'):
    
    X = pdata.layers['{}_depth{}'.format(pdata_layer, depth_limit)]
    cell_is, path_is = X.nonzero()        
    path_df = pdata[:, path_is].var.copy()
    path_df['cell_i'] = cell_is    
    
    receiver_is, targets_is = adata.layers[pdata_layer].nonzero()
    path_is = np.array(pdata.var[pdata.var.dest_i.isin(targets_is)].i)
        
    path_df = path_df[path_df.cell_i.isin(receiver_is) & \
                      path_df.dest_i.isin(targets_is) & \
                      path_df.src_i.isin(adata.uns['TF_is'])] 
        
    path_df['cell_dest'] = path_df.cell_i.astype(str) + '-' + path_df.dest_i.astype(str)
    path_df = path_df[path_df.cell_dest.isin(['{}-{}'.format(x, y) for x, y in zip(receiver_is, targets_is)])]
        
    update_X = csr_matrix(coo_matrix(X.shape))
    update_X[path_df.cell_i, path_df.i] = 1
        
    pdata.layers['{}_depth{}'.format(pdata_layer, depth_limit)] = csr_matrix(update_X)        
        
    remove_broken_RTNet_edge_aux(pdata, adata, depth_limit, pdata_layer=pdata_layer)
    
    
def remove_broken_RTNet_edge_aux(pdata, adata, depth, pdata_layer='RTNet_product_sumT'):
    if depth == 1:
        return 
    
    X = pdata.layers['{}_depth{}'.format(pdata_layer, depth)]
    cell_is, path_is = X.nonzero()
    src_is = pdata[:, path_is].var.src_i
    
    adata_X = csr_matrix(coo_matrix(adata.shape))
    adata_X[cell_is, src_is] = 1
    
    prev_X = pdata.layers['{}_depth{}'.format(pdata_layer, depth-1)]
    prev_cell_is, prev_path_is = prev_X.nonzero()    
    prev_path_df = pdata[:, prev_path_is].var.copy()
    prev_path_df['cell_i'] = prev_cell_is
    prev_dest_is = prev_path_df.dest_i

    prev_adata_X = csr_matrix(coo_matrix(adata.shape))
    prev_adata_X[prev_cell_is, prev_dest_is] = 1
    
    broken_cell_is, broken_dest_is = ((prev_adata_X - adata_X) == 1).nonzero()
    
    prev_path_df['cell_dest'] = prev_path_df.cell_i.astype(str) + '-' + prev_path_df.dest_i.astype(str)
    broken_path_df = prev_path_df[prev_path_df.cell_dest.isin(['{}-{}'.format(x, y) \
                                  for x, y in zip(broken_cell_is, broken_dest_is)])]  
    broken_path_df = broken_path_df[~(broken_path_df.src_i.isin(adata.uns['TF_is']) & \
                                    broken_path_df.dest_i.isin(adata.uns['non_TF_is']))]

        
    prev_X[broken_path_df.cell_i, broken_path_df.i] = 0    
    pdata.layers['{}_depth{}'.format(pdata_layer, depth-1)] = prev_X
    
    remove_broken_RTNet_edge_aux(pdata, adata, depth-1, pdata_layer=pdata_layer)

    
###############LRI normalization#################

    
def LRI_norm(idata, adata=None, pdata=None, method='MCMF-LR', inplace=True, 
             lr_threshold=3,
             **kwargs):
    
    if method == 'MCMF-LRT':
        prob = LRI_norm_by_MCMF_LRT(idata, adata, pdata, inplace=inplace, method='MCMF-LRT', **kwargs)
        
    if method == 'MCMF-LRTf':
        prob = LRI_norm_by_MCMF_LRT(idata, adata, pdata, inplace=inplace, method='MCMF-LRTf', **kwargs)
                        
    if method == 'MCMF-LR':
        prob = LRI_norm_by_MCMF_LRT(idata, adata, pdata, inplace=inplace, method='MCMF-LR', **kwargs)
    if method == 'MCMF-LRf':
        prob = LRI_norm_by_MCMF_LRT(idata, adata, pdata, inplace=inplace, method='MCMF-LRf', **kwargs)
            
    
    if method == 'MF-LR':
        prob = LRI_norm_by_MCMF_LRT(idata, adata, pdata, inplace=inplace, method='MF-LR', **kwargs)
        
    if method == 'LR_thresholding':
        method = 'LR_threshold {}'.format(lr_threshold)
        prob = LRI_norm_by_LR_thresholding(idata, inplace=inplace, lr_threshold=lr_threshold, tag=method, **kwargs)
        
    if method == 'LR_mean':
        prob = LRI_norm_by_LR_mean(idata, inplace=inplace, **kwargs)
            
    if adata is not None:        
        record_LRI_to_adata(idata, adata, method)
    
    return prob

    
def record_LRI_to_adata(idata, adata, method):

    sender_df = pd.DataFrame(idata.layers[method].todense(), 
                             index=idata.obs['sender_i'], columns=['s_{}'.format(x) for x in idata.var_names])
    sender_df = sender_df.groupby('sender_i').agg(sum)
    sender_df['s_total->total'] = sender_df.sum(axis=1)

    receiver_df = pd.DataFrame(idata.layers[method].todense(),
                               index=idata.obs['receiver_i'], columns=['r_{}'.format(x) for x in idata.var_names])

    receiver_df = receiver_df.groupby('receiver_i').agg(sum)
    receiver_df['r_total->total'] = receiver_df.sum(axis=1)

    receiver_df.index.name = 'sender_i'
    sender_df = pd.concat([sender_df, receiver_df], axis=1)
    sender_cols = sender_df.columns
    sender_df = adata.obs.merge(sender_df, left_on='i', right_on='sender_i', how='left')
    sender_df.index = sender_df.barcode
    sender_df = sender_df[sender_cols]
    adata.obsm[method] = sender_df

    
def LRI_norm_by_LRTNet_product(idata, adata, adata_layer='raw', inplace=True):
    
    sender_ligand_X = idata.layers['sender_ligand_x']
    receiver_RTNet_product = adata.layers['RTNet_product_sumR'][idata.obs['receiver_i'], :][:, idata.var['receptor_i']]
    
    LRTNet_product = sender_ligand_X.multiply(receiver_RTNet_product)
    idata.layers['LRTNet_product'] = LRTNet_product
    
    if inplace:
        idata.X = idata.layers['LRTNet_product']
        idata.uns['layer'] = 'LRTNet_product'    
    
    
def LRI_norm_by_LR_mean(idata, inplace=True):
    X = idata.layers['sender_ligand_x'] + idata.layers['receiver_receptor_x']
    idata.layers['LR_mean'] = X
    if inplace:
        idata.X = X
        idata.uns['layer'] = 'LR_mean'    
    
def LRI_norm_by_LR_product(idata, inplace=True):
    X = idata.layers['sender_ligand_x'].multiply(idata.layers['receiver_receptor_x'])
    X = X.sqrt()
    idata.layers['LR_product'] = X
    if inplace:
        idata.X = X
        idata.uns['layer'] = 'LR_product'
        
def LRI_norm_by_LR_thresholding(idata, inplace=True, lr_threshold=3, tag=None, **kwargs):
    l = (idata.layers['sender_ligand_x'] > lr_threshold).astype(int)
    r = (idata.layers['receiver_receptor_x'] > lr_threshold).astype(int)
    X = l.multiply(r)
    idata.layers[tag] = X
    if inplace:
        idata.X = X
        idata.uns['layer'] = tag
        
       
        
###############ILP Solver#################
      
def split_list(
    xs: list, 
    num_batches: int
):
    batch_size = len(xs) // num_batches
    ys = [xs[i:i+batch_size] for i in range(0, len(xs), batch_size)]
    return ys


def LRI_norm_by_MCMF_LRT(
    idata: ad.AnnData, 
    adata: ad.AnnData, 
    pdata: bool=None, 
    method: str='MCMF-LRT', 
    inplace: bool=True,
    thread_n: int=8, 
    batch_n: int=8, 
    verbose: bool=True):

    if method.startswith('MF'):
        prob = LpProblem(method, LpMaximize)
    elif method.startswith('MCMF'):
        prob = LpProblem(method, LpMinimize)
    else:
        return
    
    if method.endswith('LR'):
        idata_layer = 'LR_product'
    elif method.endswith('LRT'):
        idata_layer = 'LRTNet_product'
    elif method.endswith('LRTf'):
        idata_layer = 'LRTNet_product'
    else:
        return
    
    adata_layer = 'raw'
    pdata_layer = 'RTNet_product'
    
    t0 = time.time()    
    interface_lr_is = list(zip(*idata.layers[idata_layer].nonzero()))        
    sum_e_LR, e_LR_dict, n_dict = add_LR_edge_capacity(interface_lr_is, idata=idata, method=method)


    if verbose:
        print('Interface-LR edge num:', len(e_LR_dict))
    idata.uns[method] = {'Interface-LR edge num': len(e_LR_dict)}    
    t1 = time.time()
    
    if verbose:
        print('Add LR edge time:', t1 - t0)    

    if method.endswith('LRT'):
        t0 = time.time()
        cell_path_is = list(zip(*pdata.layers[pdata_layer].nonzero()))
        sum_e_path, e_path_dict, e_n_dict = add_RTNet_edge_capacity(cell_path_is, pdata, method=method, pdata_layer=pdata_layer)   
        n_dict.update(e_n_dict)
        
        t1 = time.time()
        
        if verbose:
            print('Receiver-RTNet edge num:', len(e_path_dict))
            print('Add Receiver-RTNet edge time:', t1 - t0)
        idata.uns[method] = {'Receiver-RTNet edge num:': len(e_path_dict)}

    t0 = time.time()

    
    residuals = add_node_capacity(n_dict, adata=adata, prob=prob, adata_layer=adata_layer)
    t1 = time.time()
    if verbose:
        print('Add node residual constraints:', t1 - t0)        
                
    if method.endswith('LRT'):
        if method.startswith('MCMF'):
            obj = sum_e_LR + sum_e_path + residuals
        elif method.startswith('MF'):
            obj = sum_e_LR + sum_e_path
    else:
        if method.startswith('MCMF'):
            obj = sum_e_LR + residuals
        elif method.startswith('MF'):
            obj = sum_e_LR        
                
    prob += obj

    # is the CBC solver has L1 norm?
    t0 = time.time()
    status = prob.solve(PULP_CBC_CMD(threads=thread_n, msg=0))
    
    t1 = time.time()
    if verbose:
        print('LP solver time:', t1 - t0)    
        
    val_LR_X = lil_matrix(coo_matrix(idata.shape))
    min_LR_cost = 0
    for interface_i, lr_i in zip(*idata.layers[idata_layer].nonzero()):
        e_name = "e_i{}_lr{}".format(interface_i, lr_i)
        v = value(e_LR_dict[e_name])
        if v != 0: # and v is not None:  #???
            val_LR_X[interface_i, lr_i] = v
            min_LR_cost += idata.obsm['spatial_distances'][interface_i] * v
        
    if method.endswith('LRT'):
        val_path_X = lil_matrix(coo_matrix(pdata.shape))
        min_path_cost = 0
        for cell_i, path_i in zip(*pdata.layers[pdata_layer].nonzero()):
            e_name = "e_c{}_p{}".format(cell_i, path_i) 
            v = value(e_path_dict[e_name])
            if v != 0:
                val_path_X[cell_i, path_i] = v
                min_path_cost += v
        
    if verbose:
        print("Problem status: ", LpStatus[status])
        print("Total cost:", value(prob.objective))
        print("Min LR cost:", min_LR_cost)
        if method.endswith('LRT'):
            print("Min LRNet path cost", min_path_cost)        
        print("Residual", value(residuals))

        
    idata.uns[method] = {"LP_status": LpStatus[status]}
    idata.uns[method] = {"LP_cost": value(prob.objective)}
    idata.uns[method] = {"Min LR cost": min_LR_cost}
    if method.endswith('LRT'):
        idata.uns[method] = {"Min LRNet path cost": min_path_cost}    
    idata.uns[method] = {"Residual": value(residuals)}
    
    idata.layers[method] = val_LR_X
    if inplace:
        idata.X = idata.layers[method] 
        idata.uns['layer'] = method

    if method.endswith('LRT'):        
        pdata.layers[method] = val_path_X
        if inplace:
            pdata.X = pdata.layers[method] 
            pdata.uns['layer'] = method        

    return prob   
        
    
def LRI_norm_by_MCMF_LRT_back(idata, adata, pdata=None, method='MCMF-LRT', inplace=True,
                                thread_n=1, batch_n=1, 
                                verbose=True):

    if method.startswith('MF'):
        prob = LpProblem(method, LpMaximize)
    elif method.startswith('MCMF'):
        prob = LpProblem(method, LpMinimize)
    else:
        return
    
    if method.endswith('LR'):
        idata_layer = 'LR_product'
    elif method.endswith('LRT'):
        idata_layer = 'LRTNet_product'
    elif method.endswith('LRTf'):
        idata_layer = 'LRTNet_product'
    else:
        return
    
    adata_layer = 'raw'
    pdata_layer = 'RTNet_product'
    
    t0 = time.time()    
    e_LR_dict = {}
    n_dict = {}    
    sum_e_LR, e_LR_dict, n_dict = add_LR_edge_capacity(idata, e_LR_dict, n_dict, method, idata_layer=idata_layer)

    if verbose:
        print('Interface-LR edge num:', len(e_LR_dict))
    idata.uns[method] = {'Interface-LR edge num': len(e_LR_dict)}
        
    t1 = time.time()
    if verbose:
        print('Add LR edge time:', t1 - t0)
        
    if method.endswith('LRT'):
        t0 = time.time()
        e_path_dict = {}
        sum_e_path, e_path_dict, n_dict = add_RTNet_edge_capacity(pdata, e_path_dict, n_dict, method, pdata_layer=pdata_layer)   
        t1 = time.time()
        
        if verbose:
            print('Receiver-RTNet edge num:', len(e_path_dict))
            print('Add Receiver-RTNet edge time:', t1 - t0)
        idata.uns[method] = {'Receiver-RTNet edge num:': len(e_path_dict)}
            
    t0 = time.time()
    residuals = add_node_capacity(adata, prob, n_dict, adata_layer=adata_layer)
    t1 = time.time()
    if verbose:
        print('Add node residual constraints:', t1 - t0)        
                
    if method.endswith('LRT'):
        if method.startswith('MCMF'):
            obj = sum_e_LR + sum_e_path + residuals
        elif method.startswith('MF'):
            obj = sum_e_LR + sum_e_path
    else:
        if method.startswith('MCMF'):
            obj = sum_e_LR + residuals
        elif method.startswith('MF'):
            obj = sum_e_LR        
                
    prob += obj
        
    # is the CBC solver has L1 norm?
    t0 = time.time()
    status = prob.solve(PULP_CBC_CMD(threads=thread_n, msg=0))
    
    t1 = time.time()
    if verbose:
        print('LP solver time:', t1 - t0)    
        
    val_LR_X = lil_matrix(coo_matrix(idata.shape))
    min_LR_cost = 0
    for interface_i, lr_i in zip(*idata.layers[idata_layer].nonzero()):
        e_name = "e_i{}_lr{}".format(interface_i, lr_i)
        v = value(e_LR_dict[e_name])
        if v != 0:
            val_LR_X[interface_i, lr_i] = v
            min_LR_cost += idata.obsm['spatial_distances'][interface_i] * v
        
    if method.endswith('LRT'):
        val_path_X = lil_matrix(coo_matrix(pdata.shape))
        min_path_cost = 0
        for cell_i, path_i in zip(*pdata.layers[pdata_layer].nonzero()):
            e_name = "e_c{}_p{}".format(cell_i, path_i) 
            v = value(e_path_dict[e_name])
            if v != 0:
                val_path_X[cell_i, path_i] = v
                min_path_cost += v
        
    if verbose:
        print("Problem status: ", LpStatus[status])
        print("Total cost:", value(prob.objective))
        print("Min LR cost:", min_LR_cost)
        if method.endswith('LRT'):
            print("Min LRNet path cost", min_path_cost)        
        print("Residual", value(residuals))

    idata.uns[method] = {"LP_status": LpStatus[status]}
    idata.uns[method] = {"LP_cost": value(prob.objective)}
    idata.uns[method] = {"Min LR cost": min_LR_cost}
    if method.endswith('LRT'):
        idata.uns[method] = {"Min LRNet path cost": min_path_cost}    
    idata.uns[method] = {"Residual": value(residuals)}
    
    idata.layers[method] = val_LR_X
    if inplace:
        idata.X = idata.layers[method] 
        idata.uns['layer'] = method

    if method.endswith('LRT'):        
        pdata.layers[method] = val_path_X
        if inplace:
            pdata.X = pdata.layers[method] 
            pdata.uns['layer'] = method        

    return prob           
        
def add_LR_edge_capacity(
    interface_lr_is: list,
    idata: ad.AnnData=None,
    method: str='MCMF', 
):
    sum_e = 0
    e_dict, n_dict = {}, {}
    for interface_i, lr_i in interface_lr_is:
        
        sender_i = idata.obs['sender_i'][interface_i]
        receiver_i = idata.obs['receiver_i'][interface_i]
        
        ligand_i = idata.var['ligand_i'][lr_i]
        receptor_i = idata.var['receptor_i'][lr_i]

        ligand_x = idata.layers['sender_ligand_x'][interface_i, lr_i] 
        receptor_x = idata.layers['receiver_receptor_x'][interface_i, lr_i] 
        min_x = min(ligand_x, receptor_x)
    
        e_name = "e_i{}_lr{}".format(interface_i, lr_i)
        e = LpVariable(e_name, 0, min_x)
        e_dict[e_name] = e
        
        if method.startswith('MCMF'):
            cost = idata.obsm['spatial_distances'][interface_i]
        else:
            cost = 1
            
        sum_e += e * cost

        if (sender_i, ligand_i) not in n_dict.keys():
            n_dict[(sender_i, ligand_i)] = {'src': e}
        elif 'src' not in n_dict[(sender_i, ligand_i)]:
            n_dict[(sender_i, ligand_i)]['src'] = e
        else:
            n_dict[(sender_i, ligand_i)]['src'] += e

        if (receiver_i, receptor_i) not in n_dict.keys():
            n_dict[(receiver_i, receptor_i)] = {'dest': e}
        elif 'dest' not in n_dict[(receiver_i, receptor_i)]:
            n_dict[(receiver_i, receptor_i)] = {'dest': e}
        else:
            n_dict[(receiver_i, receptor_i)]['dest'] += e

    return sum_e, e_dict, n_dict


def add_RTNet_edge_capacity(
    cell_path_is: list,
    pdata: ad.AnnData=None, 
    method: str='MCMF', 
    pdata_layer: str='RTpath_product'
):

    sum_e = None
    e_dict, n_dict = {}, {}
    for cell_i, path_i in cell_path_is:

        src_i = pdata.var['src_i'][path_i]
        dest_i = pdata.var['dest_i'][path_i]

        src_x = pdata.layers['src_x'][cell_i, path_i] 
        dest_x = pdata.layers['dest_x'][cell_i, path_i] 
        min_x = min(src_x, dest_x)
        
        e_name = "e_c{}_p{}".format(cell_i, path_i)
        e = LpVariable(e_name, 0, min_x)
        e_dict[e_name] = e
                    
        if sum_e is None:
            sum_e = e 
        else:
            sum_e += e

        if (cell_i, src_i) not in n_dict.keys():
            n_dict[(cell_i, src_i)] = {'src': e}
        else:
            if 'src' not in n_dict[(cell_i, src_i)].keys():
                n_dict[(cell_i, src_i)] = {'src': e}
            else:
                n_dict[(cell_i, src_i)]['src'] += e

        if (cell_i, dest_i) not in n_dict.keys():
            n_dict[(cell_i, dest_i)] = {'dest': e}
        else:
            if 'dest' not in n_dict[(cell_i, dest_i)].keys():
                n_dict[(cell_i, dest_i)] = {'dest': e}
            else:
                n_dict[(cell_i, dest_i)]['dest'] += e

    return sum_e, e_dict, n_dict


def add_node_capacity(
    n_dict: dict,
    adata: ad.AnnData=None, 
    prob: LpProblem=None, 
    adata_layer:str='raw'
):

    residuals = 0
    for key, v in n_dict.items():
        
        cell_i, gene_i = key
        signals = adata.layers[adata_layer][cell_i, gene_i]        
        for direction, sum_e in v.items():
            prob += sum_e <= signals
            residuals += signals - sum_e     
            
        if 'src' in v.keys() and 'dest' in v.keys():
            prob += v['src'] == v['dest']
            
    return residuals