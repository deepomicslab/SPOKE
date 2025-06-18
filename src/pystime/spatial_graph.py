import squidpy as sq
from scipy.spatial import distance_matrix
import numpy as np
import seaborn as sns
from scipy import sparse

def cell_spatial_neighbors(adata, search_method='ring', n_rings=1, n_neighs=6, distance_cutoff=150):
    
    if search_method == 'ring':
        sq.gr.spatial_neighbors(adata, n_rings=n_rings, coord_type="grid", n_neighs=n_neighs)
        adata.obsp['spatial_connectivities'].setdiag(1)
        adata.obsp['spatial_distances'].setdiag(0.001)  #?????
        
    if search_method == 'neighbor':
        sq.gr.spatial_neighbors(adata, n_neighs=n_neighs, coord_type="generic")
        adata.obsp['spatial_connectivities'].setdiag(1)
        adata.obsp['spatial_distances'].setdiag(0.001)  #?????
        
    if search_method == 'delaunay':
        sq.gr.spatial_neighbors(adata, delaunay=True, coord_type="generic")
        adata.obsp['spatial_connectivities'].setdiag(1)
        adata.obsp['spatial_distances'].setdiag(0.001)  #?????

    if search_method == 'distance':
        get_distance_matrix(adata, distance_cutoff)
        
    
def get_distance_matrix(adata, distance_cutoff):
    
    M = distance_matrix(adata.obsm["spatial"], adata.obsm["spatial"])**2
    
    distance_cutoff = distance_cutoff ** 2
    M_row, M_col = np.where(M <= distance_cutoff)
    M = sparse.csr_matrix((M[M_row,M_col], (M_row,M_col)), shape=M.shape)
    adata.obsp['spatial_distances'] = M/distance_cutoff
    adata.obsp['spatial_connectivities'] = sparse.csr_matrix(([1]*len(M_row), (M_row, M_col)), shape=M.shape)
