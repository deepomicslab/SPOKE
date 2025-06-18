import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
import argparse
import anndata
import sys
sys.path.append("..")
from src.TEC import TEC
from scipy.spatial.distance import cdist
from scipy.sparse import csr_matrix
import pdb
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from multiprocessing import Pool, cpu_count
from concurrent.futures import ProcessPoolExecutor


def calculate_cell_similarity(X,cos_sim_cell, interfaces,idata, gamma_C=0.95, K=3):
    n_cells = X.shape[0]
    S_C = np.zeros((n_cells, n_cells))


    # Calculate interface similarity
    for i in range(n_cells):
        for j in range(i + 1, n_cells):
            # Get interfaces for each cell
            index_i = idata.obs.reset_index().index[(idata.obs['sender'] == str(i))| (idata.obs['receiver']==str(i))]
            index_j = idata.obs.reset_index().index[(idata.obs['sender'] == str(j))| (idata.obs['receiver']==str(j))]
            if len(index_i)*len(index_j)!=0:
                similarities = interfaces[np.ix_(index_i,index_j)].flatten().tolist()
                # Calculate interface similarity
                top_k_similarities = sorted(similarities, reverse=True)[:K]
                avg_similarity = np.mean(top_k_similarities) if top_k_similarities else 0
                S_C[i, j] = gamma_C * cos_sim_cell[i, j] + (1 - gamma_C) * avg_similarity
                S_C[j, i] = S_C[i, j]  # Symmetric matrix
            else:
                S_C[i, j] = 0
                S_C[j, i] = S_C[i, j]

    return S_C


def calculate_cell_similarity_v1(X, cos_sim_cell, interfaces, idata, gamma_C=0.95, K=3, num_threads=8):
    n_cells = X.shape[0]
    S_C = np.zeros((n_cells, n_cells))

    # 预计算每个细胞（字符串形式）在idata.obs中的行索引
    cell_to_indices = defaultdict(set)
    sender_col = idata.obs['sender'].values
    receiver_col = idata.obs['receiver'].values

    # 遍历所有观测行，收集每个细胞出现的行索引
    for idx in range(len(idata.obs)):
        s = sender_col[idx]
        r = receiver_col[idx]
        cell_to_indices[s].add(idx)
        cell_to_indices[r].add(idx)

    # 转换为列表并过滤有效细胞（0~n_cells-1的字符串形式）
    valid_cells = {str(i) for i in range(n_cells)}
    cell_indices_dict = {}
    for cell, indices_set in cell_to_indices.items():
        if cell in valid_cells:
            cell_indices_dict[cell] = np.array(list(indices_set), dtype=int)

    # 找出有接口的细胞索引
    cells_with_interface = []
    for i in range(n_cells):
        key = str(i)
        if key in cell_indices_dict and len(cell_indices_dict[key]) > 0:
            cells_with_interface.append(i)

    # 创建线程池处理细胞对
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = []
        m = len(cells_with_interface)

        # 提交所有细胞对任务
        for idx_i in range(m):
            i = cells_with_interface[idx_i]
            indices_i = cell_indices_dict[str(i)]

            for idx_j in range(idx_i + 1, m):
                j = cells_with_interface[idx_j]
                indices_j = cell_indices_dict[str(j)]

                # 提交任务到线程池
                future = executor.submit(
                    process_cell_pair,
                    interfaces, indices_i, indices_j, K, gamma_C, cos_sim_cell[i, j], i, j
                )
                futures.append(future)

        # 收集结果并更新矩阵
        for future in as_completed(futures):
            i, j, sim_val = future.result()
            S_C[i, j] = sim_val
            S_C[j, i] = sim_val

    return S_C


def process_cell_pair(interfaces, indices_i, indices_j, K, gamma_C, cos_sim_val, i, j):
    """处理单个细胞对的函数"""
    # 提取接口子矩阵
    sub_matrix = interfaces[np.ix_(indices_i, indices_j)]
    flattened = sub_matrix.ravel()
    k_val = min(K, len(flattened))

    if k_val > 0:
        # 使用partition高效获取top-k (部分排序)
        top_k = np.partition(flattened, -k_val)[-k_val:]
        avg_similarity = np.mean(top_k)
    else:
        avg_similarity = 0

    # 计算相似度值
    sim_val = gamma_C * cos_sim_val + (1 - gamma_C) * avg_similarity
    return i, j, sim_val

# def calculate_interface_similarity(interface_i, interface_j, S_C, idata, interface_cos_sim,gamma_I=0.05):
#     # Calculate pairwise cell similarities between interfaces i and j
#     cell_il = int(idata.obs.reset_index().loc[interface_i]['sender'])
#     cell_ir = int(idata.obs.reset_index().loc[interface_i]['receiver'])
#     cell_jl = int(idata.obs.reset_index().loc[interface_j]['sender'])
#     cell_jr = int(idata.obs.reset_index().loc[interface_j]['receiver'])
#     cells_sim = (S_C[cell_il, cell_jl] + S_C[cell_ir,cell_jr])/2
#
#     return gamma_I * cells_sim + (1 - gamma_I) * interface_cos_sim[interface_i,interface_j]

def calculate_interface_similarity(i, j, S_C, cell_il, cell_ir, interface_cos_sim, gamma_I=0.05):
    cells_sim = (S_C[cell_il[i], cell_il[j]] + S_C[cell_ir[i], cell_ir[j]]) / 2
    return gamma_I * cells_sim + (1 - gamma_I) * interface_cos_sim[i, j]

def iterative_graph_construction_v0(X, interfaces, infc_adata, threshold=0.005, max_iter=5, gamma_C=0.95,gamma_I=0.05,K=3):
    cos_sim_cell = cosine_similarity(X)
    cos_sim_infc = cosine_similarity(interfaces)
    S_C = calculate_cell_similarity(X=X,cos_sim_cell=cos_sim_cell, interfaces=cos_sim_infc,idata=infc_adata, gamma_C=gamma_C, K=K)

    if np.mean(np.abs(S_C - cos_sim_cell))<threshold * np.mean(cos_sim_cell):
        return S_C, cos_sim_infc

    S_I = cos_sim_infc

    # Iterative refinement of S_C and S_I
    for t in range(max_iter):
        print("Iter " + str(t) + "..")
        S_C_prev = S_C.copy()
        S_I_prev = S_I.copy()

        # Update S_I based on S_C
        for i in range(len(interfaces)):
            for j in range(i + 1, len(interfaces)):
                S_I[i, j] = calculate_interface_similarity(interface_i=i, interface_j=j, S_C=S_C, idata=infc_adata, interface_cos_sim=cos_sim_infc,gamma_I=gamma_I)
                S_I[j, i] = S_I[i, j]  # Symmetric matrix

        # Update S_C based on S_I
        S_C = calculate_cell_similarity(X=X,cos_sim_cell=cos_sim_cell, interfaces=S_I,idata=infc_adata, gamma_C=gamma_C, K=K)  # Update S_C based on the new S_I

        # Check for convergence
        if np.mean(np.abs(S_C - S_C_prev)) < threshold * np.mean(S_C) and \
                np.mean(np.abs(S_I - S_I_prev)) < threshold * np.mean(S_I):
            break

    return S_C, S_I

def compute_similarity_s(i, j, S_C, interface_cos_sim, gamma_I):
    similarity = calculate_interface_similarity(i=i, j=j, S_C=S_C, interface_cos_sim=interface_cos_sim, gamma_I=gamma_I)
    return i, j, similarity

def iterative_graph_construction(X, interfaces, infc_adata, threshold=0.005, max_iter=5, gamma_C=0.95,gamma_I=0.05,K=3):
    cos_sim_cell = cosine_similarity(X)
    cos_sim_infc = cosine_similarity(interfaces)
    S_C = calculate_cell_similarity(X=X,cos_sim_cell=cos_sim_cell, interfaces=cos_sim_infc,idata=infc_adata, gamma_C=gamma_C, K=K)

    if np.mean(np.abs(S_C - cos_sim_cell))<threshold * np.mean(cos_sim_cell):
        return S_C, cos_sim_infc

    S_I = cos_sim_infc
    n = cos_sim_infc.shape[0]

    cell_il = infc_adata.obs.reset_index().loc[:, 'sender'].astype(int).values
    cell_ir = infc_adata.obs.reset_index().loc[:, 'receiver'].astype(int).values

    # Iterative refinement of S_C and S_I
    for t in range(max_iter):
        print("Iter "+str(t)+"..")
        S_C_prev = S_C.copy()
        S_I_prev = S_I.copy()

        with ProcessPoolExecutor() as executor:
            futures = []
            for i in range(n):
                for j in range(i + 1, n):
                    futures.append(
                        executor.submit(calculate_interface_similarity, i, j, S_C, cell_il, cell_ir, cos_sim_infc,
                                        gamma_I))

            for future in futures:
                i, j, similarity = future.result()
                S_I[i, j] = similarity

        upper_triangle = np.triu(S_I, k=1)
        S_I += upper_triangle.T

        # Update S_C based on S_I
        S_C = calculate_cell_similarity(X=X,cos_sim_cell=cos_sim_cell, interfaces=S_I,idata=infc_adata, gamma_C=gamma_C, K=K)  # Update S_C based on the new S_I

        # Check for convergence
        if np.mean(np.abs(S_C - S_C_prev)) < threshold * np.mean(S_C) and \
                np.mean(np.abs(S_I - S_I_prev)) < threshold * np.mean(S_I):
            break

    return S_C, S_I

def compute_similarity_mat(X_embed, X_celltype, spatial_pos, n_neighbor_v, neighbor_k):
    bottom_up = TEC(affinity="gaussian_kernel",
                     kernel_gamma=None,
                     sparsification="knn_neighbors_from_X",
                     n_neighbors=n_neighbor_v,
                     objective="KL",
                     strategy="bottom_up", eta=3, eta1=1,
                     eta_mode="coefficient", eta2=1,
                     merge_layers=True,
                     plot_cluster_map=True, __verbose__=True)
    bottom_up.generate_graph(X_embed)
    M1 = bottom_up.knn_m
    dist_matrix = cdist(spatial_pos, spatial_pos, 'euclidean')
    np.fill_diagonal(dist_matrix, np.inf)
    indices = np.argsort(dist_matrix, axis=1)[:, :neighbor_k]

    if not X_celltype is None:
        arr = X_celltype.values
        M3 = (arr[:, None] == arr).astype(int)
    else:
        M3 = np.zeros_like(M1)
    return M1, M3, indices


def check_and_convert_to_numpy(sparse_matrix):
    if isinstance(sparse_matrix, csr_matrix):
        numpy_matrix = sparse_matrix.toarray()
        return numpy_matrix
    else:
        return sparse_matrix

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--cell_adata", type=str, required=True, help="Adata file of cells")
    parser.add_argument("--infc_adata", type=str, required=True, help="Adata file of interface")
    parser.add_argument("--thres",type=float, required=False,default=0.005, help="Threshold of iteration")
    parser.add_argument("--max_iter", type=int, required=False, default=5, help="Max times of iteration")
    parser.add_argument("--gamma_C",type=float, required=False,default=0.95, help="Proportion parameter for cell graph")
    parser.add_argument("--gamma_I",type=float, required=False,default=0.05, help="Proportion parameter for interface graph")
    parser.add_argument("--K", type=int, required=False, default=3,
                        help="Number of nearest neighbors")
    parser.add_argument("--output", type=str, required=False, default="..",
                        help="Output directory")

    args = parser.parse_args()
    infc_adata = anndata.read_h5ad(args.infc_adata)
    cell_adata = anndata.read_h5ad(args.cell_adata)
    K = args.K
    max_iter = args.max_iter
    gamma_C = args.gamma_C
    gamma_I = args.gamma_I
    thres = args.thres
    output = args.output

    X = check_and_convert_to_numpy(cell_adata.X)
    interfaces = check_and_convert_to_numpy(infc_adata.X)  # should be embeddings

    S_C, S_I = iterative_graph_construction(X=X, interfaces=interfaces, infc_adata=infc_adata, threshold=thres, max_iter=max_iter, gamma_C=gamma_C,gamma_I=gamma_I,K=K)
    print("Finish iteration.")
    print("Saving Cell Similarity Graph")
    np.savetxt(output+"//S_C.txt", S_C)
    print("Saving Interface Similarity Graph")
    np.savetxt(output+"//S_I.txt", S_I)
