# -*- coding: utf-8 -*-
"""
    src.HE
    ~~~~~~~~~~~

    @Copyright: (c) 2021-09 by Lingxi Chen (chanlingxi@gmail.com).
    @License: LICENSE_NAME, see LICENSE for more details.
"""
import numbers
from sklearn.preprocessing import StandardScaler
from scipy.sparse import lil_matrix
import numpy as np
from math import log2
from scipy.cluster import hierarchy
from sklearn.cluster import AgglomerativeClustering
from sklearn.neighbors import kneighbors_graph
from sklearn.metrics.pairwise import pairwise_kernels
from scipy import sparse
import heapq
import itertools
import networkx as nx
from itertools import chain
import pandas as pd
from queue import Queue
# import kmeans1d
import scipy.linalg as la
from scipy.sparse import csgraph
from scipy.sparse import csr_matrix
import time
import copy
import concurrent.futures
import threading
from itertools import compress
import pdb
import random
from src import graph_metric    # commented for local debug********************
# import graph_metric            # for local debug******************************
import sys
sys.setrecursionlimit(10000)


class Node:

    def __init__(self, graph_stats, node_id, children, vs, objective_paras, node_s=None, node_v=None, parent=None,
                 is_individual=True,
                 is_leaf=True, __verbose__=False):
        self.id = node_id
        self.parent = parent
        self.children = children
        self.vs = vs
        self.dist = 1
        self.height = 0
        self.is_individual = is_individual
        self.is_leaf = is_leaf
        self.split_se = np.nan

        self.empty = False
        self.parent_tmp = self.id  # initialized as the node itself

        self.se = 0.
        self.d_square = 0.
        self.dp_i = None
        M, m, d, log_d, d_log_d, vG, log_vG, sparse_m, d_square, M_rows = graph_stats
        if objective_paras.objective == 'M':   # extra precompted data
            self.d_square = np.sum(d_square[self.vs])
        if node_s!=None and node_v!=None:
            self.V = node_v
            self.s = node_s
        else:
            self.V = graph_metric.get_v(M_rows, sparse_m, self.vs)
            self.s = graph_metric.get_s(M, sparse_m, self.vs)
        self.log_V = log2(self.V)
        self.g = self.V - self.s
        self.d_log_d = np.sum(d_log_d[self.vs])
        self.pV = None
        self.__verbose__ = __verbose__
        if parent:
            self.pV = parent.V
            self.se = graph_metric.get_node_score_v2(vG, self.g, self.V, pV=parent.V, objective_paras=objective_paras)
        else:    # root
            if objective_paras.objective == 'SE' or objective_paras.objective == 'KL':
                self.se = 0
            elif objective_paras.objective == 'M':
                if objective_paras.eta_mode == 'coefficient':
                    self.se = 1 - objective_paras.eta
                elif objective_paras.eta_mode == "exponent":
                    self.se = 0
                else:
                    print("ERROR! Wrong mode of eta input.")

    def _remove(self, node_rmv, graph_stats, objective_paras):
        """
        remove vertex (or child node) from the node, update attributes
        """
        if self.vs == node_rmv.vs:
            self.empty = True
            self.vs = []
            self.d_square = 0
            self.V = 0
            self.s = 0
            self.g = 0
            self.d_log_d = 0
            self.se = 0
        else:
            self.vs = list(set(self.vs).difference(set(node_rmv.vs)))

            M, m, d, log_d, d_log_d, vG, log_vG, sparse_m, d_square, M_rows = graph_stats
            if objective_paras.objective == 'M':  # extra precompted data
                self.d_square = np.sum(d_square[self.vs])
            self.V = graph_metric.get_v(M_rows, sparse_m, self.vs)
            self.s = graph_metric.get_s(M, sparse_m, self.vs)
            self.g = self.V - self.s
            self.d_log_d = np.sum(d_log_d[self.vs])
            if self.pV:
                self.se = graph_metric.get_node_score_v2(vG, self.g,
                                                      self.V, self.pV, objective_paras=objective_paras)
            else:  # root
                if objective_paras.objective == 'SE' or objective_paras.objective == 'KL':
                    self.se = 0
                elif objective_paras.objective == 'M':
                    if objective_paras.eta_mode == 'coefficient':
                        self.se = 1 - objective_paras.eta
                    elif objective_paras.eta_mode == "exponent":
                        self.se = 0
                    else:
                        print("ERROR! Wrong mode of eta input.")

    def _insert(self, node_inst, graph_stats, objective_paras):
        """
        insert
        """
        self.vs = self.vs + node_inst.vs

        M, m, d, log_d, d_log_d, vG, log_vG, sparse_m, d_square, M_rows = graph_stats
        if objective_paras.objective == 'M':  # extra precompted data
            self.d_square = np.sum(d_square[self.vs])
        self.V = graph_metric.get_v(M_rows, sparse_m, self.vs)
        self.s = graph_metric.get_s(M, sparse_m, self.vs)
        self.g = self.V - self.s
        self.d_log_d = np.sum(d_log_d[self.vs])
        if self.pV:
            self.se = graph_metric.get_node_score_v2(vG, self.g,
                                                  self.V, self.pV, objective_paras=objective_paras)
        else:  # root
            if objective_paras.objective == 'SE' or objective_paras.objective == 'KL':
                self.se = 0
            elif objective_paras.objective == 'M':
                if objective_paras.eta_mode == 'coefficient':
                    self.se = 1 - objective_paras.eta
                elif objective_paras.eta_mode == "exponent":
                    self.se = 0
                else:
                    print("ERROR! Wrong mode of eta input.")

    def merge(self, node_id, node1, node2, graph_stats, objective_paras,is_leaf=False):
        """
        merge node1 and node2, the children of self
        ---
        Returns
        the parent node of node1 and node2
        """
        if is_leaf:  # children are vertices
            children = node1.children + node2.children
        else:
            children = [node1, node2]
        vs = node1.vs + node2.vs
        node = Node(graph_stats, node_id, children, vs, parent=self, objective_paras=objective_paras)
        if not is_leaf:
            node.dist = max(node1.dist, node2.dist) + 1
        node.is_leaf = is_leaf
        node1.parent = node
        node2.parent = node
        self.children.append(node)  # self is root
        return node

    def __repr__(self):
        return 'id:{}'.format(self.id)


class pySETree():

    def __init__(self, aff_m, knn_m, objective_paras, min_k=2, max_k=10, auto_k=False,
                 strategy='top_down',
                 split_se_cutoff=0.05, merge_layers=False, plot_cluster_map=False, random_seed=None, __verbose__=False
                 ):
        self.strategy = strategy
        self.objective_paras = objective_paras
        self.min_k = min_k
        self.max_k = max_k
        self.auto_k = auto_k
        if self.objective_paras.objective == "M" or self.objective_paras.objective == 'KL':
            self.split_se_cutoff = 0   # for debug
        else:
            self.split_se_cutoff = split_se_cutoff  # to control the size of clubs (higher cutoff, larger clubs)

        self.vertex_num = aff_m.shape[0]
        if self.max_k > self.vertex_num:
            self.max_k = self.vertex_num - 1

        self.ks = range(self.min_k, self.max_k+1)

        if strategy == 'top_down':
            self.node_id = 2*self.vertex_num - 3  # excluding the root
        else:
            self.node_id = -2
        self.node_list = {}

        self.aff_m = aff_m
        self.knn_m = knn_m
        self.G = nx.from_numpy_array(knn_m)  # generate an undirected graph of the adjacent matrix

        self.knn_graph_stats = self.graph_stats_init(knn_m)
        self.aff_graph_stats = self.graph_stats_init(aff_m)

        self.merge_layers = merge_layers # the result of clubs is obtained from multi layers, optional: True/False
        self.__MIN = 0.0000001 # for stop merging

        self.__verbose__ = __verbose__  # for debug
        self.__long_message__ = False  # set to print Z matrix and details of merging
        self.plot_cluster_map = plot_cluster_map  # whether plot cluster map
        self.random_seed = random_seed


    # def graph_stats_init(self, sym_m):
    #     M = csr_matrix(sym_m)
    #     # 清除对角线元素并保持稀疏格式
    #     M.setdiag(0)
    #     M.eliminate_zeros()
    #
    #     # 计算度数（已排除对角线）
    #     d = M.sum(axis=1).A1  # 转换为1D数组
    #
    #     d_square = np.zeros(d.shape)
    #
    #     # 处理度数为零的节点
    #     if np.any(d == 0):
    #         # 获取度数为零的节点索引
    #         zero_deg_nodes = np.where(d == 0)[0]
    #
    #         # 转换为LIL格式便于高效修改
    #         M = M.tolil()
    #         # 添加自环边避免除零
    #         M[zero_deg_nodes, zero_deg_nodes] = 1e-3
    #         # 转回CSR格式优化计算
    #         M = M.tocsr()
    #         M.eliminate_zeros()
    #
    #         # 重新计算度数
    #         d = M.sum(axis=1).A1
    #
    #     # 计算熵相关量
    #     log_d = np.log2(d)
    #     d_log_d = d * log_d  # 逐元素相乘
    #
    #     # 计算图统计量
    #     col_M = M.sum(axis=1).A1
    #     vG = col_M.sum()  # 总度数（无向图为2m）
    #     m = vG / 2  # 边数
    #     log_vG = np.log2(vG)
    #
    #     # 根据目标函数计算d_square
    #     if self.objective_paras.objective == 'M':
    #         d_square = d ** 2
    #
    #     # 返回计算结果
    #     graph_stats = M, m, d, log_d, d_log_d, vG, log_vG, sym_m, d_square
    #     return graph_stats

    def graph_stats_init(self, sym_m):
        M = sym_m
        np.fill_diagonal(M, 0)
        d = np.sum(M, 1) - M.diagonal()
        d_square = np.zeros(d.shape)
        if np.any(d == 0):
            M += 1e-3
            np.fill_diagonal(M, 0)
            d = np.sum(M, 1) - M.diagonal()
        log_d = np.log2(d)
        d_log_d = np.multiply(d, log_d)

        # test for optimize (the space complexity)
        # sparce_m = sparse.csr_matrix(sym_m)
        # m = sparce_m.sum() / 2
        # vG = sparce_m.sum()
        # log_vG = log2(vG)
        col_M = np.sum(M, 1)
        vG = np.sum(col_M)
        m = vG/2
        log_vG = log2(vG)
        sparce_m = sym_m

        if self.objective_paras.objective == 'M':
            d_square = d ** 2  # square

        M_rows = M.sum(axis=1)
        graph_stats = M, m, d, log_d, d_log_d, vG, log_vG, sparce_m, d_square, M_rows
        return graph_stats

    def update_node_id(self, increment=True):
        if increment:
            self.node_id += 1
        else:
            self.node_id -= 1
        return self.node_id

    def get_current_node_id(self):
        return self.node_id

    def build_tree_v3(self):
        print("Building tree..")
        M, m, d, log_d, d_log_d, vG, log_vG, sparce_m, d_square, M_rows = self.knn_graph_stats
        # root: max id
        root = Node(graph_stats=self.knn_graph_stats, node_id=self.update_node_id(), children=[], vs=list(range(self.vertex_num)), is_leaf=False, objective_paras=self.objective_paras, node_v=vG, node_s=vG)
        self.node_list[root.id] = root

        if self.strategy == 'bottom_up':
            Z = self.bottom_up_v3(root)
            root = self.node_list[self.node_id]

        self.root = root
        # if self.plot_cluster_map:
        #     self.Z_ = Z[:, :4]
        #     return Z
        # else:
        #     return None
        self.Z_ = Z[:, :4]
        return Z

    def bottom_up_v3(self, root):
        # layer_nodes = []  # save the id of nodes in the current layer being computed
        # Initialize the tree with each vertex in its own community
        for i in range(self.vertex_num):
            if self.__verbose__ and self.__long_message__:
                print("id: ",i)
            node = Node(graph_stats=self.knn_graph_stats, node_id=self.update_node_id(), children=[], vs=[i], parent=root, objective_paras=self.objective_paras)
            self.node_list[node.id] = node
            root.children.append(node)
            root.vs.append(i)
            # layer_nodes.append(node.id)
        Z = self.linkage_v3(root, random_state=self.random_seed)
        # self.linkage_v2(root)
        return Z

    def _get_dividing_delta(self, node, children):
        if len(children) < 2:
            return np.nan
        M, m, d, log_d, d_log_d, vG, log_vG, sparse_m, d_square, M_rows = self.knn_graph_stats
        left, right = children
        objective, eta_mode, eta, eta1, eta2, hype_a = self.objective_paras._get_vars_()
        if objective == 'M':
            if eta_mode == "coefficient":
                delta = children[0].s/vG - eta * np.power(children[0].V/vG, 2) + children[1].s/vG - eta * np.power(children[1].V/vG, 2)
            else:
                delta = children[0].s / vG - np.power(children[0].V / vG, 2*eta) + children[1].s / vG - np.power(children[1].V / vG, 2*eta)
        elif objective == 'KL':
            left.se = graph_metric.get_node_score_v2(vG=vG, g=left.g, V=left.V, pV=node.V, objective_paras=self.objective_paras)
            right.se = graph_metric.get_node_score_v2(vG=vG, g=right.g, V=right.V, pV=node.V, objective_paras=self.objective_paras)
            if node.V - node.g <= 0 or node.pV is None:
                delta = - left.se - right.se
            else:
                delta = graph_metric.get_node_score_v2(vG=vG, g=node.g, V=node.V, pV=node.pV, objective_paras=self.objective_paras) - left.se - right.se
            delta = - delta
        elif objective == 'SE':
            if self.strategy == 'bottom_up':
                left.se = left.g/vG*log2(node.V/left.V)
                right.se = right.g/vG*log2(node.V/right.V)
            # compute delta of se
            n_V_log_V = node.V * log2(node.V)
            left_V_log_V = left.V * log2(left.V)
            right_V_log_V = right.V * log2(right.V)
            delta = (n_V_log_V - node.d_log_d)/vG \
                - (left.se + (left_V_log_V - left.d_log_d)/vG) \
                - (right.se + (right_V_log_V - right.d_log_d)/vG)
        else:
            print("ERROR! Input WRONG objective.")
        return delta


    def get_max_delta_from_heap(self, heap, row_ids, col_ids):
        """
        Return the available node pair (max_n1,max_n2) with the maximum delta
        """
        while heap:
            max_delta, max_n1, max_n2 = heapq.heappop(heap)
            if max_n1 not in row_ids or max_n2 not in col_ids:
                continue
            max_n1 = self.node_list[max_n1]
            max_n2 = self.node_list[max_n2]
            max_delta = -max_delta
            return max_n1, max_n2, max_delta
        else:
            return None, None, None

    def check_random_state(self, seed):
        """Turn seed into a np.random.RandomState instance.

        Parameters
        ----------
        seed : None | int | instance of RandomState
            If seed is None, return the RandomState singleton used by np.random.
            If seed is an int, return a new RandomState instance seeded with seed.
            If seed is already a RandomState instance, return it.
            Otherwise raise ValueError.

        """
        if seed is None or seed is np.random:
            return np.random.mtrand._rand
        if isinstance(seed, (numbers.Integral, np.integer)):
            return np.random.RandomState(seed)
        if isinstance(seed, np.random.RandomState):
            return seed
        raise ValueError("%r cannot be used to seed a numpy.random.RandomState"
                         " instance" % seed)

    def _compute_layer(self, layer_nodes, node_list):
        """
        Compute the sum of entropy of nodes in current layer
        """
        # score_value = 0
        # for node_id in layer_nodes.keys():
        #     node_score = node_list[node_id].se
        #     score_value = score_value + node_score
        score_value = sum([node_list[node_id].se for node_id in layer_nodes.keys()])
        return score_value

    def _update_layers(self, layers, layer_nodes):
        """
        Update information in 'layers' before the terminating the loop call of one-level-merge,
        for further use in binary_combination
        """
        # save the result of one-level partition
        # layer_nodes_new, node_list, ass_dict = self._renumber(layer_nodes_new,
        #                                                       node_list, root)
        # layers.append(layer_nodes_new)
        # Claim the termination of merging
        # Process the result
        ### updated: add additional layer to layers if there is
        ### node not containing a single node in the highest layer
        highest_layer = layers[-1]
        # mod_flag = False
        # for node_id in highest_layer.keys():
        #     if len(highest_layer[node_id]) == 1:
        #         continue
        #     else:
        #         mod_flag = True
        #         break
        mod_flag = any(
            len(highest_layer[node_id]) > 1 for node_id in highest_layer.keys())
        if mod_flag:
            # layer_nodes_new, node_list, ass_dict = self._renumber(layer_nodes_new,
            #                                                   node_list, root)
            # cur_graph = self._induce_graph(cur_graph, layer_nodes_new, ass_dict)
            layer_nodes = {n: [n] for n in layer_nodes.keys()}
            layers.append(layer_nodes)
        ### update over
        # layer_nodes_new, node_list, ass_dict = self._renumber(layer_nodes_new,
        #                                                       node_list, root)
        # cur_graph = self._induce_graph(cur_graph, layer_nodes_new, ass_dict)
        return layers

    def linkage_v3_v0(self, root, random_state=None):
        """
        Rewrite the procedure of merging according to community_louvain of package community
        generate random seed
        ---
        root: the root node
        (deprecated) layer_nodes: the id of nodes in one layer
        random_state: set as default
        """
        # M, m, d, log_d, d_log_d, vG, log_vG, sparse_m, d_square = self.knn_graph_stats

        N = self.vertex_num
        Z = np.zeros((N - 1, 5))
        # leaves = {n: 1 for n in range(N)}
        # individuals = {n: 1 for n in range(N)}
        layer_nodes = {n: [n] for n in range(N)}
        # cutoff value
        __PASS_MAX = -1

        # heap = []
        # heapq.heapify(heap)

        random_state = self.check_random_state(random_state)

        # special case, when there is no link
        # the best partition is everyone in its community
        if not self.G.edges(data=True):
            # save the information in Z
            for z_i in range(self.vertex_num-1):
                Z[z_i] = [z_i,z_i+1,1,1,root.id]
            if self.__verbose__:
                print('-----', self.objective_paras.objective, 'merge phase 0',
                      'start, individuals: ', self.vertex_num, ', leaves ', self.vertex_num)
                print('-----Merge ends.-----')
            return Z

        # compute the entropy of initial tree (excluding the root node and vertices)
        new_score = self._compute_layer(layer_nodes, self.node_list)
        cur_score = new_score
        layers = [layer_nodes]
        # save the current hierarchical relationship
        node_list = self.node_list
        cur_graph = self.G.copy()
        times = 0
        if self.__verbose__:
            print("Start merging..")
            time_s = time.time()
        while True:
            # one-level merging
            layer_nodes_new, node_list = self._one_level_merge(layer_nodes=layer_nodes, node_list_=node_list, cur_graph=cur_graph, root=root, randomize_=True)
            # compute the entropy the current layer
            new_score = self._compute_layer(layer_nodes_new, node_list)
            if new_score - cur_score < self.__MIN:
                if self.__verbose__ and self.__long_message__:
                    print(layer_nodes_new)
                    print("new score: ", new_score)
                layers = self._update_layers(layers=layers, layer_nodes=layer_nodes)
                if self.__verbose__:
                    print("One layer iteration terminates.")
                break
            # It means that we should renumber, induce the graph, and update 'layers'
            # process the result
            layer_nodes_new, node_list, ass_dict = self._renumber(layer_nodes_new,
                                                                      node_list, root)
            cur_graph = self._induce_graph(cur_graph, layer_nodes_new, ass_dict)
            # if times == 0:
            #     self.G = copy.deepcopy(cur_graph)
            #     Z, z_i, parent_nodes_id = self._fill_Z_v1(layer_nodes_new, self.node_list)
            #     # remapping the nodes in self.G
            #     remapping = self._remap(parent_nodes_id)
            if self.__verbose__ and self.__long_message__:
                print(layer_nodes_new)
                print("new score: ", new_score)
            # save the result of one-level partition
            layers.append(layer_nodes_new)
            if not self.merge_layers:
                # Considering merge_flag = False, the iteration should be conducted only once.
                # We need to suppose it as merging terminates in the second time instead of the first one.
                layer_nodes = {n: [n] for n in layer_nodes_new.keys()}
                layers.append(layer_nodes)
                break
            # renumber and induce the graph
            layer_nodes = {n: [n] for n in layer_nodes_new.keys()}
            cur_score = new_score
            times += 1
        if self.__verbose__:
            time_e = time.time()
            print("Merging used: ", str(time_e-time_s), "s")
        self.G = cur_graph
        Z, z_i, i = self._binary_combine_v2(root=root, layers=layers)
        if self.__verbose__ and self.__long_message__:
            print(Z)
        # get clubs
        if len(layers) == 1:
            layer_tmp = layers[0]
        else:
            if self.merge_layers is True:
                # the process of binary combination needs one-layer tree
                layer_tmp = self._merge_into_one_layer(layers)
            else:
                # merge the multi-layer result into one-layer
                layer_tmp = layers[1]
        self.clubs = self._get_clubs_v2(layer_tmp)
        if self.auto_k:
            # set max_k as the number of clubs
            self.max_k = len(set(self.clubs))
            self.ks = range(self.min_k, self.max_k + 1)
        return Z

    def linkage_v3(self, root, random_state=None):
        """
        Rewrite the procedure of merging according to community_louvain of package community
        generate random seed
        ---
        root: the root node
        (deprecated) layer_nodes: the id of nodes in one layer
        random_state: set as default
        """
        # M, m, d, log_d, d_log_d, vG, log_vG, sparse_m, d_square = self.knn_graph_stats

        N = self.vertex_num
        Z = np.zeros((N - 1, 5))
        # leaves = {n: 1 for n in range(N)}
        # individuals = {n: 1 for n in range(N)}
        layer_nodes = {n: [n] for n in range(N)}
        # cutoff value
        __PASS_MAX = -1

        # heap = []
        # heapq.heapify(heap)

        random_state = self.check_random_state(random_state)

        # special case, when there is no link
        # the best partition is everyone in its community
        if not self.G.edges(data=True):
            # save the information in Z
            n = N - 1
            Z[:n, 0] = np.arange(n)
            Z[:n, 1] = np.arange(1, n + 1)
            Z[:n, 2:4] = 1
            Z[:n, 4] = root.id
            if self.__verbose__:
                print('-----', self.objective_paras.objective, 'merge phase 0',
                      'start, individuals: ', self.vertex_num, ', leaves ', self.vertex_num)
                print('-----Merge ends.-----')
            return Z

        # compute the entropy of initial tree (excluding the root node and vertices)
        new_score = self._compute_layer(layer_nodes, self.node_list)
        cur_score = new_score
        layers = [layer_nodes]
        # save the current hierarchical relationship
        node_list = self.node_list
        cur_graph = self.G.copy()
        times = 0
        if self.__verbose__:
            print("Start merging..")
            time_s = time.time()
        while True:
            # one-level merging
            layer_nodes_new, node_list = self._one_level_merge(layer_nodes=layer_nodes, node_list_=node_list, cur_graph=cur_graph, root=root, randomize_=True)
            # compute the entropy the current layer
            new_score = self._compute_layer(layer_nodes_new, node_list)
            if new_score - cur_score < self.__MIN:
                if self.__verbose__ and self.__long_message__:
                    print(layer_nodes_new)
                    print("new score: ", new_score)
                layers = self._update_layers(layers=layers, layer_nodes=layer_nodes)
                if self.__verbose__:
                    print("One layer iteration terminates.")
                break
            # It means that we should renumber, induce the graph, and update 'layers'
            # process the result
            layer_nodes_new, node_list, ass_dict = self._renumber(layer_nodes_new,
                                                                      node_list, root)
            cur_graph = self._induce_graph(cur_graph, layer_nodes_new, ass_dict)
            # if times == 0:
            #     self.G = copy.deepcopy(cur_graph)
            #     Z, z_i, parent_nodes_id = self._fill_Z_v1(layer_nodes_new, self.node_list)
            #     # remapping the nodes in self.G
            #     remapping = self._remap(parent_nodes_id)
            if self.__verbose__ and self.__long_message__:
                print(layer_nodes_new)
                print("new score: ", new_score)
            # save the result of one-level partition
            layers.append(layer_nodes_new)
            if not self.merge_layers:
                # Considering merge_flag = False, the iteration should be conducted only once.
                # We need to suppose it as merging terminates in the second time instead of the first one.
                layer_nodes = {n: [n] for n in layer_nodes_new.keys()}
                layers.append(layer_nodes)
                break
            # renumber and induce the graph
            layer_nodes = {n: [n] for n in layer_nodes_new.keys()}
            cur_score = new_score
            times += 1
        if self.__verbose__:
            time_e = time.time()
            print("Merging used: ", str(time_e-time_s), "s")
        self.G = cur_graph
        Z, z_i, i = self._binary_combine_v2(root=root, layers=layers)
        if self.__verbose__ and self.__long_message__:
            print(Z)
        # get clubs
        if len(layers) == 1:
            layer_tmp = layers[0]
        else:
            if self.merge_layers is True:
                # the process of binary combination needs one-layer tree
                layer_tmp = self._merge_into_one_layer(layers)
            else:
                # merge the multi-layer result into one-layer
                layer_tmp = layers[1]
        self.clubs = self._get_clubs_v2(layer_tmp)
        if self.auto_k:
            # set max_k as the number of clubs
            self.max_k = len(set(self.clubs))
            self.ks = range(self.min_k, self.max_k + 1)
        return Z


    def _merge_into_one_layer(self, layers):
        """
        merge the multi-layer tree into a one-layer tree
        """
        if len(layers) <= 2:
            return layers[-1]
        cur_layer = layers[1]
        for i in range(1, len(layers)-1):
            higher_layer = layers[i+1]
            for com_id in higher_layer.keys():
                # com_new = [-1]
                # for node_id in higher_layer[com_id]:
                #     com_new.extend(cur_layer[node_id])
                # com_new.remove(-1)
                # higher_layer[com_id] = com_new.copy()
                nodes = higher_layer[com_id]
                com_new = [cur_layer[node_id] for node_id in nodes]
                higher_layer[com_id] = [node for sublist in com_new for node in sublist]
            cur_layer = higher_layer
        return cur_layer

    def _merge_into_one_layer_v0(self, layers):
        """
        合并多层树结构为一层
        tried to optimize
        can lead to errors when getting clubs
        2025/03/07
        """
        if len(layers) <= 2:
            return layers[-1]

        # 预先构建父层到子层的映射字典
        parent_map = {}
        for i in range(len(layers) - 1):
            current_layer = layers[i]
            next_layer = layers[i + 1]
            for com_id, nodes in next_layer.items():
                for node in nodes:
                    if node not in parent_map:
                        parent_map[node] = []
                    parent_map[node].append(com_id)

        # 将每一层转换为使用集合存储node_id
        for i in range(len(layers)):
            layers[i] = {com_id: set(nodes) for com_id, nodes in layers[i].items()}

        cur_layer = layers[1]
        for i in range(1, len(layers) - 1):
            higher_layer = layers[i + 1]
            new_higher_layer = {}
            for com_id in higher_layer:
                # 获取当前层的节点列表
                nodes = higher_layer[com_id]
                # 合并所有子节点
                merged_nodes = []
                for node in nodes:
                    if node in parent_map:
                        merged_nodes.extend(parent_map[node])
                new_higher_layer[com_id] = set(merged_nodes)
            cur_layer = new_higher_layer

        return cur_layer

    # def _merge_into_one_layer(self, layers):
    #     """
    #     merge the multi-layer tree into a one-layer tree
    #     """
    #     if len(layers) == 1:
    #         #  should be avoided
    #         # pdb.set_trace()
    #         return layers[0]
    #     if len(layers) == 2:
    #         return layers[1]
    #     cur_layer = layers[1].copy()
    #     for i in range(1, len(layers) - 1):
    #         higher_layer = layers[i + 1].copy()
    #         for com_id in higher_layer.keys():
    #             com_new = [-1]
    #             for node_id in higher_layer[com_id]:
    #                 com_new.extend(cur_layer[node_id])
    #             com_new.remove(-1)
    #             higher_layer[com_id] = com_new.copy()
    #         cur_layer = higher_layer
    #     return cur_layer

    def _get_clubs_v2(self, layer_nodes):
        """
        Return the clubs as a list
        """
        clubs = np.zeros((self.vertex_num),dtype=int)
        flag = 0
        for i in layer_nodes.keys():
            children = layer_nodes[i]
            clubs[children] = flag
            flag += 1
        clubs = list(clubs)
        return clubs

    def _randomize(self, items, random_state):
        """Returns a List containing a random permutation of items"""
        randomized_items = list(items)
        random_state.shuffle(randomized_items)
        return randomized_items

    def _one_level_merge(self, layer_nodes, node_list_, root, cur_graph, is_leaf=True, randomize_=True, random_state=None):
        """
        Compute one level of nodes (communities)
        """
        M, m, d, log_d, d_log_d, vG, log_vG, sparse_m, d_square, M_rows = self.knn_graph_stats
        layer_nodes_new = copy.deepcopy(layer_nodes)
        modified = True
        cur_score = self._compute_layer(layer_nodes, node_list_)
        new_score = cur_score
        node_list = copy.deepcopy(node_list_)
        # node_list_org = copy.deepcopy(node_list_)
        node_list_org = node_list_
        random_state = self.check_random_state(random_state)
        while modified:
            cur_score = new_score
            modified = False
            if randomize_:
                visit_layer_nodes = self._randomize(layer_nodes_new.keys(), random_state)
            else:
                visit_layer_nodes = layer_nodes_new.keys()
            for node_id in visit_layer_nodes:
                node_ = node_list_org[node_id]
                best_delta = 0
                par_com = node_.parent_tmp
                best_com = par_com
                # get the list of all the neighbors
                neighbor_nodes = cur_graph.neighbors(node_id)
                node1 = node_list[par_com]
                parent_node1 = node_list[par_com]
                # compute the value of delta
                visited = set()

                if randomize_:
                    visit_neighs = self._randomize(neighbor_nodes,random_state)
                else:
                    visit_neighs = neighbor_nodes
                for neigh_id in visit_neighs:
                    neigh_node = node_list[neigh_id]
                    if neigh_node.parent_tmp in visited or neigh_node.parent_tmp == par_com:
                        continue
                    neigh_com_node = node_list[neigh_node.parent_tmp]
                    visited.add(neigh_node.parent_tmp)
                    delta_tmp = graph_metric.get_delta_merge_score_v3(M=M, sparse_m=sparse_m, vG=vG, d=d, node1=node1, node0=node_, node2=neigh_com_node, objective_paras=self.objective_paras, parent_status=False, parent_node1=parent_node1, parent_node2=neigh_com_node)
                    if self.__verbose__ and self.__long_message__:
                        print("           delta: ", delta_tmp)
                    if delta_tmp > best_delta:
                        best_delta = delta_tmp
                        best_com = neigh_node.parent_tmp
                if self.__verbose__ and self.__long_message__:
                    print("node id: ", node_id, "best com: ", best_com, "best delta: ", best_delta)
                if best_com != par_com:
                    # modified = True
                    # # max_n1 and max_n2 are merged, their parent is new_node
                    # # the parent of new_node is root
                    # new_node = root.merge(self.update_node_id(), node_, self.node_list[best_com],
                    #                       is_leaf=is_leaf)
                    # # Z[z_i] = [max_n1.id, max_n2.id, new_node.dist, len(new_node.vs),
                    # #           new_node.id]  # id of node1, id of node2, dist(height) of node merged with node1 and node2, # of vs, id of new node
                    # self.node_list[new_node.id] = new_node
                    # # update the layer_nodes_new
                    # layer_nodes_newer.remove(node_id)
                    # layer_nodes_newer.remove(best_com)
                    # layer_nodes_newer.append(new_node.id)
                    # # update the graph G
                    # self.G.add_node(new_node.id)
                    # for x in set(chain(self.G.neighbors(node_id),
                    #                    self.G.neighbors(best_com))):
                    #     if x == node_id or x == best_com:
                    #         continue
                    #     # node = self.node_list[x]
                    #     self.G.add_edge(x, new_node.id, weight=1)
                    # self.G.remove_node(node_id)
                    # self.G.remove_node(best_com)
                    modified = True
                    node_list[node_id].parent_tmp = best_com
                    node_list_org[node_id].parent_tmp = best_com
                    node_list[best_com]._insert(node_,self.knn_graph_stats, self.objective_paras)
                    layer_nodes_new[best_com].append(node_id)
                    layer_nodes_new[par_com].remove(node_id)
                    node_list[par_com]._remove(node_, self.knn_graph_stats, self.objective_paras)


            new_score = self._compute_layer(layer_nodes_new, node_list)
            if new_score - cur_score < self.__MIN:
                break
            # node_list_org = copy.deepcopy(node_list)
        layer_nodes_new = self._clean_layer(layer_nodes_new)
        return layer_nodes_new, node_list


    def _clean_layer(self, layer_nodes):
        """
        clean the layer, remove empty nodes without any child
        """
        # layer_nodes_new = {}
        # for node_id in layer_nodes.keys():
        #     list_tmp = layer_nodes[node_id]
        #     if len(list_tmp) == 0:
        #         continue
        #     else:
        #         layer_nodes_new[node_id] = copy.deepcopy(list_tmp)
        layer_nodes_new = {node_id: nodes for node_id, nodes in layer_nodes.items() if
                           nodes}
        return layer_nodes_new

    def _renumber(self, layer_nodes, node_list, root):
        """
        Renumber the nodes and save the hierarchical relationship
        """
        layer_nodes_new = {}
        node_list_new = {}
        ass_dict = {}
        # for i in layer_nodes.keys():
        #     if layer_nodes[i]:
        #         com_node = node_list[i]
        #         new_node = Node(self.knn_graph_stats, self.update_node_id(), [], com_node.vs, parent=root, objective_paras=self.objective_paras)
        #         new_node_id = new_node.id
        #         node_list_new[new_node_id] = new_node
        #         self.node_list[new_node_id] = new_node
        #         layer_nodes_new[new_node_id] = layer_nodes[i]
        #         for node_id in layer_nodes[i]:
        #             ass_dict[node_id] = new_node_id
        #             self.node_list[node_id].parent_tmp = new_node_id
        for i, nodes in layer_nodes.items():
            if nodes:
                com_node = node_list[i]
                new_node = Node(self.knn_graph_stats, self.update_node_id(), [], com_node.vs, parent=root, objective_paras=self.objective_paras)
                new_node_id = new_node.id
                node_list_new[new_node_id] = new_node
                self.node_list[new_node_id] = new_node
                layer_nodes_new[new_node_id] = nodes
                for node_id in nodes:
                    ass_dict[node_id] = new_node_id
                    self.node_list[node_id].parent_tmp = new_node_id
        return layer_nodes_new, node_list_new, ass_dict

    def _renumber_v0(self, layer_nodes, node_list, root):
        """
        Renumber the nodes and save the hierarchical relationship
        tried but failed
        2024/03/08
        """
        layer_nodes_new = {}
        node_list_new = {}
        ass_dict = {}
        # 缓存实例变量到局部变量以加速访问
        node_list_inst = self.node_list
        update_node_id = self.update_node_id  # 缓存方法引用

        for i, nodes in layer_nodes.items():
            if not nodes:
                continue

            # 创建新节点并获取必要属性
            new_node = Node(
                self.knn_graph_stats,
                update_node_id(),  # 使用缓存的方法引用
                [],
                node_list[i].vs,  # 直接访问属性，避免中间变量
                parent=root,
                objective_paras=self.objective_paras
            )
            new_node_id = new_node.id
            new_parent = new_node_id  # 缓存到局部变量

            # 更新数据结构
            node_list_new[new_node_id] = node_list_inst[new_node_id] = new_node
            layer_nodes_new[new_node_id] = nodes

            # 批量更新ass_dict
            ass_dict.update(zip(nodes, [new_parent] * len(nodes)))

            # 预取节点对象并批量更新parent_tmp
            nodes_to_update = (node_list_inst[node_id] for node_id in nodes)
            for node in nodes_to_update:
                node.parent_tmp = new_parent

        return layer_nodes_new, node_list_new, ass_dict


    def _fill_Z_tmp(self, layers):
        """
        fill in the data Z (multi level)
        only according to the result of merging operation
        (instead of using binary combination randomly)
        """
        Z = np.zeros((self.vertex_num - 1, 5))
        z_i = 0
        height = len(layers)
        # if height == 1:
        #     print("height: 1")
        if height == 1:
            mapping = dict(zip(list(range(self.vertex_num)), list(range(self.vertex_num))))
            # fill Z correctly
            parent_id = self.vertex_num
            Z[z_i] = [0,1,1,2,parent_id]
            z_i += 1
            vertices_num = 2
            while z_i < self.vertex_num-1:
                Z[z_i] = [z_i+1, parent_id,1,vertices_num+1,parent_id+1]
                parent_id += 1
                z_i += 1
                vertices_num += 1
        else:   # height>=2
            cur_layer = copy.deepcopy(layers[1])
            parent_node_id = self.vertex_num - 1
            # single_children = []
            parent_nodes_fin = []
            # first one layer
            for new_node_id in cur_layer.keys():
                children = cur_layer[new_node_id]
                if len(children) == 1:
                    # single_children.append(children[0])
                    parent_nodes_fin.append(children[0])
                else:
                    n1 = children[0]
                    n2 = children[1]
                    len_vs = len(self.node_list[n1].vs) + len(self.node_list[n2].vs)
                    parent_node_id += 1
                    Z[z_i] = [n1, n2, 1, len_vs, parent_node_id]
                    z_i += 1
                    for index in range(2, len(children)):
                        n1 = children[index]
                        n2 = parent_node_id
                        len_vs += len(self.node_list[n1].vs)
                        parent_node_id += 1
                        Z[z_i] = [n1, n2, 1, len_vs, parent_node_id]
                        z_i += 1
                    parent_nodes_fin.append(parent_node_id) # the node id not in the first two column of Z
            mapping = {}  # map the node id in the node list to the node_id in Z
            # save mapping from singleton nodes
            if height > 2:
                p_i = 0
                for org_node_id in layers[1].keys():
                    # map the indices in self.node_list to indices in parent_nodes_fin
                    mapping[org_node_id] = parent_nodes_fin[p_i]
                    p_i += 1
                for h in range(2, height):
                    cur_layer = copy.deepcopy(layers[h])
                    for new_node_id in cur_layer.keys():
                        children = cur_layer[new_node_id]
                        if len(children) == 1:
                            # single_children.append(children[0])
                            # cur_node_id = mapping[children[0]]
                            # parent_nodes_fin.append(cur_node_id)
                            if children[0] in mapping:
                                mapping[new_node_id] = mapping[children[0]]
                            else:
                                mapping[new_node_id] = children[0]
                            continue
                        else:
                            if self.__verbose__ and self.__long_message__:
                                print(mapping)
                            if children[0] not in mapping.keys() or children[1] not in mapping.keys():
                                print("mapping: ", mapping)
                                print("layers: ", layers)
                                print("children: ", children)
                                pdb.set_trace()
                            n1 = mapping[children[0]]
                            n2 = mapping[children[1]]
                            # i1_find = np.where(Z[:, 4] == n1)[0][0]
                            # i2_find = np.where(Z[:, 4] == n2)[0][0]
                            # len_vs = Z[i1_find, 3] + Z[i2_find, 3]
                            if n1 not in Z[:z_i, 4]:
                                len_vs_n1 = 1
                            else:
                                i1_find = np.where(Z[:z_i, 4] == n1)[0][0]
                                len_vs_n1 = Z[i1_find, 3]
                            if n2 not in Z[:z_i, 4]:
                                len_vs_n2 = 1
                            else:
                                i2_find = np.where(Z[:z_i, 4] == n2)[0][0]
                                len_vs_n2 = Z[i2_find, 3]
                            len_vs = len_vs_n1 + len_vs_n2
                            parent_node_id += 1
                            Z[z_i] = [n1, n2, h, len_vs, parent_node_id]
                            z_i += 1
                            del mapping[children[0]]
                            del mapping[children[1]]
                            for index in range(2, len(children)):
                                if not (children[index]) in mapping:
                                    print("mapping: ", mapping)
                                    print("children: ", children)
                                    print("index: ", index)
                                    print(children[index])
                                n1 = mapping[children[index]]
                                n2 = parent_node_id
                                if not n1 in Z[:z_i,4]:
                                    len_vs += 1
                                else:
                                    i1_find = np.where(Z[:z_i, 4] == n1)[0][0]
                                    len_vs += Z[i1_find, 3]
                                parent_node_id += 1
                                Z[z_i] = [n1, n2, h, len_vs, parent_node_id]
                                z_i += 1
                                del mapping[children[index]]
                            mapping[new_node_id] = parent_node_id
                parent_nodes_fin = list(mapping.values())
                # additional process of matrix Z
                nodes_id = list(layers[-1].keys())
                for node_id in nodes_id:
                    if len(layers[height - 1][node_id]) > 1:
                        print("not single child in the highest layer")
                        continue
                    node_id_prev = layers[height - 1][node_id][0]
                    for h in range(height - 2, -1, -1):
                        if len(layers[h][node_id_prev]) > 1:
                            mapping[node_id] = mapping[node_id_prev]
                            break
                        if h == 0:
                            mapping[node_id] = layers[h][node_id_prev][0]
                            break
                        node_id_prev = layers[h][node_id_prev][0]
        # print(Z)
        return Z, z_i, mapping

    def _order(self, layer):
        """
        output the cell order
        """
        first_flag = True
        for parent in layer.keys():
            children = layer[parent]
            if first_flag:
                order_list = children
            else:
                order_list.extend(children)
        self.order_list = order_list

    def _order_v0(self, layer):
        """
        output the cell order
        tried to optimize
        2025/03/07
        """
        self.order_list = []
        for children in layer.values():
            self.order_list.extend(children)

    def _induce_graph(self, cur_graph, layer_nodes_new, ass_dict, weight="weight"):
        """
        Produce the graph where nodes are the communities
        """
        ret = nx.Graph()
        ret.add_nodes_from(layer_nodes_new.keys())
        for node1, node2, datas in cur_graph.edges(data=True):
            edge_weight = datas.get(weight, 1)
            com1 = ass_dict[node1]
            com2 = ass_dict[node2]
            w_prec = ret.get_edge_data(com1, com2, {weight: 0}).get(weight, 1)
            ret.add_edge(com1, com2, **{weight: w_prec + edge_weight})
        return ret


    def _update_edges(self, n1_id, n2_id, new_node_id, weight="weight"):
        nodes_set = {n1_id, n2_id}
        for node1_id, node2_id, datas in list(self.G.edges(data=True)):
            if node1_id not in nodes_set and node2_id not in nodes_set:
                continue
            if node1_id == n1_id and node2_id == n2_id or node1_id == n2_id and node2_id == n1_id:
                # n1 and n2 are combined, the edges are inside and should be ignored
                continue
            if node1_id == node2_id:
                # self connected edges are ignored
                continue
            if node1_id in nodes_set:
                node_c_id = node2_id
            elif node2_id in nodes_set:
                node_c_id = node1_id
            else:
                print("Wrong node!")
                pdb.set_trace()
                continue
            edge_weight = datas.get(weight, 1)
            w_prec = self.G.get_edge_data(new_node_id, node_c_id, {weight: 0}).get(weight, 1)
            self.G.add_edge(new_node_id, node_c_id, **{weight: w_prec + edge_weight})

    def _update_edges_v0(self, n1_id, n2_id, new_node_id, weight="weight"):
        """
        tried to optimize
        2025/03/08
        """
        # 获取与n1_id和n2_id相关的所有边
        edges_n1 = self.G.edges(n1_id, data=True)
        edges_n2 = self.G.edges(n2_id, data=True)

        # 避免重复处理相同的边，特别是n1- n2之间的边
        seen_edges = set()

        for node1, node2, datas in edges_n1:
            if (node1 == n1_id and node2 == n2_id) or (node1 == n2_id and node2 == n1_id):
                continue  # 跳过n1和n2之间的边
            if node1 == node2:
                continue  # 跳过自环边

            edge_weight = datas.get(weight, 1)

            # 确定连接的节点
            connected_node = node2 if node1 == n1_id else node1

            key = (new_node_id, connected_node)
            if key not in seen_edges:
                seen_edges.add(key)

                # 获取现有权重或默认为0
                existing_weight = self.G.get_edge_data(new_node_id, connected_node, default={}).get(weight, 0)
                total_weight = existing_weight + edge_weight

                # 更新边的权重
                if total_weight > 0:
                    self.G.add_edge(new_node_id, connected_node, **{weight: total_weight})

        for node1, node2, datas in edges_n2:
            if (node1 == n1_id and node2 == n2_id) or (node1 == n2_id and node2 == n1_id):
                continue  # 跳过n1和n2之间的边
            if node1 == node2:
                continue  # 跳过自环边

            edge_weight = datas.get(weight, 1)

            # 确定连接的节点
            connected_node = node2 if node1 == n2_id else node1

            key = (new_node_id, connected_node)
            if key not in seen_edges:
                seen_edges.add(key)

                # 获取现有权重或默认为0
                existing_weight = self.G.get_edge_data(new_node_id, connected_node, default={}).get(weight, 0)
                total_weight = existing_weight + edge_weight

                # 更新边的权重
                if total_weight > 0:
                    self.G.add_edge(new_node_id, connected_node, **{weight: total_weight})


    def _binary_combine_v2(self, root, layers, only_positive_delta=False):
        """
        combine to obtain the binary tree
        """
        if self.__verbose__:
            print("Computing Z..")
            time_s = time.time()
        Z, z_i, mapping = self._fill_Z_tmp(layers)
        if self.__verbose__:
            time_e = time.time()
            print("Building Z matrix used: ", str(time_e-time_s), "s")
        M, m, d, log_d, d_log_d, vG, log_vG, sparse_m, d_square, M_rows = self.knn_graph_stats
        # obtain leaves
        if len(layers) == 1 or self.merge_layers:
            layer = layers[-1]
        else:
            layer = layers[1]
        self.leaves = [l for l in layer.keys()]
        # leaves = {}
        # for l in layer.keys():
        #     leaves[l] = 1
        leaves = set(layer.keys())
        heap = []
        heapq.heapify(heap)
        i = 0
        ns = [(n1, n2) for n1, n2, _ in self.G.edges(data=True) if
              n1 != n2]  # node pairs (excluding self-self)
        if not ns:
            ns = itertools.combinations(leaves, 2)
        # obtain the heap
        for n1, n2 in ns:  # for all possible pairs
            node1, node2 = self.node_list[n1], self.node_list[n2]
            # try to add a parent of n1 and n2
            if self.objective_paras.objective == 'SE':
                delta = graph_metric.get_delta_combine_se_v2(M, sparse_m, vG, root, node1,
                                                          node2)
            elif self.objective_paras.objective == 'M':
                delta = graph_metric.get_delta_combine_nm(M=M, sparse_m=sparse_m, vG=vG,
                                                          parent=root, node1=node1,
                                                          node2=node2,
                                                          objective_paras=self.objective_paras)
                # print("delta: ", delta)
            elif self.objective_paras.objective == 'KL':
                delta = graph_metric.get_delta_combine_kl_v2(M=M, sparse_m=sparse_m, vG=vG,
                                                          parent=root, node1=node1,
                                                          node2=node2, objective_paras=self.objective_paras)
            if only_positive_delta and delta < 0:
                continue
            heapq.heappush(heap, (-delta, n1, n2))
            i += 1
        if self.__verbose__:
            time_2 = time.time()
            print("Building heap used: ", str(time_2 - time_e), "s")
        while z_i < self.vertex_num - 1:
            max_n1, max_n2, max_delta = self.get_max_delta_from_heap(heap, leaves, leaves)

            if max_n1 is None:
                break

            new_node = root.merge(self.update_node_id(), max_n1, max_n2, self.knn_graph_stats, objective_paras=self.objective_paras,is_leaf=False)
            self.node_list[new_node.id] = new_node
            new_node_id = Z[:,4].max() + 1
            mapping[new_node.id] = new_node_id

            max_n1_id = mapping[max_n1.id]
            max_n2_id = mapping[max_n2.id]
            Z[z_i] = [max_n1_id, max_n2_id, new_node.dist, len(new_node.vs), new_node_id]
            if (Z[z_i,3] == Z[z_i-1,3]) and (Z[z_i,0]==0 or Z[z_i,1]==0):
                print(Z)
                print("line 1570: There may something wrong with matrix Z!")
                # pdb.set_trace()

            # update
            # del leaves[max_n1.id]
            # del leaves[max_n2.id]
            leaves.remove(max_n1.id)
            leaves.remove(max_n2.id)
            self.G.add_node(new_node.id)  # O(k)

            xs = set(chain(self.G.neighbors(max_n1.id), self.G.neighbors(max_n2.id)))  # all neighbors of max_n
            # xs = []  # not solving HC problem
            if not xs:
                xs = leaves
            for x in xs:
                node = self.node_list[x]
                if self.objective_paras.objective == 'SE':
                    delta = graph_metric.get_delta_combine_se_v2(M, sparse_m, vG, root, node, new_node)
                elif self.objective_paras.objective == 'M':
                    delta = graph_metric.get_delta_combine_nm(M=M, sparse_m=sparse_m, vG=vG, parent=root, node1=node, node2=new_node, objective_paras=self.objective_paras)
                elif self.objective_paras.objective == 'KL':
                    delta = graph_metric.get_delta_combine_kl_v2(M=M, sparse_m=sparse_m, vG=vG, parent=root, node1=node, node2=new_node, objective_paras=self.objective_paras)
                else:
                    print("ERROR! Input wrong objective.")
                if only_positive_delta and delta < 0:
                    continue
                heapq.heappush(heap, (-delta, x, new_node.id))
                i += 1
            # update information of edges
            self._update_edges(max_n1.id, max_n2.id, new_node.id)
            self.G.remove_node(max_n1.id)
            self.G.remove_node(max_n2.id)
            # leaves[new_node.id] = 1
            leaves.add(new_node.id)

            z_i += 1
        if self.__verbose__:
            time_3 = time.time()
            print("Meaningful combination used: ", str(time_3 - time_2), "s")
        if z_i < self.vertex_num-1:
            if self.__verbose__:
                print("No more meaningful combination")
            nodes_id = list(self.G.nodes)
            n1_id = nodes_id[0]
            n2_id = nodes_id[1]
            n1 = self.node_list[n1_id]
            n2 = self.node_list[n2_id]
            new_node = root.merge(self.update_node_id(), n1, n2, self.knn_graph_stats,objective_paras=self.objective_paras,is_leaf=False)
            self.node_list[new_node.id] = new_node
            new_node_id = Z[:, 4].max() + 1
            mapping[new_node.id] = new_node_id
            Z[z_i] = [mapping[n1_id], mapping[n2_id], new_node.dist, len(new_node.vs),
                      new_node_id]
            if Z[z_i,3] == Z[z_i-1,3]:
                print(Z)
                print("line 1614: corrupt Z warning")
            z_i += 1
            for node_id in nodes_id[2:]:
                n2 = new_node
                n1 = self.node_list[node_id]
                # find the key of the given value
                new_node = root.merge(self.update_node_id(), n1, n2, self.knn_graph_stats,objective_paras=self.objective_paras, is_leaf=False)
                self.node_list[new_node.id] = new_node
                new_node_id = Z[:, 4].max() + 1
                mapping[new_node.id] = new_node_id
                n1_id = mapping[node_id]
                Z[z_i] = [n1_id, Z[:, 4].max(), new_node.dist, len(new_node.vs), new_node_id]
                if Z[z_i,3] == Z[z_i-1,3]:
                    print(Z)
                    print("may be error of Z")
                z_i += 1
        if self.__verbose__:
            time_4 = time.time()
            print("Non-meaningful combination used: ", str(time_4 - time_3), "s")
        return Z, z_i, i   # i: the times of heappush operation

    def order_tree(self):
        M, m, d, log_d, d_log_d, vG, log_vG, sparse_m, d_square, M_rows = self.knn_graph_stats
        node = self.root
        A = M[np.ix_(node.vs, node.vs)]
        L = csgraph.laplacian(A, normed=True)
        eig_values, eig_vectors = la.eigh(L)

        unique_eig_values = np.sort(list(set(eig_values.real)))
        fiedler_pos = np.where(eig_values.real == unique_eig_values[1])[0][0]
        fiedler_vector = np.transpose(eig_vectors)[fiedler_pos]
        self._order_tree_aux(node, fiedler_vector)

    def _order_tree_aux(self, node, vector): #switch the left and right to ensure median of vs of left is no larger than that of right
        if len(node.children) == 0:
            return
        left_v = np.median(vector[node.children[0].vs])
        right_v = np.median(vector[node.children[1].vs])
        if left_v > right_v:  # swith
            node.chidlren = node.children[::-1]
            row = self.Z[node.id - self.vertex_num]
            self.Z[node.id - self.vertex_num] = [row[1], row[0], row[2], row[3]]

        for child in node.children:
            self._order_tree_aux(child, vector)


    def contract_tree_v2(self, Z, n_clusters):
        # update node distance
        se_scores, ks_subpopulation_node_ids, ks_clusters, optimal_k = self._contract_tree_dp_v2(self.root)
        tmp = pd.DataFrame(np.matrix(se_scores), columns=self.ks).T
        tmp['K'] = tmp.index
        tmp.columns = ['SE Score', 'K']
        self.se_scores = tmp
        delta_se_scores = se_scores[1:] - se_scores[:-1]
        tmp = pd.DataFrame(np.matrix(delta_se_scores), columns=self.ks[1:]).T
        tmp['K'] = tmp.index
        tmp.columns = ['Delta SE Score', 'K']
        self.delta_se_scores = tmp
        self.ks_clusters = pd.DataFrame(np.matrix(ks_clusters).T, columns=['K={}'.format(k) for k in self.ks])
        Z_clusters = hierarchy.cut_tree(Z[:, :4], n_clusters=n_clusters)
        self.Z_clusters = pd.DataFrame(np.matrix(Z_clusters), columns=['K={}'.format(k) for k in n_clusters])

        self.optimal_k = optimal_k
        self.optimal_clusters = self.ks_clusters['K={}'.format(self.optimal_k)].tolist()
        self.optimal_subpopulation_node_ids = ks_subpopulation_node_ids[optimal_k-1]
        return


    def _contract_tree_dp_v2(self, root):
        if self.strategy == 'bottom_up':
            node_ids = self.leaves + list(range(self.leaves[-1]+1, len(self.node_list)-1))
        else:
            node_ids = self.leaves[::-1] + self.beyond_leaves[::-1]
        # print(node_ids)  # the nodes ordered from leaf to parent

        cost_m = np.zeros((len(node_ids), self.max_k+1))
        cost_m.fill(np.nan)
        cutoff_m = np.zeros((len(node_ids), self.max_k+1))
        cutoff_m.fill(-1)

        self._dp_compute_cost(cost_m, cutoff_m, node_ids)
        if self.__verbose__:
            print('se cost')
            print(cost_m[-1, :])
        # print(cost_m)
        # print(cost_m.shape)
        # print(cutoff_m)
        ks_clusters = []  # list of lists, the result of clustering
        ks_subpopulation_node_ids = []
        # if self.min_k == 1:
        #     ks_clusters.append([0]*self.vertex_num)
        #     ks_subpopulation_node_ids.append([self.root.id])
        # k_start = max(self.min_k,2)
        # for k in range(k_start, self.max_k):
        for k in self.ks:
            if k == 1:
                ks_clusters.append([0]*self.vertex_num)
                ks_subpopulation_node_ids.append([self.root.id])
                continue
            subpopulation_node_ids = []
            self._trace_back(root, cost_m, cutoff_m, subpopulation_node_ids, k)
            clusters = [(v, i) for i, c in enumerate(subpopulation_node_ids) for v in self.node_list[c].vs]
            clusters = sorted(clusters)
            clusters = [c for v, c in clusters]
            if len(clusters) != self.vertex_num:  # happens in bottom up node if k larger than number of leaves
                ks_clusters.append([0]*self.vertex_num)
                ks_subpopulation_node_ids.append([])
                continue

            ks_clusters.append(clusters)
            ks_subpopulation_node_ids.append(subpopulation_node_ids)

        optimal_k = self.max_k - np.nanargmin(cost_m[-1, 2:][::-1])
        if self.__verbose__:
            print(cost_m[-1, :][optimal_k])
        return cost_m[-1, 1:], ks_subpopulation_node_ids, ks_clusters, optimal_k

    def _dp_compute_cost(self, cost_m, cutoff_m, node_ids, contract_type='2D'):

        M, m, d, log_d, d_log_d, vG, log_vG, sparse_m, d_square, M_rows = self.knn_graph_stats
        aff_M, aff_m, aff_d, aff_log_d, aff_d_log_d, aff_vG, aff_log_vG, aff_sparse_m, aff_d_square, aff_M_rows = self.aff_graph_stats
        for dp_i, node_id in enumerate(node_ids):
            node = self.node_list[node_id]
            node.dp_i = dp_i
            node_g = graph_metric.get_g(aff_M, sparse_m, node.vs,M_rows)  # sparse_m is not used
            node_V = graph_metric.get_v(aff_M_rows, sparse_m, node.vs)
            if node.parent and contract_type == 'high-dimensional':
                parent_V = graph_metric.get_v(aff_M_rows, sparse_m, node.parent.vs)
            else:  # root node or 2D hirerachy
                parent_V = aff_vG
            # node_se = node_g/aff_vG*log2(parent_V/node_V)
            node_se = graph_metric.get_node_score_v2(vG=aff_vG, g=node_g, V=node_V, pV=parent_V,
                                                     objective_paras=self.objective_paras)
            node_se = - node_se
            if self.min_k == 1:
                node_d_log_d = np.sum(aff_d_log_d[node.vs])
                node_d_square = np.sum(aff_d_square[node.vs])
                cost_m[node.dp_i, 1] = node_se - graph_metric.get_node_vertices_score(vG=aff_vG, d_log_d=node_d_log_d, V=node_V, d_square=node_d_square, objective_paras=self.objective_paras)
            k_start = max(self.min_k, 2)
            # for k in self.ks:
            for k in range(k_start, self.max_k+1):
                # if node_id>self.vertex_num:
                #     pdb.set_trace()
                # if k == 1:
                    # node_d_log_d = np.sum(aff_d_log_d[node.vs])
                    # node_d_square = np.sum(aff_d_square[node.vs])
                    # cost_m[node.dp_i, k] = node_se - graph_metric.get_node_vertices_score(vG=aff_vG, d_log_d=node_d_log_d, V=node_V, d_square=node_d_square, objective_paras=self.objective_paras)
                        # pdb.set_trace()
                    # - (node_d_log_d - node_V*log2(node_V))/aff_vG
                    # print(node.id, k, cost_m[node.dp_i, k], len(node.vs))
                    # continue

                # the number of vertices less than k or node is a leaf node
                # if not node.children and k>1, the loop can be broken
                if len(node.vs) < k or not node.children:
                    cost_m[node.dp_i, k] = np.inf
                    continue

                l_id = node.children[0].dp_i
                r_id = node.children[1].dp_i
                min_i = None
                min_cost = np.inf
                cost_m_l = cost_m[l_id, 1:k]
                cost_m_r = cost_m[r_id, k-1:0:-1]

                if contract_type == 'high_dimensional':
                    cost_m_l += node_se
                cost_m_sum = cost_m_l + cost_m_r
                min_cost = np.nanmin(cost_m_sum)
                min_i = np.nanargmin(cost_m_sum) + 1

                if node.dp_i == None:
                    pdb.set_trace()
                cost_m[node.dp_i, k] = min_cost
                cutoff_m[node.dp_i, k] = min_i
                # print("node_id: ",node_id, "k: ", k)
                # print("cost_m_l: ",cost_m_l)
                # print("cost_m_r: ", cost_m_r)
                # print("cost_m_sum: ", cost_m_sum)
                # pdb.set_trace()
                # for i in range(1, k):
                #     cost = cost_m[l_id, i] + cost_m[r_id, k-i]
                #     if contract_type == 'high_dimensional':
                #         cost += node_se
                #     if cost < min_cost:
                #         min_cost = cost
                #         min_i = i
                # cost_m[node.dp_i, k] = min_cost
                # # pdb.set_trace()
                # cutoff_m[node.dp_i, k] = min_i
                # # print(node.id, k, cost_m[node.dp_i, k], len(node.vs))
        # if self.__verbose__:
        #     print(cost_m)
        # pdb.set_trace()

    def _trace_back(self, node, cost_m, cutoff_m, clusters, k_hat):
        k_prime = cutoff_m[node.dp_i, k_hat]
        if np.isnan(k_prime) or k_prime == -1:
            return
        k_prime = int(k_prime)
        left_node = node.children[0]
        right_node = node.children[-1]

        if k_prime > 1:
            self._trace_back(left_node, cost_m, cutoff_m, clusters, k_prime)
        else:
            clusters.append(left_node.id)
        if k_prime < k_hat-1:  # k_prime<=1
            self._trace_back(right_node, cost_m, cutoff_m, clusters, k_hat-k_prime)
        else:
            clusters.append(right_node.id)

    def get_tree_se(self):
        return 0
        tree_se = 0
        self.get_tree_se_aux(self.root, tree_se)
        print(tree_se)

    def get_tree_se_aux(self, node, tree_se):
        tree_se += node.se
        if len(node.children) != 0:
            tree_se += node.se

    def to_newick(self):
        return '({});'.format(self._to_newick_aux(self.root, is_root=True))

    def _to_newick_aux(self, node, is_root=False):
        if len(node.vs) == 1:
            return 'n{}:{}'.format(node.id, 1)

        if node.is_leaf:
            if self.strategy == 'bottom_up':
                res = self._to_newick_leaf_bottom_up(node)
            else:
                res = self._to_newick_leaf_top_down(node)
        else:
            res = ','.join([self._to_newick_aux(c) for c in node.children])

        return '({})n{}:{}'.format(res, node.id, 1)

    def _to_newick_leaf_bottom_up(self, node):
        if len(node.vs) == 1:
            return 'n{}'.format(node.id, 1)
        else:
            return ','.join([self._to_newick_leaf_bottom_up(self.node_list[v]) for v in node.vs])

    def _to_newick_leaf_top_down(self, node):
        try:
            return ','.join([self._to_newick_leaf_top_down(v) for v in node.vs])
        except Exception:
            return 'n{}:{}'.format(node, 1)

    def to_split_se(self):
        split_dict_list = []
        self._to_split_se_aux(self.root, split_dict_list)
        df = pd.DataFrame(split_dict_list)
        return df

    def _to_split_se_aux(self, node, split_dict_list, subpopulation=False, club=False):
        if node.id in self.optimal_subpopulation_node_ids:
            subpopulation = True
        if node.id in self.leaves:
            club = True
        if self.strategy == 'bottom_up':
            split_se = self._get_dividing_delta(node, node.children)
        else:
            split_se = node.split_se
        split_dict_list.append({
            'node_id': node.id,
            'split_se': split_se,
            'vertex_num': len(node.vs),
            'subpopulation': subpopulation,
            'club': club,
        })
        for child in node.children:
            self._to_split_se_aux(child, split_dict_list, subpopulation, club)


class SEATParas:
    def __init__(self, objective='SE',
                 eta_mode="coefficient", eta=1,
                 eta1=1,
                 eta2=1, hype_a=0):
        """
        save the paras
        """
        self.objective = objective
        self.eta_mode = eta_mode
        self.eta = eta
        self.eta1 = eta1
        self.eta2 = eta2
        self.hype_a = hype_a

    def _get_vars_(self):
        return self.objective, self.eta_mode, self.eta, self.eta1, self.eta2, self.hype_a


class TEC(AgglomerativeClustering):

    def __init__(self, min_k=1, max_k=10, auto_k=False,
                 a=None,
                 affinity='precomputed',
                 sparsification='knn_neighbors',
                 knn_m=None,
                 strategy='bottom_up',
                 objective='SE',
                 n_neighbors=10,
                 dist_topk=5,
                 split_se_cutoff=0.05,
                 kernel_gamma=None,
                 eta=2,
                 eta_mode="coefficient",
                 eta1=1,
                 eta2=1, 
                 hype_a=0, merge_layers=False, plot_cluster_map=False, random_seed=None,
                 __verbose__=False
                 ):
        self.min_k = min_k
        self.max_k = max_k
        self.auto_k = auto_k
        self.ks = range(min_k, max_k+1)
        self.a = a
        self.affinity = affinity
        self.sparsification = sparsification
        self.strategy = strategy
        self.objective_paras = SEATParas(objective=objective,eta_mode=eta_mode,eta=eta,eta1=eta1,eta2=eta2, hype_a=hype_a)
        self.n_neighbors = n_neighbors
        self.dist_topk = dist_topk
        self.split_se_cutoff = split_se_cutoff

        self.kernel_gamma = kernel_gamma

        self.knn_m = knn_m

        self.merge_layers = merge_layers
        self.__MIN = 0.0000001  # for stop merging
        self.__verbose__ = __verbose__
        self.plot_cluster_map = plot_cluster_map
        self.random_seed = random_seed  # random seed, int or None

    def construct_affinity(self, X, spatial=None, spatial_weight=None):
        # https://scikit-learn.org/stable/modules/metrics.html

        if self.affinity == 'precomputed':
            aff_m = X

        elif self.affinity == 'cosine_similarity':
            aff_m = pairwise_kernels(X, metric='cosine')

        elif self.affinity == 'linear_kernel':
            aff_m = pairwise_kernels(X, metric='linear')

        elif self.affinity == 'gaussian_kernel':
            if self.kernel_gamma:
                sigma = self.kernel_gamma
            else:
                sigma = X.std()
            print('sigma', sigma)
            gamma = 1/(sigma*sigma)
            aff_m = pairwise_kernels(X, metric='rbf', gamma=gamma)

        elif self.affinity == 'gaussian_kernel_with_spatial':
            if self.kernel_gamma:
                sigma = self.kernel_gamma
            else:
                sigma = X.std()
            print('sigma', sigma)
            gamma = 1/(sigma*sigma)
            aff_m = pairwise_kernels(X, metric='rbf', gamma=gamma)
            spa_count = pairwise_kernels(spatial, metric='rbf', gamma=gamma)
            aff_m = aff_m + spa_count * spatial_weight


        elif self.affinity == 'gaussian_kernel_topk':
            if self.kernel_gamma:
                sigma = self.kernel_gamma
            else:
                sigma = X.std()
            if self.dist_topk > X.shape[1]:
                dist_topk = X.shape[1]
            else:
                dist_topk = self.dist_topk
            n = X.shape[0]
            aff_m = np.zeros((n, n))
            for i in range(n):
                for j in range(i, n):
                    dist = (X[i] - X[j])**2
                    dist.sort()
                    dist = dist[-dist_topk:]
                    v = np.exp(-np.sum(dist)/sigma/sigma)
                    aff_m[i][j] = v
                    aff_m[j][i] = v

        elif self.affinity == 'laplacian_kernel':
            if self.kernel_gamma:
                gamma = self.kernel_gamma
            else:
                gamma = 0.1
            aff_m = pairwise_kernels(X, metric='laplacian', gamma=gamma)

        elif self.affinity == 'knn_neighbors_from_X':
            aff_m = kneighbors_graph(X, self.n_neighbors).toarray()
            aff_m = (aff_m + aff_m.T)/2
            aff_m[np.nonzero(aff_m)] = 1

        if (aff_m < 0).any():
            aff_m = aff_m + np.abs(np.min(aff_m))

        self.aff_m = aff_m
        print("Successfully construct affinity matrix.")

    def graph_sparsification(self, X):
        knn_m = None
        if self.sparsification == 'affinity':
            knn_m = copy.deepcopy(self.aff_m)
        elif self.sparsification == 'precomputed':
            knn_m = self.knn_m
        if self.sparsification == 'knn_neighbors':
            k = self.n_neighbors
            knn_m = np.zeros((X.shape[0], X.shape[0]))
            for i in range(X.shape[0]):
                ids = np.argpartition(self.aff_m[i], -k)[-k:]
                top_set = set(self.aff_m[i, ids])
                if len(top_set) == 1:
                    b = self.aff_m[i] == top_set.pop()
                    ids = []
                    offset = 1
                    left = True
                    while len(ids) < k:
                        if left:
                            idx = i + offset
                        else:
                            idx = i - offset
                        if idx < 0 or idx > len(b)-1:
                            offset += 1
                            left = not left
                            continue
                        if b[idx]:
                            ids.append(idx)
                        offset += 1
                        left = not left
                knn_m[i, ids] = 1

            knn_m = (knn_m + knn_m.T)/2
            knn_m[np.nonzero(knn_m)] = 1

        if self.sparsification == 'knn_neighbors_from_X':
            knn_m = kneighbors_graph(X, self.n_neighbors).toarray()
            # knn_m = (knn_m + knn_m.T)/2
            # knn_m[np.nonzero(knn_m)] = 1
            X = StandardScaler().fit_transform(X)
            # connectivity = 0.5 * (connectivity + connectivity.T)
            # connectivity = connectivity.toarray()
            knn_m = kneighbors_graph(X, self.n_neighbors, mode='connectivity', include_self=False)
            # knn_m = knn_m.maximum(knn_m.T)  # org mode
            # knn_m = lil_matrix(knn_m, dtype='int8')
            knn_m[knn_m.nonzero()] = 1
            knn_m = 0.5 * (knn_m + knn_m.T)
            knn_m = knn_m.toarray()

        self.knn_m = knn_m
        print("sparse done")


    def generate_graph(self, X, spatial=None, spatial_weight=None):
        """
        construct affinity graph and sparse graph
        """
        X = self._validate_data(X, ensure_min_samples=2, estimator=self)

        if self.min_k is not None and self.min_k <= 0:
            raise ValueError("min_k should be an integer greater than 0."
                             " %s was provided." % str(self.min_k))

        if self.max_k is not None and self.max_k <= 2:
            raise ValueError("max_k should be an integer greater than 2."
                             " %s was provided." % str(self.max_k))

        if self.affinity not in ['precomputed', 'gaussian_kernel',
                                 'gaussian_kernel_topk', 'linear_kernel',
                                 'cosine_similarity', 'knn_neighbors_from_X',
                                 'laplacian_kernel']:
            raise ValueError(
                "affinity should be precomputed, gaussian_kernel, linear_kernel, cosine_similarity, knn_neighbors_from_X, "
                "laplacian_kernel. "
                "%s was provided." % str(self.affinity))

        if self.sparsification not in ['affinity', 'precomputed', 'knn_neighbors',
                                       'knn_neighbors_from_X']:
            raise ValueError(
                "sparsification should be affinity, precomputed, knn_neighbors, knn_neighbors_from_X."
                " %s was provided." % str(self.sparsification))

        if self.strategy not in ['bottom_up', 'top_down']:
            raise ValueError("affinity should be bottom_up, top_down."
                             " %s was provided." % str(self.strategy))

        print('fit', self.strategy)
        self.construct_affinity(X,spatial=spatial,spatial_weight=spatial_weight)
        self.graph_sparsification(X)

    def fit_v3(self, X,spatial=None,spatial_weight=None, y=None):

        X = self._validate_data(X, ensure_min_samples=2, estimator=self)

        if self.min_k is not None and self.min_k <= 0:
            raise ValueError("min_k should be an integer greater than 0."
                             " %s was provided." % str(self.min_k))

        if self.max_k is not None and self.max_k <= 2:
            raise ValueError("max_k should be an integer greater than 2."
                             " %s was provided." % str(self.max_k))

        if self.affinity not in ['precomputed', 'gaussian_kernel', 'gaussian_kernel_topk', 'linear_kernel', 'cosine_similarity', 'knn_neighbors_from_X', 'laplacian_kernel']:
            raise ValueError("affinity should be precomputed, gaussian_kernel, linear_kernel, cosine_similarity, knn_neighbors_from_X, "
                             "laplacian_kernel. "
                             "%s was provided." % str(self.affinity))

        if self.sparsification not in ['affinity', 'precomputed', 'knn_neighbors', 'knn_neighbors_from_X']:
            raise ValueError("sparsification should be affinity, precomputed, knn_neighbors, knn_neighbors_from_X."
                             " %s was provided." % str(self.sparsification))

        if self.strategy not in ['bottom_up', 'top_down']:
            raise ValueError("affinity should be bottom_up, top_down."
                             " %s was provided." % str(self.strategy))

        print('fit', self.strategy)
        self.construct_affinity(X,spatial=spatial,spatial_weight=spatial_weight)
        self.graph_sparsification(X)

        # build the tree

        setree_class = pySETree
        # setree_class = seat_wrapper.SETree

        se_tree = setree_class(aff_m=self.aff_m, knn_m=self.knn_m,
                               objective_paras=self.objective_paras, min_k=self.min_k,
                               max_k=self.max_k, auto_k=self.auto_k,
                               strategy=self.strategy,
                               split_se_cutoff=self.split_se_cutoff, merge_layers=self.merge_layers,
                               plot_cluster_map=self.plot_cluster_map, random_seed=self.random_seed, __verbose__=self.__verbose__)
        self.se_tree = se_tree
        #  1. build tree
        t1 = time.time()
        # se_tree.build_tree_v2()
        self.Z = se_tree.build_tree_v3()
        self.Z_ = se_tree.Z_
        t2 = time.time()
        print('build tree time', t2 - t1)
        self.clubs = se_tree.clubs
        self.aff_m = se_tree.aff_m

        # se_tree.order_tree()
        self.tree_se = se_tree.get_tree_se()

        #  2. contract tree
        if self.__verbose__:
            time1 = time.time()
        se_tree.contract_tree_v2(self.Z, self.ks)
        if self.__verbose__:
            time2 = time.time()
            print('contract tree time', time2 - time1)
        self.vertex_num = se_tree.vertex_num
        self.ks = list(se_tree.ks)
        self.se_scores = se_tree.se_scores
        self.delta_se_scores = se_tree.delta_se_scores
        self.labels_ = se_tree.optimal_clusters
        self.optimal_k = self._amend_optimal_k(
            se_tree.optimal_k)  # amend the value of optimal k
        self.Z_ = self.Z[:, :4]
        self.leaves_list = hierarchy.leaves_list(self.Z_)
        self.order = self._order()
        self.ks_clusters = se_tree.ks_clusters
        self.Z_clusters = se_tree.Z_clusters
        self.clubs = self._get_clubs()

        self.club_k = len(se_tree.leaves)

        self.newick = se_tree.to_newick()

        self.split_se = se_tree.to_split_se()

    def _amend_optimal_k(self, op_k_org):
        """
        if the value of optimal k is larger than one, there should be k clusters,
        instead of only one
        """
        if op_k_org <= 1:
            return op_k_org
        ac_k = len(set(self.labels_))
        if ac_k > 1:
            return op_k_org
        else:
            return 1

    def _order(self):
        # hierarchy.leaves_list(self.Z_)
        order = [(l, i) for i, l in enumerate(self.leaves_list)]
        order.sort()
        return [i for l, i in order]

    def _get_clubs(self):
        leaves = sorted([(self.order[self.se_tree.node_list[l].vs[0]], l) for l in self.se_tree.leaves])
        order = [(v, i) for i, l in enumerate(leaves) for v in self.se_tree.node_list[l[1]].vs]
        order.sort()
        return [i for n, i in order]

    def oval_embedding(self, a=3, b=2, k=0.2):
        # http://www.mathematische-basteleien.de/eggcurves.htm
        angle = np.array([self.order])*(2*np.pi/len(self.order))
        xcor = a*np.cos(angle)
        ycor = b*np.sin(angle)*1/np.sqrt(np.exp(k*np.cos(angle)))
        plane_coordinate = np.concatenate((xcor, ycor), axis=0).T
        return plane_coordinate
