# -*- coding: utf-8 -*-
"""
    src.HE
    ~~~~~~~~~~~

    @Copyright: (c) 2022-07 by Lingxi Chen (chanlingxi@gmail.com).
    @License: LICENSE_NAME, see LICENSE for more details.
"""
import pdb

import numpy as np
from math import log2
from scipy.sparse import csr_matrix


def get_v_sparse(M, sparse_m, vs):
    v = M[vs, :].sum()
    # 对角线修正（如果需要）
    # v -= M_sparse_csr[vs, vs].sum()  # 更高效的实现
    return v if v != 0 else 1


def get_v_v0_1(M, sparse_m, vs):
    # intra- and inter- affinity
    # return np.sum(M[vs, :]) - np.sum([M[i, i] for i in vs])
    # test for reduce memory cost 2023/12/25
    # v = np.sum(M[vs, :])
    _max_read_ = 8000
    if len(vs) < _max_read_:
        v = M[vs, :].sum()
    else:
        _max_read_ = 8000
        len_vs = len(vs)
        times = (len_vs + _max_read_ - 1) // _max_read_  # 更高效地计算分块次数
        v = 0
        for t in range(times):
            start = t * _max_read_
            end = start + _max_read_
            vs_tmp = vs[start:end]
            v += np.sum(M[vs_tmp, :])
    if v == 0:
        v += 1
    return v


def get_v(row_sums, sparse_m, vs):
    v = row_sums[vs].sum()
    return v if v != 0 else 1


def get_v_4_large(row_sums, sparse_m, vs):
    _max_read_ = 8000
    len_vs = len(vs)
    if len_vs < _max_read_:
        v = row_sums[vs].sum()
    else:
        _max_read_ = 8000
        times = (len_vs + _max_read_ - 1) // _max_read_  # 更高效地计算分块次数
        v = 0
        for t in range(times):
            start = t * _max_read_
            end = start + _max_read_
            vs_tmp = vs[start:end]
            v += row_sums[vs_tmp].sum()
    return v if v != 0 else 1


def get_v_v0(M, sparse_m, vs):
    # intra- and inter- affinity
    # return np.sum(M[vs, :]) - np.sum([M[i, i] for i in vs])
    # test for reduce memory cost 2023/12/25
    # v = np.sum(M[vs, :])
    _max_read_ = 8000
    if len(vs) < _max_read_:
        v = np.sum(M[vs, :])
    else:
        len_vs = len(vs)
        times = int(len_vs/_max_read_) + 1
        v = 0
        for t in range(times):
            pos_end = min((t+1)*_max_read_,len_vs)
            vs_tmp = vs[t*_max_read_:pos_end]
            v += np.sum(M[vs_tmp, :])
    if v == 0:
        v += 1
    return v
    # return sparse_m[vs, :].sum()  # slow


def get_s_sparse(M, sparse_m, vs):
    s = M[np.ix_(vs, vs)].sum()
    return s


def get_s(M, sparse_m, vs):
    # intra-affinity
    _max_read_ = 8000
    len_vs = len(vs)
    if len_vs < _max_read_:
        s = M[np.ix_(vs, vs)].sum()
    else:
        s = 0
        times = (len_vs + _max_read_ - 1) // _max_read_  # 计算块数
        for t in range(times):
            start = t * _max_read_
            end = min((t + 1) * _max_read_, len_vs)
            vs_tmp = vs[start:end]
            s += M[np.ix_(vs_tmp, vs)].sum()
    return s


def get_s_v0(M, sparse_m, vs):
    # intra-affinity
    _max_read_ = 8000
    if len(vs) < _max_read_:
        s = np.sum(M[np.ix_(vs, vs)])
    else:
        len_vs = len(vs)
        times = int(len_vs / _max_read_) + 1
        s = 0
        for t in range(times):
            pos_end = min((t+1)*_max_read_,len_vs)
            vs_tmp = vs[t*_max_read_:pos_end]
            s += np.sum(M[np.ix_(vs_tmp,vs)])
    return s
    # return np.sum(M[np.ix_(vs, vs)])  # test for reduce memory cost
    ''' slow
    c = 0
    for i in vs:
        for j in vs:
            c += M[i, j]
    return c
    '''


def get_g(M, sparse_m, vs, M_rows):
    # inter-affinity
    return get_v(M_rows, sparse_m, vs) - get_s(M, sparse_m, vs)


# structure entropy
def get_node_se(vG, g, V, pV):
    return float(g) / float(vG) * log2(float(pV) / float(V))


def get_node_vertices_se(vG, d_log_d, V):
    return - (d_log_d - V * log2(V)) / vG

# structure entropy
def get_node_se_v2(vG, g, V, pV):
    # pV = vG
    return - float(g) / float(vG) * log2(float(pV) / float(V))


def get_node_vertices_se_v2(vG, d_log_d, V):
    return (d_log_d - V * log2(V)) / vG


# modularity
def get_node_nm(vG, g, V, eta, eta_mode):
    if eta_mode is "coefficient":
        return (float(V - g) / float(vG) - eta * pow((float(V) / float(vG)),
                                                     2))  # take negative for minimization
    elif eta_mode is "exponent":
        return (float(V - g) / float(vG) - pow((float(V) / float(vG)),
                                               2 * eta))  # take negative for minimization


def get_node_vertices_nm(vG, V, d_square, eta, eta_mode):
    if eta_mode is "coefficient":
        if d_square ==0 or float(vG==0) or pow(float(vG),2)==0:
            return 0
        return -(eta * float(d_square) / pow(float(vG),
                                              2))  # take negative for minimization
    elif eta_mode is "exponent":
        return -(pow((float(d_square) / float(vG ** 2)),
                      eta))  # take negative for minimization

# KL divergence
def get_node_kl_v2(vG, g, V, pV, eta, eta1, eta2):
    """
    use pV as denominator instead of vG in v1
    """
    s = V - g
    if s == 0:
        return 0
    return np.power((float(V) / float(pV)), eta2) * (
                eta1 * log2(float(s) / float(pV)) - eta * log2(float(V) / float(pV)))

# KL divergence
def get_node_kl_v3(vG, g, V, pV, eta, eta1, eta2, hype_a=0):
    """
    subtract a hype-parameter "a" to the equation
    """
    s = V - g
    if s == 0:
        return 0
    pV = vG
    return np.power((float(s) / float(vG)), eta2) * (
                eta1 * log2(float(s) / float(pV)) - hype_a - eta * log2(float(V) / float(pV)))


def get_node_vertices_kl():
    return 0


def get_delta_merge_se(M, sparse_m, vG, d, parent, node1, node2):
    # merge node1 and node2 into a new one
    # leaf merge
    new_V = node1.V + node2.V
    new_s = node1.s + node2.s
    new_s_tmp = 0
    for b1 in node1.vs:
        for b2 in node2.vs:
            new_s_tmp += M[b1, b2]
    # new_s_tmp = np.sum(M[np.ix_(node1.vs, node2.vs)])   # slow
    new_s += 2 * new_s_tmp
    new_se = (new_V * log2(new_V) - node1.d_log_d - node2.d_log_d) / vG  # vertex
    if parent is None:
        new_se += -(new_V - new_s) / vG * log2(new_V / vG)  # node
    else:
        new_se += -(new_V - new_s) / vG * log2(new_V / parent.V)  # node

    old_se = node1.se + node2.se \
             + (
                         node1.V * node1.log_V - node1.d_log_d + node2.V * node2.log_V - node2.d_log_d) / vG
    delta = old_se - new_se
    return delta  # minimise se

def get_delta_merge_se_v2(M, sparse_m, vG, d, parent, node1, node2):
    # merge node1 and node2 into a new one
    # leaf merge
    new_V = node1.V + node2.V
    new_s = node1.s + node2.s
    new_s_tmp = 0
    for b1 in node1.vs:
        for b2 in node2.vs:
            new_s_tmp += M[b1, b2]
    # new_s_tmp = np.sum(M[np.ix_(node1.vs, node2.vs)])   # slow
    new_s += 2 * new_s_tmp
    new_se = (new_V * log2(new_V) - node1.d_log_d - node2.d_log_d) / vG  # vertex
    if parent is None:
        new_se += -(new_V - new_s) / vG * log2(new_V / vG)  # node
    else:
        new_se += -(new_V - new_s) / vG * log2(new_V / parent.V)  # node

    old_se = node1.se + node2.se \
             + (
                         node1.V * node1.log_V - node1.d_log_d + node2.V * node2.log_V - node2.d_log_d) / vG
    delta = old_se - new_se
    return - delta  # minimise se

def get_delta_merge_se_v3(M, sparse_m, vG, d, parent, node1, node2, node0):
    """
    remove child node in node1 and insert into node2
    """
    node1_p_V = node1.V - node0.V
    node2_p_V = node2.V + node0.V
    if node1_p_V == 0:
        node1_p_s = 0
    else:
        s_tmp = 0
        tmp_vs = list(set(node1.vs) - set(node0.vs))
        for b1 in tmp_vs:
            for b2 in node0.vs:
                s_tmp += M[b1, b2]
        node1_p_s = node1.s - 2 * s_tmp - node0.s
    node2_p_s = node0.s + node2.s
    new_s_tmp = 0
    for b1 in node0.vs:
        for b2 in node2.vs:
            new_s_tmp += M[b1, b2]
    # new_s_tmp = np.sum(M[np.ix_(node1.vs, node2.vs)])   # slow
    node2_p_s += 2 * new_s_tmp
    if node1_p_V <= 0:
        node1_p_vlogv = 0
        node1_p_se = 0
    else:
        if node1_p_s<=0:
            node1_p_s = 0
        node1_p_vlogv = node1_p_V * log2(node1_p_V)
        node1_p_se = -(node1_p_V - node1_p_s) / vG * log2(node1_p_V / vG)
    new_se = (node1_p_vlogv + node2_p_V * log2(
        node2_p_V) - node1.d_log_d - node2.d_log_d) / vG  # vertex
    parent = None
    if parent is None:
        new_se += -(node2_p_V - node2_p_s) / vG * log2(
            node2_p_V / vG) + node1_p_se  # node
    else:
        new_se += -(node2_p_V - node2_p_s) / vG * log2(
            node2_p_V / node2.pV) + node1_p_se  # node
    old_se = -node1.se + -node2.se \
             + (
                     node1.V * node1.log_V - node1.d_log_d + node2.V * node2.log_V - node2.d_log_d) / vG
    delta = old_se - new_se
    return delta  # minimise se


# modularity
def get_delta_merge_nm(M, sparse_m, vG, parent, node1, node2, eta_mode="coefficient",
                       eta=1):
    """
    eta_mode: the position of eta, possible flags: None, exponent
    eta: the value of eta if eta mode is not None, eta>=0
    """
    if eta_mode is "coefficient":  # coefficient, not exponent
        if eta < 0:
            print("ERROR! Wrong eta input.")
        if len(node1.vs) == 1:
            node1_q = - eta * np.power(node1.V / vG, 2)
            # node1_q = node1.g/vG - r*np.power(node1.V/vG, 2)
        else:
            node1_q = node1.s / vG - eta * np.power(node1.V / vG, 2)
        if len(node2.vs) == 1:
            node2_q = - eta * np.power(node2.V / vG, 2)
            # node2_q = node2.g/vG - r*np.power(node2.V/vG, 2)
        else:
            node2_q = node2.s / vG - eta * np.power(node2.V / vG, 2)
        new_node_V = node1.V + node2.V
        new_node_s = node1.s + node2.s
        new_s_tmp = 0
        for b1 in node1.vs:
            for b2 in node2.vs:
                new_s_tmp += M[b1, b2]
        new_node_s += 2 * new_s_tmp
        # new_node_g = new_node_V - new_node_s

        after = new_node_s / vG - eta * np.power(new_node_V / vG, 2)
    else:
        if eta < 0:
            print("ERROR! Wrong eta input.")
        else:
            if len(node1.vs) == 1:  # only one vertex
                node1_q = - np.power(node1.V / vG, 2 * eta)
                # node1_q = node1.g/vG - r*np.power(node1.V/vG, 2)
            else:
                node1_q = node1.s / vG - np.power(node1.V / vG, 2 * eta)
            if len(node2.vs) == 1:
                node2_q = - np.power(node2.V / vG, 2 * eta)
                # node2_q = node2.g/vG - r*np.power(node2.V/vG, 2)
            else:
                node2_q = node2.s / vG - np.power(node2.V / vG, 2 * eta)
            new_node_V = node1.V + node2.V
            new_node_s = node1.s + node2.s
            new_s_tmp = 0
            for b1 in node1.vs:
                for b2 in node2.vs:
                    new_s_tmp += M[b1, b2]
            new_node_s += 2 * new_s_tmp
            # new_node_g = new_node_V - new_node_s

            after = new_node_s / vG - np.power(new_node_V / vG, 2 * eta)
    return after - (node1_q + node2_q)  # maximize modularity, take negative value


def get_delta_merge_nm_v2(M, sparse_m, vG, parent, node1, node2, node0,
                          eta_mode="coefficient", eta=1):
    """
    eta_mode: the position of eta, possible flags: None, exponent
    eta: the value of eta if eta mode is not None, eta>=0
    """
    node1_p_V = node1.V - node0.V
    node2_p_V = node2.V + node0.V
    if node1_p_V == 0:
        node1_p_s = 0
    else:
        s_tmp = 0
        tmp_vs = list(set(node1.vs) - set(node0.vs))
        for b1 in tmp_vs:
            for b2 in node0.vs:
                s_tmp += M[b1, b2]
        node1_p_s = node1.s - 2 * s_tmp - node0.s
    node2_p_s = node0.s + node2.s
    new_s_tmp = 0
    for b1 in node0.vs:
        for b2 in node2.vs:
            new_s_tmp += M[b1, b2]
    # new_s_tmp = np.sum(M[np.ix_(node1.vs, node2.vs)])   # slow
    node2_p_s += 2 * new_s_tmp
    if eta_mode is "coefficient":  # coefficient, not exponent
        if eta < 0:
            print("ERROR! Wrong eta input.")
        if len(node1.vs) == 1:
            node1_q = - eta * np.power(node1.V / vG, 2)
            # node1_q = node1.g/vG - r*np.power(node1.V/vG, 2)
        else:
            node1_q = node1.s / vG - eta * np.power(node1.V / vG, 2)
        if len(node2.vs) == 1:
            node2_q = - eta * np.power(node2.V / vG, 2)
            # node2_q = node2.g/vG - r*np.power(node2.V/vG, 2)
        else:
            node2_q = node2.s / vG - eta * np.power(node2.V / vG, 2)
        # new_node_g = new_node_V - new_node_s

        after = node1_p_s / vG - eta * np.power(node1_p_V / vG,
                                                2) + node2_p_s / vG - eta * np.power(
            node2_p_V / vG, 2)
    else:
        if eta < 0:
            print("ERROR! Wrong eta input.")
        else:
            if len(node1.vs) == 1:  # only one vertex
                node1_q = - np.power(node1.V / vG, 2 * eta)
                # node1_q = node1.g/vG - r*np.power(node1.V/vG, 2)
            else:
                node1_q = node1.s / vG - np.power(node1.V / vG, 2 * eta)
            if len(node2.vs) == 1:
                node2_q = - np.power(node2.V / vG, 2 * eta)
                # node2_q = node2.g/vG - r*np.power(node2.V/vG, 2)
            else:
                node2_q = node2.s / vG - np.power(node2.V / vG, 2 * eta)
            # new_node_g = new_node_V - new_node_s

            after = node1_p_s / vG - np.power(node1_p_V / vG,
                                              2 * eta) + node2_p_s / vG - np.power(
                node2_p_V / vG, 2 * eta)
    return after - (node1_q + node2_q)  # maximize modularity, take negative value



def compute_kl(s, V, vG, parent_V, eta, eta1, eta2, hype_a, is_root):
    """
    计算合并后的节点kl值的辅助函数
    """
    if s <= 0 or V <= 0:
        return 0.0
    if is_root:
        ratio_s = s / vG
        ratio_V = V / vG
    else:
        if parent_V <= 0:
            return 0.0
        ratio_s = s / parent_V
        ratio_V = V / parent_V
    log_s = np.log2(ratio_s) if ratio_s > 0 else 0.0
    log_V = np.log2(ratio_V) if ratio_V > 0 else 0.0
    term = eta1 * log_s - hype_a - eta * log_V
    se = (ratio_s ** eta2) * term if not is_root else (s / vG) ** eta2 * term
    return se if se > 0 else 0.0

def get_delta_merge_kl_v5_v0(M, vG, node1, node2, node0, parent_status=False,
                                   parent_node1=None, parent_node2=None, eta=0.5, eta1=1,
                                   eta2=1, hype_a=0):
    """
    tried to optimize
    2025/03/07
    """
    # 转换为NumPy数组（假设node.vs原本为列表，此处需在数据结构中预先处理）
    node1_vs = np.array(node1.vs)
    node0_vs = np.array(node0.vs)
    node2_vs = np.array(node2.vs)

    # 计算合并后的V值
    node1_p_V = node1.V - node0.V
    node2_p_V = node2.V + node0.V
    vG_float = float(vG)

    # 计算node1_p_s
    tmp_vs = np.setdiff1d(node1_vs, node0_vs, assume_unique=True)
    if tmp_vs.size == 0 or node0_vs.size == 0:
        s_tmp = 0
    else:
        s_tmp = M[np.ix_(tmp_vs, node0_vs)].sum()
    node1_p_s = node1.s - 2 * s_tmp - node0.s

    # 计算node2_p_s
    if node0_vs.size == 0 or node2_vs.size == 0:
        new_s_tmp = 0
    else:
        new_s_tmp = M[np.ix_(node0_vs, node2_vs)].sum()
    node2_p_s = node0.s + node2.s + 2 * new_s_tmp

    # 计算SE值
    if parent_status:
        parent_node1_V = parent_node1.V - node0.V
        parent_node2_V = parent_node2.V + node0.V
        node1_p_se = compute_kl(node1_p_s, node1_p_V, vG_float, parent_node1_V, eta, eta1, eta2, hype_a, False)
        node2_p_se = compute_kl(node2_p_s, node2_p_V, vG_float, parent_node2_V, eta, eta1, eta2, hype_a, False)
    else:
        node1_p_se = compute_kl(node1_p_s, node1_p_V, vG_float, None, eta, eta1, eta2, hype_a, True)
        node2_p_se = compute_kl(node2_p_s, node2_p_V, vG_float, None, eta, eta1, eta2, hype_a, True)

    delta = (node1_p_se + node2_p_se) - (node1.se + node2.se)
    return delta




def get_delta_merge_kl_v5(M, vG, node1, node2, node0, parent_status=False,
                          parent_node1=None, parent_node2=None, eta=0.5, eta1=1,
                          eta2=1, hype_a=0):
    """
    remove child node in node1 and insert into node2
    derived from get_delta_merge_kl_v3: add - log {\eta} to the equation
    """
    node1_p_V = node1.V - node0.V
    node2_p_V = node2.V + node0.V
    if node1_p_V == 0:
        node1_p_s = 0
    else:
        # s_tmp = 0
        tmp_vs = list(set(node1.vs) - set(node0.vs))
        # for b1 in tmp_vs:
        #     for b2 in node0.vs:
        #         s_tmp += M[b1, b2]
        s_tmp = M[np.ix_(tmp_vs, node0.vs)].sum()  # slow
        node1_p_s = node1.s - 2 * s_tmp - node0.s
    node2_p_s = node0.s + node2.s
    # new_s_tmp = 0
    # for b1 in node0.vs:
    #     for b2 in node2.vs:
    #         new_s_tmp += M[b1, b2]
    new_s_tmp = M[np.ix_(node0.vs, node2.vs)].sum()   # slow
    node2_p_s += 2 * new_s_tmp
    #     tmp_vs = list(set(node1.vs) - set(node0.vs))
    #     node1_p_s = np.sum(M[np.ix_(tmp_vs, tmp_vs)])
    # tmp_vs = list(set(node2.vs).union(set(node0.vs)))
    # node2_p_s = np.sum(M[np.ix_(tmp_vs, tmp_vs)])
    if parent_status is False:  # root
        if node1_p_s <= 0:
            node1_p_se = 0
        else:
            if node1_p_V <= 0:
                logV = 0
            else:
                logV = log2(float(node1_p_V) / vG)

            node1_p_se = np.power(float(node1_p_s) / vG, eta2) * (eta1 * log2(
                float(node1_p_s) / vG) - hype_a - eta * logV)  # the parent node is root
        if node2_p_s == 0:
            node2_p_se = 0
        else:
            if node2_p_V == 0:
                logV = 0
            else:
                logV = log2(float(node2_p_V) / vG)
            node2_p_se = np.power(float(node2_p_s) / vG, eta2) * (eta1 * log2(
                float(node2_p_s) / vG) - hype_a - eta * logV)  # the parent node is root
    else:
        if node1_p_s == 0:
            node1_p_se = 0
        else:
            if node1_p_V == 0:
                node1_p_se = 0
            else:
                parent_node1_V = parent_node1.V - node0.V
                node1_p_se = np.power(float(node1_p_s) / vG, eta2) * (
                        eta1 * log2(float(node1_p_s) / parent_node1_V) - hype_a - eta * log2(
                    float(node1_p_V) / parent_node1_V))  # the parent node is not a root
        if node2_p_s == 0:
            node2_p_se = 0
        else:
            if node2_p_V == 0:
                node2_p_se = 0
            else:
                parent_node2_V = parent_node2.V + node0.V
                node2_p_se = np.power(float(node2_p_s) / vG, eta2) * (
                            eta1 * log2(float(node2_p_s) / parent_node2_V) - hype_a -eta * log2(
                        float(node2_p_V) / parent_node2_V))  # the parent node is root
        # if node1_p_s * node2_p_s != 0:
        #     print("*********prev: node1: ", node1.V," parent node 1: ", parent_node1.V, " parent_node_2: ",
        #           parent_node2.V, " node0: ", node0.V)
        #     print("*********post: node1: ", node1_p_V," parent node 1: ", parent_node1_V, " parent_node_2: ",
        #           parent_node2_V)
    after = node1_p_se + node2_p_se
    delta = after - (node1.se + node2.se)
    return delta


def get_delta_combine_se(M, sparse_m, vG, parent, node1, node2):
    # add a parent node of node1 and node2
    new_node_V = node1.V + node2.V
    new_node_s = node1.s + node2.s
    new_s_tmp = 0
    for b1 in node1.vs:
        for b2 in node2.vs:
            new_s_tmp += M[b1, b2]
    # new_s_tmp = np.sum(M[np.ix_(node1.vs, node2.vs)])  # slow
    new_node_s += 2 * new_s_tmp
    new_node_g = new_node_V - new_node_s
    new_node_se = get_node_se(vG, new_node_g, new_node_V, parent.V)
    new_node1_se = get_node_se(vG, node1.g, node1.V, new_node_V)
    new_node2_se = get_node_se(vG, node2.g, node2.V, new_node_V)
    return node1.se + node2.se - new_node_se - new_node1_se - new_node2_se


def get_delta_combine_se_v2(M, sparse_m, vG, parent, node1, node2):
    # add a parent node of node1 and node2
    new_node_V = node1.V + node2.V
    new_node_s = node1.s + node2.s
    new_s_tmp = 0
    for b1 in node1.vs:
        for b2 in node2.vs:
            new_s_tmp += M[b1, b2]
    # new_s_tmp = np.sum(M[np.ix_(node1.vs, node2.vs)])  # slow
    new_node_s += 2 * new_s_tmp
    new_node_g = new_node_V - new_node_s
    new_node_se = get_node_se_v2(vG, new_node_g, new_node_V, parent.V)
    new_node1_se = get_node_se_v2(vG, node1.g, node1.V, new_node_V)
    new_node2_se = get_node_se_v2(vG, node2.g, node2.V, new_node_V)
    return (node1.se + node2.se - new_node_se - new_node1_se - new_node2_se)


# KL divergence
def get_delta_combine_kl_v2(M, sparse_m, vG, parent, node1, node2, objective_paras):
    objective, eta_mode, eta, eta1, eta2, hype_a = objective_paras._get_vars_()
    # add a parent node of node1 and node2
    new_node_V = node1.V + node2.V
    new_node_s = node1.s + node2.s
    new_s_tmp = 0
    # for b1 in node1.vs:
    #     for b2 in node2.vs:
    #         new_s_tmp += M[b1, b2]
    new_s_tmp = np.sum(M[np.ix_(node1.vs, node2.vs)])  # slow
    new_node_s += 2 * new_s_tmp
    new_node_g = new_node_V - new_node_s
    new_node_kl = get_node_kl_v3(vG=vG, g=new_node_g, V=new_node_V, pV=parent.V, eta=eta,
                              eta1=eta1, eta2=eta2, hype_a=hype_a)  # (vG, g, V, pV, eta)
    new_node1_kl = get_node_kl_v3(vG=vG, g=node1.g, V=node1.V, pV=new_node_V, eta=eta,
                               eta1=eta1, eta2=eta2, hype_a=hype_a)
    new_node2_kl = get_node_kl_v3(vG=vG, g=node2.g, V=node2.V, pV=new_node_V, eta=eta,
                               eta1=eta1, eta2=eta2, hype_a=hype_a)
    return (node1.se + node2.se - new_node_kl - new_node1_kl - new_node2_kl)


# topology entropy
def get_node_te(vG, g, V, pV, eta=1):
    vG, g, V, pV = float(vG), float(g), float(V), float(pV)
    s = V - g
    if s == 0:
        return 0
    diff = log2(s / pV) - 2 * eta * log2(V / pV)
    return -s / vG * diff


def get_node_vertices_te(vG, d_log_d, V):
    return 0


# topology entropy louvain like
def get_node_lte(vG, g, V, pV, eta=1):
    vG, g, V, pV = float(vG), float(g), float(V), float(pV)
    s = V - g
    if s == 0:
        return 0
    diff = log2(s / pV) - log2(eta) - 2 * log2(V / pV)
    return -s / vG * diff


def get_node_vertices_lte(vG, d_log_d, V):
    return 0


def get_node_score_v2(vG, g, V, pV, objective_paras):
    objective, eta_mode, eta, eta1, eta2, hype_a = objective_paras._get_vars_()
    if objective == 'SE':
        return get_node_se_v2(vG, g, V, pV)
    elif objective == 'TE':
        return get_node_te(vG, g, V, pV, eta)
    elif objective == 'LTE':
        return get_node_lte(vG, g, V, pV, eta)
    elif objective == 'M':
        return get_node_nm(vG, g, V, eta, eta_mode)
    elif objective == 'KL':
        return get_node_kl_v3(vG, g, V, pV, eta, eta1, eta2, hype_a)


def get_node_vertices_score(vG, d_log_d, V, d_square, objective_paras):
    objective, eta_mode, eta, eta1, eta2, hype_a = objective_paras._get_vars_()
    if objective == 'SE':
        return get_node_vertices_se_v2(vG, d_log_d, V)
    elif objective == 'TE':
        return get_node_vertices_te(vG, d_log_d, V)
    elif objective == 'LTE':
        return get_node_vertices_lte(vG, d_log_d, V)
    elif objective == 'M':
        return get_node_vertices_nm(vG, V, d_square, eta, eta_mode)
    elif objective == 'KL':
        return get_node_vertices_kl()


def get_delta_merge_score_v3(M, sparse_m, vG, d, node1, node2, node0, objective_paras,
                             parent=None, parent_status=False, parent_node1=None,
                             parent_node2=None):
    objective, eta_mode, eta, eta1, eta2, hype_a = objective_paras._get_vars_()
    if objective == 'SE':
        return get_delta_merge_se_v3(M, sparse_m, vG, d, parent, node1, node2,
                                     node0=node0)
    elif objective == 'M':
        return get_delta_merge_nm_v2(M, sparse_m, vG, parent, node1, node2, node0,
                                     eta_mode, eta)
    elif objective == 'KL':
        return get_delta_merge_kl_v5(M=M, vG=vG, node1=node1, node2=node2, node0=node0,
                                     parent_status=parent_status,
                                     parent_node1=parent_node1,
                                     parent_node2=parent_node2, eta=eta, eta1=eta1,
                                     eta2=eta2, hype_a=hype_a)
    else:
        print("ERROR! No appropriate objective.")
        return


def get_delta_combine_nm(M, sparse_m, vG, parent, node1, node2, objective_paras):
    objective, eta_mode, eta, eta1, eta2, hype_a = objective_paras._get_vars_()
    # add a parent node of node1 and node2
    new_node_V = node1.V + node2.V
    new_node_s = node1.s + node2.s
    # new_s_tmp = 0
    # for b1 in node1.vs:
    #     for b2 in node2.vs:
    #         new_s_tmp += M[b1, b2]
    new_s_tmp = np.sum(M[np.ix_(node1.vs, node2.vs)])  # slow
    new_node_s += 2 * new_s_tmp
    new_node_g = new_node_V - new_node_s
    new_node_nm = get_node_nm(vG, new_node_g, new_node_V, eta, eta_mode)
    new_node1_nm = get_node_nm(vG, node1.g, node1.V, eta, eta_mode)
    new_node2_nm = get_node_nm(vG, node2.g, node2.V, eta, eta_mode)
    return -(node1.se + node2.se - new_node_nm - new_node1_nm - new_node2_nm)
