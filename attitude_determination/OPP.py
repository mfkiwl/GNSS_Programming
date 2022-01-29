# -*- coding: utf-8 -*-
"""

@title:	expirement
@author: iDeal0103
@status:	Active
@type:	Process
@created:	25-Jan-2022
@post-History:	25-Jan-2022

comment：
    1.OPP
    2.
    3.

"""

# iport module
import numpy as np
from scipy.spatial.transform.rotation import Rotation as rota
import random

def diagonalize_number_and_squarematrix(a, b):
    # 构造合并后大小的矩阵
    b_row = b.shape[0]
    b_col = b.shape[1]
    c = np.zeros((1 + b_row, 1 + b_col))
    c[:1, :1] = a
    c[1:, 1:] = b
    return c

def get_matrix_from_vectors(list_of_vectors):
    """
    param
    list_of_vectors : list , 向量所组成的列表如:[[x1,y1,z1],[x2,y2,z2]...[xn,yn,zn]]
    return
    vector_matrix : np.darray , 向量所组成的矩阵如:np.array(b1,b2,b3,...bn) bn=[xn,yn,zn]T(其中通过转置变为列向量)
    """
    matrix = np.array(list_of_vectors)
    vector_matrix = matrix.T
    return vector_matrix

def add_perturbation_to_vector(vector, sigmas):
    """
    vector : list , 三维向量 [x, y, z]
    sigmas : list , x、y、z对应的最大最大扰动量 [dx_max, dy_max, dz_max]
    """
    x, y, z = vector.tolist()
    x += 2 * (random.random() - 0.5) * sigmas[0]
    y += 2 * (random.random() - 0.5) * sigmas[1]
    z += 2 * (random.random() - 0.5) * sigmas[2]
    return np.array([x, y, z])


def solve_OPP_withSVD(F, B, W):
    """
    F: np.darray, 由列向量组成的参考坐标系框架下的向量矩阵 [f1, f2, f3....]
    B_hat:  np.darray, 由列向量组成的参考坐标系框架下的向量矩阵 [b1, b2, b3....]
    W: list, 各向量的权重矩阵
    return R_check: np.darray, 旋转矩阵

    note: 本算法与rotation.align_vectors()相同
    """
    the_product = F @ np.diag(W) @ B.T
    if abs(np.linalg.det(the_product)) > 1e-8:
        U, S, VT = np.linalg.svd(the_product)  # 有大坑！注意svd分解得到的结果为V.T
        V = VT.T
        VTU = np.linalg.det(V.T @ U)
        M = diagonalize_number_and_squarematrix(VTU, np.eye(2))
        R_check = V @ M @ U.T
    else:
        print("the product of [F, B_hat, W] is singular!")
        raise SystemExit
    return R_check




if __name__ == "__main__":
    # 在参考框架坐标系下
    f1 = np.array([-2.2775, 2.4688, 5.1465])    # CUCC -> CUBB
    f2 = np.array([3.9066, 1.7174, 0.1474])    # CUT0 -> CUBB
    f3 = np.array([-2.0202, 4.1602, 7.0336])    # CUAA -> CUT0
    f4 = np.array([-6.1844, 0.7520, 4.9993])    # CUCC -> CUT0

    # 在空间指教坐标系下
    r = rota.from_euler('zyx', [5, 7, 10], degrees=True)
    print(r.as_matrix())
    # b1 = r.apply(f1)    # CUCC -> CUBB
    # b2 = r.apply(f2)    # CUT0 -> CUBB
    # b3 = r.apply(f3)    # CUAA -> CUT0
    # b4 = r.apply(f4)    # CUCC -> CUT0
    # b1 = add_perturbation_to_vector(r.apply(f1), [0.1, 0.1, 0.1])    # CUCC -> CUBB
    # b2 = add_perturbation_to_vector(r.apply(f2), [0.1, 0.1, 0.1])   # CUT0 -> CUBB
    # b3 = add_perturbation_to_vector(r.apply(f3), [0.1, 0.1, 0.1])   # CUAA -> CUT0
    # b4 = add_perturbation_to_vector(r.apply(f4), [0.1, 0.1, 0.1])    # CUCC -> CUT0
    # b1 = add_perturbation_to_vector(r.apply(f1), [0.02, 0.02, 0.02])    # CUCC -> CUBB
    # b2 = add_perturbation_to_vector(r.apply(f2), [0.02, 0.02, 0.02])   # CUT0 -> CUBB
    # b3 = add_perturbation_to_vector(r.apply(f3), [0.02, 0.02, 0.02])   # CUAA -> CUT0
    # b4 = add_perturbation_to_vector(r.apply(f4), [0.02, 0.02, 0.02])    # CUCC -> CUT0
    b1 = add_perturbation_to_vector(r.apply(f1), [0.002, 0.002, 0.002])  # CUCC -> CUBB
    b2 = add_perturbation_to_vector(r.apply(f2), [0.002, 0.002, 0.002])  # CUT0 -> CUBB
    b3 = add_perturbation_to_vector(r.apply(f3), [0.002, 0.002, 0.002])  # CUAA -> CUT0
    b4 = add_perturbation_to_vector(r.apply(f4), [0.002, 0.002, 0.002])    # CUCC -> CUT0

    # 组成向量矩阵
    F = get_matrix_from_vectors([f1.tolist(), f2.tolist(), f3.tolist(), f4.tolist()])
    B = get_matrix_from_vectors([b1.tolist(), b2.tolist(), b3.tolist(), b4.tolist()])

    # 解旋转
    # 用rota.align_vectors
    r_check, loss = rota.align_vectors(B.T, F.T, [1, 1, 1, 1])
    # 用自己编写的函数
    r_check_SVD = solve_OPP_withSVD(F, B, [1, 1, 1, 1])
    print(r_check_SVD)
    print(r.as_euler('zyx', degrees=True))
    print(r_check.as_euler('zyx', degrees=True))
    print(rota.from_matrix(r_check_SVD).as_euler('zyx', degrees=True))







