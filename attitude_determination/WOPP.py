# -*- coding: utf-8 -*-
"""

@title:	expirement
@author: iDeal0103
@status:	Active
@type:	Process
@created:	27-Jan-2022
@post-History:	27-Jan-2022

comment：
    1.WOPP
    2.
    3.

"""

# iport module
import numpy as np
from scipy.spatial.transform.rotation import Rotation as rota
import random
from attitude_determination.OPP import *
import sympy


def diagonalize_squarematrix(a, b):
    # 构造合并后大小的矩阵
    a_row = a.shape[0]
    a_col = a.shape[1]
    b_row = b.shape[0]
    b_col = b.shape[1]
    c = np.zeros((a_row + b_row, a_col + b_col))
    c[:a_row, :a_col] = a
    c[a_row:, a_col:] = b
    return c

def solve_WOPP_withLagrangianMultipliers(B, F, QBB):
    """
    [27.33] GNSS Carrier Phase-based Attitude Determination  p51
    A solution based on the Lagrangian multipliers method
    --------
    B: np.darray, 由列向量组成的参考坐标系框架下的向量矩阵 [b1, b2, b3....]
    F: np.darray, 由列向量组成的参考坐标系框架下的向量矩阵 [f1, f2, f3....]
    QBB: 基线解算的vc矩阵
    """
    def get_lamb_matrix(lamb_list):
        lamb1, lamb2, lamb3, lamb4, lamb5, lamb6 = lamb_list
        lamb_matrix = np.array([[lamb1, lamb4/2, lamb5/2], [lamb4/2, lamb2, lamb6/2], [lamb5/2, lamb6/2, lamb3]])
        return lamb_matrix

    def vec_joint(vec1, vec2):
        new_vec = np.array(vec1.tolist() + vec2.tolist())
        return new_vec

    def get_jocabian_matrix(Qbb, R, lamb_matrix):
        # 获得R中各元素
        x1 = R[0, 0]; x2 = R[1, 0]; x3 = R[2, 0]
        x4 = R[0, 1]; x5 = R[1, 1]; x6 = R[2, 1]
        x7 = R[0, 2]; x8 = R[1, 2]; x9 = R[2, 2]
        # 构造雅可比矩阵的各分块矩阵
        part11 = np.linalg.inv(Q) - np.kron(lamb_matrix, np.eye(3))
        part12 = np.array([
            [-x1, 0, 0, -0.5*x4, -0.5*x7, 0],
            [-x2, 0, 0, -0.5*x5, -0.5*x8, 0],
            [-x3, 0, 0, -0.5*x6, -0.5*x9, 0],
            [0, -x4, 0, -0.5*x1, 0, -0.5*x7],
            [0, -x5, 0, -0.5*x2, 0, -0.5*x8],
            [0, -x6, 0, -0.5*x3, 0, -0.5*x9],
            [0, 0, -x7, 0, -0.5*x1, -0.5*x4],
            [0, 0, -x8, 0, -0.5*x2, -0.5*x5],
            [0, 0, -x9, 0, -0.5*x3, -0.5*x6]
        ])
        part21 = np.array([
            [2*x1, 2*x2, 2*x3, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 2*x4, 2*x5, 2*x6, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 2*x7, 2*x8, 2*x9],
            [x4, x5, x6, x1, x2, x3, 0, 0, 0],
            [x7, x8, x9, 0, 0, 0, x1, x2, x3],
            [0, 0, 0, x7, x8, x9, x4, x5, x6]
        ])
        part22 = np.zeros((6, 6))
        # 构造雅可比矩阵
        jacobian_matrix = np.block([[part11, part12], [part21, part22]])
        return jacobian_matrix

    def from_y_get_Rlamb(y):
        R = np.array([y[0:3], y[3:6], y[6:9]]).T
        lamb_list= y[9:].tolist()
        return R, lamb_list

    def from_Rlamb_get_y(R, lamb_list):
        y = []
        y += np.ravel(R, order='F').tolist()
        y += lamb_list
        y = np.array(y)
        return y

    def calculate_QvecR_check(F, QBB):
        F_matrix_part = np.linalg.inv(F @ F.T) @ F
        F_matrix = diagonalize_squarematrix(diagonalize_squarematrix(F_matrix_part, F_matrix_part),F_matrix_part)
        Qvec_B = QBB
        for i in range(B.shape[1]-1):
            Qvec_B = diagonalize_squarematrix(Qvec_B, QBB)
        Qr_check = F_matrix @ Qvec_B @ F_matrix.T
        return Qr_check


    """
    以下为计算过程
    """
    # 根据最小二乘计算R_check
    R_check = B @ F.T @ np.linalg.inv(F @ F.T)
    Q = calculate_QvecR_check(F, QBB)
    # Q = np.eye(9)
    # X = np.random.rand(9**2).reshape(9, 9)
    # X = np.triu(X)
    # X += X.T - np.diag(X.diagonal())
    # Q = X

    # 开始迭代
    lamb_list = [1, 1, 1, 1, 1, 1]
    R = R_check
    while True:
        lamb_matrix = get_lamb_matrix(lamb_list)
        g1 = (np.linalg.inv(Q) - np.kron(lamb_matrix, np.eye(3))) @ np.ravel(R, order='F') - np.linalg.inv(Q) @ np.ravel(R_check, order='F')
        g2_R = R.T @ R - np.eye(3)
        g2 = np.array([g2_R[0, 0], g2_R[1, 1], g2_R[2, 2], g2_R[1, 0], g2_R[2, 0], g2_R[2, 1]])
        g = vec_joint(g1, g2)
        y = from_Rlamb_get_y(R, lamb_list)
        J = get_jocabian_matrix(Q, R, lamb_matrix)
        dy = -np.linalg.inv(J) @ g
        # print(dy)
        y += dy
        R, lamb_list = from_y_get_Rlamb(y)
        if np.all(abs(dy) < 1e-8):
            # 计算vc阵
            Qrr = J[:9, :9] @ Q @ J[:9, :9].T
            break
    return R, Qrr



if __name__ == "__main__":
    # 在参考框架坐标系下
    f1 = np.array([-2.2775, 2.4688, 5.1465])  # CUCC -> CUBB
    f2 = np.array([3.9066, 1.7174, 0.1474])  # CUT0 -> CUBB
    f3 = np.array([-2.0202, 4.1602, 7.0336])  # CUAA -> CUT0
    f4 = np.array([-6.1844, 0.7520, 4.9993])  # CUCC -> CUT0

    # 在空间指教坐标系下
    r = rota.from_euler('zyx', [5, 7, 10], degrees=True)
    print(r.as_matrix())
    # b1 = r.apply(f1)    # CUCC -> CUBB
    # b2 = r.apply(f2)    # CUT0 -> CUBB
    # b3 = r.apply(f3)    # CUAA -> CUT0
    # b4 = r.apply(f4)    # CUCC -> CUT0
    # b1 = add_perturbation_to_vector(r.apply(f1), [0.05, 0.05, 0.05])  # CUCC -> CUBB
    # b2 = add_perturbation_to_vector(r.apply(f2), [0.05, 0.05, 0.05])  # CUT0 -> CUBB
    # b3 = add_perturbation_to_vector(r.apply(f3), [0.05, 0.05, 0.05])  # CUAA -> CUT0
    # b4 = add_perturbation_to_vector(r.apply(f4), [0.05, 0.05, 0.05])  # CUCC -> CUT0
    b1 = add_perturbation_to_vector(r.apply(f1), [0.002, 0.002, 0.002])  # CUCC -> CUBB
    b2 = add_perturbation_to_vector(r.apply(f2), [0.002, 0.002, 0.002])  # CUT0 -> CUBB
    b3 = add_perturbation_to_vector(r.apply(f3), [0.002, 0.002, 0.002])  # CUAA -> CUT0
    b4 = add_perturbation_to_vector(r.apply(f4), [0.002, 0.002, 0.002])  # CUCC -> CUT0


    # 组成向量矩阵
    F = get_matrix_from_vectors([f1.tolist(), f2.tolist(), f3.tolist(), f4.tolist()])
    B = get_matrix_from_vectors([b1.tolist(), b2.tolist(), b3.tolist(), b4.tolist()])

    # 解算WOPP问题
    Qbb = np.array([[0.02, 0, 0], [0, 0.02, 0], [0, 0, 0.02]])
    R, Qrr= solve_WOPP_withLagrangianMultipliers(B, F, Qbb)
    print(R)
    print(Qrr)

    print(r.as_euler('zyx', degrees=True))
    print(rota.from_matrix(R).as_euler('zyx', degrees=True))









