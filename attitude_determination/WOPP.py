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
    QBB: vec(B)的vc矩阵
    """
    def get_lamb_matrix_p3(lamb_list):
        lamb1, lamb2, lamb3, lamb4, lamb5, lamb6 = lamb_list
        lamb_matrix = np.array([[lamb1, lamb4/2, lamb5/2], [lamb4/2, lamb2, lamb6/2], [lamb5/2, lamb6/2, lamb3]])
        return lamb_matrix

    def get_lamb_matrix_p2(lamb_list):
        lamb1, lamb2, lamb3 = lamb_list
        lamb_matrix = np.array([[lamb1, lamb3 / 2], [lamb3 / 2, lamb2]])
        return lamb_matrix

    def vec_joint(vec1, vec2):
        new_vec = np.array(vec1.tolist() + vec2.tolist())
        return new_vec

    def get_jocabian_matrix_p3(Qbb, R, lamb_matrix):
        # 获得R中各元素
        x1 = R[0, 0]; x2 = R[1, 0]; x3 = R[2, 0]
        x4 = R[0, 1]; x5 = R[1, 1]; x6 = R[2, 1]
        x7 = R[0, 2]; x8 = R[1, 2]; x9 = R[2, 2]
        # 构造雅可比矩阵的各分块矩阵
        part11 = np.linalg.inv(Qbb) - np.kron(lamb_matrix, np.eye(3))
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

    def get_jocabian_matrix_p2(Qbb, R, lamb_matrix):
        # 获得R中各元素
        x1 = R[0, 0]; x2 = R[1, 0]; x3 = R[2, 0]
        x4 = R[0, 1]; x5 = R[1, 1]; x6 = R[2, 1]
        x7 = R[0, 2]; x8 = R[1, 2]; x9 = R[2, 2]
        # 构造雅可比矩阵的各分块矩阵
        part11 = np.linalg.inv(Qbb) - np.kron(lamb_matrix, np.eye(3))
        part12 = np.array([
            [-x1, 0, 0],
            [-x2, 0, 0],
            [-x3, 0, 0],
            [0, -x4, 0],
            [0, -x5, 0],
            [0, -x6, 0]
        ])
        part21 = np.array([
            [2*x1, 2*x2, 2*x3, 0, 0, 0],
            [0, 0, 0, 2*x4, 2*x5, 2*x6],
            [0, 0, 0, 0, 0, 0]
        ])
        part22 = np.zeros((6, 6))
        # 构造雅可比矩阵
        jacobian_matrix = np.block([[part11, part12], [part21, part22]])
        return jacobian_matrix

    def from_y_get_Rlamb(y):
        R = np.array([y[0:3], y[3:6], y[6:9]]).T
        # R = np.array([y[0:3], y[3:6], y[6:9]])
        lamb_list= y[9:].tolist()
        return R, lamb_list

    def from_Rlamb_get_y(R, lamb_list):
        y = []
        y += np.ravel(R, order='F').tolist()
        y += lamb_list
        y = np.array(y)
        return y

    def calculate_QvecR_check(F, QBB):    # QBB为vec(B)对应的vc阵
        F_matrix_part = np.linalg.inv(F @ F.T) @ F
        # F_matrix_part =F    # 错误
        F_matrix = diagonalize_squarematrix(diagonalize_squarematrix(F_matrix_part, F_matrix_part), F_matrix_part)
        # Qvec_B = QBB
        # for i in range(B.shape[1]-1):
        #     Qvec_B = diagonalize_squarematrix(Qvec_B, QBB)
        Qr_check = F_matrix @ QBB @ F_matrix.T
        return Qr_check

    def from_Rvec_get_Rmatrix(R_vec):
        R_matrix = np.array([R_vec[:3].tolist(), R_vec[3:6].tolist(), R_vec[6:9].tolist()]).T
        return R_matrix


    """
    以下为计算过程
    """
    # 根据F(与几何构造有关)的秩得到O空间的维数p
    p = np.linalg.matrix_rank(F)

    # 根据最小二乘计算R_check
    # [27.33]解法，要求基线多于两条
    # R_check = B @ F.T @ np.linalg.inv(F @ F.T)
    # Q = calculate_QvecR_check(F, QBB)

    # Springer Handbook上的解法
    Q = np.linalg.inv(np.kron(F.T, np.eye(3)).T @ np.linalg.inv(QBB) @ np.kron(F.T, np.eye(3)))
    R_check = Q @ np.kron(F.T, np.eye(3)).T @ np.linalg.inv(QBB) @ np.ravel(B, order='F')
    R_check = from_Rvec_get_Rmatrix(R_check)


    # 开始迭代
    R=R_check
    if p == 3:
        lamb_list = [1, 1, 1, 1, 1, 1]
        n = 0
        while True:
            lamb_matrix = get_lamb_matrix_p3(lamb_list)
            g1 = (np.linalg.inv(Q) - np.kron(lamb_matrix, np.eye(3))) @ np.ravel(R, order='F') - np.linalg.inv(Q) @ np.ravel(R_check, order='F')
            g2_R = R.T @ R - np.eye(3)
            g2 = np.array([g2_R[0, 0], g2_R[1, 1], g2_R[2, 2], g2_R[1, 0], g2_R[2, 0], g2_R[2, 1]])
            g = vec_joint(g1, g2)
            y = from_Rlamb_get_y(R, lamb_list)
            J = get_jocabian_matrix_p3(Q, R, lamb_matrix)
            dy = -np.linalg.inv(J) @ g
            # print(dy)
            y += dy
            R, lamb_list = from_y_get_Rlamb(y)
            # Q = J[:9, :9] @ Q @ J[:9, :9].T
            if np.all(abs(dy) < 1e-8):
                # 计算vc阵
                Qrr = J[:9, :9] @ Q @ J[:9, :9].T
                break
            n+=1
            if n>10:
                Qrr = []
                break

    elif p == 2:
        lamb_list = [1, 1, 1]
        n=0
        while True:
            lamb_matrix = get_lamb_matrix_p2(lamb_list)
            g1 = (np.linalg.inv(Q) - np.kron(lamb_matrix, np.eye(3))) @ np.ravel(R, order='F') - np.linalg.inv(Q) @ np.ravel(R_check, order='F')
            g2_R = R.T @ R - np.eye(3)
            g2 = np.array([g2_R[0, 0], g2_R[1, 1], g2_R[2, 2], g2_R[1, 0], g2_R[2, 0], g2_R[2, 1]])
            g = vec_joint(g1, g2)
            y = from_Rlamb_get_y(R, lamb_list)
            J = get_jocabian_matrix_p2(Q, R, lamb_matrix)
            dy = -np.linalg.inv(J) @ g
            y += dy
            R, lamb_list = from_y_get_Rlamb(y)
            # Q = J[:9, :9] @ Q @ J[:9, :9].T
            if np.all(abs(dy) < 1e-8):
                # 计算vc阵
                Qrr = J[:9, :9] @ Q @ J[:9, :9].T
                break
            n+=1
            if n>10:
                Qrr=[]
                break

    return R, Qrr, R_check



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
    b1 = add_perturbation_to_vector(r.apply(f1), [0.02, 0.02, 0.02])  # CUCC -> CUBB
    b2 = add_perturbation_to_vector(r.apply(f2), [0.02, 0.02, 0.02])  # CUT0 -> CUBB
    b3 = add_perturbation_to_vector(r.apply(f3), [0.02, 0.02, 0.02])  # CUAA -> CUT0
    b4 = add_perturbation_to_vector(r.apply(f4), [0.02, 0.02, 0.02])  # CUCC -> CUT0
    # b1 = add_perturbation_to_vector(r.apply(f1), [0.001, 0.001, 0.001])    # CUCC -> CUBB
    # b2 = add_perturbation_to_vector(r.apply(f2), [0.001, 0.001, 0.001])   # CUT0 -> CUBB
    # b3 = add_perturbation_to_vector(r.apply(f3), [0.001, 0.001, 0.001])   # CUAA -> CUT0
    # b4 = add_perturbation_to_vector(r.apply(f4), [0.001, 0.001, 0.001])    # CUCC -> CUT0


    # 组成向量矩阵
    F = get_matrix_from_vectors([f1.tolist(), f2.tolist(), f3.tolist(), f4.tolist()])
    B = get_matrix_from_vectors([b1.tolist(), b2.tolist(), b3.tolist(), b4.tolist()])
    # F = get_matrix_from_vectors([f1.tolist(), f2.tolist()])
    # B = get_matrix_from_vectors([b1.tolist(), b2.tolist()])

    # 解算WOPP问题
    Qbb = np.array([[0.01, 0.02, 0.01], [0.02, 0.01, 0.03], [0.01, 0.01, 0.01]])
    QBB = Qbb
    for i in range(F.shape[1] - 1):
        QBB = diagonalize_squarematrix(QBB, Qbb)
    R, Qrr, R_check = solve_WOPP_withLagrangianMultipliers(B, F, QBB)

    # print(R-R_check)
    print(R_check - r.as_matrix())
    print(R - r.as_matrix())
    # print(R)
    # print(Qrr)

    # print(r.as_euler('zyx', degrees=True))
    print(rota.from_matrix(R_check).as_euler('zyx', degrees=True))
    print(rota.from_matrix(R).as_euler('zyx', degrees=True))









