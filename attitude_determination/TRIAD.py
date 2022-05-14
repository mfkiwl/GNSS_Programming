"""

@title:	TRIAD
@author: iDeal0103
@status:	Active
@type:	Process
@created:	25-Jan-2022
@post-History:	25-Jan-2022

comment：
    1.TRIAD

"""


import math
from attitude_determination.OPP import *
import numpy as np
from scipy.spatial.transform.rotation import Rotation as rota


def get_unit_vector(vector):
    """
    vector : np.darray, 向量坐标
    """
    vector_len = math.sqrt(vector.T @ vector)
    unit_vector = vector / vector_len
    return unit_vector, vector_len

def get_vector_length(vector):
    """
    vector : np.darray, 向量坐标
    """
    vector_len = math.sqrt(vector.T @ vector)


def solve_TRIAD(F, B):
    """
    F: np.darray, 由列向量组成的参考坐标系框架下的向量矩阵 [f1, f2]
    B_hat: np.darray, 由列向量组成的参考坐标系框架下的向量矩阵 [b1, b2]
    """
    b1 = B[:, 0]
    b2 = B[:, 1]
    f1 = F[:, 0]
    f2 = F[:, 1]
    # 取第一组对应基线计算
    v1 = get_unit_vector(b1)[0]
    u1 = get_unit_vector(f1)[0]
    # 获取第二组正交基
    v2 = get_unit_vector(np.cross(b1, b2))[0]
    u2 = get_unit_vector(np.cross(f1, f2))[0]
    # 获取第三组正交基
    v3 = np.cross(v1, v2)
    u3 = np.cross(u1, u2)
    # 计算旋转矩阵
    R_TRIAD = np.mat(v1).T @ np.mat(u1) + np.mat(v2).T @ np.mat(u2) + np.mat(v3).T @ np.mat(u3)
    return R_TRIAD

if __name__ == "__main__":
    # 在参考框架坐标系下
    f1 = np.array([-2.2775, 2.4688, 5.1465])  # CUCC -> CUBB
    f2 = np.array([3.9066, 1.7174, 0.1474])  # CUT0 -> CUBB

    # 在空间指教坐标系下
    r = rota.from_euler('zyx', [5, 7, 10], degrees=True)
    print(r.as_matrix())
    # b1 = r.apply(f1)    # CUCC -> CUBB
    # b2 = r.apply(f2)    # CUT0 -> CUBB
    b1 = add_perturbation_to_vector(r.apply(f1), [0.01, 0.01, 0.01])  # CUCC -> CUBB
    b2 = add_perturbation_to_vector(r.apply(f2), [0.01, 0.01, 0.01])  # CUT0 -> CUBB

    F = get_matrix_from_vectors([f1.tolist(), f2.tolist()])
    B = get_matrix_from_vectors([b1.tolist(), b2.tolist()])

    R = solve_TRIAD(F, B)
    print(rota.from_matrix(R).as_euler('zyx', degrees=True))

