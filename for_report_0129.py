# -*- coding: utf-8 -*-
"""

@title:	for_report
@author: iDeal0103
@status:	Active
@type:	Process
@created:	29-Jan-2022
@post-History:	29-Jan-2022

comment：
    进展汇报
    1.单点定位
    2.差分定位
      1）载波相位双差双历元
      2）载波相位+伪距双差单历元
    3.单

    选取 Curting University的观测站CUAA CUBB CUCC CUT0 的观测数据

"""



# import module
import datetime

import numpy as np

import SinglePointPosition as SPP
import RTK
import RTK2
import utils.DoFile as DoFile
from attitude_determination.WOPP import *
from attitude_determination.TRIAD import *
import matplotlib.pyplot as plt
import utils.ResultAnalyse as ResultAnalyse

# 设置文件
# broadcast_file = r"edata\attitude_project\rinex2\brdc0030.21n"
# CUAA_observation_file = r"edata\attitude_project\rinex2\cuaa0030.21o"
# CUBB_observation_file = r"edata\attitude_project\rinex2\cubb0030.21o"
# CUCC_observation_file = r"edata\attitude_project\rinex2\cucc0030.21o"
# CUT0_observation_file = r"edata\attitude_project\rinex2\cut00030.21o"
#
# # 读取文件
# br_records = DoFile.read_GPS_nFile(broadcast_file)
# CUAA_ob_records = DoFile.read_Rinex2_oFile(CUAA_observation_file)
# CUBB_ob_records = DoFile.read_Rinex2_oFile(CUBB_observation_file)
# # CUCC_ob_records = DoFile.read_Rinex2_oFile(CUCC_observation_file)
# # CUT0_ob_records = DoFile.read_Rinex2_oFile(CUT0_observation_file)
# # print('文件读取完毕!')
#
# # 设置观测站参考坐标
# CUAA_coordinate = [-2364336.1554, 4870280.8223, -3360815.9725]
# CUBB_coordinate = [-2364334.0912, 4870285.4649, -3360807.7523]
# CUCC_coordinate = [-2364332.0167, 4870284.0337, -3360814.1380]
# CUT0_coordinate = [-2364338.2011, 4870284.7857, -3360809.1387]
#
# # 设置起始时间
# Tr_strat = datetime.datetime(2021, 1, 3, 0, 1, 0)
# Tr_end = datetime.datetime(2021, 1, 3, 0, 59, 00)

# 单点定位
# true_coors = []
# cal_coors = []
# vs = []
# Tr = Tr_strat
# print("开始解算各历元单点定位!")
# while Tr < Tr_end:
#     Xk, Yk, Zk, Q, v = SPP.SPP_on_broadcastrecords(CUBB_ob_records, br_records, Tr, cutoff=0,
#                                                    init_coor=CUBB_coordinate, recalP=True, doTDC=True, doIDC=True)
#     print(Xk, Yk, Zk)
#     cal_coors.append([Xk, Yk, Zk])
#     true_coors.append(CUBB_coordinate)
#     vs.append(v)
#     Tr += datetime.timedelta(seconds=30)
# SPP.cal_NEUerrors(true_coors, cal_coors)
# SPP.cal_XYZerrors(true_coors, cal_coors)
# SPP.cal_Coorerrors(true_coors, cal_coors)
# print("neu各方向RMSE:", ResultAnalyse.get_NEU_rmse(true_coors, cal_coors))
# print("坐标RMSE:", ResultAnalyse.get_coor_rmse(true_coors, cal_coors))

# # 差分定位
# true_coors = []
# cal_coors = []
# Tr = Tr_strat
# print("开始解算各历元差分定位!")
# while Tr < Tr_end:
#     # Tr2 = Tr + datetime.timedelta(seconds=30*60)
#     print(Tr.hour, Tr.minute, Tr.second)
#     CoorXYZ, Q = RTK2.DD_onCarrierPhase_and_Pseudorange_1known(CUAA_ob_records, CUBB_ob_records, br_records, Tr,
#                                   CUAA_coordinate, CUBB_coordinate, cutoff=5, ambi_fix=True)
#     Xk, Yk, Zk = CoorXYZ
#     cal_coors.append([Xk, Yk, Zk])
#     true_coors.append(CUBB_coordinate)
#     Tr += datetime.timedelta(seconds=30)
# SPP.cal_NEUerrors(true_coors, cal_coors)
# SPP.cal_XYZerrors(true_coors, cal_coors)
# SPP.cal_Coorerrors(true_coors, cal_coors)
# print("neu各方向RMSE:", ResultAnalyse.get_NEU_rmse(true_coors, cal_coors))
# print("坐标RMSE:", ResultAnalyse.get_coor_rmse(true_coors, cal_coors))

# true_coors = []
# cal_coors = []
# Tr = Tr_strat
# print("开始解算各历元差分定位!")
# while Tr < Tr_end:
#     Tr2 = Tr + datetime.timedelta(seconds=30*60)
#     print(Tr.hour, Tr.minute, Tr.second)
#     CoorXYZ, Q = RTK.DD_onCarrierPhase_1known(CUT0_ob_records, CUBB_ob_records, br_records, Tr, Tr2,
#                                   CUT0_coordinate, CUBB_coordinate, cutoff=15, ambi_fix=False)
#     Xk, Yk, Zk = CoorXYZ
#     cal_coors.append([Xk, Yk, Zk])
#     true_coors.append(CUBB_coordinate)
#     Tr += datetime.timedelta(seconds=30)
# SPP.cal_NEUerrors(true_coors, cal_coors)
# SPP.cal_XYZerrors(true_coors, cal_coors)


# OPP和WOPP定姿
# 在参考框架坐标系下
f1 = np.array([-2.2775, 2.4688, 5.1465])  # CUCC -> CUBB
f2 = np.array([3.9066, 1.7174, 0.1474])  # CUT0 -> CUBB
f3 = np.array([-2.0202, 4.1602, 7.0336])  # CUAA -> CUT0
f4 = np.array([-6.1844, 0.7520, 4.9993])  # CUCC -> CUT0
F = get_matrix_from_vectors([f1.tolist(), f2.tolist(), f3.tolist(), f4.tolist()])
# F = get_matrix_from_vectors([f1.tolist(), f2.tolist()])
#
r = rota.from_euler('zyx', [5, 7, 10], degrees=True)

WOPP_residual = []
# OPP_residual = []
# TRIAD_residual = []
dr = 0.01
for i in range(1000):
    b1 = add_perturbation_to_vector(r.apply(f1), [dr, dr, dr])    # CUCC -> CUBB
    b2 = add_perturbation_to_vector(r.apply(f2), [dr, dr, dr])   # CUT0 -> CUBB
    b3 = add_perturbation_to_vector(r.apply(f3), [dr, dr, dr])   # CUAA -> CUT0
    b4 = add_perturbation_to_vector(r.apply(f4), [dr, dr, dr])    # CUCC -> CUT0
    B = get_matrix_from_vectors([b1.tolist(), b2.tolist(), b3.tolist(), b4.tolist()])
    # B = get_matrix_from_vectors([b1.tolist(), b2.tolist()])
#
#     # OPP
#     r_check_SVD = solve_OPP_withSVD(F, B, [1, 1, 1, 1])
#     OPP_residual.append((np.array(rota.from_matrix(r_check_SVD).as_euler('zyx', degrees=True)) - np.array([5, 7, 10])).tolist())
#
    # WOPP
    Qbb = np.array([[0.02, 0.01, 0.01], [0.01, 0.02, 0.01], [0.01, 0.01, 0.02]])
    R, Qrr= solve_WOPP_withLagrangianMultipliers(B, F, Qbb)
    WOPP_residual.append((np.array(rota.from_matrix(R).as_euler('zyx', degrees=True)) - np.array([5, 7, 10])).tolist())

    # # TRIAD
    # R_TRIAD = solve_TRIAD(F, B)
    # TRIAD_residual.append((np.array(rota.from_matrix(R_TRIAD).as_euler('zyx', degrees=True)) - np.array([5, 7, 10])).tolist())

#
# # 画图
x = list(range(1000))
# plt.rcParams['font.sans-serif'] = ['SimHei']
# plt.rcParams['axes.unicode_minus'] = False
# # plt.plot(np.array(OPP_residual)[:, 0], color="r", label="delta z / 度")
# # plt.plot(np.array(OPP_residual)[:, 1], color="g", label="delta x / 度")
# # plt.plot(np.array(OPP_residual)[:, 2], color="b", label="delta y / 度")
# plt.scatter(x, np.array(OPP_residual)[:, 0], color="r", label="delta z / 度")
# plt.scatter(x, np.array(OPP_residual)[:, 1], color="g", label="delta x / 度")
# plt.scatter(x, np.array(OPP_residual)[:, 2], color="b", label="delta y / 度")
# plt.legend(loc='upper right')
# plt.title("OPP姿态角度解算模拟(dr=%sm)"%str(dr))
# plt.show()
#
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False
# plt.plot(np.array(WOPP_residual)[:, 0], color="r", label="delta z / 度")
# plt.plot(np.array(WOPP_residual)[:, 1], color="g", label="delta x / 度")
# plt.plot(np.array(WOPP_residual)[:, 2], color="b", label="delta y / 度")
plt.scatter(x, np.array(WOPP_residual)[:, 0], color="r", label="delta z / 度")
plt.scatter(x, np.array(WOPP_residual)[:, 1], color="g", label="delta x / 度")
plt.scatter(x, np.array(WOPP_residual)[:, 2], color="b", label="delta y / 度")
plt.legend(loc='upper right')
plt.title("WOPP姿态角度解算模拟(dr=%sm)"%str(dr))
plt.show()

# plt.rcParams['font.sans-serif'] = ['SimHei']
# plt.rcParams['axes.unicode_minus'] = False
# plt.plot(np.array(TRIAD_residual)[:, 0], color="r", label="delta z / 度")
# plt.plot(np.array(TRIAD_residual)[:, 1], color="g", label="delta x / 度")
# plt.plot(np.array(TRIAD_residual)[:, 2], color="b", label="delta y / 度")
# plt.scatter(x, np.array(TRIAD_residual)[:, 0], color="r", label="delta z / 度")
# plt.scatter(x, np.array(TRIAD_residual)[:, 1], color="g", label="delta x / 度")
# plt.scatter(x, np.array(TRIAD_residual)[:, 2], color="b", label="delta y / 度")
# plt.legend(loc='upper right')
# plt.title("TRIAD姿态角度解算模拟(dr=%sm)"%str(dr))
# plt.show()














