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
import SinglePointPosition as SPP
import RTK
import RTK2
import utils.DoFile as DoFile


# 设置文件
broadcast_file = r"edata\attitude_project\rinex2\brdc0030.21n"
CUAA_observation_file = r"edata\attitude_project\rinex2\cuaa0030.21o"
CUBB_observation_file = r"edata\attitude_project\rinex2\cubb0030.21o"
CUCC_observation_file = r"edata\attitude_project\rinex2\cucc0030.21o"
CUT0_observation_file = r"edata\attitude_project\rinex2\cut00030.21o"

# 读取文件
br_records = DoFile.read_GPS_nFile(broadcast_file)
CUAA_ob_records = DoFile.read_Rinex2_oFile(CUAA_observation_file)
CUBB_ob_records = DoFile.read_Rinex2_oFile(CUBB_observation_file)
CUCC_ob_records = DoFile.read_Rinex2_oFile(CUCC_observation_file)
CUT0_ob_records = DoFile.read_Rinex2_oFile(CUT0_observation_file)
print('文件读取完毕!')

# 设置观测站参考坐标
CUAA_coordinate = [-2364336.1554, 4870280.8223, -3360815.9725]
CUBB_coordinate = [-2364334.0912, 4870285.4649, -3360807.7523]
CUCC_coordinate = [-2364332.0167, 4870284.0337, -3360814.1380]
CUT0_coordinate = [-2364338.2011, 4870284.7857, -3360809.1387]

# 设置起始时间
Tr_strat = datetime.datetime(2021, 1, 3, 0, 1, 0)
Tr_end = datetime.datetime(2021, 1, 3, 0, 59, 00)

# 单点定位
# true_coors = []
# cal_coors = []
# vs = []
# Tr = Tr_strat
# print("开始解算各历元单点定位!")
# while Tr < Tr_end:
#     Xk, Yk, Zk, Q, v = SPP.SPP_on_broadcastrecords(CUAA_ob_records, br_records, Tr, cutoff=15,
#                                                    init_coor=CUAA_coordinate, recalP=True, doTDC=True, doIDC=True)
#     print(Xk, Yk, Zk)
#     cal_coors.append([Xk, Yk, Zk])
#     true_coors.append(CUAA_coordinate)
#     vs.append(v)
#     Tr += datetime.timedelta(seconds=30)
# SPP.cal_NEUerrors(true_coors, cal_coors, ylimit=10)
# SPP.cal_XYZerrors(true_coors, cal_coors, ylimit=10)

# 差分定位
true_coors = []
cal_coors = []
Tr = Tr_strat
print("开始解算各历元差分定位!")
while Tr < Tr_end:
    Tr2 = Tr + datetime.timedelta(seconds=30*60)
    print(Tr.hour, Tr.minute, Tr.second)
    CoorXYZ, Q = RTK.DD_onCarrierPhase_1known(CUCC_ob_records, CUBB_ob_records, br_records, Tr, Tr2,
                                  CUCC_coordinate, CUBB_coordinate, cutoff=15, ambi_fix=True)
    Xk, Yk, Zk = CoorXYZ
    cal_coors.append([Xk, Yk, Zk])
    true_coors.append(CUBB_coordinate)
    Tr += datetime.timedelta(seconds=30)
SPP.cal_NEUerrors(true_coors, cal_coors)








