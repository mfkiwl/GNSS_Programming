# -*- coding: utf-8 -*-
"""

@title:	dynamic baseline resolution
@author: iDeal0103
@status:	Active
@type:	Process
@created:	25-Apr-2022
@post-History:	25-Apr-2022

comment：
     1. 动态基线的解算

"""
import math

import RTK_withfilter
import RTK_withfilter_baselineconstrain
import utils.DoFile as DoFile
import datetime
import matplotlib.pyplot as plt



antenna1_observation_file = r"D:\Desktop\毕业设计\数据\origin\main.22O"
antenna2_observation_file = r"D:\Desktop\毕业设计\数据\origin\AUX1.22O"
# broadcast_file = r"D:\Desktop\毕业设计\数据\origin\base.22p"
# broadcast_file = r"D:\Desktop\毕业设计\数据\origin\BRDC00IGS_R_20220870000_01D_MN.rnx"    # Rinex3.04
broadcast_file = r"D:\Desktop\毕业设计\数据\origin\brdc0870.22n"        # Rinex 2

# 读取文件
knownAntenna_ob_records = DoFile.read_Rinex3_oFile(antenna1_observation_file)
unknownAntenna_ob_records = DoFile.read_Rinex3_oFile(antenna2_observation_file)
# br_records = DoFile.read_Renix304_nFile(broadcast_file)
br_records = DoFile.read_GPS_nFile(broadcast_file, False)
print("数据读取完毕！")

antenna1_coor = [-2850374.3006, 4651839.9140, 3292945.7368]
antenna2_coor = [-2850374.3006, 4651839.9140, 3292945.7368]

strat_time = datetime.datetime(2022, 3, 28, 13, 30, 16)
# end_time = datetime.datetime(2022, 3, 28, 13, 40, 40)
end_time = datetime.datetime(2022, 3, 28, 13, 44, 40)
cal_coors = RTK_withfilter.constaneously_RTK_withfilter(strat_time, end_time, knownAntenna_ob_records, unknownAntenna_ob_records,
                                         br_records, antenna1_coor, antenna2_coor, 'L1C', 'C1C', 1.0, ele_limit=7)

bl_len = 0.735
# cal_coors = RTK_withfilter_baselineconstrain.constaneously_RTK_withfilter_baseline_length_constrained(strat_time, end_time, knownAntenna_ob_records, unknownAntenna_ob_records,
#                                          br_records, antenna1_coor, antenna2_coor, 'L1C', 'C1C', 1.0, bl_len, 0.1, ele_limit=7)


# cal_coors = RTK_withfilter_baselineconstrain.constaneously_RTK_withfilter_SPP_baseline_length_constrained(strat_time, end_time, knownAntenna_ob_records, unknownAntenna_ob_records,
#                                          br_records, antenna1_coor, 'L1C', 'C1C', 1.0, bl_len, ele_limit=7)

baseline_vectors = []
for coor in cal_coors:
    baseline_vector = [coor[0]-antenna1_coor[0], coor[1]-antenna1_coor[1], coor[2]-antenna1_coor[2]]
    baseline_vectors.append(baseline_vector)

bl_len = 0.735
bl_len_errors = []
for bl in baseline_vectors:
    bl_len_error = math.sqrt(bl[0]**2 + bl[1]**2 + bl[2]**2) - bl_len
    bl_len_errors.append(bl_len_error)

# 画图
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False
plt.plot(bl_len_errors, color="r", label="delta l / m")
plt.legend(loc='upper right')
plt.title("基线长误差序列图")
plt.show()


