# -*- coding: utf-8 -*-
"""

@title:	ResultAnalyse
@author: iDeal0103
@status:	Active
@type:	Process
@created:	5-Mar-2022
@post-History:	5-Mar-2022

comment：
    1.计算结果的均方根误差


"""
#import module
import math
import utils.CoorTransform as CoorTransform



# 计算平均值
def get_average(records):
    return sum(records) / len(records)

# 计算方差
def get_variance(records):
    """
    方差 反映一个数据集的离散程度
    """
    average = get_average(records)
    return sum([(x - average) ** 2 for x in records]) / len(records)

# 计算均方差
def get_standard_deviation(records):
    """
    标准差 == 均方差 反映一个数据集的离散程度
    """
    variance = get_variance(records)
    return math.sqrt(variance)

#计算均方根值
def get_rms(records):
    """
    均方根值 反映的是有效值而不是平均值
    """
    return math.sqrt(sum([x ** 2 for x in records]) / len(records))

# 计算均方误差
def get_mse(records_real, records_predict):
    """
    均方误差 估计值与真值 偏差
    """
    if len(records_real) == len(records_predict):
        return sum([(x - y) ** 2 for x, y in zip(records_real, records_predict)]) / len(records_real)
    else:
        print('数据列表长度不相符')

def get_rmse(records_real, records_predict):
    """
    均方根误差：是均方误差的算术平方根
    """
    mse = get_mse(records_real, records_predict)
    if mse:
        return math.sqrt(mse)
    # else:
    #     print('数据列表长度不相符')

def get_mae(records_real, records_predict):
    """
    平均绝对误差
    """
    if len(records_real) == len(records_predict):
        return sum([abs(x - y) for x, y in zip(records_real, records_predict)]) / len(records_real)
    else:
        print('数据列表长度不相符')


# 计算NEU方向上的均方根误差
def get_NEU_rmse(true_coors, cal_coors):
    delta_N = []
    delta_E = []
    delta_U = []
    for i in range(len(true_coors)):
        n, e, u = CoorTransform.cal_NEU(true_coors[i], cal_coors[i])
        delta_N.append(n)
        delta_E.append(e)
        delta_U.append(u)
    rms_N = get_rms(delta_N)
    rms_E = get_rms(delta_E)
    rms_U = get_rms(delta_U)
    return rms_N, rms_E, rms_U

# 计算XYZ方向上的均方根误差
def get_XYZ_rmse(true_coors, cal_coors):
    delta_X = []
    delta_Y = []
    delta_Z = []
    for i in range(len(true_coors)):
        dx = cal_coors[i][0] - true_coors[i][0]
        dy = cal_coors[i][1] - true_coors[i][1]
        dz = cal_coors[i][2] - true_coors[i][2]
        delta_X.append(dx)
        delta_Y.append(dy)
        delta_Z.append(dz)
    rms_X = get_rms(delta_X)
    rms_Y = get_rms(delta_Y)
    rms_Z = get_rms(delta_Z)
    return rms_X, rms_Y, rms_Z

# 计算坐标总的均方根误差
def get_coor_rmse(true_coors, cal_coors):
    delta_d = []
    for i in range(len(true_coors)):
        dx = cal_coors[i][0] - true_coors[i][0]
        dy = cal_coors[i][1] - true_coors[i][1]
        dz = cal_coors[i][2] - true_coors[i][2]
        d = math.sqrt(dx**2 + dy**2 + dz**2)
        delta_d.append(d)
    rms_d = get_rms(delta_d)
    return rms_d













