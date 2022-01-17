# -*- coding: utf-8 -*-
"""

@title:	Style Guide for Python Code
@author: iDeal0103
@status:	Active
@type:	Process
@created:	14-Apr-2021
@post-History:	14-Apr-2021

comment：
    1.单点定位
    2.

"""

# import module
import numpy as np
import math
import datetime
import utils.DoFile as DoFile
import utils.SatellitePosition as SatellitePosition
import utils.TimeSystem as TimeSystem
import utils.CoorTransform as CoorTransform
import matplotlib.pyplot as plt


# 电离层改正
def Ionospheric_Delay_Correction(record):
    """
    Parameters
    record : GPS_observation_record , 一条观测记录

    Returns
    P : 电离层改正后的伪距
    """
    # 如果两个频道有缺失,则返回“”
    if record.data['P1']['observation'] == "" or record.data['P2']['observation'] == "":
        print("数据缺失")
        P = ""
    # 如果两个频道都有数据,则返回电离层改正后的伪距观测值
    else:
        P1 = int(record.data['P1']['observation'])
        P2 = int(record.data['P2']['observation'])
        f1 = 1575.42  # Hz
        f2 = 1227.60  # Hz
        P = f1 ** 2 / (f1 ** 2 - f2 ** 2) * P1 - f2 ** 2 / (f1 ** 2 - f2 ** 2) * P2
    return P


# 对流层改正计算过程中的函数
# Saastamoinen对流层改正
def Tropospheric_Delay_Correction_Saastamoinen(P, XYZ_pos, sat_pos, Ps=1013.75, es=8, T=276):
    """
    P : m , 伪距值
    XYZ_pos : float list [X,Y,Z] , 地面点直角坐标
    sat_pos : float list [X,Y,Z] , 卫星直角坐标
    """
    # 计算天顶对流层延迟模型
    #  fai:float,弧度表示的纬度     h:float,大地高
    fai, L, h = CoorTransform.cal_XYZ2BLH(XYZ_pos[0], XYZ_pos[1], XYZ_pos[2])
    ele, A = CoorTransform.cal_ele_and_A(XYZ_pos, sat_pos)
    delta_Phz = 0.0022767 * Ps / (1 - 0.00266 * math.cos(2 * fai) - 0.00028 * h)
    delta_Pwz = 0.0022767 * (1255 / T + 0.05) * es / (1 - 0.00266 * math.cos(2 * fai) - 0.00028 * h)
    # 计算投影函数
    mh = cal_mh_Ifadi(ele, Ps, T, es)
    mw = cal_mw_Ifadi(ele, Ps, T, es)
    # 计算对流层改正数值
    deltaP = delta_Phz * mh + delta_Pwz * mw
    # 计算对流层改正之后的数值
    P_corrected = P - deltaP
    return P_corrected


# Hopfield对流层改正
def Tropospheric_Delay_Correction_Hopfield(P, XYZ_pos, sat_pos, Ps=1012.75, es=75, Ts=298.15):
    """
    P : m , 伪距值
    XYZ_pos : float list [X,Y,Z] , 地面点直角坐标
    sat_pos : float list [X,Y,Z] , 卫星直角坐标
    """
    # 计算天顶对流层延迟模型
    #  fai:float,弧度表示的纬度     h:float,大地高
    fai, L, h = CoorTransform.cal_XYZ2BLH(XYZ_pos[0], XYZ_pos[1], XYZ_pos[2])
    ele, A = CoorTransform.cal_ele_and_A(XYZ_pos, sat_pos)
    # 计算函数中的分量
    hd = 40136 + 148.72 * (Ts - 273.16)
    hw = 11000
    Kd = 155.2e-7 * Ps / Ts * (hd - h)
    Kw = 155.2e-7 * 4810 / Ts ** 2 * es * (hw - h)
    # 计算投影函数
    # mh=cal_mh(ele,Ps,Ts,es)
    # mw=cal_mw(ele,Ps,Ts,es)
    mh = 1 / math.sin((ele ** 2 + du2angle(2.5) ** 2) ** 0.5)
    mw = 1 / math.sin((ele ** 2 + du2angle(1.5) ** 2) ** 0.5)
    # 计算对流层改正数值
    deltaP = Kd * mh + Kw * mw
    # 计算对流层改正之后的数值
    P_corrected = P - deltaP
    return P_corrected


# UNB3对流层改正
def Tropospheric_Delay_Correction_UNB3(P, t, XYZ_pos, sat_pos):
    """
    P : m , 伪距值
    t : datetime.datetime , 日期
    XYZ_pos : float list [X,Y,Z] , 地面点直角坐标
    sat_pos : float list [X,Y,Z] , 卫星直角坐标
    """
    # 计算天顶对流层延迟模型
    #  fai:float,弧度表示的纬度     h:float,大地高
    fai, L, h = CoorTransform.cal_XYZ2BLH(XYZ_pos[0], XYZ_pos[1], XYZ_pos[2])
    # #如果纬度<0,则不进行对流层改正
    # if fai<0:
    #     return P
    ele, A = CoorTransform.cal_ele_and_A(XYZ_pos, sat_pos)
    # 计算doy
    doy = cal_doy(t)
    # 计算插值项
    beita = cal_beita(angle2du(fai), doy)
    T = cal_T(angle2du(fai), doy)
    lamb = cal_lambda(angle2du(fai), doy)
    e = cal_e(angle2du(fai), doy)
    Pa = cal_Pa(angle2du(fai), doy)
    # 计算 t_hz
    t_hz = cal_t_hz(fai, h)
    # 计算Kh
    Kh = cal_Kh(beita, h, T)
    # 计算deltaP_hz
    deltaP_hz = t_hz * Kh * Pa
    # 计算t_wz
    t_wz = cal_t_wz(lamb, beita, fai, h, T)
    # 计算kw
    Kw = cal_Kw(beita, h, T, lamb)
    # 计算deltaP_wz
    deltaP_wz = t_wz * Kw * e / T
    # 计算投影函数
    mh = cal_mh_Neil(fai, doy, ele, h)
    mw = cal_mw_Neil(fai, ele, h)
    # 计算对流层延迟
    # 不加入大气梯度
    deltaP = mh * deltaP_hz + mw * deltaP_wz
    # 加入大气梯度
    deltaP = mh * deltaP_hz + mw * deltaP_wz + mh * (math.cos(A) + math.sin(A)) + mw * (math.cos(A) + math.sin(A))
    P_corrected = P - deltaP
    return P_corrected


# 计算gm
def cal_gm(fai, h):
    # fai 弧度
    gm = 9.784 * (1 - 2.66e-3 * math.cos(2 * fai) - 2.8e-7 * h)
    return gm


# 计算Kh
def cal_Kh(beita, h, T):
    g = 9.80665  # m/s^2
    Rh = 287.054  # J/kg/K
    Kh = math.pow((1 - beita * h / T), g / (Rh * beita))
    return Kh


# 计算Kw
def cal_Kw(beita, h, T, lamb):
    g = 9.80665  # m/s^2
    Rh = 287.054  # J/kg/K
    Kw = math.pow((1 - beita * h / T), (1 + lamb) * g / (Rh * beita) - 1)
    return Kw


# 计算t_hz
def cal_t_hz(fai, h):
    # fai 弧度
    k1 = 77.60  # K/hPa
    Rh = 287.054  # J/kg/K
    gm = 9.784 * (1 - 2.66e-3 * math.cos(2 * fai) - 2.8e-7 * h)
    t_hz = (1e-6 * k1 * Rh) / gm
    return t_hz


# 计算t_wz
def cal_t_wz(lamb, beita, fai, h, T):
    # fai 弧度
    k3 = 377600  # K^2/hpa
    k2 = 16.6  # K/hpa
    Rh = 287.054  # J/kg/K
    gm = 9.784 * (1 - 2.66e-3 * math.cos(2 * fai) - 2.8e-7 * h)
    Tm = T * (1 - beita * Rh / (gm * (1 + lamb)))
    t_wz = (1e-6 * (k3 + k2 * Tm) * Rh) / (gm * (1 + lamb) - Rh * beita)
    return t_wz


# UNB3模型插值函数计算
# 计算年积日
def cal_doy(time):
    """
    time : datetime.datetime类型
    """
    year = time.year
    month = time.month
    day = time.day
    # 闰年各月天数
    leap_year = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    # 平年各月天数
    no_leap_year = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    if (year % 4 == 0 and year % 100 != 0) or year % 400 == 0:
        doy = sum(leap_year[:month - 1]) + day
    else:
        doy = sum(no_leap_year[:month - 1]) + day
    return doy


# 完成角度与弧度之间的转换
def angle2du(angle):
    du = angle / math.pi * 180
    return du


def du2angle(du):
    angle = du / 180 * math.pi
    return angle


# 插值计算气压
def cal_Pa(fai, doy):
    # fai 角度
    Ps_av = [1013.25, 1017.25, 1015.75, 1011.75, 1013.00]
    Ps_am = [0, -3.75, -2.25, -1.75, -0.50]
    fais = [15, 30, 45, 60, 75]
    abs_fai = abs(fai)
    if abs_fai < 15:
        P_av = Ps_av[0]
        P_am = Ps_am[0]
    elif abs_fai > 15 and abs_fai < 30:
        P_av = Ps_av[0] + (Ps_av[1] - Ps_av[0]) / 15 * (abs_fai - fais[0])
        P_am = Ps_am[0] + (Ps_am[1] - Ps_am[0]) / 15 * (abs_fai - fais[0])
    elif abs_fai > 30 and abs_fai < 45:
        P_av = Ps_av[1] + (Ps_av[2] - Ps_av[1]) / 15 * (abs_fai - fais[1])
        P_am = Ps_am[1] + (Ps_am[2] - Ps_am[1]) / 15 * (abs_fai - fais[1])
    elif abs_fai > 45 and abs_fai < 60:
        P_av = Ps_av[2] + (Ps_av[3] - Ps_av[2]) / 15 * (abs_fai - fais[2])
        P_am = Ps_am[2] + (Ps_am[3] - Ps_am[2]) / 15 * (abs_fai - fais[2])
    elif abs_fai > 60 and abs_fai < 75:
        P_av = Ps_av[3] + (Ps_av[4] - Ps_av[3]) / 15 * (abs_fai - fais[3])
        P_am = Ps_am[3] + (Ps_am[4] - Ps_am[3]) / 15 * (abs_fai - fais[3])
    else:
        P_av = Ps_av[4]
        P_am = Ps_am[4]
    if fai > 0:
        Pa = P_av - P_am * math.cos(2 * math.pi * (doy - 28) / 365.25)
    else:
        Pa = P_av - P_am * math.cos(2 * math.pi * (doy - 211) / 365.25)
    return Pa


# 插值计算温度
def cal_T(fai, doy):
    # fai 角度
    Ts_av = [299.65, 294.15, 283.15, 272.15, 263.65]
    Ts_am = [0.00, 7.00, 11.00, 15.00, 14.50]
    fais = [15, 30, 45, 60, 75]
    abs_fai = abs(fai)
    if abs_fai<15:
        T_av = Ts_av[0]
        T_am = Ts_am[0]
    elif abs_fai>15 and abs_fai<30:
        T_av = Ts_av[0] + (Ts_av[1] - Ts_av[0]) / 15 * (abs_fai - fais[0])
        T_am = Ts_am[0] + (Ts_am[1] - Ts_am[0]) / 15 * (abs_fai - fais[0])
    elif abs_fai > 30 and abs_fai < 45:
        T_av = Ts_av[1] + (Ts_av[2] - Ts_av[1]) / 15 * (abs_fai - fais[1])
        T_am = Ts_am[1] + (Ts_am[2] - Ts_am[1]) / 15 * (abs_fai - fais[1])
    elif abs_fai > 45 and abs_fai < 60:
        T_av = Ts_av[2] + (Ts_av[3] - Ts_av[2]) / 15 * (abs_fai - fais[2])
        T_am = Ts_am[2] + (Ts_am[3] - Ts_am[2]) / 15 * (abs_fai - fais[2])
    elif abs_fai > 60 and abs_fai < 75:
        T_av = Ts_av[3] + (Ts_av[4] - Ts_av[3]) / 15 * (abs_fai - fais[3])
        T_am = Ts_am[3] + (Ts_am[4] - Ts_am[3]) / 15 * (abs_fai - fais[3])
    else:
        T_av = Ts_av[4]
        T_am = Ts_am[4]
    if fai>0:
        T = T_av - T_am * math.cos(2 * math.pi * (doy - 28) / 365.25)
    else:
        T = T_av - T_am * math.cos(2 * math.pi * (doy - 211) / 365.25)
    return T


# 插值计算湿度
def cal_e(fai, doy):
    # fai 角度
    es_av = [26.31, 21.79, 11.66, 6.78, 4.11]
    es_am = [0.00, 8.85, 7.24, 5.36, 3.39]
    fais = [15, 30, 45, 60, 75]
    abs_fai = abs(fai)
    if abs_fai<15:
        e_av = es_av[0]
        e_am = es_am[0]
    elif abs_fai>15 and abs_fai<30:
        e_av = es_av[0] + (es_av[1] - es_av[0]) / 15 * (abs_fai - fais[0])
        e_am = es_am[0] + (es_am[1] - es_am[0]) / 15 * (abs_fai - fais[0])
    elif abs_fai > 30 and abs_fai < 45:
        e_av = es_av[1] + (es_av[2] - es_av[1]) / 15 * (abs_fai - fais[1])
        e_am = es_am[1] + (es_am[2] - es_am[1]) / 15 * (abs_fai - fais[1])
    elif abs_fai > 45 and abs_fai < 60:
        e_av = es_av[2] + (es_av[3] - es_av[2]) / 15 * (abs_fai - fais[2])
        e_am = es_am[2] + (es_am[3] - es_am[2]) / 15 * (abs_fai - fais[2])
    elif abs_fai > 60 and abs_fai < 75:
        e_av = es_av[3] + (es_av[4] - es_av[3]) / 15 * (abs_fai - fais[3])
        e_am = es_am[3] + (es_am[4] - es_am[3]) / 15 * (abs_fai - fais[3])
    else:
        e_av = es_av[4]
        e_am = es_am[4]
    if fai>0:
        e = e_av - e_am * math.cos(2 * math.pi * (doy - 28) / 365.25)
    else:
        e = e_av - e_am * math.cos(2 * math.pi * (doy - 211) / 365.25)
    return e


# 插值计算温度变化率
def cal_beita(fai, doy, lowerindex=None, upperindex=None):
    # fai 角度
    beitas_av = [6.30e-3, 6.05e-3, 5.58e-3, 5.39e-3, 4.53e-3]
    beitas_am = [0.00e-3, 0.25e-3, 0.32e-3, 0.81e-3, 0.62e-3]
    # beitas_av=[6.30,6.05,5.58,5.39,4.53]
    # beitas_am=[0.00,0.25,0.32,0.81,0.62]
    fais = [15, 30, 45, 60, 75]
    abs_fai = abs(fai)
    if abs_fai < 15:
        beita_av = beitas_av[0]
        beita_am = beitas_am[0]
    elif abs_fai > 15 and abs_fai < 30:
        beita_av = beitas_av[0] + (beitas_av[1] - beitas_av[0]) / 15 * (abs_fai - fais[0])
        beita_am = beitas_am[0] + (beitas_am[1] - beitas_am[0]) / 15 * (abs_fai - fais[0])
    elif abs_fai > 30 and abs_fai < 45:
        beita_av = beitas_av[1] + (beitas_av[2] - beitas_av[1]) / 15 * (abs_fai - fais[1])
        beita_am = beitas_am[1] + (beitas_am[2] - beitas_am[1]) / 15 * (abs_fai - fais[1])
    elif abs_fai > 45 and abs_fai < 60:
        beita_av = beitas_av[2] + (beitas_av[3] - beitas_av[2]) / 15 * (abs_fai - fais[2])
        beita_am = beitas_am[2] + (beitas_am[3] - beitas_am[2]) / 15 * (abs_fai - fais[2])
    elif abs_fai > 60 and abs_fai < 75:
        beita_av = beitas_av[3] + (beitas_av[4] - beitas_av[3]) / 15 * (abs_fai - fais[3])
        beita_am = beitas_am[3] + (beitas_am[4] - beitas_am[3]) / 15 * (abs_fai - fais[3])
    else:
        beita_av = beitas_av[4]
        beita_am = beitas_am[4]
    if fai > 0:
        beita = beita_av - beita_am * math.cos(2 * math.pi * (doy - 28) / 365.25)
    else:
        beita = beita_av - beita_am * math.cos(2 * math.pi * (doy - 211) / 365.25)
    return beita


# 插值计算变量lambda
def cal_lambda(fai, doy, lowerindex=None, upperindex=None):
    # fai 角度
    lambdas_av = [2.77, 3.15, 2.57, 1.81, 1.55]
    lambdas_am = [0.00, 0.33, 0.46, 0.74, 0.30]
    fais = [15, 30, 45, 60, 75]
    abs_fai = abs(fai)
    if abs_fai < 15:
        lambda_av = lambdas_av[0]
        lambda_am = lambdas_am[0]
    elif abs_fai > 15 and abs_fai < 30:
        lambda_av = lambdas_av[0] + (lambdas_av[1] - lambdas_av[0]) / 15 * (abs_fai - fais[0])
        lambda_am = lambdas_am[0] + (lambdas_am[1] - lambdas_am[0]) / 15 * (abs_fai - fais[0])
    elif abs_fai > 30 and abs_fai < 45:
        lambda_av = lambdas_av[1] + (lambdas_av[2] - lambdas_av[1]) / 15 * (abs_fai - fais[1])
        lambda_am = lambdas_am[1] + (lambdas_am[2] - lambdas_am[1]) / 15 * (abs_fai - fais[1])
    elif abs_fai > 45 and abs_fai < 60:
        lambda_av = lambdas_av[2] + (lambdas_av[3] - lambdas_av[2]) / 15 * (abs_fai - fais[2])
        lambda_am = lambdas_am[2] + (lambdas_am[3] - lambdas_am[2]) / 15 * (abs_fai - fais[2])
    elif abs_fai > 60 and abs_fai < 75:
        lambda_av = lambdas_av[3] + (lambdas_av[4] - lambdas_av[3]) / 15 * (abs_fai - fais[3])
        lambda_am = lambdas_am[3] + (lambdas_am[4] - lambdas_am[3]) / 15 * (abs_fai - fais[3])
    else:
        lambda_av = lambdas_av[4]
        lambda_am = lambdas_am[4]
    if fai > 0:
        lamb = lambda_av - lambda_am * math.cos(2 * math.pi * (doy - 28) / 365.25)
    else:
        lamb = lambda_av - lambda_am * math.cos(2 * math.pi * (doy - 211) / 365.25)
    return lamb



# 定义投影函数
# 计算静力学投影系数
def cal_mh_Ifadi(ele, Ps, Ts, es):
    # ele 弧度
    ah = 0.1273e-2 + 0.1316e-6 * (Ps - 1000) + 0.1378e-5 * (Ts - 15) + 0.8057e-5 * math.sqrt(es)
    bh = 0.3333e-2 + 0.1946e-6 * (Ps - 1000) + 0.1040e-6 * (Ts - 15) + 0.1747e-6 * math.sqrt(es)
    ch = 0.078
    mh = 1 / (math.sin(ele) + (ah / (math.sin(ele) + bh / (math.sin(ele) + ch))))
    return mh


# 计算湿延迟投影系数
def cal_mw_Ifadi(ele, Ps, Ts, es):
    aw = 0.5236e-2 + 0.2471e-6 * (Ps - 1000) - 0.1724e-6 * (Ts - 15) + 0.1328e-4 * math.sqrt(es)
    bw = 0.1705e-2 + 0.7384e-6 * (Ps - 1000) + 0.3767e-6 * (Ts - 15) + 0.2147e-4 * math.sqrt(es)
    cw = 0.05917
    mw = 1 / (math.sin(ele) + (aw / (math.sin(ele) + bw / (math.sin(ele) + cw))))
    return mw


# 计算Neil投影函数
def cal_mh_Neil(fai, doy, ele, h):
    # fai 角度
    # ele 弧度
    fais = [15, 30, 45, 60, 75]
    a_avs = [1.2769934e-3, 1.2683230e-3, 1.2465397e-3, 1.2196049e-3, 1.2045996e-3]
    a_ams = [0.00, 1.2709626e-3, 2.6523662e-5, 3.4000452e-5, 4.1202191e-5]
    b_avs = [2.9153695e-3, 2.9152299e-3, 2.9288445e-3, 2.9022565e-3, 2.9024912e-3]
    b_ams = [0.00, 2.1414979e-5, 3.0160779e-5, 7.2562722e-5, 11.723375e-5]
    c_avs = [62.610505e-3, 62.837393e-3, 63.721774e-3, 63.824265e-3, 62.258455e-3]
    c_ams = [0.00, 9.0128400e-5, 4.3497037e-5, 84.795348e-5, 170.37206e-5]
    abs_fai = abs(fai)
    if abs_fai < 15:
        a_av = a_avs[0]
        a_am = a_ams[0]
        b_av = b_avs[0]
        b_am = b_ams[0]
        c_av = c_avs[0]
        c_am = c_ams[0]
    elif abs_fai>15 and abs_fai<30:
        a_av = a_avs[0] + (a_avs[1] - a_avs[0]) / 15 * (fai - fais[0])
        a_am = a_ams[0] + (a_ams[1] - a_ams[0]) / 15 * (fai - fais[0])
        b_av = b_avs[0] + (b_avs[1] - b_avs[0]) / 15 * (fai - fais[0])
        b_am = b_ams[0] + (b_ams[1] - b_ams[0]) / 15 * (fai - fais[0])
        c_av = c_avs[0] + (c_avs[1] - c_avs[0]) / 15 * (fai - fais[0])
        c_am = c_ams[0] + (c_ams[1] - c_ams[0]) / 15 * (fai - fais[0])
    elif abs_fai>30 and abs_fai<45:
        a_av = a_avs[1] + (a_avs[2] - a_avs[1]) / 15 * (fai - fais[1])
        a_am = a_ams[1] + (a_ams[2] - a_ams[1]) / 15 * (fai - fais[1])
        b_av = b_avs[1] + (b_avs[2] - b_avs[1]) / 15 * (fai - fais[1])
        b_am = b_ams[1] + (b_ams[2] - b_ams[1]) / 15 * (fai - fais[1])
        c_av = c_avs[1] + (c_avs[2] - c_avs[1]) / 15 * (fai - fais[1])
        c_am = c_ams[1] + (c_ams[2] - c_ams[1]) / 15 * (fai - fais[1])
    elif abs_fai > 45 and abs_fai < 60:
        a_av = a_avs[2] + (a_avs[3] - a_avs[2]) / 15 * (fai - fais[2])
        a_am = a_ams[2] + (a_ams[3] - a_ams[2]) / 15 * (fai - fais[2])
        b_av = b_avs[2] + (b_avs[3] - b_avs[2]) / 15 * (fai - fais[2])
        b_am = b_ams[2] + (b_ams[3] - b_ams[2]) / 15 * (fai - fais[2])
        c_av = c_avs[2] + (c_avs[3] - c_avs[2]) / 15 * (fai - fais[2])
        c_am = c_ams[2] + (c_ams[3] - c_ams[2]) / 15 * (fai - fais[2])
    elif abs_fai>60 and abs_fai<75:
        a_av = a_avs[3] + (a_avs[4] - a_avs[3]) / 15 * (fai - fais[3])
        a_am = a_ams[3] + (a_ams[4] - a_ams[3]) / 15 * (fai - fais[3])
        b_av = b_avs[3] + (b_avs[4] - b_avs[3]) / 15 * (fai - fais[3])
        b_am = b_ams[3] + (b_ams[4] - b_ams[3]) / 15 * (fai - fais[3])
        c_av = c_avs[3] + (c_avs[4] - c_avs[3]) / 15 * (fai - fais[3])
        c_am = c_ams[3] + (c_ams[4] - c_ams[3]) / 15 * (fai - fais[3])
    else:
        a_av = a_avs[4]
        a_am = a_ams[4]
        b_av = b_avs[4]
        b_am = b_ams[4]
        c_av = c_avs[4]
        c_am = c_ams[4]
    if fai>0:
        a = a_av - a_am * math.cos(2 * math.pi * (doy - 28) / 365.25)
        b = b_av - b_am * math.cos(2 * math.pi * (doy - 28) / 365.25)
        c = c_av - c_am * math.cos(2 * math.pi * (doy - 28) / 365.25)
    else:
        a = a_av - a_am * math.cos(2 * math.pi * (doy - 211) / 365.25)
        b = b_av - b_am * math.cos(2 * math.pi * (doy - 211) / 365.25)
        c = c_av - c_am * math.cos(2 * math.pi * (doy - 211) / 365.25)
    # 计算高程为0处的投影函数
    mh = (1 + a / (1 + b / (1 + c))) * 1 / (math.sin(ele) + (a / (math.sin(ele) + b / (math.sin(ele) + c))))
    # 计算高程改正项
    ah = 2.53e-5
    bh = 5.49e-3
    ch = 1.14e-3
    f = (1 + ah / (1 + bh / (1 + ch))) * 1 / (math.sin(ele) + (ah / (math.sin(ele) + bh / (math.sin(ele) + ch))))
    delta_mh = (1 / math.sin(ele) - f) * (h / 1000)
    # 对mh进行高程改正
    mh += delta_mh
    return mh


def cal_mw_Neil(fai, ele, h):
    fais = [15, 30, 45, 60, 75]
    a_avs = [5.8021897e-4, 5.6794847e-4, 5.8118019e-4, 5.9727542e-4, 6.1641693e-4]
    b_avs = [1.4275268e-3, 1.5138625e-3, 1.4572752e-3, 1.5007428e-3, 1.7599082e-3]
    c_avs = [4.3472961e-2, 4.6729510e-2, 4.3908931e-2, 4.4626982e-2, 5.4736038e-2]
    abs_fai = abs(fai)
    if abs_fai < 15:
        a_av = a_avs[0]
        b_av = b_avs[0]
        c_av = c_avs[0]
    elif abs_fai>15 and abs_fai<30:
        a_av = a_avs[0] + (a_avs[1] - a_avs[0]) / 15 * (fai - fais[0])
        b_av = b_avs[0] + (b_avs[1] - b_avs[0]) / 15 * (fai - fais[0])
        c_av = c_avs[0] + (c_avs[1] - c_avs[0]) / 15 * (fai - fais[0])
    elif abs_fai>30 and abs_fai<45:
        a_av = a_avs[1] + (a_avs[2] - a_avs[1]) / 15 * (fai - fais[1])
        b_av = b_avs[1] + (b_avs[2] - b_avs[1]) / 15 * (fai - fais[1])
        c_av = c_avs[1] + (c_avs[2] - c_avs[1]) / 15 * (fai - fais[1])
    elif abs_fai > 45 and abs_fai < 60:
        a_av = a_avs[2] + (a_avs[3] - a_avs[2]) / 15 * (fai - fais[2])
        b_av = b_avs[2] + (b_avs[3] - b_avs[2]) / 15 * (fai - fais[2])
        c_av = c_avs[2] + (c_avs[3] - c_avs[2]) / 15 * (fai - fais[2])
    elif abs_fai>60 and abs_fai<75:
        a_av = a_avs[3] + (a_avs[4] - a_avs[3]) / 15 * (fai - fais[3])
        b_av = b_avs[3] + (b_avs[4] - b_avs[3]) / 15 * (fai - fais[3])
        c_av = c_avs[3] + (c_avs[4] - c_avs[3]) / 15 * (fai - fais[3])
    else:
        a_av = a_avs[4]
        b_av = b_avs[4]
        c_av = c_avs[4]
    a = a_av
    b = b_av
    c = c_av
    mw = (1 + a / (1 + b / (1 + c))) * 1 / (math.sin(ele) + (a / (math.sin(ele) + b / (math.sin(ele) + c))))
    return mw


# 蒙冠龙UNB3
def UNB3(B,H,D):
    k1=77.604;k2=382000;Rd=284.054;gm=9.784;g=9.80665
    d=math.cos(2*math.pi*(D-28)/365.25)
    if B<=15/180*math.pi:
        P=1013.25;T=299.65;e=26.3;bt=6.3*10**(-3);ld=2.77
    elif 15/180*math.pi<B<30/180*math.pi:
        m=(B-15/180*math.pi)/(15/180*math.pi)
        P=1013.25+(1017.25-1013.25)*m+(-3.75*m)*d
        T=299.65+(294.15-299.65)*m+7*m*d
        e=26.31+(21.79-26.31)*m+8.85*m*d
        bt=6.3*10**(-3)+(6.05-6.3)*10**(-3)*m+0.25*10**(-3)*m*d
        ld=2.77+(3.15-2.77)*m+0.33*m*d
    elif 30/180*math.pi<=B<45/180*math.pi:
        m=(B-30/180*math.pi)/(15/180*math.pi)
        P=1017.25+(1015.25-1017.25)*m+(-3.75+(-2.25+3.75)*m)*d
        T=294.15+(283.15-294.15)*m+(7+(11-7*m)*d)
        e=21.79+(11.66-21.79)*m+(8.85+(7.24-8.85)*m)*d
        bt=(6.05*+(5.58-6.05)*m+(0.25+(0.32-0.25)*m)*d)*10**(-3)
        ld=3.15+(2.57-3.15)*m+(0.33+(0.46-0.33)*m)*d
    elif 45*180*math.pi<=B<60/180*math.pi:
        m=(B-45/180*math.pi)/(15/180*math.pi)
        P=1015.75+(1011.75-1015.75)*m+(-2.25+(-1.75+2.25)*m)*d
        T=283.15+(272.15-283.15)*m+(11+(15-11)*m)*d
        e=(11.66+(6.78-11.66)*m+(7.24+(5.36-7.24)*m)*d)
        bt=(5.58+(5.39-5.58)*m+(0.32+(0.81-0.32)*m)*d)*10**(-3)
        ld=2.57+(1.81-2.57)*m+(0.46+(0.74-0.46)*m)*d
    elif 60/180*math.pi<=B<75/180*math.pi:
        m=(B-60/180*math.pi)/(15/180*math.pi)
        P=1011.75+(1013.00-1011.75)*m+(-1.75+(-0.5+1.75)*m)*d
        T=272.15+(263.65-272.15)*m+(15+(14.5-15)*m)*d
        e=6.78+(4.11-6.78)*m+(5.36+(3.39-5.36)*m)*d
        bt=(5.39+(4.53-5.39)*m+(0.81+(0.62-0.81)*m)*d)*10**(-3)
        ld=1.81+(1.55-1.81)*m+(0.74+(0.3-0.74)*m)*d
    elif B>=75/180*math.pi:
        P=1013-0.5*d
        T=263.65+14.5*d
        e=4.11+3.39*d
        bt=(4.53+0.62*d)*10**(-3)
        ld=1.55+0.3*d
    Zdry=10**(-6)*k1*Rd*P/(gm)
    Zwet=10**(-6)*k2*Rd/(gm*(ld+1)-bt*Rd)*e/T
    Ddry=(1-bt*H/T)**(g/(Rd*bt))*Zdry
    Dwet=(1-bt*H/T)**((ld+1)*g/(Rd*bt)-1)*Zwet
    return([Ddry,Dwet])

def Niell(B,H,E,D):
    d=math.cos(2*math.pi*(D-28)/365.25)
    if B<=15/180*math.pi :
        ahyd=1.2769934*10**(-3);bhyd=2.9153695*10**(-3);chyd=62.610505*10**(-3)
        awet=5.8021897*10**(-4);bwet=1.4275268*10**(-3);cwet=4.3472961*10**(-2)
    elif 15/180*math.pi <B<30/180*math.pi:
        m=(B-15/180*math.pi)/(15/180*math.pi)
        ahyd=1.2769934*10**(-3)+(1.268323-1.2769934)*10**(-3)*m+(1.2709626*10**(-3)*m)*d
        bhyd=2.9153695*10**(-3)+(2.9152299-2.9153695)*10**(-3)*m+(2.1414979*10**(-3)*m)*d
        chyd=62.610505*10**(-3)+(62.837393-62.610505)*10**(-3)*m+(9.01284*10**(-3)*m)*d
        awet=5.8021897*10**(-4)+(5.6494847-5.8021897)*10**(-4)*m
        bwet=1.4275268*10**(-3)+(1.5138625-1.4275268)*10**(-3)*m
        cwet=4.3472961*10**(-2)+(4.6729510-4.3472961)*10**(-2)*m
    elif 30/180*math.pi<=B<45/180*math.pi:
        m=(B-30/180*math.pi)/(15/180*math.pi)
        ahyd=1.268323*10**(-3)+(1.2465397-1.268323)*10**(-3)*m+(1.2709626*10**(-3)+(2.653662*10**(-3)-1.2709626*10**(-3))*m)*d
        bhyd=2.9152299*10**(-3)+(2.9288445-2.9152299)*10**(-3)*m+(2.1414979*10**(-3)+(3.0160779*10**(-3)-2.1414979*10**(-3))*m)*d
        chyd=62.837393*10**(-3)+(63.721774-62.837393)*10**(-3)*m+(9.01284*10**(-3)+(4.3497037*10**(-3)-9.01284*10**(-3))*m)*d
        awet=5.6794847*10**(-4)+(5.8118019-5.6794847)*10**(-4)*m
        bwet=1.5138625*10**(-3)+(1.4572752-1.5138625)*10**(-3)*m
        cwet=4.6729510*10**(-2)+(4.3908931-4.6729510)*10**(-2)*m
    elif 45*180*math.pi<=B<60/180*math.pi:
        m=(B-45/180*math.pi)/(15/180*math.pi)
        ahyd=1.2465397*10**(-3)+(1.2196049-1.2465397)*10**(-3)*m+(2.6523662*10**(-3)+(3.4000452*10**(-3)-2.653662*10**(-3))*m)*d
        bhyd=2.9288445*10**(-3)+(2.9022565-2.9288445)*10**(-3)*m+(3.0160779*10**(-3)+(7.2562722*10**(-3)-3.0160779*10**(-3))*m)*d
        chyd=63.721774*10**(-3)+(63.824265-63.721774)*10**(-3)*m+(4.3497037*10**(-3)+(84.795348*10**(-3)-4.3497037*10**(-3))*m)*d
        awet=5.8118019*10**(-4)+(5.9727542-5.8118019)*10**(-4)*m
        bwet=1.4572752*10**(-3)+(1.5007428-1.4572752)*10**(-3)*m
        cwet=4.3908931*10**(-2)+(4.4626982-4.3908931)*10**(-2)*m
    elif 60/180*math.pi<=B<75/180*math.pi:
        m=(B-60/180*math.pi)/(15/180*math.pi)
        ahyd=1.2196049*10**(-3)+(1.2045996-1.2196049)*10**(-3)*m+(3.4000452*10**(-3)+(4.1202191*10**(-3)-3.4000452*10**(-3))*m)*d
        bhyd=2.9022565*10**(-3)+(2.9024912-2.9022565)*10**(-3)*m+(7.2562722*10**(-3)+(11.723375*10**(-3)-7.2562722*10**(-3))*m)*d
        chyd=63.824265*10**(-3)+(64.258455-63.824265)*10**(-3)*m+(84.795348*10**(-3)+(170.37206*10**(-3)-84.795348*10**(-3))*m)*d
        awet=5.9757542*10**(-4)+(6.1641693-5.9727542)*10**(-4)*m
        bwet=1.5007428*10**(-3)+(1.7599082-1.5007428)*10**(-3)*m
        cwet=4.4626982*10**(-2)+(5.4736038-4.4626982)*10**(-2)*m
    elif B>=75/180*math.pi:
        ahyd=1.2045996*10**(-3)+4.1202191*10**(-3)*d
        bhyd=2.9024912*10**(-3)+11.723375*10**(-3)*d
        chyd=64.258455*10**(-3)+170.37206*10**(-3)*d
        awet=6.1644693*10**(-4);bwet=1.7599082*10**(-3);cwet=5.4736038*10**(-2)
    ahgt=2.53*10**(-5);bhgt=5.49*10**(-3);chgt=1.14*10**(-3)
    Mhyd=(1+ahyd/(1+bhyd/(1+chyd)))/(math.sin(E)+ahyd/(math.sin(E)+bhyd/(math.sin(E)+chyd)))+(1/math.sin(E)-(1+ahgt/(1+bhgt/(1+chgt)))/(math.sin(E)+ahgt/(math.sin(E)+bhgt/(math.sin(E)+chgt))))*(H/1000)
    Mwet=(1+awet/(1+bwet/(1+cwet)))/(math.sin(E)+awet/(math.sin(E)+bwet/(math.sin(E)+cwet)))
    return ([Mhyd,Mwet])