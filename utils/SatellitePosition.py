# -*- coding: utf-8 -*-
"""

@title:	Style Guide for Python Code
@author: iDeal0103
@status:	Active
@type:	Process
@created:	13-Apr-2021
@post-History:	13-Apr-2021

comment:
    1.计算某时刻某卫星坐标
    2.计算某时刻某卫星钟差
    
"""


#import
import utils.const as const
import utils.TimeSystem as TimeSystem
import math
import utils.DoFile as DoFile
import datetime
import utils.RecordFilter as RecordFilter
import numpy as np


# 以某条记录为基础,计算卫星位置
def cal_SatellitePosition_GPS_datetime(time, serial_no, brs):
    """
    Parameters
    ----------
        UTCtime : datetime.datetime,所求时刻的datetime格式的GPS时间
        brs : list[GPS_brdc_record class],所依据的广播星历记录
    Returns
    -------
        X,Y,Z ：卫星的三维坐标(单位为m)
    """ 
    #筛选出最接近的一条记录
    br = RecordFilter.find_closest_record(brs, time, serial_no)
    # (1)Time from ephemeris epoch toe
    tk = TimeSystem.from_datetime_cal_GPSws(time)[1] - br.toe
    # (2)Calculate the semimajor axis
    a = br.sqrt_a**2
    # (3)Compute mean motion – rad/s
    n0 = math.sqrt(const.miu/a**3)
    # (4)Correct mean motion
    n = n0+br.delta_n
    # (5)Mean anomaly Mk at epoch t
    Mk = br.M0+n*tk
    # (6)Eccentricity anomaly Ek at epoch t
    Ek = cal_Ek(Mk,br.e)
    # (7)True anomaly vk at epoch t
    vk = 2*math.atan(math.sqrt(1-br.e**2)/(1-br.e)*math.tan(Ek/2))
    # (8)argument of latitude
    uk = vk+br.w
    # (9)Corrections for second harmonic perturbations
    delta_uk = br.Cuc*math.cos(2*uk)+br.Cus*math.sin(2*uk)
    delta_rk = br.Crc*math.cos(2*uk)+br.Crs*math.sin(2*uk)
    delta_ik = br.Cic*math.cos(2*uk)+br.Cis*math.sin(2*uk)
    # (10) Corrected argument of latitude, radius and inclination
    u = uk+delta_uk
    r = a*(1-br.e*math.cos(Ek))+delta_rk
    i = br.i0+delta_ik+br.i_dot*tk
    # (11) Satellites’ positions in orbital plane
    x = r*math.cos(u)
    y = r*math.sin(u)
    # (12) Corrected longitude of ascending node at epoch t
    lamb = br.omega0+(br.omega_dot-const.we)*tk-const.we*br.toe
    # (13) Satellites coordinates at ECEF system
    poscoor_ECEF_X = r*(math.cos(u)*math.cos(lamb)-math.sin(u)*math.cos(i)*math.sin(lamb))
    poscoor_ECEF_Y = r*(math.cos(u)*math.sin(lamb)+math.sin(u)*math.cos(i)*math.cos(lamb))
    poscoor_ECEF_Z = r*math.sin(u)*math.sin(i)
    #输出坐标单位为m
    return poscoor_ECEF_X,poscoor_ECEF_Y,poscoor_ECEF_Z

# 以某条记录为基础,计算卫星钟差
def cal_ClockError_GPS_datetime(time,serial_no,brs):
    """
    Parameters
    ----------
        time : datetime.datetime,所求时刻的datetime格式的GPS时间
        brs : list[GPS_brdc_record class],所依据的广播星历记录
    Returns
    -------
        clockerror ：卫星的钟差,单位s
    """ 
    #筛选出最接近的一条记录
    br=RecordFilter.find_closest_record(brs,time,serial_no)   
    print(time,br.toc)
    tc= TimeSystem.from_datetime_cal_GPSws(time)[1] - br.toe
    print(tc, TimeSystem.from_datetime_cal_GPSws(time)[1], br.toe)
    clockerror=br.a0+br.a1*tc+br.a2*tc**2
    return clockerror



#卫星位置计算中求解Ek所用函数
def cal_Ek(Mk,e):
    Ek1 = Mk
    Ek0 = 0.0
    while (abs(Ek1 - Ek0)>1.0e-12):
        Ek0 = Ek1
        Ek1 = Mk + e*math.sin(Ek0)
    Ek = Ek1
    return Ek

# 以某条记录为基础,计算卫星位置
def cal_SatellitePosition_GPS_GPSws(time,serial_no,brs):
    """
    Parameters
    ----------
        time : GPSws,所求时刻的GPSws类的GPS时间
        brs : lsit[GPS_brdc_record class],所依据的广播星历记录
    Returns
    -------
        X,Y,Z ：卫星的三维坐标(单位为m)
        dt : 卫星的钟差
    """ 
    # 筛选出最接近的一条记录
    datetime_time = TimeSystem.from_GPSws_cal_datetime_2(time)
    br = RecordFilter.find_closest_record(brs, datetime_time, serial_no)
    # (1)Time from ephemeris epoch toe
    tk = time.GpsSecond-br.toe
    if abs(tk) > 7200:
        print('效果可能较差!')
    if tk > 302400:
        tk = tk - 604800
    elif tk < -302400:
        tk = tk + 604800
    # (2)Calculate the semimajor axis
    a = br.sqrt_a**2
    # (3)Compute mean motion – rad/s
    n0 = math.sqrt(const.miu/a**3)
    # (4)Correct mean motion
    n = n0 + br.delta_n
    # (5)Mean anomaly Mk at epoch t
    Mk = br.M0+n*tk
    # (6)Eccentricity anomaly Ek at epoch t
    Ek = cal_Ek(Mk, br.e)
    # (7)True anomaly vk at epoch t
    vk = 2*math.atan(math.sqrt((1+br.e)/(1-br.e))*math.tan(Ek/2))
    # (8)argument of latitude
    uk = vk+br.w
    # (9)Corrections for second harmonic perturbations
    delta_uk = br.Cuc*math.cos(2*uk)+br.Cus*math.sin(2*uk)
    delta_rk = br.Crc*math.cos(2*uk)+br.Crs*math.sin(2*uk)
    delta_ik = br.Cic*math.cos(2*uk)+br.Cis*math.sin(2*uk)
    # (10) Corrected argument of latitude, radius and inclination
    u = uk + delta_uk
    r = a*(1-br.e*math.cos(Ek))+delta_rk
    i = br.i0 + delta_ik + br.i_dot * tk
    # (11) Satellites’ positions in orbital plane
    x = r*math.cos(u)
    y = r*math.sin(u)
    # (12) Corrected longitude of ascending node at epoch t
    lamb = br.omega0+(br.omega_dot-const.we)*tk-const.we*br.toe
    # (13) Satellites coordinates at ECEF system
    poscoor_ECEF_X = r*(math.cos(u)*math.cos(lamb)-math.sin(u)*math.cos(i)*math.sin(lamb))
    poscoor_ECEF_Y = r*(math.cos(u)*math.sin(lamb)+math.sin(u)*math.cos(i)*math.cos(lamb))
    poscoor_ECEF_Z = r*math.sin(u)*math.sin(i)
    # 坐标单位化为m
    return poscoor_ECEF_X, poscoor_ECEF_Y, poscoor_ECEF_Z

# 以某条记录为基础,计算BDS卫星位置
def cal_SatellitePosition_BDS_GPSws(time,serial_no,brs):
    """
    Parameters
    ----------
        time : GPSws,所求时刻的GPSws类的GPS时间
        brs : list[Renix304_navigation_record_BDS class],所依据的广播星历记录
    Returns
    -------
        X,Y,Z ：卫星在地心地固坐标系下的的三维坐标(单位为m)
    """
    # 筛选出最接近的一条记录
    datetime_time = TimeSystem.from_GPSws_cal_datetime_2(time)
    BDSws = TimeSystem.from_GPSws_get_BDSws(time)
    br = RecordFilter.find_closest_record(brs, datetime_time, serial_no)
    # (1)Time from ephemeris epoch toe
    tk = BDSws.BDSSecond-br.toe
    if abs(tk) > 7200:
        print('效果可能较差!')
    if tk > 302400:
        tk = tk - 604800
    elif tk < -302400:
        tk = tk + 604800
    # (2)Calculate the semimajor axis
    a = br.sqrt_a**2
    # (3)Compute mean motion – rad/s
    n0 = math.sqrt(const.miu/a**3)
    # (4)Correct mean motion
    n = n0 + br.delta_n
    # (5)Mean anomaly Mk at epoch t
    Mk = br.M0+n*tk
    # (6)Eccentricity anomaly Ek at epoch t
    Ek = cal_Ek(Mk, br.e)
    # (7)True anomaly vk at epoch t
    vk = 2*math.atan(math.sqrt((1+br.e)/(1-br.e))*math.tan(Ek/2))
    # (8)argument of latitude
    uk = vk+br.w
    # (9)Corrections for second harmonic perturbations
    delta_uk = br.Cuc*math.cos(2*uk)+br.Cus*math.sin(2*uk)
    delta_rk = br.Crc*math.cos(2*uk)+br.Crs*math.sin(2*uk)
    delta_ik = br.Cic*math.cos(2*uk)+br.Cis*math.sin(2*uk)
    # (10) Corrected argument of latitude, radius and inclination
    u = uk + delta_uk
    r = a*(1-br.e*math.cos(Ek))+delta_rk
    i = br.i0 + delta_ik + br.i_dot * tk
    # (11) Satellites’ positions in orbital plane
    x = r*math.cos(u)
    y = r*math.sin(u)
    if br.SVN in ["C01", "C02", "C03", "C04", "C05", "C59", "C60"]:  # GEO卫星的计算
        ff1 = const.we * tk
        ff2 = -5 / 180 * math.pi
        lamb = br.omega0+br.omega_dot*tk-const.we*br.toe
        Xgk = (x * math.cos(lamb) - y * math.cos(i) * math.sin(lamb))
        Ygk = (x * math.sin(lamb) + y * math.cos(i) * math.cos(lamb))
        Zgk = y * math.sin(i)
        Rx = ([1, 0, 0], [0, math.cos(ff2), math.sin(ff2)], [0, -math.sin(ff2), math.cos(ff2)])
        Rz = ([math.cos(ff1), math.sin(ff1), 0], [-math.sin(ff1), math.cos(ff1), 0], [0, 0, 1])
        Rx = np.mat(Rx)
        Rz = np.mat(Rz)
        a = Rz * Rx * np.mat([Xgk, Ygk, Zgk]).T
        poscoor_ECEF_X = float(a[0])
        poscoor_ECEF_Y = float(a[1])
        poscoor_ECEF_Z = float(a[2])
    else:
        # (12) Corrected longitude of ascending node at epoch t
        lamb = br.omega0+(br.omega_dot-const.we)*tk-const.we*br.toe
        # (13) Satellites coordinates at ECEF system
        poscoor_ECEF_X = r*(math.cos(u)*math.cos(lamb)-math.sin(u)*math.cos(i)*math.sin(lamb))
        poscoor_ECEF_Y = r*(math.cos(u)*math.sin(lamb)+math.sin(u)*math.cos(i)*math.cos(lamb))
        poscoor_ECEF_Z = r*math.sin(u)*math.sin(i)
    # 坐标单位化为m
    return poscoor_ECEF_X, poscoor_ECEF_Y, poscoor_ECEF_Z

# 以某条记录(GPSweek和GPSsecond)为基础,计算卫星钟差
def cal_ClockError_GPSws(time, SVN, brs):
    """
    Parameters
    ----------
        time : GPSws,所求时刻的GPSws类的GPS时间
        SVN : str,卫星的SVN
        brs : list[GPS_brdc_record class],所依据的广播星历记录
    Returns
    -------
        clockerror : 卫星钟差,单位s
    """ 
    # 筛选出最接近的一条记录
    datetime_time = TimeSystem.from_GPSws_cal_datetime_2(time)
    br = RecordFilter.find_closest_record(brs, datetime_time, SVN)
    tc = time.GpsSecond-br.toe
    clockerror = br.a0+br.a1*tc+br.a2*tc**2
    return clockerror

# 以某条记录(GPSweek和GPSsecond)为基础,计算卫星钟差(包含相对论效应)
def cal_ClockError_GPSws_withRelativisticEffect(time, SVN, brs):
    """
    Parameters
    ----------
        time : GPSws,所求时刻的GPSws类的GPS时间
        SVN : str,卫星的SVN号
        brs : list[GPS_brdc_record class],所依据的广播星历记录
    Returns
    -------
        clockerror : 卫星钟差,单位s
    """ 
    # 筛选出最接近的一条记录
    # 由GPSws类数据得到datetime.datetime类型时间,便于筛选数据
    datetime_time = TimeSystem.from_GPSws_cal_datetime_2(time)
    br = RecordFilter.find_closest_record(brs, datetime_time, SVN)
    # (1)Time from ephemeris epoch toe
    tk = time.GpsSecond-br.toe
    # (2)Calculate the semimajor axis
    a = br.sqrt_a**2
    # (3)Compute mean motion – rad/s
    n0 = math.sqrt(const.miu/a**3)
    # (4)Correct mean motion
    n = n0+br.delta_n
    # (5)Mean anomaly Mk at epoch t
    Mk = br.M0+n*tk
    # (6)Eccentricity anomaly Ek at epoch t
    Ek = cal_Ek(Mk, br.e)
    # 计算钟差偏移和漂移部分
    clockerror_biasdrift = br.a0+br.a1*tk+br.a2*tk**2
    # 计算钟差相对论效应部分
    clockerror_re = -4.443e-10*br.e*math.sqrt(a)*math.sin(Ek)
    # 合并钟差
    clockerror = clockerror_biasdrift+clockerror_re
    return clockerror

if __name__ == "__main__":
    pass



