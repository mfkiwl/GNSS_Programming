# -*- coding: utf-8 -*-
"""

@name:EphemerisInterpolation
@title:	Style Guide for Python Code
@author: iDeal0103
@status:	Active
@type:	Process
@created:	14-Apr-2021
@post-History:	14-Apr-2021

comment：
    1.精密星历插值(拉格朗日插值计算坐标,线性插值计算钟差)
    2.

"""

# import 
import datetime


# 便于插值中的时间差表示
def cal_delta_time(time1,time2):
    if time1>time2:
        delta_time=(time1-time2).seconds
    else:
        delta_time=-(time2-time1).seconds
    return delta_time


#拉格朗日插值计算X、Y、Z
def Lagrange_interpolation(time_tobe_cal,cal_based_records):
    X_result=0
    Y_result=0
    Z_result=0
    for i in range(len(cal_based_records)):
        part_result=1
        for j in list(range(i))+list(range(i+1,len(cal_based_records))):
            part_result*=cal_delta_time(time_tobe_cal,cal_based_records[j].time)/cal_delta_time(cal_based_records[i].time,cal_based_records[j].time)
        X_result+=part_result*cal_based_records[i].X
        Y_result+=part_result*cal_based_records[i].Y
        Z_result+=part_result*cal_based_records[i].Z
    return X_result,Y_result,Z_result



#线性插值计算时间
def linear_interpolation_emp(time_tobe_cal,cal_based_records):
    a=cal_delta_time(cal_based_records[0].time,time_tobe_cal)
    b=cal_delta_time(time_tobe_cal,cal_based_records[1].time)
    c=cal_delta_time(cal_based_records[0].time, cal_based_records[1].time)
    result=(a*cal_based_records[0].dT+b*cal_based_records[1].dT)/c
    return result


# 基于某精密卫星星历的数据部分记录,插值计算卫星位置及钟差
def cal_XYZ_basedon_preciseephemeris(pe_records,the_time,the_prn,lagran_n=9):
    """
    Parameters
    ----------
        pe_records: satelliteobit_XYZdT class , 从精密星历中解析得到的所有satelliteobit_XYZdT类型的记录
        the_time : datetime.datetime , 所要计算位置及钟差的时间
        the_prn : str , 卫星编号
        lagran_n : int , 拉格朗日插值阶数.
    Returns
    -------
        X : float
        Y : float
        Z : float
            依次对应卫星X、Y、Z坐标
    """
    #筛选对应prn编号的卫星记录
    prn_records=list(filter(lambda o:o.PRN==the_prn,pe_records))
    #获得拉格朗日插值所需要的最临近n个记录
    prn_records.sort(key=lambda o:abs(cal_delta_time(the_time,o.time)))
    the_records=prn_records[:lagran_n]
    # print(the_records)
    #进行插值计算
    X,Y,Z=Lagrange_interpolation(the_time,the_records)
    return X,Y,Z


# 基于某精密卫星星历的数据部分记录,插值计算卫星钟差
def cal_dT_basedon_preciseephemeris(pe_records,the_time,the_prn):
    """
    Parameters
    ----------
        pe_records: satelliteobit_XYZdT class , 从精密星历中解析得到的所有satelliteobit_XYZdT类型的记录
        the_time : datetime.datetime , 所要计算位置及钟差的时间
        the_prn : str , 卫星编号
    Returns
    -------
        dT : float,卫星钟差
    """
    #筛选对应prn编号的卫星记录
    prn_records=list(filter(lambda o:o.PRN==the_prn,pe_records))
    #获得所需要的最临近2个记录
    prn_records.sort(key=lambda o:abs(cal_delta_time(the_time,o.time)))
    the_records=prn_records[:2]
    #进行插值计算
    dT=linear_interpolation_emp(the_time,the_records)
    return dT


#用精密钟差线性插值计算时间
def linear_interpolation_clk(time_tobe_cal,cal_based_records):
    a=cal_delta_time(cal_based_records[0].time,time_tobe_cal)
    b=cal_delta_time(time_tobe_cal,cal_based_records[1].time)
    c=cal_delta_time(cal_based_records[0].time, cal_based_records[1].time)
    result=(a*cal_based_records[0].data[0]+b*cal_based_records[1].data[0])/c
    return result


# 基于某精密钟差的数据部分记录,插值计算卫星钟差
def cal_dT_basedon_clk(clk_satellite_records,the_time,the_prn):
    """
    Parameters
    ----------
        clk_satellite_records:clk_satellite_record class , 从精密钟差中解析得到的所有clk_satellite_record类型的记录
        the_time : datetime.datetime , 所要计算位置及钟差的时间
        the_prn : str , 卫星编号
    Returns
    -------
        dT : float,卫星钟差
    """
    #筛选对应prn编号的卫星记录
    prn_records=list(filter(lambda o:o.PRN==the_prn,clk_satellite_records))
    #获得所需要的最临近2个记录
    prn_records.sort(key=lambda o:abs(cal_delta_time(the_time,o.time)))
    the_records=prn_records[:2]
    #进行插值计算
    dT=linear_interpolation_clk(the_time,the_records)
    return dT



# import utils.DoFile as DoFile
# igs_file=r"E:\大三下\卫星与导航定位\代码集合\Satellite_Navigation_and_Positioning\data\sat_obit\igs21304.sp3"
# pe_records=DoFile.read_GPS_sp3File(igs_file)
# the_time=datetime.datetime(2020,11,5 ,12,45,0)
# cal_XYZdT_basedon_preciseephemeris(pe_records,the_time,"01")





