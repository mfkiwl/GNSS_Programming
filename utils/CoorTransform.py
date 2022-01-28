# -*- coding: utf-8 -*-
"""

@title:	Style Guide for Python Code
@author: iDeal0103
@status:	Active
@type:	Process
@created:	13-Apr-2021
@post-History:	13-Apr-2021

comment:
    1.XYZ与BLH转换
    2.BLH与NEU转换
    
"""


#import module
import math
from math import sin
from math import cos
import numpy as np

#def cal_XYZ2BLH(X,Y,Z,e=0.08181333402,a=6378245.0):
def cal_XYZ2BLH(X,Y,Z,e=0.08181919084,a=6378137.0): 
    L=math.atan2(Y, X)
    dB=100
    B=math.atan2(Z, math.sqrt(X**2+Y**2))
    while abs(dB)>4.848e-11:
        N=a/math.sqrt(1-e**2*math.sin(B)**2)
        B2=math.atan2(Z+N*e**2*math.sin(B), math.sqrt(X**2+Y**2))
        dB=B2-B
        B=B2
    N=a/math.sqrt(1-e**2*math.sin(B)**2)
    H=Z/math.sin(B)-N*(1-e**2)
    #H=math.sqrt(X**2+Y**2)/math.cos(B)-N
    return B, L, H


'''
由设站点XYZ坐标和目标点XYZ坐标
计算目标点在设站点的站心地平坐标系中的坐标(NEU)
'''
def cal_NEU(stationcenter_coor,object_coor):
    Xs,Ys,Zs=stationcenter_coor
    B,L,H=cal_XYZ2BLH(Xs,Ys,Zs)
    R=np.array([[-sin(B)*cos(L),-sin(L),cos(B)*cos(L)],[-sin(B)*sin(L),cos(L),cos(B)*sin(L)],[cos(B),0,sin(B)]])
    dxyz=np.array(object_coor)-np.array(stationcenter_coor)
    NEU=np.linalg.inv(R)@dxyz
    N=NEU[0]
    E=NEU[1]
    U=NEU[2]
    return N,E,U


'''
由设站点XYZ坐标和目标点XYZ坐标
计算目标点在设站点的站心地平坐标系中的高度角
返回的角度值均为弧度
'''
def cal_ele_and_A(stationcenter_coor,object_coor):
    Xr, Yr, Zr = stationcenter_coor
    B, L, H = cal_XYZ2BLH(Xr, Yr, Zr)
    R = np.array([[-sin(B)*cos(L), -sin(L), cos(B)*cos(L)], [-sin(B)*sin(L),cos(L),cos(B)*sin(L)], [cos(B),0,sin(B)]]).astype(float)
    dxyz = np.array(object_coor)-np.array(stationcenter_coor)
    NEU = np.linalg.inv(R)@dxyz
    N = NEU[0]
    E = NEU[1]
    U = NEU[2]
    ele = math.atan2(U, math.sqrt(N**2+E**2))
    # 根据所在象限计算方位角(依次为一、二、三、四象限)
    if N > 0 and E > 0:
        A = math.atan2(E, N)
    elif N > 0 and E < 0:
        A = math.atan2(E, N)  # 此处A为负
        A = 2 * math.pi + A
    elif N < 0 and E < 0:   # 此处A为正
        A = math.atan2(E, N)
        A = math.pi - A
    elif N < 0 and E > 0 :
        A = math.atan2(E, N)   # 此处A为负
        A=math.pi + A
    return ele, A


#计算地球自转改正项
def earth_rotation_correction(coor, dt, w=7.2921151467e-5):
    """
    Parameters
    ----------
        coor : [Xecef,Yecef,Zecef]
        w : 地球自转的角速度
        dt : 信号从卫星传播到观测站的时间
    Returns
    -------
        Xeci,Yeci,Zeci : 发射信号时刻卫星,在观测站接收到信号时刻的ECEF坐标系中的坐标

    """
    coor_ecef = np.array(coor)
    R = np.array([[math.cos(w*dt), math.sin(w*dt), 0], [-math.sin(w*dt), math.cos(w*dt), 0], [0, 0, 1]])
    corrected_coor = R@coor_ecef
    Xeci = corrected_coor[0]
    Yeci = corrected_coor[1]
    Zeci = corrected_coor[2]
    return Xeci, Yeci, Zeci


# 计算两点之间距离
def cal_distance(coor1,coor2):
    """
    Parameters
    ----------
        coor1/coor2 : coordinate_list , [X,Y,Z]
    Returns
    -------
        distance : float
    """
    X1,Y1,Z1=coor1
    X2,Y2,Z2=coor2
    distance=math.sqrt((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2)
    return distance


#计算某大地直角坐标系坐标处,站心地平坐标系到站心赤道坐标系的旋转矩阵
def cal_R_THS2TES(stationcenter_coor):
    """
    stationcenter_coor : [X,Y,Z]
    """
    Xs,Ys,Zs=stationcenter_coor
    B,L,H=cal_XYZ2BLH(Xs,Ys,Zs)
    R=np.array([[-sin(B)*cos(L),-sin(L),cos(B)*cos(L)],[-sin(B)*sin(L),cos(L),cos(B)*sin(L)],[cos(B),0,sin(B)]])
    return R




