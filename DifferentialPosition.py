# -*- coding: utf-8 -*-
"""

@title:	Style Guide for Python Code
@author: iDeal0103
@status:	Active
@type:	Process
@created:	31-Jan-2022
@post-History:	09-Jan-2021

comment：
    1.单差伪距定位 (Single Differential Position)
    2.双差载波相位定位 (Double Differential Poisition)

"""

# import module
import numpy as np
import SinglePointPosition as SPP
import math
import datetime
import utils.DoFile as DoFile
import utils.SatellitePosition as SatellitePosition
import utils.TimeSystem as TimeSystem
import utils.CoorTransform as CoorTransform
from utils.ErrorReduction import *
import matplotlib.pyplot as plt

# 基于伪距的站间单差
def SD_onPseudorange(knownStation_ob_records, unknownStation_ob_records, br_records,
                      Tr, knownStation_coor, CRC=True, c=299792458, init_coor=[2000, 2000, 2000]):
    """
    knownStation_ob_records : list[GPS_observation_record class] , 所使用的已知站观测文件记录
    unknownStation_ob_records : list[GPS_observation_record class] , 所使用的待测站观测文件记录
    br_records :  list[GPS_brdc_record class] , 所使用的卫星广播星历记录
    Tr : 接收机接收到信号的时刻,GPS时刻
    knownStation_coor : list ,已知观测站坐标
    CRC : bool , 是否进行相对论钟差改正
    c : const , 光速(单位为m/s)
    init_coor : list , 观测站坐标初值
    """
    # 筛选出某时刻的记录
    the_knownStation_ob_record = list(filter(lambda o:o.system=="G" and o.time==Tr and o.data!="",knownStation_ob_records))
    the_unknownStation_ob_record = list(filter(lambda o:o.system=="G" and o.time==Tr and o.data!="",unknownStation_ob_records))
    # 开始迭代前knownStation与unknownStation接收机钟差dtr_k和dtr_u为0
    # dtr_k=0
    # dtr_u=0
    dtr_uk = 0
    # 初始地面点坐标
    Xu,Yu,Zu = init_coor
    Xk,Yk,Zk = knownStation_coor
    # 初始对流层改正和电离层改正
    # dtTD_uk=0
    # dtID_uk=0
    # 电离层改正系数
    f1 = 1575.42    #Hz
    f2 = 1227.60    #Hz
    miuf = f1**2/f2**2
    # 平差迭代次数计数
    no = 0
    # 迭代平差求解差分定位模型
    while True:
        # 如果超出平差迭代求解超出8次则跳出
        if no > 8:
            break
        no += 1
        A = []
        l = []
        Ps = []
        for unknownStation_record in the_unknownStation_ob_record:
            # 对the_knownStation_ob_record进行同一颗卫星的筛选
            knownStation_record_base = list(filter(lambda o:o.PRN==unknownStation_record.PRN,the_knownStation_ob_record))
            if len(knownStation_record_base) == 1:
                knownStation_record = knownStation_record_base[0]
            else:
                knownStation_record = ""
            # 如果所选观测值数据为空则跳过
            if unknownStation_record.data['P1']['observation'] == "" \
                    or unknownStation_record.data['P2']['observation'] == "" \
                    or knownStation_record == "" or knownStation_record.data['P1']['observation'] == "" \
                    or knownStation_record.data['P2']['observation'] == "":
                continue
            # 不加电离层改正
            # P_u=unknownStation_record.data['P1']['observation']
            # P_k=knownStation_record.data['P1']['observation']
            #加入电离层改正
            P_u = Ionospheric_Delay_Correction(unknownStation_record)
            P_k = Ionospheric_Delay_Correction(knownStation_record)
            the_prn = unknownStation_record.PRN
            ## 根据接收时间,计算信号发射时间及此时卫星所在位置
            # 将datetime格式时间转换为精度更高的自定义GPSws时间类
            w, s = TimeSystem.from_datetime_cal_GPStime_2(Tr)
            Tr_GPSws = TimeSystem.GPSws(w, s)
            # tr=Tr_GPSws.minus(dtr_u)     # 此处有对接收机钟差的改正
            tr=Tr_GPSws      # 不考虑对接收机钟差的改正
            Ts=tr.cal_minus_result(P_u/c)
            # Ts = tr.minus(P_u/c)
            ts = Ts
            t_del = 10
            dts1 = 0
            # 迭代计算发射时间
            while abs(t_del) > 1e-8:
                # 计算卫星的位置和钟差(单位为m)
                coorX,coorY,coorZ=SatellitePosition.cal_SatellitePosition_GPS_GPSws(ts, the_prn, br_records)
                if CRC:
                    # 加入相对论钟差改正
                    dts2 = SatellitePosition.cal_ClockError_GPS_GPSws_withRelativisticEffect(ts, the_prn, br_records)
                else:
                    dts2 = SatellitePosition.cal_ClockError_GPS_GPSws(ts,the_prn,br_records)
                t_del = dts2-dts1
                dts1 = dts2
                # 此处计算消除卫星钟差后的时间
                ts = Ts.cal_minus_result(dts2)
                # ts = Ts.minus(dts2)
            # 修正地球自转影响
            dt = TimeSystem.cal_deltatime_second_GPSws(tr, ts)
            Xeci, Yeci, Zeci=CoorTransform.earth_rotation_correction([coorX, coorY, coorZ], dt)
            # 对流层延迟改正
            P_u = Tropospheric_Delay_Correction_Hopfield(P_u, [Xu, Yu, Zu], [Xeci, Yeci, Zeci])
            P_k = Tropospheric_Delay_Correction_Hopfield(P_k, [Xk, Yk, Zk], [Xeci, Yeci, Zeci])
            # 计算系数阵成分
            lou_u = CoorTransform.cal_distance([Xu, Yu, Zu], [Xeci, Yeci, Zeci])
            lou_k = CoorTransform.cal_distance([Xk, Yk, Zk], [Xeci, Yeci, Zeci])
            axki = (Xu-Xeci)/lou_u
            ayki = (Yu-Yeci)/lou_u
            azki = (Zu-Zeci)/lou_u
            # A.append([-axki,-ayki,-azki,1,-1,1,miuf])
            A.append([-axki, -ayki, -azki, 1])
            # 计算常数阵成分
            # li=(P_k-P_u)-(lou_k-lou_u)+c*(dtr_k-dtr_u)+dtTD_uk+miuf*dtID_uk
            li = (P_k-P_u)-(lou_k-lou_u)-c*dtr_uk
            l.append(li)
            # 更新权阵成分
            Ps.append(1)
        # 求解
        if Ps == []:
            continue
        Pz = np.diag(Ps)
        A = np.array(A)
        l = np.array(l)
        # 改正数发散太过严重则不再继续平差
        if abs(max(l.tolist())) > 1e10:
            break
        x = np.linalg.inv(A.T@Pz@A)@(A.T@Pz@l)
        # print(x)
        # 更新参数
        dXu = x[0]
        dYu = x[1]
        dZu = x[2]
        ddtr_uk = x[3]/c
        # ddtr_k = x[3]/c
        # ddtr_u = x[4]/c
        # ddtTD_uk = x[5]
        # ddtID_uk = x[6]
        Xu += dXu
        Yu += dYu
        Zu += dZu
        dtr_uk += ddtr_uk
        # dtr_k += ddtr_k
        # dtr_u += ddtr_u
        # dtTD_uk += ddtTD_uk
        # dtID_uk += ddtID_uk
        # 判断迭代停止条件
        if abs(dXu) < 1e-4 and abs(dYu) < 1e-4 and abs(dZu) < 1e-4:
            break
    return Xu, Yu, Zu
    

    
# 主程序
if __name__ == "__main__":
    unknownStation_observation_file = r"edata\obs\warn3100.20o"
    knownStation_observation_file = r"edata\obs\leij3100.20o"
    broadcast_file = r"edata\sat_obit\brdc3100.20n"
    # 读入观测文件内容,得到类型对象列表
    unknownStation_ob_records = DoFile.read_GPS_oFile(unknownStation_observation_file)
    knownStation_ob_records = DoFile.read_GPS_oFile(knownStation_observation_file)
    br_records = DoFile.read_GPS_nFile(broadcast_file)
    Tr = datetime.datetime(2020, 11, 5, 0, 0, 0)
    init_coor = [3658785.6000, 784471.1000, 5147870.7000]
    knownStation_coor = [0.389873613453103E+07, 0.855345521080705E+06, 0.495837257579542E+07]
    true_coors = []
    cal_coors = []
    while Tr < datetime.datetime(2020, 11, 5, 23, 45, 16):
        Xk, Yk, Zk=SD_onPseudorange(knownStation_ob_records, unknownStation_ob_records, br_records, Tr, knownStation_coor)
        cal_coors.append([Xk, Yk, Zk])
        true_coors.append([0.365878555276965E+07, 0.784471127238666E+06, 0.514787071062059E+07])  #warn
        #true_coors.append([0.389873613453103E+07,0.855345521080705E+06,0.495837257579542E+07])   #leij
        #true_coors.append([-0.267442768572702E+07,0.375714305701559E+07,0.439152148514515E+07])  #chan
        Tr+=datetime.timedelta(seconds=30)
    SPP.cal_NEUerrors(true_coors, cal_coors)
    




