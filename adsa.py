
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
from utils.const import *
import utils.LAMBDA as LAMBDA


def observation_isnot_null(station_record, FreqBandList=['L1', 'P1']):
    '''
    station_record :  list[GPS_observation_record class] , 所使用观测文件读取的各条记录
    FreqBandList : list[str] , 检查数据是否为空的频率波段 (如 'C1','L2','P1' )
    '''
    isnot_null_flag = True
    for band in FreqBandList:
        # 判断L1波段的C1数据是否齐全
        if station_record.data[band]['observation'] == "" or station_record.data[band]['observation'] == "":
            isnot_null_flag = False
        else:
            isnot_null_flag = True
    return isnot_null_flag

def cal_EmitTime_from_datetime(Tr, dtr, the_prn, P, br_records, doCRC, c=299792458):
    '''

    :param Tr: datetime.datetime , 观测站观测的时刻
    :param the_prn: str , 观测到卫星的PRN , "xx"
    :param dtr:  float , 接收机钟差
    :param P: 观测所得到的伪距(也可以是改正过的)
    :param br_records: list[GPS_brdc_record class] , 所使用的卫星广播星历记录
    :param doCRC: bool , 是否进行相对论效应改正
    :param c: float , 光速(单位 m/s)
    :return: tr: GPSws class , 得到接收机的信号接收时刻(已考虑接收机钟差)
    :return: ts: GPSws class , 得到卫星的信号发射时刻(已考虑卫星钟差)
    :return: dts: float , 得到卫星信号发射时刻对应的钟差(单位：s)
    '''
    # 将datetime格式时间转换为精度更高的自定义GPSws时间类
    w, s = TimeSystem.from_datetime_cal_GPStime_2(Tr)
    Tr_GPSws = TimeSystem.GPSws(w, s)
    # 减去接收机钟差
    tr = Tr_GPSws.cal_minus_result(dtr)
    # tr=Tr_GPSws   #不考虑接收机钟差
    Ts = tr.cal_minus_result(float(P) / c)
    # ts为考虑钟差的卫星发射时间
    ts = Ts
    t_del = 10
    dts1 = 0
    # 迭代计算发射时间
    while abs(t_del) > 1e-8:
        # 计算卫星的钟差(单位为s)
        if doCRC:
            # 加入相对论钟差改正
            dts2 = SatellitePosition.cal_ClockError_GPS_GPSws(ts, the_prn, br_records)
        else:
            dts2 = SatellitePosition.cal_ClockError_GPS_GPSws_withRelativisticEffect(ts, the_prn, br_records)
        t_del = dts2 - dts1
        dts1 = dts2
        # 此处计算消除卫星钟差后的时间
        ts = Ts.cal_minus_result(dts2)
    return tr, ts, dts2

def make_Ambiguity_coefficient_matrix_row(no, total_conut, lamb):
    efficient_row = []
    for i in range(total_conut):
        if i == (no-1) :
            efficient_row.append(1 * lamb)
        else:
            efficient_row.append(0)
    return efficient_row

def diagonalize_squarematrix(a, b):
    # 构造合并后大小的矩阵
    a_row = a.shape[0]
    a_col = a.shape[1]
    b_row = b.shape[0]
    b_col = b.shape[1]
    c = np.zeros((a_row + b_row, a_col + b_col))
    c[:a_row, :a_col] = a
    c[a_row:, a_col:] = b
    return c

def get_DD_Pmatrix(nDD, sigma=1):
    Pmatrix = np.full((nDD, nDD), -1).astype(float)
    for i in range(nDD):
        Pmatrix[i, i] = nDD
    constparam = 1/ ((2*sigma**2) * (nDD+1))
    Pmatrix *= constparam
    return Pmatrix

def MAPmethod(b_hat, a_hat, Qaa, Qba, a_check):
    b_check = b_hat - Qba @ np.linalg.inv(Qaa) @ (a_hat - a_check)
    return b_check



# 基于载波相位的双差(两个站皆为未知站点)
def DD_onCarrierPhase_2unknown(station1_ob_records, station2_ob_records, br_records,
                      Tr1, Tr2, CRC=True, c=299792458,
                      station1_init_coor=[2000, 2000, 2000], station2_init_coor=[2000, 2000, 2000]):
    """
    station1_ob_records : list[GPS_observation_record class] , 所使用的观测站1观测文件记录
    station2_ob_records : list[GPS_observation_record class] , 所使用的观测站2观测文件记录
    br_records :  list[GPS_brdc_record class] , 所使用的卫星广播星历记录
    Tr1 : 接收机接收到信号的时刻1,GPS时刻
    Tr2 : 接收机接收到信号的时刻2,GPS时刻
    CRC : bool , 是否进行相对论钟差改正
    c : const , 光速(单位为m/s)
    station1_init_coor : list , 观测站1坐标初值  (如果是在1为已知站的情况下，可直接给入真值)
    station2_init_coor : list , 观测站2坐标初值
    """
    # 筛选出某时刻的记录
    station1_Tr1_ob_records = list(filter(lambda o: o.system == "G" and o.time == Tr1 and o.data != "", station1_ob_records))
    station2_Tr1_ob_records = list(filter(lambda o: o.system == "G" and o.time == Tr1 and o.data != "", station2_ob_records))
    station1_Tr2_ob_records = list(filter(lambda o: o.system == "G" and o.time == Tr2 and o.data != "", station1_ob_records))
    station2_Tr2_ob_records = list(filter(lambda o: o.system == "G" and o.time == Tr2 and o.data != "", station2_ob_records))

    # 选择相邻两个历元均存在的四颗卫星
    available_PRNs = []
    for station2_Tr1_record in station2_Tr1_ob_records:
        # 对station1_ob_record进行同一颗卫星的筛选
        station1_Tr1_record_base = list(filter(lambda o: o.PRN == station2_Tr1_record.PRN, station1_Tr1_ob_records))
        station2_Tr2_record_base = list(filter(lambda o: o.PRN == station2_Tr1_record.PRN, station2_Tr2_ob_records))
        station1_Tr2_record_base = list(filter(lambda o: o.PRN == station2_Tr1_record.PRN, station1_Tr2_ob_records))
        if len(station1_Tr1_record_base) != 1 or len(station2_Tr2_record_base) != 1 or len(station1_Tr2_record_base) != 1:
            continue
        else:
            if (observation_isnot_null(station1_Tr1_record_base[0]) and observation_isnot_null(station2_Tr1_record)
                    and observation_isnot_null(station1_Tr2_record_base[0]) and observation_isnot_null(station2_Tr2_record_base[0])):
                available_PRNs.append(station2_Tr1_record.PRN)
            else:
                continue

    # 判断卫星数是否足够
    num_flag = True
    if len(available_PRNs) >= 6:
        num_flag = True
    elif len(available_PRNs)<6:
        num_flag = False

    # 卫星数足够，则开始迭代进行双差观测的平差求解
    if num_flag:
        # 初始地面点坐标
        X1, Y1, Z1 = station1_init_coor
        X2, Y2, Z2 = station2_init_coor
        # 载波的波长
        lamb =  lamb_L1
        # 平差迭代次数计数
        no = 0
        while True:
            # 如果超出平差迭代求解超出8次则跳出
            if no > 8:
                break
            no += 1
            A = []
            l = []
            Ps = []

            the_prn = available_PRNs[0]    # 卫星1
            # 根据所筛选的PRN构造第一个历元的双差观测方程
            station1_Tr1_base_record = list(filter(lambda o: o.PRN == the_prn, station1_Tr1_ob_records))
            station2_Tr1_base_record = list(filter(lambda o: o.PRN == the_prn, station2_Tr1_ob_records))
            station1_Tr2_base_record = list(filter(lambda o: o.PRN == the_prn, station1_Tr2_ob_records))
            station2_Tr2_base_record = list(filter(lambda o: o.PRN == the_prn, station2_Tr2_ob_records))
            for available_PRN in available_PRNs[1:]:    # 卫星2

                """
                根据PRN对第一个历元两个站观测记录的的筛选
                """

                station1_Tr1_record = list(filter(lambda o: o.PRN == available_PRN, station1_Tr1_ob_records))[0]
                station2_Tr1_record = list(filter(lambda o: o.PRN == available_PRN, station2_Tr1_ob_records))[0]
                # 构造双差方程
                C1obs_sta1sat1_Tr1 = station1_Tr1_base_record[0].data['L1']['observation']
                C1obs_sta2sat1_Tr1 = station2_Tr1_base_record[0].data['L1']['observation']
                C1obs_sta1sat2_Tr1 = station1_Tr1_record.data['L1']['observation']
                C1obs_sta2sat2_Tr1 = station2_Tr1_record.data['L1']['observation']

                # 计算卫星发出信号时刻及发出信号时刻在ECEF坐标系中的位置，以及信号发射时刻站星距离
                # 站1到卫星1
                tr_sta1sat1_Tr1, ts_sta1sat1_Tr1, dts_sta1_Tr1 = cal_EmitTime_from_datetime(Tr1, 0, the_prn, station1_Tr1_base_record[0].data['P1']['observation'], br_records, doCRC=True)
                coorX_sta1sat1_Tr1, coorY_sta1sat1_Tr1, coorZ_sta1sat1_Tr1 = SatellitePosition.cal_SatellitePosition_GPS_GPSws(ts_sta1sat1_Tr1, the_prn, br_records)
                dt_sta1sat1_Tr1 = TimeSystem.cal_deltatime_second_GPSws(tr_sta1sat1_Tr1, ts_sta1sat1_Tr1)
                Xeci_sta1sat1_Tr1, Yeci_sta1sat1_Tr1, Zeci_sta1sat1_Tr1 = CoorTransform.earth_rotation_correction([coorX_sta1sat1_Tr1, coorY_sta1sat1_Tr1, coorZ_sta1sat1_Tr1], dt_sta1sat1_Tr1)
                lou_sta1sat1_Tr10= CoorTransform.cal_distance([X1, Y1, Z1], [Xeci_sta1sat1_Tr1, Yeci_sta1sat1_Tr1, Zeci_sta1sat1_Tr1])

                # 站2到卫星1
                tr_sta2sat1_Tr1, ts_sta2sat1_Tr1, dts_sta2_Tr1 = cal_EmitTime_from_datetime(Tr1, 0, the_prn, station2_Tr1_base_record[0].data['P1']['observation'], br_records, doCRC=True)
                coorX_sta2sat1_Tr1, coorY_sta2sat1_Tr1, coorZ_sta2sat1_Tr1 = SatellitePosition.cal_SatellitePosition_GPS_GPSws(ts_sta2sat1_Tr1, the_prn, br_records)
                dt_sta2sat1_Tr1 = TimeSystem.cal_deltatime_second_GPSws(tr_sta2sat1_Tr1, ts_sta2sat1_Tr1)
                Xeci_sta2sat1_Tr1, Yeci_sta2sat1_Tr1, Zeci_sta2sat1_Tr1 = CoorTransform.earth_rotation_correction([coorX_sta2sat1_Tr1, coorY_sta2sat1_Tr1, coorZ_sta2sat1_Tr1], dt_sta2sat1_Tr1)
                lou_sta2sat1_Tr10= CoorTransform.cal_distance([X2, Y2, Z2], [Xeci_sta2sat1_Tr1, Yeci_sta2sat1_Tr1, Zeci_sta2sat1_Tr1])

                # 站1到卫星2
                tr_sta1sat2_Tr1, ts_sta1sat2_Tr1, dts_sta1_Tr1 = cal_EmitTime_from_datetime(Tr1, 0, available_PRN, station1_Tr1_record.data['P1']['observation'], br_records, doCRC=True)
                coorX_sta1sat2_Tr1, coorY_sta1sat2_Tr1, coorZ_sta1sat2_Tr1 = SatellitePosition.cal_SatellitePosition_GPS_GPSws(ts_sta1sat2_Tr1, available_PRN, br_records)
                dt_sta1sat2_Tr1 = TimeSystem.cal_deltatime_second_GPSws(tr_sta1sat2_Tr1, ts_sta1sat2_Tr1)
                Xeci_sta1sat2_Tr1, Yeci_sta1sat2_Tr1, Zeci_sta1sat2_Tr1 = CoorTransform.earth_rotation_correction([coorX_sta1sat2_Tr1, coorY_sta1sat2_Tr1, coorZ_sta1sat2_Tr1], dt_sta1sat2_Tr1)
                lou_sta1sat2_Tr10 = CoorTransform.cal_distance([X1, Y1, Z1], [Xeci_sta1sat2_Tr1, Yeci_sta1sat2_Tr1, Zeci_sta1sat2_Tr1])

                # 站2到卫星2
                tr_sta2sat2_Tr1, ts_sta2sat2_Tr1, dts_sta2_Tr1 = cal_EmitTime_from_datetime(Tr1, 0, available_PRN, station2_Tr1_record.data['P1']['observation'], br_records, doCRC=True)
                coorX_sta2sat2_Tr1, coorY_sta2sat2_Tr1, coorZ_sta2sat2_Tr1 = SatellitePosition.cal_SatellitePosition_GPS_GPSws(ts_sta2sat2_Tr1, available_PRN, br_records)
                dt_sta2sat2_Tr1 = TimeSystem.cal_deltatime_second_GPSws(tr_sta2sat2_Tr1, ts_sta2sat2_Tr1)
                Xeci_sta2sat2_Tr1, Yeci_sta2sat2_Tr1, Zeci_sta2sat2_Tr1 = CoorTransform.earth_rotation_correction([coorX_sta2sat2_Tr1, coorY_sta2sat2_Tr1, coorZ_sta2sat2_Tr1], dt_sta2sat2_Tr1)
                lou_sta2sat2_Tr10 = CoorTransform.cal_distance([X2, Y2, Z2], [Xeci_sta2sat2_Tr1, Yeci_sta2sat2_Tr1, Zeci_sta2sat2_Tr1])

                # 构造系数阵
                a_sta1_X = - (X1 - Xeci_sta1sat2_Tr1) / lou_sta1sat2_Tr10 + (X1 - Xeci_sta1sat1_Tr1) / lou_sta1sat1_Tr10
                a_sta1_Y = - (Y1 - Yeci_sta1sat2_Tr1) / lou_sta1sat2_Tr10 + (Y1 - Yeci_sta1sat1_Tr1) / lou_sta1sat1_Tr10
                a_sta1_Z = - (Z1 - Zeci_sta1sat2_Tr1) / lou_sta1sat2_Tr10 + (Z1 - Zeci_sta1sat1_Tr1) / lou_sta1sat1_Tr10
                a_sta2_X = - (X2 - Xeci_sta2sat2_Tr1) / lou_sta2sat2_Tr10 + (X2 - Xeci_sta2sat1_Tr1) / lou_sta2sat1_Tr10
                a_sta2_Y = - (Y2 - Yeci_sta2sat2_Tr1) / lou_sta2sat2_Tr10 + (Y2 - Yeci_sta2sat1_Tr1) / lou_sta2sat1_Tr10
                a_sta2_Z = - (Z2 - Zeci_sta1sat2_Tr1) / lou_sta2sat2_Tr10 + (Z2 - Zeci_sta2sat1_Tr1) / lou_sta2sat1_Tr10
                N_DD = make_Ambiguity_coefficient_matrix_row(available_PRNs.index(available_PRN), len(available_PRNs)-1, lamb)
                A_Tr1 = [a_sta1_X, a_sta1_Y, a_sta1_Z, a_sta2_X, a_sta2_Y, a_sta2_Z] + N_DD
                # 构造常数阵
                l_Tr1 = lamb * (C1obs_sta2sat2_Tr1 - C1obs_sta1sat2_Tr1 - C1obs_sta2sat1_Tr1 + C1obs_sta1sat1_Tr1) - lou_sta2sat2_Tr10 + lou_sta2sat1_Tr10 + lou_sta1sat2_Tr10 - lou_sta1sat1_Tr10

                # 加入总矩阵中
                A.append(A_Tr1)
                l.append(l_Tr1)
                Ps.append(1)    # todo: 此处后面再考虑加权的问题

                """
                根据PRN对第二个历元两个站观测记录的的筛选
                """

                station1_Tr2_record = list(filter(lambda o: o.PRN == available_PRN, station1_Tr2_ob_records))[0]
                station2_Tr2_record = list(filter(lambda o: o.PRN == available_PRN, station2_Tr2_ob_records))[0]
                # 构造双差方程
                C1obs_sta1sat1_Tr2 = station1_Tr2_base_record[0].data['L1']['observation']
                C1obs_sta2sat1_Tr2 = station2_Tr2_base_record[0].data['L1']['observation']
                C1obs_sta1sat2_Tr2 = station1_Tr2_record.data['L1']['observation']
                C1obs_sta2sat2_Tr2 = station2_Tr2_record.data['L1']['observation']

                # 计算卫星发出信号时刻及发出信号时刻在ECEF坐标系中的位置，以及信号发射时刻站星距离
                # 站1到卫星1
                tr_sta1sat1_Tr2, ts_sta1sat1_Tr2, dts_sta1_Tr2 = cal_EmitTime_from_datetime(Tr2, 0, the_prn, station1_Tr2_base_record[0].data['P1']['observation'], br_records, doCRC=True)
                coorX_sta1sat1_Tr2, coorY_sta1sat1_Tr2, coorZ_sta1sat1_Tr2 = SatellitePosition.cal_SatellitePosition_GPS_GPSws(ts_sta1sat1_Tr2, the_prn, br_records)
                dt_sta1sat1_Tr2 = TimeSystem.cal_deltatime_second_GPSws(tr_sta1sat1_Tr2, ts_sta1sat1_Tr2)
                Xeci_sta1sat1_Tr2, Yeci_sta1sat1_Tr2, Zeci_sta1sat1_Tr2 = CoorTransform.earth_rotation_correction([coorX_sta1sat1_Tr2, coorY_sta1sat1_Tr2, coorZ_sta1sat1_Tr2], dt_sta1sat1_Tr2)
                lou_sta1sat1_Tr20 = CoorTransform.cal_distance([X1, Y1, Z1], [Xeci_sta1sat1_Tr2, Yeci_sta1sat1_Tr2, Zeci_sta1sat1_Tr2])

                # 站2到卫星1
                tr_sta2sat1_Tr2, ts_sta2sat1_Tr2, dts_sta2_Tr2 = cal_EmitTime_from_datetime(Tr2, 0, the_prn, station2_Tr2_base_record[0].data['P1']['observation'], br_records, doCRC=True)
                coorX_sta2sat1_Tr2, coorY_sta2sat1_Tr2, coorZ_sta2sat1_Tr2 = SatellitePosition.cal_SatellitePosition_GPS_GPSws(ts_sta2sat1_Tr2, the_prn, br_records)
                dt_sta2sat1_Tr2 = TimeSystem.cal_deltatime_second_GPSws(tr_sta2sat1_Tr2, ts_sta2sat1_Tr2)
                Xeci_sta2sat1_Tr2, Yeci_sta2sat1_Tr2, Zeci_sta2sat1_Tr2 = CoorTransform.earth_rotation_correction([coorX_sta2sat1_Tr2, coorY_sta2sat1_Tr2, coorZ_sta2sat1_Tr2], dt_sta2sat1_Tr2)
                lou_sta2sat1_Tr20 = CoorTransform.cal_distance([X2, Y2, Z2], [Xeci_sta2sat1_Tr2, Yeci_sta2sat1_Tr2, Zeci_sta2sat1_Tr2])

                # 站1到卫星2
                tr_sta1sat2_Tr2, ts_sta1sat2_Tr2, dts_sta1_Tr2 = cal_EmitTime_from_datetime(Tr2, 0, available_PRN,station1_Tr2_record.data['P1']['observation'], br_records, doCRC=True)
                coorX_sta1sat2_Tr2, coorY_sta1sat2_Tr2, coorZ_sta1sat2_Tr2 = SatellitePosition.cal_SatellitePosition_GPS_GPSws(ts_sta1sat2_Tr2, available_PRN, br_records)
                dt_sta1sat2_Tr2 = TimeSystem.cal_deltatime_second_GPSws(tr_sta1sat2_Tr2, ts_sta1sat2_Tr2)
                Xeci_sta1sat2_Tr2, Yeci_sta1sat2_Tr2, Zeci_sta1sat2_Tr2 = CoorTransform.earth_rotation_correction([coorX_sta1sat2_Tr2, coorY_sta1sat2_Tr2, coorZ_sta1sat2_Tr2], dt_sta1sat2_Tr2)
                lou_sta1sat2_Tr20 = CoorTransform.cal_distance([X1, Y1, Z1], [Xeci_sta1sat2_Tr2, Yeci_sta1sat2_Tr2, Zeci_sta1sat2_Tr2])

                # 站2到卫星2
                tr_sta2sat2_Tr2, ts_sta2sat2_Tr2, dts_sta2_Tr2 = cal_EmitTime_from_datetime(Tr2, 0, available_PRN,station2_Tr2_record.data['P1']['observation'], br_records, doCRC=True)
                coorX_sta2sat2_Tr2, coorY_sta2sat2_Tr2, coorZ_sta2sat2_Tr2 = SatellitePosition.cal_SatellitePosition_GPS_GPSws(ts_sta2sat2_Tr2, available_PRN, br_records)
                dt_sta2sat2_Tr2 = TimeSystem.cal_deltatime_second_GPSws(tr_sta2sat2_Tr2, ts_sta2sat2_Tr2)
                Xeci_sta2sat2_Tr2, Yeci_sta2sat2_Tr2, Zeci_sta2sat2_Tr2 = CoorTransform.earth_rotation_correction([coorX_sta2sat2_Tr2, coorY_sta2sat2_Tr2, coorZ_sta2sat2_Tr2], dt_sta2sat2_Tr2)
                lou_sta2sat2_Tr20 = CoorTransform.cal_distance([X2, Y2, Z2], [Xeci_sta2sat2_Tr2, Yeci_sta2sat2_Tr2, Zeci_sta2sat2_Tr2])

                # 构造系数阵
                a_sta1_X = - (X1 - Xeci_sta1sat2_Tr2) / lou_sta1sat2_Tr20 + (X1 - Xeci_sta1sat1_Tr2) / lou_sta1sat1_Tr20
                a_sta1_Y = - (Y1 - Yeci_sta1sat2_Tr2) / lou_sta1sat2_Tr20 + (Y1 - Yeci_sta1sat1_Tr2) / lou_sta1sat1_Tr20
                a_sta1_Z = - (Z1 - Zeci_sta1sat2_Tr2) / lou_sta1sat2_Tr20 + (Z1 - Zeci_sta1sat1_Tr2) / lou_sta1sat1_Tr20
                a_sta2_X = (X2 - Xeci_sta2sat2_Tr2) / lou_sta2sat2_Tr20 - (X2 - Xeci_sta2sat1_Tr2) / lou_sta2sat1_Tr20
                a_sta2_Y = (Y2 - Yeci_sta2sat2_Tr2) / lou_sta2sat2_Tr20 - (Y2 - Yeci_sta2sat1_Tr2) / lou_sta2sat1_Tr20
                a_sta2_Z = (Z2 - Zeci_sta1sat2_Tr2) / lou_sta2sat2_Tr20 - (Z2 - Zeci_sta2sat1_Tr2) / lou_sta2sat1_Tr20
                N_DD = make_Ambiguity_coefficient_matrix_row(available_PRNs.index(available_PRN), len(available_PRNs) - 1, lamb)
                A_Tr2 = [a_sta1_X, a_sta1_Y, a_sta1_Z, a_sta2_X, a_sta2_Y, a_sta2_Z] + N_DD
                # 构造常数阵
                l_Tr2 = lamb * (C1obs_sta2sat2_Tr2 - C1obs_sta1sat2_Tr2 - C1obs_sta2sat1_Tr2 + C1obs_sta1sat1_Tr2) - lou_sta2sat2_Tr20 + lou_sta2sat1_Tr20 + lou_sta1sat2_Tr20 - lou_sta1sat1_Tr20

                # 加入总矩阵中
                A.append(A_Tr2)
                l.append(l_Tr2)
                Ps.append(1)  # todo: 此处后面再考虑加权的问题

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

            # 更新参数
            dX1 = x[0]
            dY1 = x[1]
            dZ1 = x[2]
            dX2 = x[3]
            dY2 = x[4]
            dZ2 = x[5]
            N = x[6:]
            X1 += dX1
            Y1 += dY1
            Z1 += dZ1
            X2 += dX2
            Y2 += dY2
            Z2 += dZ2
            # 判断迭代停止条件
            if abs(dX1) < 1e-4 and abs(dY1) < 1e-4 and abs(dZ1) < 1e-4 and abs(dX2) < 1e-4 and abs(dY2) < 1e-4 and abs(dZ2) < 1e-4:
                break
    return [X1, Y1, Z1], [X2, Y2, Z2]


# 基于载波相位的双差 (其中一个观测站为已知坐标站点)
def DD_onCarrierPhase_1known(station1_ob_records, station2_ob_records, br_records,
                      Tr1, Tr2, station1_init_coor, CRC=True, c=299792458, ambi_fix=True,
                      station2_init_coor=[2000, 2000, 2000]):
    """
    station1_ob_records : list[GPS_observation_record class] , 所使用的观测站1观测文件记录
    station2_ob_records : list[GPS_observation_record class] , 所使用的观测站2观测文件记录
    br_records :  list[GPS_brdc_record class] , 所使用的卫星广播星历记录
    Tr1 : 接收机接收到信号的时刻1,GPS时刻
    Tr2 : 接收机接收到信号的时刻2,GPS时刻
    station1_init_coor : list , 观测站1坐标已知值
    CRC : bool , 是否进行相对论钟差改正
    c : const , 光速(单位为m/s)
    ambi_fix : bool , 是否进行模糊度固定
    station2_init_coor : list , 观测站2坐标初值
    """
    # 筛选出某时刻的记录
    station1_Tr1_ob_records = list(
        filter(lambda o: o.system == "G" and o.time == Tr1 and o.data != "", station1_ob_records))
    station2_Tr1_ob_records = list(
        filter(lambda o: o.system == "G" and o.time == Tr1 and o.data != "", station2_ob_records))
    station1_Tr2_ob_records = list(
        filter(lambda o: o.system == "G" and o.time == Tr2 and o.data != "", station1_ob_records))
    station2_Tr2_ob_records = list(
        filter(lambda o: o.system == "G" and o.time == Tr2 and o.data != "", station2_ob_records))

    # 选择相邻两个历元均存在的四颗卫星
    available_PRNs = []
    for station2_Tr1_record in station2_Tr1_ob_records:
        # 对station1_ob_record进行同一颗卫星的筛选
        station1_Tr1_record_base = list(filter(lambda o: o.PRN == station2_Tr1_record.PRN, station1_Tr1_ob_records))
        station2_Tr2_record_base = list(filter(lambda o: o.PRN == station2_Tr1_record.PRN, station2_Tr2_ob_records))
        station1_Tr2_record_base = list(filter(lambda o: o.PRN == station2_Tr1_record.PRN, station1_Tr2_ob_records))
        if len(station1_Tr1_record_base) != 1 or len(station2_Tr2_record_base) != 1 or len(
                station1_Tr2_record_base) != 1:
            continue
        else:
            if (observation_isnot_null(station1_Tr1_record_base[0]) and observation_isnot_null(
                    station2_Tr1_record)
                    and observation_isnot_null(station1_Tr2_record_base[0]) and observation_isnot_null(
                        station2_Tr2_record_base[0])):
                available_PRNs.append(station2_Tr1_record.PRN)
            else:
                continue

    # 判断卫星数是否足够
    num_flag = True
    print(available_PRNs)
    ob_num = len(available_PRNs)
    if ob_num < 4:
        num_flag = False
    elif ob_num >= 4:
        num_flag = True

    # 卫星数足够，则开始迭代进行双差观测的平差求解
    # if num_flag == False:
    #     return
    if num_flag:
        # 初始地面点坐标
        X1, Y1, Z1 = station1_init_coor
        X2, Y2, Z2 = station2_init_coor
        # 载波的波长
        lamb = lamb_L1
        # 平差迭代次数计数
        no = 0
        while True:
            # 如果超出平差迭代求解超出8次则跳出
            if no > 8:
                break
            no += 1

            # 初始化各历元矩阵
            A1 = []
            l1 = []
            A2 = []
            l2 = []

            the_prn = available_PRNs[0]  # 卫星1

            # 获取卫星1的两站观测记录
            station1_Tr1_base_record = list(filter(lambda o: o.PRN == the_prn, station1_Tr1_ob_records))[0]
            station2_Tr1_base_record = list(filter(lambda o: o.PRN == the_prn, station2_Tr1_ob_records))[0]
            station1_Tr2_base_record = list(filter(lambda o: o.PRN == the_prn, station1_Tr2_ob_records))[0]
            station2_Tr2_base_record = list(filter(lambda o: o.PRN == the_prn, station2_Tr2_ob_records))[0]

            for available_PRN in available_PRNs[1:]:  # 卫星2

                """
                根据PRN对第一个历元两个站观测记录的的筛选
                """

                station1_Tr1_record = list(filter(lambda o: o.PRN == available_PRN, station1_Tr1_ob_records))[0]
                station2_Tr1_record = list(filter(lambda o: o.PRN == available_PRN, station2_Tr1_ob_records))[0]

                # 构造双差方程
                L1obs_sta1sat1_Tr1 = station1_Tr1_base_record.data['L1']['observation']
                L1obs_sta2sat1_Tr1 = station2_Tr1_base_record.data['L1']['observation']
                L1obs_sta1sat2_Tr1 = station1_Tr1_record.data['L1']['observation']
                L1obs_sta2sat2_Tr1 = station2_Tr1_record.data['L1']['observation']

                # 计算卫星发出信号时刻及发出信号时刻在ECEF坐标系中的位置，以及信号发射时刻站星距离
                # 站1到卫星1
                tr_sta1sat1_Tr1, ts_sta1sat1_Tr1, dts_sta1_Tr1 = cal_EmitTime_from_datetime(Tr1, 0, the_prn, station1_Tr1_base_record.data['P1']['observation'], br_records, doCRC=True)
                coorX_sta1sat1_Tr1, coorY_sta1sat1_Tr1, coorZ_sta1sat1_Tr1 = SatellitePosition.cal_SatellitePosition_GPS_GPSws(ts_sta1sat1_Tr1, the_prn, br_records)
                dt_sta1sat1_Tr1 = TimeSystem.cal_deltatime_second_GPSws(tr_sta1sat1_Tr1, ts_sta1sat1_Tr1)
                Xeci_sta1sat1_Tr1, Yeci_sta1sat1_Tr1, Zeci_sta1sat1_Tr1 = CoorTransform.earth_rotation_correction([coorX_sta1sat1_Tr1, coorY_sta1sat1_Tr1, coorZ_sta1sat1_Tr1], dt_sta1sat1_Tr1)
                lou_sta1sat1_Tr10 = CoorTransform.cal_distance([X1, Y1, Z1], [Xeci_sta1sat1_Tr1, Yeci_sta1sat1_Tr1,Zeci_sta1sat1_Tr1])

                # 站2到卫星1
                tr_sta2sat1_Tr1, ts_sta2sat1_Tr1, dts_sta2_Tr1 = cal_EmitTime_from_datetime(Tr1, 0, the_prn, station2_Tr1_base_record.data['P1']['observation'], br_records, doCRC=True)
                coorX_sta2sat1_Tr1, coorY_sta2sat1_Tr1, coorZ_sta2sat1_Tr1 = SatellitePosition.cal_SatellitePosition_GPS_GPSws(ts_sta2sat1_Tr1, the_prn, br_records)
                dt_sta2sat1_Tr1 = TimeSystem.cal_deltatime_second_GPSws(tr_sta2sat1_Tr1, ts_sta2sat1_Tr1)
                Xeci_sta2sat1_Tr1, Yeci_sta2sat1_Tr1, Zeci_sta2sat1_Tr1 = CoorTransform.earth_rotation_correction([coorX_sta2sat1_Tr1, coorY_sta2sat1_Tr1, coorZ_sta2sat1_Tr1], dt_sta2sat1_Tr1)
                lou_sta2sat1_Tr10 = CoorTransform.cal_distance([X2, Y2, Z2], [Xeci_sta2sat1_Tr1, Yeci_sta2sat1_Tr1,Zeci_sta2sat1_Tr1])

                # 站1到卫星2
                tr_sta1sat2_Tr1, ts_sta1sat2_Tr1, dts_sta1_Tr1 = cal_EmitTime_from_datetime(Tr1, 0, available_PRN, station1_Tr1_record.data['P1']['observation'],br_records, doCRC=True)
                coorX_sta1sat2_Tr1, coorY_sta1sat2_Tr1, coorZ_sta1sat2_Tr1 = SatellitePosition.cal_SatellitePosition_GPS_GPSws(ts_sta1sat2_Tr1, available_PRN, br_records)
                dt_sta1sat2_Tr1 = TimeSystem.cal_deltatime_second_GPSws(tr_sta1sat2_Tr1, ts_sta1sat2_Tr1)
                Xeci_sta1sat2_Tr1, Yeci_sta1sat2_Tr1, Zeci_sta1sat2_Tr1 = CoorTransform.earth_rotation_correction([coorX_sta1sat2_Tr1, coorY_sta1sat2_Tr1, coorZ_sta1sat2_Tr1], dt_sta1sat2_Tr1)
                lou_sta1sat2_Tr10 = CoorTransform.cal_distance([X1, Y1, Z1], [Xeci_sta1sat2_Tr1, Yeci_sta1sat2_Tr1,Zeci_sta1sat2_Tr1])

                # 站2到卫星2
                tr_sta2sat2_Tr1, ts_sta2sat2_Tr1, dts_sta2_Tr1 = cal_EmitTime_from_datetime(Tr1, 0, available_PRN, station2_Tr1_record.data['P1']['observation'], br_records, doCRC=True)
                coorX_sta2sat2_Tr1, coorY_sta2sat2_Tr1, coorZ_sta2sat2_Tr1 = SatellitePosition.cal_SatellitePosition_GPS_GPSws(
                    ts_sta2sat2_Tr1, available_PRN, br_records)
                dt_sta2sat2_Tr1 = TimeSystem.cal_deltatime_second_GPSws(tr_sta2sat2_Tr1, ts_sta2sat2_Tr1)
                Xeci_sta2sat2_Tr1, Yeci_sta2sat2_Tr1, Zeci_sta2sat2_Tr1 = CoorTransform.earth_rotation_correction(
                    [coorX_sta2sat2_Tr1, coorY_sta2sat2_Tr1, coorZ_sta2sat2_Tr1], dt_sta2sat2_Tr1)
                lou_sta2sat2_Tr10 = CoorTransform.cal_distance([X2, Y2, Z2], [Xeci_sta2sat2_Tr1, Yeci_sta2sat2_Tr1,
                                                                              Zeci_sta2sat2_Tr1])

                # 构造系数阵
                a_sta2_X =  (X2 - Xeci_sta2sat2_Tr1) / lou_sta2sat2_Tr10 - (X2 - Xeci_sta2sat1_Tr1) / lou_sta2sat1_Tr10
                a_sta2_Y =  (Y2 - Yeci_sta2sat2_Tr1) / lou_sta2sat2_Tr10 - (Y2 - Yeci_sta2sat1_Tr1) / lou_sta2sat1_Tr10
                a_sta2_Z =  (Z2 - Zeci_sta1sat2_Tr1) / lou_sta2sat2_Tr10 - (Z2 - Zeci_sta2sat1_Tr1) / lou_sta2sat1_Tr10

                N_DD = make_Ambiguity_coefficient_matrix_row(available_PRNs.index(available_PRN), len(available_PRNs) - 1, lamb)
                A_Tr1 = [a_sta2_X, a_sta2_Y, a_sta2_Z] + N_DD
                # 构造常数阵
                l_Tr1 = lamb * (L1obs_sta2sat2_Tr1 - L1obs_sta1sat2_Tr1 - L1obs_sta2sat1_Tr1 + L1obs_sta1sat1_Tr1) - lou_sta2sat2_Tr10 + lou_sta2sat1_Tr10 + lou_sta1sat2_Tr10 - lou_sta1sat1_Tr10

                # 加入历元1矩阵中
                A1.append(A_Tr1)
                l1.append(l_Tr1)

                """
                根据PRN对第一个历元两个站观测记录的的筛选
                """

                # 根据PRN对第一个历元两个站观测记录的的筛选
                station1_Tr2_record = list(filter(lambda o: o.PRN == available_PRN, station1_Tr2_ob_records))[0]
                station2_Tr2_record = list(filter(lambda o: o.PRN == available_PRN, station2_Tr2_ob_records))[0]
                # 构造双差方程
                L1obs_sta1sat1_Tr2 = station1_Tr2_base_record.data['L1']['observation']
                L1obs_sta2sat1_Tr2 = station2_Tr2_base_record.data['L1']['observation']
                L1obs_sta1sat2_Tr2 = station1_Tr2_record.data['L1']['observation']
                L1obs_sta2sat2_Tr2 = station2_Tr2_record.data['L1']['observation']

                # 计算卫星发出信号时刻及发出信号时刻在ECEF坐标系中的位置，以及信号发射时刻站星距离
                # 站1到卫星1
                tr_sta1sat1_Tr2, ts_sta1sat1_Tr2, dts_sta1_Tr2 = cal_EmitTime_from_datetime(Tr2, 0, the_prn,station1_Tr2_base_record.data['P1']['observation'],br_records, doCRC=True)
                coorX_sta1sat1_Tr2, coorY_sta1sat1_Tr2, coorZ_sta1sat1_Tr2 = SatellitePosition.cal_SatellitePosition_GPS_GPSws(ts_sta1sat1_Tr2, the_prn, br_records)
                dt_sta1sat1_Tr2 = TimeSystem.cal_deltatime_second_GPSws(tr_sta1sat1_Tr2, ts_sta1sat1_Tr2)
                Xeci_sta1sat1_Tr2, Yeci_sta1sat1_Tr2, Zeci_sta1sat1_Tr2 = CoorTransform.earth_rotation_correction([coorX_sta1sat1_Tr2, coorY_sta1sat1_Tr2, coorZ_sta1sat1_Tr2], dt_sta1sat1_Tr2)
                lou_sta1sat1_Tr20 = CoorTransform.cal_distance([X1, Y1, Z1], [Xeci_sta1sat1_Tr2, Yeci_sta1sat1_Tr2, Zeci_sta1sat1_Tr2])

                # 站2到卫星1
                tr_sta2sat1_Tr2, ts_sta2sat1_Tr2, dts_sta2_Tr2 = cal_EmitTime_from_datetime(Tr2, 0, the_prn, station2_Tr2_base_record.data['P1']['observation'], br_records, doCRC=True)
                coorX_sta2sat1_Tr2, coorY_sta2sat1_Tr2, coorZ_sta2sat1_Tr2 = SatellitePosition.cal_SatellitePosition_GPS_GPSws(ts_sta2sat1_Tr2, the_prn, br_records)
                dt_sta2sat1_Tr2 = TimeSystem.cal_deltatime_second_GPSws(tr_sta2sat1_Tr2, ts_sta2sat1_Tr2)
                Xeci_sta2sat1_Tr2, Yeci_sta2sat1_Tr2, Zeci_sta2sat1_Tr2 = CoorTransform.earth_rotation_correction([coorX_sta2sat1_Tr2, coorY_sta2sat1_Tr2, coorZ_sta2sat1_Tr2], dt_sta2sat1_Tr2)
                lou_sta2sat1_Tr20 = CoorTransform.cal_distance([X2, Y2, Z2], [Xeci_sta2sat1_Tr2, Yeci_sta2sat1_Tr2,Zeci_sta2sat1_Tr2])

                # 站1到卫星2
                tr_sta1sat2_Tr2, ts_sta1sat2_Tr2, dts_sta1_Tr2 = cal_EmitTime_from_datetime(Tr2, 0, available_PRN, station1_Tr2_record.data['P1']['observation'], br_records, doCRC=True)
                coorX_sta1sat2_Tr2, coorY_sta1sat2_Tr2, coorZ_sta1sat2_Tr2 = SatellitePosition.cal_SatellitePosition_GPS_GPSws(ts_sta1sat2_Tr2, available_PRN, br_records)
                dt_sta1sat2_Tr2 = TimeSystem.cal_deltatime_second_GPSws(tr_sta1sat2_Tr2, ts_sta1sat2_Tr2)
                Xeci_sta1sat2_Tr2, Yeci_sta1sat2_Tr2, Zeci_sta1sat2_Tr2 = CoorTransform.earth_rotation_correction([coorX_sta1sat2_Tr2, coorY_sta1sat2_Tr2, coorZ_sta1sat2_Tr2], dt_sta1sat2_Tr2)
                lou_sta1sat2_Tr20 = CoorTransform.cal_distance([X1, Y1, Z1], [Xeci_sta1sat2_Tr2, Yeci_sta1sat2_Tr2, Zeci_sta1sat2_Tr2])

                # 站2到卫星2
                tr_sta2sat2_Tr2, ts_sta2sat2_Tr2, dts_sta2_Tr2 = cal_EmitTime_from_datetime(Tr2, 0, available_PRN, station2_Tr2_record.data['P1']['observation'], br_records, doCRC=True)
                coorX_sta2sat2_Tr2, coorY_sta2sat2_Tr2, coorZ_sta2sat2_Tr2 = SatellitePosition.cal_SatellitePosition_GPS_GPSws(ts_sta2sat2_Tr2, available_PRN, br_records)
                dt_sta2sat2_Tr2 = TimeSystem.cal_deltatime_second_GPSws(tr_sta2sat2_Tr2, ts_sta2sat2_Tr2)
                Xeci_sta2sat2_Tr2, Yeci_sta2sat2_Tr2, Zeci_sta2sat2_Tr2 = CoorTransform.earth_rotation_correction([coorX_sta2sat2_Tr2, coorY_sta2sat2_Tr2, coorZ_sta2sat2_Tr2], dt_sta2sat2_Tr2)
                lou_sta2sat2_Tr20 = CoorTransform.cal_distance([X2, Y2, Z2], [Xeci_sta2sat2_Tr2, Yeci_sta2sat2_Tr2,Zeci_sta2sat2_Tr2])

                # 构造系数阵

                a_sta2_X = (X2 - Xeci_sta2sat2_Tr2) / lou_sta2sat2_Tr20 - (X2 - Xeci_sta2sat1_Tr2) / lou_sta2sat1_Tr20
                a_sta2_Y = (Y2 - Yeci_sta2sat2_Tr2) / lou_sta2sat2_Tr20 - (Y2 - Yeci_sta2sat1_Tr2) / lou_sta2sat1_Tr20
                a_sta2_Z = (Z2 - Zeci_sta1sat2_Tr2) / lou_sta2sat2_Tr20 - (Z2 - Zeci_sta2sat1_Tr2) / lou_sta2sat1_Tr20
                N_DD = make_Ambiguity_coefficient_matrix_row(available_PRNs.index(available_PRN),len(available_PRNs) - 1, lamb)
                A_Tr2 = [a_sta2_X, a_sta2_Y, a_sta2_Z] + N_DD
                # 构造常数阵
                l_Tr2 = lamb * (
                        L1obs_sta2sat2_Tr2 - L1obs_sta1sat2_Tr2 - L1obs_sta2sat1_Tr2 + L1obs_sta1sat1_Tr2) - lou_sta2sat2_Tr20 + lou_sta2sat1_Tr20 + lou_sta1sat2_Tr20 - lou_sta1sat1_Tr20

                # 加入总矩阵中
                A2.append(A_Tr2)
                l2.append(l_Tr2)

            # 构造权阵并求解
            Ps1 = get_DD_Pmatrix(ob_num-1)
            Ps2 = get_DD_Pmatrix(ob_num-1)
            Pz = diagonalize_squarematrix(Ps1, Ps2)
            A = np.array(A1 + A2)
            l = np.array(l1 + l2)

            # 改正数发散太过严重则不再继续平差
            if abs(max(l.tolist())) > 1e10:
                break
            x = np.linalg.inv(A.T @ Pz @ A) @ (A.T @ Pz @ l)
            Q = np.linalg.inv(A.T @ Pz @ A)
            Qcoor = Q[:3, :3]

            # 更新参数
            dX2 = x[0]
            dY2 = x[1]
            dZ2 = x[2]
            N_float = x[3:]
            X2 += dX2
            Y2 += dY2
            Z2 += dZ2
            print(dX2, dY2, dZ2)
            # 判断迭代停止条件
            if abs(dX2) < 1e-4 and abs(dY2) < 1e-4 and abs(dZ2) < 1e-4:
                break

        # 进行模糊度固定
        if ambi_fix:
            # 调用LAMBDA方法进行整数估计
            Qaa = Q[3:, 3:]
            Qba = Q[:3, 3:]
            N_fixed, sqnorm, Ps, Qzhat, Z, nfixed, mu = LAMBDA.main(N_float, Qaa)
            # 更新参数估计
            b_hat = np.array([X2, Y2, Z2])
            a_hat = N_float
            Coor = MAPmethod(b_hat, a_hat, Qaa, Qba, N_fixed[:, 1])
            # todo 计算MAP计算后的坐标方差
        else:
            Coor = [X2, Y2, Z2]

    return Coor, Qcoor


if __name__ == "__main__":
    station2_observation_file = r"edata\obs\warn3100.20o"
    station1_observation_file = r"edata\obs\leij3100.20o"    # 已知站点
    broadcast_file = r"edata\sat_obit\brdc3100.20n"
    # 读入观测文件内容,得到类型对象列表
    knownStation_ob_records = DoFile.read_GPS_oFile(station1_observation_file)
    unknownStation_ob_records = DoFile.read_GPS_oFile(station2_observation_file)
    br_records = DoFile.read_GPS_nFile(broadcast_file)
    print("数据读取完毕！")
    Tr = datetime.datetime(2020, 11, 5, 0, 0, 0)
    init_coor = [3658785.6000, 784471.1000, 5147870.7000]
    # init_coor = [0, 0, 0]
    knownStation_coor = [0.389873613453103E+07, 0.855345521080705E+06, 0.495837257579542E+07]
    true_coors = []
    cal_coors = []
    while Tr < datetime.datetime(2020, 11, 5, 1, 45, 16):
        Tr2 = Tr + datetime.timedelta(seconds=30)
        CoorXYZ, Q = DD_onCarrierPhase_1known(knownStation_ob_records, unknownStation_ob_records, br_records, Tr, Tr2,
                                      knownStation_coor, station2_init_coor=init_coor, ambi_fix=False)
        Xk, Yk, Zk = CoorXYZ
        cal_coors.append([Xk, Yk, Zk])
        true_coors.append([0.365878555276965E+07, 0.784471127238666E+06, 0.514787071062059E+07])  # warn
        # true_coors.append([0.389873613453103E+07,0.855345521080705E+06,0.495837257579542E+07])   #leij
        # true_coors.append([-0.267442768572702E+07,0.375714305701559E+07,0.439152148514515E+07])  #chan
        Tr += datetime.timedelta(seconds=30)
    SPP.cal_NEUerrors(true_coors, cal_coors)