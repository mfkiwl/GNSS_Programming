# -*- coding: utf-8 -*-
"""

@title:	multi-antenna RTK
@author: iDeal0103
@status:	Active
@type:	Process
@created:	23-May-2022
@post-History:	23-May-2022

comment:
    1. multi-antenna RTK
    2.

"""


# import
import datetime
import math

import numpy as np

import RTK2
from utils.const import *
import SinglePointPosition as SPP
import utils.SatellitePosition as SatellitePosition
import utils.CoorTransform as CoorTransform
import RTK
import utils.LAMBDA as LAMBDA
import utils.CycleslipsDetection as CycleslipsDetection
import utils.DoFile as DoFile
import utils.ResultAnalyse as ResultAnalyse
import utils.TimeSystem as TimeSystem
import utils.PictureResults as PictureResults
import RTK_manage

#
# class DDamb_difference_storage:
#     def __init__(self, station1_name, station2_name):
#         self.station1_name = station1_name
#         self.station2_name = station2_name
#
#     def update_amb_transfer(self, ambiguity_transfer1, ambiguity_transfer2):
#         self.ambiguity_transfer1 = ambiguity_transfer1
#         self.ambiguity_transfer2 = ambiguity_transfer2
#         if self.ambiguity_transfer1.base_svn != self.ambiguity_transfer2.base_svn:
#             self.ambiguity_transfer2.shift_to_new_basesats(self.ambiguity_transfer1.base_svn)
#         self.base_svn = self.ambiguity_transfer1.base_svn
#         self.svns = []
#         self.DDamb_difference = []
#         self.get_DDamb_difference()
#
#     def get_DDamb_difference(self):
#         for svn in self.ambiguity_transfer1.satellites:
#             if svn in self.ambiguity_transfer2.satellites:
#                 self.svns.append(svn)
#                 at1_svn_index = self.ambiguity_transfer1.satellites.index(svn)
#                 at2_svn_index = self.ambiguity_transfer2.satellites.index(svn)
#                 at1_DDamb = self.ambiguity_transfer1.ambiguitys[at1_svn_index]
#                 at2_DDamb = self.ambiguity_transfer2.ambiguitys[at2_svn_index]
#                 self.DDamb_difference.append(at2_DDamb - at1_DDamb)
#
#     def get_station2_DDamb(self, station1_DDamb, svn):
#         station2_DDamb = station1_DDamb + self.DDamb_difference[svn]
#         return station2_DDamb



# 储存各站间各波段模糊度差值类
class DD_ambiguilty_between_stations:
    def __init__(self, main_station, ambiguity_transfer1):
        self.main_station = main_station
        self.main_ambiguity_transfer = ambiguity_transfer1
        self.difference_stations = []
        self.DDamb_transfers = []
        self.DDamb_difference_storages = []

    def update_main_antenna_transfer(self, ambiguity_transfer1):
        self.main_ambiguity_transfer = ambiguity_transfer1
        for ambtransfer in self.DDamb_transfers:
            ambtransfer.imitate_ambtransfer(ambiguity_transfer1)

    def add_station_DDambiguity_transfer(self, station, DD_ambiguity_transfer):
        self.difference_stations.append(station)
        self.DDamb_transfers.append(DD_ambiguity_transfer)

    def update_DDamb_difference_storages(self, the_station=""):
        if the_station:    # 更新某个站
            pass
        else:     # 全部更新
            self.DDamb_difference_storages = []
            for station in self.difference_stations:
                station_index = self.difference_stations.index(station)
                DDamb_difference = []
                # for svn in self.main_ambiguity_transfer.satellites:
                #     at1_svn_index = self.main_ambiguity_transfer.satellites.index(svn)
                #     at2_svn_index = self.DDamb_transfers[station_index].satellites.index(svn)
                #     at1_DDamb = self.main_ambiguity_transfer.ambiguitys[at1_svn_index]
                #     at2_DDamb = self.DDamb_transfers[station_index].ambiguitys[at2_svn_index]
                #     DDamb_difference.append(at2_DDamb - at1_DDamb)
                # 认为卫星顺序已经排列好了
                for i in range(len(self.main_ambiguity_transfer.satellites)):
                    at1_DDamb = self.main_ambiguity_transfer.ambiguitys[i]
                    at2_DDamb = self.DDamb_transfers[station_index].ambiguitys[i]
                    DDamb_difference.append(at2_DDamb - at1_DDamb)
                self.DDamb_difference_storages.append(DDamb_difference)


# 阵列天线网解的卡尔曼滤波
class Kalman_Filter_forMultiAntenna:
    def __init__(self, stations_DD_records_atTr, DD_ambiguilty_relation, cp_band, pr_band, mode_int):
        """
        stations_DD_records_atTr : list[DD_records_atTr class], 所有天线双差观测集合， 第一个是主天线
        DD_ambiguilty_relation : DD_ambiguilty_between_stations
        cp_band : 载波相位
        pr_band : 伪距
        mode_int : 滤波解模式，0为静态，1为动态
        """
        self.stations_DD_records_atTr = stations_DD_records_atTr
        self.cp_band = cp_band
        self.pr_band = pr_band
        self.mode_int = mode_int
        self.DD_ambiguilty_relation = DD_ambiguilty_relation
        self.sta_num = len(self.stations_DD_records_atTr)

    def getF(self):
        DD_ambiguilty_num = len(self.stations_DD_records_atTr[0].DD_observations_data[self.cp_band])   # 双差模糊度个数
        n = DD_ambiguilty_num + 3 * self.sta_num  # 双差位置模糊度参数+位置参数
        Fmatrix = np.diag([1 for i in range(n)])
        return Fmatrix

    def getH(self, nav_data, baseref_coor, stations_coor, x_array):
        """
        baseref_coor : 设置固定参考站
        stations_coor : 各天线坐标（第一个是主天线位置）
        """
        H11 = []
        H21 = []
        hx = []
        for index in range(self.sta_num):
            DD_records = self.stations_DD_records_atTr[index]
            sta2_coor = stations_coor[index]
            H_cp_part = []   # 单个站载波相位对应系数矩阵
            H_pr_part = []   # 单个站伪距对应系数矩阵
            for DD_data in DD_records.DD_observations_data[self.cp_band]:
                # 计算站1星1元素
                ts_sta1sat1, dts_sta1sat1 = SPP.cal_EmitTime_from_datetime(DD_data.T, DD_data.sat1, DD_data.obs1_sat1.data[self.pr_band]['observation'], nav_data, doCRC=True)
                coorX_sta1sat1, coorY_sta1sat1, coorZ_sta1sat1 = SatellitePosition.cal_SatellitePosition_GPS_GPSws(ts_sta1sat1, DD_data.sat1, nav_data)
                dt_sta1sat1 = DD_data.obs1_sat1.data[self.pr_band]['observation'] / c
                Xeci_sta1sat1, Yeci_sta1sat1, Zeci_sta1sat1 = CoorTransform.earth_rotation_correction([coorX_sta1sat1, coorY_sta1sat1, coorZ_sta1sat1], dt_sta1sat1)
                lou_sta1sat1 = CoorTransform.cal_distance(baseref_coor, [Xeci_sta1sat1, Yeci_sta1sat1, Zeci_sta1sat1])

                # 计算站1星2元素
                ts_sta1sat2, dts_sta1sat2 = SPP.cal_EmitTime_from_datetime(DD_data.T, DD_data.sat2, DD_data.obs1_sat2.data[self.pr_band]['observation'], nav_data, doCRC=True)
                coorX_sta1sat2, coorY_sta1sat2, coorZ_sta1sat2 = SatellitePosition.cal_SatellitePosition_GPS_GPSws(ts_sta1sat2, DD_data.sat2, nav_data)
                dt_sta1sat2 = DD_data.obs1_sat2.data[self.pr_band]['observation'] / c
                Xeci_sta1sat2, Yeci_sta1sat2, Zeci_sta1sat2 = CoorTransform.earth_rotation_correction([coorX_sta1sat2, coorY_sta1sat2, coorZ_sta1sat2], dt_sta1sat2)
                lou_sta1sat2 = CoorTransform.cal_distance(baseref_coor, [Xeci_sta1sat2, Yeci_sta1sat2, Zeci_sta1sat2])

                # 计算站2星1元素
                ts_sta2sat1, dts_sta2sat1 = SPP.cal_EmitTime_from_datetime(DD_data.T, DD_data.sat1, DD_data.obs2_sat1.data[self.pr_band]['observation'], nav_data, doCRC=True)
                coorX_sta2sat1, coorY_sta2sat1, coorZ_sta2sat1 = SatellitePosition.cal_SatellitePosition_GPS_GPSws(ts_sta2sat1, DD_data.sat1, nav_data)
                dt_sta2sat1 = DD_data.obs2_sat1.data[self.pr_band]['observation'] / c
                Xeci_sta2sat1, Yeci_sta2sat1, Zeci_sta2sat1 = CoorTransform.earth_rotation_correction([coorX_sta2sat1, coorY_sta2sat1, coorZ_sta2sat1], dt_sta2sat1)
                lou_sta2sat1 = CoorTransform.cal_distance(sta2_coor, [Xeci_sta2sat1, Yeci_sta2sat1, Zeci_sta2sat1])

                # 计算站2星2元素
                ts_sta2sat2, dts_sta2sat2 = SPP.cal_EmitTime_from_datetime(DD_data.T, DD_data.sat2, DD_data.obs2_sat2.data[self.pr_band]['observation'], nav_data, doCRC=True)
                coorX_sta2sat2, coorY_sta2sat2, coorZ_sta2sat2 = SatellitePosition.cal_SatellitePosition_GPS_GPSws(ts_sta2sat2, DD_data.sat2, nav_data)
                dt_sta2sat2 = DD_data.obs2_sat2.data[self.pr_band]['observation'] / c
                Xeci_sta2sat2, Yeci_sta2sat2, Zeci_sta2sat2 = CoorTransform.earth_rotation_correction([coorX_sta2sat2, coorY_sta2sat2, coorZ_sta2sat2], dt_sta2sat2)
                lou_sta2sat2 = CoorTransform.cal_distance(sta2_coor, [Xeci_sta2sat2, Yeci_sta2sat2, Zeci_sta2sat2])

                # 构造几何系数阵
                X2, Y2, Z2 = sta2_coor
                a_sta2_X = (X2 - Xeci_sta2sat2) / lou_sta2sat2 - (X2 - Xeci_sta2sat1) / lou_sta2sat1
                a_sta2_Y = (Y2 - Yeci_sta2sat2) / lou_sta2sat2 - (Y2 - Yeci_sta2sat1) / lou_sta2sat1
                a_sta2_Z = (Z2 - Zeci_sta2sat2) / lou_sta2sat2 - (Z2 - Zeci_sta2sat1) / lou_sta2sat1

                H_cp_part.append([a_sta2_X, a_sta2_Y, a_sta2_Z])
                H_pr_part.append([a_sta2_X, a_sta2_Y, a_sta2_Z])
                hx.append(lou_sta2sat2 - lou_sta1sat2 - lou_sta2sat1 + lou_sta1sat1)
            H11.append(H_cp_part)
            H21.append(H_pr_part)

        # 构造系数矩阵的左半部分
        H11 = RTK.diagonalize_several_squarematrix(H11)
        H21 = RTK.diagonalize_several_squarematrix(H21)
        # 构造系数矩阵的右半部分
        n = H11.shape[0]  # 得到行数
        H12 = np.eye(n) * get_lamb_from_band(self.cp_band)
        H22 = np.zeros((n, n))
        H = np.block([[np.array(H11), H12], [np.array(H21), H22]])

        # 由参数值得到的观测值
        hx1 = np.array(2*hx)  # 空间距离项
        DD_ambiguitys = x_array[3 * self.sta_num:].tolist()   # 获得主天线的双差模糊度
        DD_ambiguitys_difference = [0 for i in range(len(DD_ambiguitys))]
        DD_ambiguitys += 2 * [0 for i in range(len(DD_ambiguitys))]  # 补齐
        for i in range(self.sta_num - 1):  # 将各站的模糊度读进来
            # for svn in self.DD_ambiguilty_relation.main_ambiguity_transfer.satellites:
            #     svn_index = self.DD_ambiguilty_relation.main_ambiguity_transfer.satellites.index(svn)
            #     DD_ambiguitys.append(self.DD_ambiguilty_relation.DDamb_difference_storages[i][svn_index])    # 将该站该卫星的模糊度加进来
            DD_ambiguitys_difference += self.DD_ambiguilty_relation.DDamb_difference_storages[i]   # 将该站该卫星的模糊度加进来
        hx21 = np.array(DD_ambiguitys) + np.array(DD_ambiguitys_difference)
        hx22 = [0 for i in range(H21.shape[0])]   # 伪距部分不需要额外加项
        hx2 = get_lamb_from_band(self.cp_band) * np.array(hx21+hx22)
        hx = hx1+hx2
        return H, hx

    def gety(self):
        y1 = []
        y2 = []
        # 加入各天线数据
        for DD_records in self.stations_DD_records_atTr:
            for DD_data in DD_records.DD_observations_data[self.cp_band]:
                y1.append(get_lamb_from_band(self.cp_band) * DD_data.DD_obs)
            for DD_data in DD_records.DD_observations_data[self.pr_band]:
                y2.append(DD_data.DD_obs)
        y = np.array(y1+y2)
        return y

    def getx(self, stations_coor):
        # 转换类型
        sta_coors = []
        for sta_coor in stations_coor:
            if isinstance(sta_coor, np.ndarray):
                sta_coor = sta_coor.tolist()
            sta_coors.append(sta_coor)
        if isinstance(self.stations_DD_records_atTr[0].DD_ambiguitys[self.cp_band], np.ndarray):
            DD_ambiguitys = self.stations_DD_records_atTr[0].DD_ambiguitys[self.cp_band].tolist()
        else:
            DD_ambiguitys = self.stations_DD_records_atTr[0].DD_ambiguitys[self.cp_band]
        x = np.array(sta_coors + DD_ambiguitys)
        return x

    def getQ(self):
        Qs = self.stations_DD_records_atTr[0].DD_ambiguitys_noise[self.cp_band]
        # 静态基线模式
        if self.mode == 0:
            Q = np.diag([1e-1, 1e-1, 1e-1] * self.sta_num + Qs)
        # 动态基线模式
        elif self.mode == 1:
            Q = np.diag([10, 10, 10] * self.sta_num + Qs)
        return Q

    def getR(self, sat_n, sigma1=0.02, sigma2=3):
        vc_matrix_uint = np.array(self.get_DDobs_vcmatrix(sat_n, 1))
        # 构造矩阵两种观测值的vc元阵
        vc_matrix_part = []
        vc_matrix_row = []
        for i in range(self.sta_num):
            vc_matrix_row.append(vc_matrix_uint)
        for i in range(self.sta_num):
            vc_matrix_part.append(vc_matrix_row)
        for i in range(self.sta_num):
            vc_matrix_part[i][i] *= 2
        # 构造R阵
        cpobs_vc_matrix = np.block(vc_matrix_part) * sigma1 ** 2
        probs_vc_matrix = np.block(vc_matrix_part) * sigma2 ** 2
        zeros = np.zeros((sat_n*self.sta_num, sat_n*self.sta_num))
        R = np.block([[cpobs_vc_matrix, zeros], [zeros, probs_vc_matrix]])
        return R

    def get_DDobs_vcmatrix(self, n, sigma):
        vc_matrix = np.full((n, n), 1).astype(float)
        for i in range(n):
            vc_matrix[i, i] = 2
        vc_matrix *= sigma**2
        return vc_matrix


    def ekf_estimation(self, P_before, sta1_coor, stations_coor, nav_records):
        # 预测
        F = self.getF()
        Q = self.getQ()
        x_pri = self.getx(stations_coor)
        P_pri = F @ P_before @ F.T + Q
        # 更新
        H, hx = self.getH(nav_records, sta1_coor, stations_coor, x_pri)
        R = self.getR()
        y = self.gety()
        K = P_pri @ H.T @ np.linalg.inv(H @ P_pri @ H.T + R)
        x_est = x_pri + K @ (y - hx)
        # 输出残差
        H, hx = self.getH(nav_records, sta1_coor, x_pri.tolist()[:3], x_pri)
        P_est = (np.eye(len(x_est)) - K @ H) @ P_pri
        return x_est, P_est


def get_common_element(list1, list2):
    common_elements = []
    for element in list1:
        if element in list2:
            common_elements.append(element)
    return common_elements



def three_antenna_RTK_Kalmanfilter(start_time, end_time, basestation_ob_records, antennas_ob_records, br_records,
                                 basetation_coor, antenna_init_coor, cp_band, pr_band,
                                 interval_time, mode_int, cycle_slip_detect=False, baseline_length=0, l_sigma=0, ele_limit=13):
    """
    start_time : datetime.datetime , 时间段开始时刻
    end_time : datetime.datetime , 时间段结束时刻
    basestation_ob_records : 基准站多时刻观测数据
    antennas_ob_records : list[ob_records] , 各个天线的多时刻观测数据
    br_records : 卫星星历数据
    basetation_coor : 基准站坐标
    antenna_init_coor : 未知的天线坐标（多个合用一个，因为互相之间比较接近）
    cp_band : 载波波段
    pr_band : 伪距波段
    interval_time : 间隔时间, s
    mode_int : 定位模式, 0为静态,1为动态
    baseline_length : 基线长度, m
    l_sigma : 基线长度标准差, m
    ele_limit : 截止高度角, 度
    """
    Tr = start_time

    # 先对几个天线与基准站形成的基线进行单基线解，确定模糊度
    preprocess_time = start_time + datetime.timedelta(seconds=600)
    antenna1_coordinates, antenna1_P, antenna1_amb_trans = RTK_manage.constaneously_RTK_withfilter(start_time, preprocess_time,
                                 basestation_ob_records, antennas_ob_records[0], br_records, basetation_coor,
                                 antenna_init_coor, cp_band, pr_band, interval_time, cycle_slip_detect, ele_limit)
    antenna2_coordinates, antenna2_P, antenna2_amb_trans = RTK_manage.constaneously_RTK_withfilter(start_time, preprocess_time,
                                 basestation_ob_records, antennas_ob_records[1], br_records, basetation_coor,
                                 antenna_init_coor, cp_band, pr_band, interval_time, cycle_slip_detect, ele_limit)
    antenna3_coordinates, antenna3_P, antenna3_amb_trans = RTK_manage.constaneously_RTK_withfilter(start_time, preprocess_time,
                                 basestation_ob_records, antennas_ob_records[2], br_records, basetation_coor,
                                 antenna_init_coor, cp_band, pr_band, interval_time, cycle_slip_detect, ele_limit)
    # 统一amb_transfer的卫星的组成和排列顺序
    antenna2_amb_trans.imitate_ambtransfer(antenna1_amb_trans)
    antenna3_amb_trans.imitate_ambtransfer(antenna1_amb_trans)
    antenna1_coor_P = antenna1_P[:3, :3]
    antenna2_coor_P = antenna2_P[:3, :3]
    antenna3_coor_P = antenna3_P[:3, :3]
    DD_P = antenna1_P[3:, 3:]
    Tr = preprocess_time

    # 构造双差模糊度差值存储对象
    DDamb_between_stations = DD_ambiguilty_between_stations("A", antenna1_amb_trans)
    # antenna12_DDamb_difference = DDamb_difference_storage("A", "B")
    # antenna12_DDamb_difference.update_amb_transfer(antenna1_amb_trans, antenna2_amb_trans)
    # antenna13_DDamb_difference = DDamb_difference_storage("A", "C")
    # antenna13_DDamb_difference.update_amb_transfer(antenna1_amb_trans, antenna3_amb_trans)
    DDamb_between_stations.add_station_DDambiguity_transfer("B", antenna2_amb_trans)
    DDamb_between_stations.add_station_DDambiguity_transfer("C", antenna3_amb_trans)
    DDamb_between_stations.update_DDamb_difference_storages()

    # 阵列天线网解
    # 参数继承
    DD_GF_collection = ''
    antenna1_coor = antenna1_coordinates[-1]
    antenna2_coor = antenna2_coordinates[-1]
    antenna3_coor = antenna3_coordinates[-1]
    while Tr < end_time:
        # Tr时刻基准站和各天线的观测记录
        basestation_ob_records_atTr = list(filter(lambda o:o.time == Tr and o.SVN[0] == "G" and o.managed_data_flag[cp_band] and o.managed_data_flag[pr_band] and o.managed_data_flag['L2_L'], basestation_ob_records))
        antenna1_ob_records_atTr = list(filter(lambda o:o.time == Tr and o.SVN[0] == "G" and o.managed_data_flag[cp_band] and o.managed_data_flag[pr_band] and o.managed_data_flag['L2_L'], antennas_ob_records[0]))
        antenna2_ob_records_atTr = list(filter(lambda o: o.time == Tr and o.SVN[0] == "G" and o.managed_data_flag[cp_band] and o.managed_data_flag[pr_band] and o.managed_data_flag['L2_L'], antennas_ob_records[1]))
        antenna3_ob_records_atTr = list(filter(lambda o: o.time == Tr and o.SVN[0] == "G" and o.managed_data_flag[cp_band] and o.managed_data_flag[pr_band] and o.managed_data_flag['L2_L'], antennas_ob_records[2]))

        # 获取符合条件的共有svn
        basestation_svns = [x.SVN for x in basestation_ob_records_atTr]
        antenna1_svns = [x.SVN for x in antenna1_ob_records_atTr]
        antenna2_svns = [x.SVN for x in antenna2_ob_records_atTr]
        antenna3_svns = [x.SVN for x in antenna3_ob_records_atTr]
        svn_Tr = get_common_element(get_common_element(get_common_element(basestation_svns, antenna1_svns), antenna2_svns), antenna3_svns)

        # 计算各个卫星的高度角, 并筛选符合高度角条件的数据
        sats_ele = {}
        svn_ele_qualitified = []
        w, s = TimeSystem.from_datetime_cal_GPSws(Tr)
        t = TimeSystem.GPSws(w, s)
        for svn in svn_Tr:
            sat_coor = SatellitePosition.cal_SatellitePosition_GPS_GPSws(t, svn, br_records)
            sats_ele[svn] = CoorTransform.cal_ele_and_A(antenna1_coor, sat_coor)[0] * 180 / math.pi
            if sats_ele[svn] > ele_limit:
                svn_ele_qualitified.append(svn)

        # 提取符合高度角要求的观测记录
        basestation_ob_records_atTr = list(filter(lambda o: o.SVN in svn_ele_qualitified, basestation_ob_records_atTr))
        antenna1_ob_records_atTr = list(filter(lambda o: o.SVN in svn_ele_qualitified, antenna1_ob_records_atTr))
        antenna2_ob_records_atTr = list(filter(lambda o: o.SVN in svn_ele_qualitified, antenna2_ob_records_atTr))
        antenna3_ob_records_atTr = list(filter(lambda o: o.SVN in svn_ele_qualitified, antenna3_ob_records_atTr))

        # 构造双差观测值
        # 基站-天线1
        base_antenna1_DD_records_collection = RTK_manage.DD_records_atTr(Tr, [cp_band, pr_band])
        trans_matrix = base_antenna1_DD_records_collection.set_elebest_base_satellite_from_ambiguity_transfer(antenna1_amb_trans, basestation_ob_records_atTr, antenna1_ob_records_atTr)
        # 若发生参考星变化则改变其他天线的参考星
        if not np.all(trans_matrix == 0):
            antenna2_amb_trans.shift_to_new_basesats(antenna1_amb_trans.base_svn)
            antenna3_amb_trans.shift_to_new_basesats(antenna1_amb_trans.base_svn)
        svn_removed_index, new_in, new_in_svn = base_antenna1_DD_records_collection.add_satellites(svn_ele_qualitified, basestation_ob_records_atTr, antenna1_ob_records_atTr)
        # 判断是否发生了参考星的变化
        if np.all(trans_matrix == 0):  # 参考星不变
            # 调整较上个历元被去除的卫星对应的vc阵
            if svn_removed_index:  # 上个历元的卫星有缺失
                DD_P = RTK_manage.remove_cross_from_matrix(DD_P, svn_removed_index, svn_removed_index)
        else:  # 参考星变化
            # 调整较上个历元被去除的卫星对应的vc阵
            if svn_removed_index:  # 上个历元的卫星要进行删减的
                DD_P = RTK_manage.remove_cross_from_matrix(DD_P, svn_removed_index, svn_removed_index)
                # trans_matrix = np.delete(trans_matrix, svn_removed_index, axis=1)
                trans_matrix = RTK_manage.remove_cross_from_matrix(trans_matrix, svn_removed_index, svn_removed_index)
            DD_P = trans_matrix @ DD_P @ trans_matrix.T
        # 更新由卫星数量增加导致的模糊度参数vc阵的变化
        if new_in != 0:  # 此历元有新增卫星
            eyes = 10000 * np.eye(new_in)  # 前面的常数为新增卫星的方差
            DD_P = RTK.diagonalize_squarematrix(DD_P, eyes)

        # 基站-天线2
        base_antenna2_DD_records_collection = RTK_manage.DD_records_atTr(Tr, [cp_band, pr_band])
        trans_matrix = base_antenna2_DD_records_collection.set_base_satellite_from_ambiguity_transfer(antenna2_amb_trans, basestation_ob_records_atTr, antenna2_ob_records_atTr)
        svn_removed_index, new_in, new_in_svn = base_antenna2_DD_records_collection.add_satellites(svn_ele_qualitified, basestation_ob_records_atTr, antenna2_ob_records_atTr)
        # 不考虑参考星的变化
        if not np.all(trans_matrix == 0):  # 参考星不变
            print(Tr, 'antenna2参考星想要发生变化！')
            # 调整较上个历元被去除的卫星对应的vc阵
        if svn_removed_index:  # 上个历元的卫星有缺失
            print(Tr, 'antenna2卫星数据缺失：', svn_removed_index)

        # 基站-天线3
        base_antenna3_DD_records_collection = RTK_manage.DD_records_atTr(Tr, [cp_band, pr_band])
        trans_matrix = base_antenna3_DD_records_collection.set_base_satellite_from_ambiguity_transfer(antenna3_amb_trans, basestation_ob_records_atTr, antenna3_ob_records_atTr)
        svn_removed_index, new_in, new_in_svn = base_antenna3_DD_records_collection.add_satellites(svn_ele_qualitified, basestation_ob_records_atTr, antenna3_ob_records_atTr)
        # 不考虑参考星的变化
        if not np.all(trans_matrix == 0):  # 参考星不变
            print(Tr, 'antenna3参考星想要发生变化！')
            # 调整较上个历元被去除的卫星对应的vc阵
        if svn_removed_index:  # 上个历元的卫星有缺失
            print(Tr, 'antenna3卫星数据缺失：', svn_removed_index)

        # 准备滤波解需要的数据
        P = RTK.diagonalize_several_squarematrix([antenna1_coor_P, antenna2_coor_P, antenna3_coor_P, DD_P])
        antennas_coors = [antenna1_coordinates, antenna2_coordinates, antenna3_coordinates]
        KM_estimator = Kalman_Filter_forMultiAntenna(basestation_ob_records, DD_ambiguilty_between_stations, cp_band, pr_band, mode_int)
        x, P = KM_estimator.ekf_estimation(P, sta1_coor=basetation_coor, sta2_coor=antennas_coors, nav_records=br_records)

        # 更新参数和协方差阵
        antenna1_coordinates = x[:3]
        antenna2_coordinates = x[3:6]
        antenna3_coordinates = x[6:9]
        base_antenna1_DD_records_collection.DD_ambiguitys[cp_band] = x[9:].tolist()
        antenna1_coor_P = P[:3, :3]
        antenna2_coor_P = P[3:6, 3:6]
        antenna3_coor_P = P[6:9, 6:9]
        DD_P = P[9:, 9:]

        # 构造amb_transfer
        ambiguitys = base_antenna1_DD_records_collection.DD_ambiguitys[cp_band]
        base_svn = base_antenna1_DD_records_collection.base_svn
        satellites = base_antenna1_DD_records_collection.qualitified_satellites
        amb_nos = [0.000001 for i in range(len(P) - 3)]
        antenna1_amb_trans = RTK_manage.ambiguity_transfer(cp_band, ambiguitys, base_svn, satellites, amb_nos)



























