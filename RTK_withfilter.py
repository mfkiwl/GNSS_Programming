# -*- coding: utf-8 -*-
"""

@title:	RTK with filter
@author: iDeal0103
@status:	Active
@type:	Process
@created:	11-March-2022
@post-History:	11-March-2022

comment:
    1. 带滤波器的基线解算

"""

# import
import datetime

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

# 双差观测值类
class DD_recordclass():
    def __init__(self, Tr, sat1, sat2, band):
        self.sat1 = sat1
        self.sat2 = sat2
        self.band = band
        self.T = Tr

    def obs_isnot_null(self, band, obs1_sat1, obs1_sat2, obs2_sat1, obs2_sat2):
        # 判断波段数据是否完整
        if not (obs1_sat1.data.__contains__(band) and obs2_sat1.data.__contains__(band)):
            isnot_null_flag = False
        elif obs1_sat1.data[band]['observation'] == '' or obs1_sat2.data[band]['observation'] == ''\
                or obs2_sat1.data[band]['observation'] == '' or obs2_sat2.data[band]['observation'] == '':
            isnot_null_flag = False
        else:
            isnot_null_flag = True
        return isnot_null_flag

    def get_DD_obs(self, obs1_sat1, obs1_sat2, obs2_sat1, obs2_sat2):
        if self.obs_isnot_null(self.band, obs1_sat1, obs1_sat2, obs2_sat1, obs2_sat2):
            # 构造双差观测值
            obs_sta1sat1 = float(obs1_sat1.data[self.band]['observation'])
            obs_sta1sat2 = float(obs1_sat2.data[self.band]['observation'])
            obs_sta2sat1 = float(obs2_sat1.data[self.band]['observation'])
            obs_sta2sat2 = float(obs2_sat2.data[self.band]['observation'])
            self.DD_obs = obs_sta2sat2 - obs_sta1sat2 - obs_sta2sat1 + obs_sta1sat1
            self.obs1_sat1 = obs1_sat1
            self.obs1_sat2 = obs1_sat2
            self.obs2_sat1 = obs2_sat1
            self.obs2_sat2 = obs2_sat2
        else:
            self.DD_obs = ''

    def get_DD_GF(self, band2, obs1_sat1, obs1_sat2, obs2_sat1, obs2_sat2):
        band1 = self.band
        self.GF_band = [band1, band2]
        if self.obs_isnot_null(band1, obs1_sat1, obs1_sat2, obs2_sat1, obs2_sat2) \
                and self.obs_isnot_null(band2, obs1_sat1, obs1_sat2, obs2_sat1, obs2_sat2):
            # 分别进行双差观测值构造
            obs_sta1sat1_band1 = float(obs1_sat1.data[band1]['observation'])
            obs_sta1sat2_band1 = float(obs1_sat2.data[band1]['observation'])
            obs_sta2sat1_band1 = float(obs2_sat1.data[band1]['observation'])
            obs_sta2sat2_band1 = float(obs2_sat2.data[band1]['observation'])
            DD_obs_band1 = obs_sta2sat2_band1 - obs_sta1sat2_band1 -obs_sta2sat1_band1 + obs_sta1sat1_band1
            obs_sta1sat1_band2 = float(obs1_sat1.data[band2]['observation'])
            obs_sta1sat2_band2 = float(obs1_sat2.data[band2]['observation'])
            obs_sta2sat1_band2 = float(obs2_sat1.data[band2]['observation'])
            obs_sta2sat2_band2 = float(obs2_sat2.data[band2]['observation'])
            DD_obs_band2 = obs_sta2sat2_band2 - obs_sta1sat2_band2 -obs_sta2sat1_band2 +obs_sta1sat1_band2
            self.DD_GF = abs(DD_obs_band1 - DD_obs_band2)
        else:
            self.DD_GF = 100



class DD_GF_record():
    def __init__(self, base_svn):
        self.base_svn = base_svn
        self.GF_records = {}

    def add_DD_GF(self, svn, value):
        self.GF_records[svn] = value


# 参考星变换问题
class DD_records_atTr:
    def __init__(self, Tr, bands, sta1='1', sta2='2'):
        self.Tr = Tr     # GPS时间
        self.sta1 = sta1
        self.sta2 = sta2
        self.qualitified_satellites = []      # 待差卫星
        self.unqualitified_satellites = []      # 不适合双差的卫星
        self.bands = bands
        self.DD_observations_data = {}  # 双差观测值记录
        self.DD_ambiguitys = {}   # 双差模糊度记录
        self.DD_ambiguitys_noise = {}    # 模糊度噪声
        for band in self.bands:
            self.DD_observations_data[band] = []
            if band in GPS_carrier_phase_list:
                self.DD_ambiguitys[band] = []
                self.DD_ambiguitys_noise[band] = []
        self.amb_noise = 10000
        self.DD_GF = {}

    def set_base_satellite(self, base_SVN, sta1obs_records_Tr, sta2obs_records_Tr):
        self.base_svn = base_SVN
        self.base_sta1obs_record = list(filter(lambda o: o.SVN == base_SVN, sta1obs_records_Tr))[0]
        self.base_sta2obs_record = list(filter(lambda o: o.SVN == base_SVN, sta2obs_records_Tr))[0]

    # 增加新卫星进行双差
    def add_satellite(self, the_svn, sta1obs_records_Tr, sta2obs_records_Tr):
        if self.base_svn:
            other_sta1obs_record = list(filter(lambda o: o.SVN == the_svn, sta1obs_records_Tr))[0]
            other_sta2obs_record = list(filter(lambda o: o.SVN == the_svn, sta2obs_records_Tr))[0]
            # 对卫星数据进行检查
            qualified_flag = True
            DD_obses = []
            for band in self.bands:
                DD_record = DD_recordclass(self.Tr, self.base_svn, the_svn, band)
                DD_record.get_DD_obs(self.base_sta1obs_record, other_sta1obs_record, self.base_sta2obs_record, other_sta2obs_record)
                if DD_record.DD_obs == '':
                    qualified_flag = False
                    break
                DD_obses.append(DD_record)
            # 根据卫星数据进行是否加入待差卫星的判断, 并增加数据
            if qualified_flag:
                self.qualitified_satellites.append(the_svn)
                for band in self.bands:
                    self.DD_observations_data[band].append(DD_obses[self.bands.index(band)])
                    if band in GPS_carrier_phase_list:
                        lamb = get_lamb_from_band(band)
                        self.DD_ambiguitys[band].append(DD_obses[self.bands.index(band)].DD_obs/lamb)
                        self.DD_ambiguitys_noise[band].append(self.amb_noise)
            else:
                self.unqualitified_satellites.append(the_svn)
                print(the_svn + '卫星数据不符合要求！')
        else:
            qualified_flag = False
            print('请设置参考星！')
        return qualified_flag

    # def get_DD_GF_record(self):
    #     DD_GFs = DD_GF_record(self.base_svn)
    #     for svn in self.qualitified_satellites:
    #         L1_DD = list(filter(lambda o: o.sat2 == svn, self.DD_observations_data['L1']))[0]
    #         L2_DD = list(filter(lambda o: o.sat2 == svn, self.DD_observations_data['L2']))[0]
    #         DD_GF = L1_DD - L2_DD
    #         DD_GFs.add_DD_GF(svn, DD_GF)
    #     return DD_GFs
    def get_DD_GF_record(self, band):
        for data in self.DD_observations_data[band]:
            self.DD_GF[data.sat2] = data.DD_GF


    def add_satellites(self, svns, sta1obs_records_Tr, sta2obs_records_Tr, DD_GF_before):
        qualitified_satellites_copy = self.qualitified_satellites.copy()
        if self.qualitified_satellites:
            # 先获得所有数据缺失卫星的原始索引
            svn_satellite_miss = []
            for svn in self.qualitified_satellites:
                if svn not in svns:
                    svn_satellite_miss.append(svn)
            # 将所有数据缺失的上历元合格卫星删除
            self.delete_from_qualitified_svn(svn_satellite_miss)
            # for svn in self.qualitified_satellites:
            #     if svn not in svns:
            #         svn_index = self.qualitified_satellites.index(svn)
            #         self.qualitified_satellites.remove(svn)
            #         for band in self.bands:
            #             if band in GPS_carrier_phase_list:
            #                 self.DD_ambiguitys[band].pop(svn_index)
            #                 self.DD_ambiguitys_noise[band].pop(svn_index)
            # 对本历元有记录的上历元合格卫星进行数据完整性检查, 记录数据不完整的
            svn_data_incomplete = []
            for svn in self.qualitified_satellites:
                qualified_flag = True
                DD_obses = []      # 某卫星所有波段的双差观测数据
                sta1obs_records_Tr_svn = list(filter(lambda o: o.SVN == svn, sta1obs_records_Tr))[0]
                sta2obs_records_Tr_svn = list(filter(lambda o: o.SVN == svn, sta2obs_records_Tr))[0]
                for band in self.bands:
                    DD_record = DD_recordclass(self.Tr, self.base_svn, svn, band)
                    DD_record.get_DD_obs(self.base_sta1obs_record, sta1obs_records_Tr_svn, self.base_sta2obs_record,
                                         sta2obs_records_Tr_svn)
                    if band in GPS_carrier_phase_list:
                        DD_record.get_DD_GF('L2', self.base_sta1obs_record, sta1obs_records_Tr_svn, self.base_sta2obs_record,
                                             sta2obs_records_Tr_svn)    # 计算GF双差值
                    if DD_record.DD_obs == '':
                        qualified_flag = False
                        break
                    DD_obses.append(DD_record)
                # 根据上述卫星数据检查判断是否加入待差卫星
                if qualified_flag:
                    for band in self.bands:
                        self.DD_observations_data[band].append(DD_obses[self.bands.index(band)])
                else:
                    svn_data_incomplete.append(svn)
                    print(svn + '卫星数据不符合要求！')
            # 将数据不完整的对象删除
            if svn_data_incomplete:
                self.delete_from_qualitified_svn(svn_data_incomplete)
            # 进行周跳探测, 将认为有周跳的卫星删除
            cycle_slip_svn = []
            if DD_GF_before:
                self.get_DD_GF_record('L1')
                for svn in self.qualitified_satellites:
                    GF_before = DD_GF_before[svn]
                    GF_now = self.DD_GF[svn]
                    if CycleslipsDetection.GF_detector(GF_before, GF_now):
                        cycle_slip_svn.append(svn)
                        print(svn + '卫星发生周跳！')
                if cycle_slip_svn:
                    self.delete_from_qualitified_svn(cycle_slip_svn)
                # for svn in cycle_slip_svn:
                #     svn_index = self.qualitified_satellites.index(svn)
                #     self.qualitified_satellites.remove(svn)
                #     for band in self.bands:
                #         if band in GPS_carrier_phase_list:
                #             self.DD_ambiguitys[band].pop(svn_index)
                #             self.DD_ambiguitys_noise[band].pop(svn_index)
            # 获得所有被移除的上一历元合格卫星的索引
            svn_removed = svn_satellite_miss + svn_data_incomplete + cycle_slip_svn
            svn_removed_index = []
            for svn in svn_removed:
                svn_index = qualitified_satellites_copy.index(svn)
                svn_removed_index.append(svn_index)
            # 将新卫星加入
            new_in = 0
            for svn in svns:
                if svn not in self.qualitified_satellites and svn not in self.unqualitified_satellites and svn != self.base_svn:
                    new_in_flag = self.add_satellite(svn, sta1obs_records_Tr, sta2obs_records_Tr)
                    if new_in_flag:
                        new_in += 1
        else:
            new_in = 0
            for svn in svns:
                if svn != self.base_svn:
                    new_in_flag = self.add_satellite(svn, sta1obs_records_Tr, sta2obs_records_Tr)
                    if new_in_flag:
                        new_in += 1
        return svn_removed_index, new_in

    def delete_from_qualitified_svn(self, svns):
        if svns:
            for svn in svns:
                svn_index = self.qualitified_satellites.index(svn)
                self.qualitified_satellites.remove(svn)
                for band in self.bands:
                    # # 删除双差观测记录
                    # self.DD_observations_data[band].pop(svn_index)
                    if band in GPS_carrier_phase_list:
                        # 删除模糊度和模糊度参数噪声
                        self.DD_ambiguitys[band].pop(svn_index)
                        self.DD_ambiguitys_noise[band].pop(svn_index)

    def set_base_satellite_from_ambiguity_transfer(self, ambiguity_transfer, sta1obs_records_Tr, sta2obs_records_Tr):
        base_sta1obs_record = list(filter(lambda o: o.SVN == ambiguity_transfer.base_svn, sta1obs_records_Tr))
        base_sta2obs_record = list(filter(lambda o: o.SVN == ambiguity_transfer.base_svn, sta2obs_records_Tr))
        transfer_matrix = np.array([0])
        # 判断上一历元继承的参考星存在及本历元中数据完整, 否则需要更换参考星
        if not self.obs_is_qualitified(base_sta1obs_record, base_sta2obs_record):
            # 更换参考星
            print('更换参考星!原参考星为：' + ambiguity_transfer.base_svn)
            for svn in ambiguity_transfer.satellites:
                sta1obs_record = list(filter(lambda o: o.SVN == svn, sta1obs_records_Tr))
                sta2obs_record = list(filter(lambda o: o.SVN == svn, sta2obs_records_Tr))
                if self.obs_is_qualitified(sta1obs_record, sta2obs_record):
                    self.set_base_satellite(svn, sta1obs_record, sta2obs_record)
                    transfer_matrix = ambiguity_transfer.shift_to_new_basesats(svn)
                    self.qualitified_satellites = ambiguity_transfer.satellites
                    self.DD_ambiguitys[ambiguity_transfer.band] = ambiguity_transfer.ambiguitys
                    self.DD_ambiguitys_noise[ambiguity_transfer.band] = ambiguity_transfer.amb_noises
                    print('更换参考星!现参考星为：' + self.base_svn)
                    break
        else:
            self.set_base_satellite(ambiguity_transfer.base_svn, base_sta1obs_record, base_sta2obs_record)
            self.qualitified_satellites = ambiguity_transfer.satellites
            self.DD_ambiguitys[ambiguity_transfer.band] = ambiguity_transfer.ambiguitys
            self.DD_ambiguitys_noise[ambiguity_transfer.band] = ambiguity_transfer.amb_noises
        # set_base_satellite_from_ambiguity_transfer 与 add_satellites 配合使用
        return transfer_matrix

    # 判断svn对应的观测数据是否符合要求
    def obs_is_qualitified(self, obs1_list, obs2_list):
        # 判断波段数据是否完整
        if not (obs1_list and obs2_list):
            is_qualitified = [False]
        else:
            obs1 = obs1_list[0]
            obs2 = obs2_list[0]
            is_qualitified = []
            for band in self.bands:
                if obs1.data[band]['observation'] == '' or obs2.data[band]['observation'] == '':
                    isnot_null_flag = False
                else:
                    isnot_null_flag = True
                is_qualitified.append(isnot_null_flag)
        return False not in is_qualitified



class ambiguity_transfer():
    def __init__(self, band, ambiguitys, base_svn, satellites, amb_noises):
        self.band = band
        self.base_svn = base_svn
        self.ambiguitys = ambiguitys
        self.satellites = satellites
        self.amb_noises = amb_noises
        self.init_ambiguity_noise = 10000

    def shift_to_new_basesats(self, new_base_svn):
        self.new_base_svn = new_base_svn
        new_base_svn_index = self.satellites.index(new_base_svn)
        transfer_matrix_pre = np.eye(len(self.satellites)-1)
        transfer_matrix = np.insert(transfer_matrix_pre, new_base_svn_index, np.full((len(self.satellites)-1), -1), axis=1)
        self.satellites.remove(new_base_svn)
        self.ambiguitys = (transfer_matrix @ np.array(self.ambiguitys)).tolist()
        # F1.直接给参考星变换后的模糊度的方差一个大值
        # self.amb_noises = [100000 for i in range(len(self.ambiguitys))]
        # F2.参考星变换后模糊度的方差进行传递求得
        # am_vc = transfer_matrix @ self.amb_noises @ transfer_matrix.T
        am_vc = transfer_matrix @ np.diag(self.amb_noises) @ transfer_matrix.T
        # self.amb_noises = am_vc
        self.amb_noises = get_diaglist_from_matrix(am_vc)
        return transfer_matrix



# 卡尔曼滤波器
class Kalman_Filter():
    def __init__(self, DD_records_atTr, band1, band2):
        """
        DD_records_atTr : DD_records_atTr class, 双差观测集合
        band1 : 载波相位
        band2 : 伪距
        """
        self.DD_records = DD_records_atTr
        self.band1 = band1
        self.band2 = band2
        # print(len(DD_records_atTr.DD_observations_data[band1]), len(DD_records_atTr.DD_ambiguitys_noise[band1]))

    def getF(self, band):
        n = len(self.DD_records.DD_observations_data[band])+3
        Fmatrix = np.diag([1 for i in range(n)])
        return Fmatrix

    def getH(self, nav_data, sta1_coor, sta2_coor, band):
        H11 = []
        H21 = []
        hx = []
        for DD_data in self.DD_records.DD_observations_data[band]:
            # 计算站1星1元素
            ts_sta1sat1, dts_sta1sat1 = SPP.cal_EmitTime_from_datetime(DD_data.T, DD_data.sat1, DD_data.obs1_sat1.data['P2']['observation'], nav_data, doCRC=True)
            coorX_sta1sat1, coorY_sta1sat1, coorZ_sta1sat1 = SatellitePosition.cal_SatellitePosition_GPS_GPSws(ts_sta1sat1, DD_data.sat1, nav_data)
            dt_sta1sat1 = DD_data.obs1_sat1.data['P2']['observation'] / c
            Xeci_sta1sat1, Yeci_sta1sat1, Zeci_sta1sat1 = CoorTransform.earth_rotation_correction([coorX_sta1sat1, coorY_sta1sat1, coorZ_sta1sat1], dt_sta1sat1)
            lou_sta1sat1 = CoorTransform.cal_distance(sta1_coor, [Xeci_sta1sat1, Yeci_sta1sat1, Zeci_sta1sat1])

            # 计算站1星2元素
            ts_sta1sat2, dts_sta1sat2 = SPP.cal_EmitTime_from_datetime(DD_data.T, DD_data.sat2, DD_data.obs1_sat2.data['P2']['observation'], nav_data, doCRC=True)
            coorX_sta1sat2, coorY_sta1sat2, coorZ_sta1sat2 = SatellitePosition.cal_SatellitePosition_GPS_GPSws(ts_sta1sat2, DD_data.sat2, nav_data)
            dt_sta1sat2 = DD_data.obs1_sat2.data['P2']['observation'] / c
            Xeci_sta1sat2, Yeci_sta1sat2, Zeci_sta1sat2 = CoorTransform.earth_rotation_correction([coorX_sta1sat2, coorY_sta1sat2, coorZ_sta1sat2], dt_sta1sat2)
            lou_sta1sat2 = CoorTransform.cal_distance(sta1_coor, [Xeci_sta1sat2, Yeci_sta1sat2, Zeci_sta1sat2])

            # 计算站2星1元素
            ts_sta2sat1, dts_sta2sat1 = SPP.cal_EmitTime_from_datetime(DD_data.T, DD_data.sat1, DD_data.obs2_sat1.data['P2']['observation'], nav_data, doCRC=True)
            coorX_sta2sat1, coorY_sta2sat1, coorZ_sta2sat1 = SatellitePosition.cal_SatellitePosition_GPS_GPSws(ts_sta2sat1, DD_data.sat1, nav_data)
            dt_sta2sat1 = DD_data.obs2_sat1.data['P2']['observation'] / c
            Xeci_sta2sat1, Yeci_sta2sat1, Zeci_sta2sat1 = CoorTransform.earth_rotation_correction([coorX_sta2sat1, coorY_sta2sat1, coorZ_sta2sat1], dt_sta2sat1)
            lou_sta2sat1 = CoorTransform.cal_distance(sta2_coor, [Xeci_sta2sat1, Yeci_sta2sat1, Zeci_sta2sat1])

            # 计算站2星2元素
            ts_sta2sat2, dts_sta2sat2 = SPP.cal_EmitTime_from_datetime(DD_data.T, DD_data.sat2, DD_data.obs2_sat2.data['P2']['observation'], nav_data, doCRC=True)
            coorX_sta2sat2, coorY_sta2sat2, coorZ_sta2sat2 = SatellitePosition.cal_SatellitePosition_GPS_GPSws(ts_sta2sat2, DD_data.sat2, nav_data)
            dt_sta2sat2 = DD_data.obs2_sat2.data['P2']['observation'] / c
            Xeci_sta2sat2, Yeci_sta2sat2, Zeci_sta2sat2 = CoorTransform.earth_rotation_correction([coorX_sta2sat2, coorY_sta2sat2, coorZ_sta2sat2], dt_sta2sat2)
            lou_sta2sat2 = CoorTransform.cal_distance(sta2_coor, [Xeci_sta2sat2, Yeci_sta2sat2, Zeci_sta2sat2])

            # 构造几何系数阵
            X2, Y2, Z2 = sta2_coor
            a_sta2_X = (X2 - Xeci_sta2sat2) / lou_sta2sat2 - (X2 - Xeci_sta2sat1) / lou_sta2sat1
            a_sta2_Y = (Y2 - Yeci_sta2sat2) / lou_sta2sat2 - (Y2 - Yeci_sta2sat1) / lou_sta2sat1
            a_sta2_Z = (Z2 - Zeci_sta2sat2) / lou_sta2sat2 - (Z2 - Zeci_sta2sat1) / lou_sta2sat1

            H11.append([a_sta2_X, a_sta2_Y, a_sta2_Z])
            H21.append([a_sta2_X, a_sta2_Y, a_sta2_Z])
            hx.append(lou_sta2sat2 - lou_sta1sat2 - lou_sta2sat1 + lou_sta1sat1)

        n = len(self.DD_records.DD_observations_data[band])
        H12 = np.eye(n) * get_lamb_from_band(band)
        H22 = np.zeros((n, n))
        H = np.block([[np.array(H11), H12], [np.array(H21), H22]])
        hx = np.array(2*hx)
        return H, hx

    def gety(self, band1, band2):
        y1 = []
        y2 = []
        for DD_data in self.DD_records.DD_observations_data[band1]:
            y1.append(get_lamb_from_band(band1)*DD_data.DD_obs)
        for DD_data in self.DD_records.DD_observations_data[band2]:
            y2.append(DD_data.DD_obs)
        y = np.array(y1+y2)
        return y

    def getx(self, band, sta2_coor):
        # 转换类型
        if isinstance(sta2_coor, np.ndarray):
            sta2_coor = sta2_coor.tolist()
        if isinstance(self.DD_records.DD_ambiguitys[band], np.ndarray):
            DD_ambiguitys = self.DD_records.DD_ambiguitys[band].tolist()
        else:
            DD_ambiguitys = self.DD_records.DD_ambiguitys[band]
        x = np.array(sta2_coor + DD_ambiguitys)
        return x

    def getQ(self, band):
        Qs = self.DD_records.DD_ambiguitys_noise[band]
        # Q = np.diag([0.3, 0.5, 0.6] + Qs)
        Q = np.diag([10, 10, 10] + Qs)
        return Q

    def getR(self, band1, band2, sigma1=0.002, sigma2=5):
        nDD = len(self.DD_records.DD_observations_data[band1])
        covDD = np.full((nDD, nDD), 1).astype(float)
        for i in range(nDD):
            covDD[i, i] = 2
        covDD1 = 2 * sigma1**2 * covDD
        covDD2 = 2 * sigma2**2 * covDD
        R = RTK.diagonalize_squarematrix(covDD1, covDD2)
        return R

    def ekf_estimation(self, P_before, sta1_coor, sta2_coor, nav_records):
        # 预测
        F = self.getF(self.band1)
        Q = self.getQ(self.band1)
        x_pri = self.getx(self.band1, sta2_coor)
        P_pri = F @ P_before @ F.T + Q

        # 更新
        H, hx1 = self.getH(nav_records, sta1_coor, sta2_coor, self.band1)
        hx2 = get_lamb_from_band(self.band1) * np.array(x_pri[3:].tolist()+[0 for i in range(len(x_pri[3:].tolist()))])
        hx = hx1+hx2
        R = self.getR(self.band1, self.band2)
        y = self.gety(self.band1, self.band2)
        K = P_pri @ H.T @ np.linalg.inv(H @ P_pri @ H.T + R)
        x_est = x_pri + K @ (y - hx)
        # x_est = x_pri + K @ (y - H @ x_pri)
        P_est = (np.eye(len(x_est)) - K @ H) @ P_pri

        return x_est, P_est


def ambiguity_resolution_solution(Q, x):
    Qaa = RTK.get_symmetric_matrix(Q[3:, 3:])
    Qba = Q[:3, 3:]
    N_fixed, sqnorm, Ps, Qzhat, Z, nfixed, mu = LAMBDA.main(x[3:], Qaa)
    # 更新参数估计
    b_hat = x[:3]
    a_hat = x[3:]
    b_check = RTK.MAPmethod(b_hat, a_hat, Qaa, Qba, N_fixed[:, 0])
    return b_check

def get_diaglist_from_matrix(matrix):
    """
    获取方阵对角元素并组成列表
    """
    le = len(matrix)
    diag_list = []
    for i in range(le):
        diag_list.append(matrix[i, i])
    return diag_list


# def obs_isnot_null(band, obs1_sat1, obs1_sat2, obs2_sat1, obs2_sat2):
#     # 判断波段数据是否完整
#     if not (obs1_sat1.data.__contains__(band) and obs2_sat1.data.__contains__(band)):
#         isnot_null_flag = False
#     elif obs1_sat1.data[band]['observation'] == '' or obs1_sat2.data[band]['observation'] == ''\
#             or obs2_sat1.data[band]['observation'] == '' or obs2_sat2.data[band]['observation'] == '':
#         isnot_null_flag = False
#     else:
#         isnot_null_flag = True
#     return isnot_null_flag
#
# def get_DD_GF_from_observationdata(svns, band1, band2, obs1_sat1, obs1_sat2, obs2_sat1, obs2_sat2):
#     GF_band = [band1, band2]
#     DD_GF = {}
#     for svn in svns:
#         if obs_isnot_null(band1) and obs_isnot_null(band2):
#             # 分别进行双差观测值构造
#             obs_sta1sat1_band1 = float(obs1_sat1.data[band1]['observation'])
#             obs_sta1sat2_band1 = float(obs1_sat2.data[band1]['observation'])
#             obs_sta2sat1_band1 = float(obs2_sat1.data[band1]['observation'])
#             obs_sta2sat2_band1 = float(obs2_sat2.data[band1]['observation'])
#             DD_obs_band1 = obs_sta2sat2_band1 - obs_sta1sat2_band1 -obs_sta2sat1_band1 +obs_sta1sat1_band1
#             obs_sta1sat1_band2 = float(obs1_sat1.data[band2]['observation'])
#             obs_sta1sat2_band2 = float(obs1_sat2.data[band2]['observation'])
#             obs_sta2sat1_band2 = float(obs2_sat1.data[band2]['observation'])
#             obs_sta2sat2_band2 = float(obs2_sat2.data[band2]['observation'])
#             DD_obs_band2 = obs_sta2sat2_band2 - obs_sta1sat2_band2 -obs_sta2sat1_band2 +obs_sta1sat1_band2
#             DD_GF[svn] = DD_obs_band1 - DD_obs_band2
#         else:
#             DD_GF[svn] = 100


def constaneously_RTK_withfilter(start_time, end_time, knownStation_ob_records, unknownStation_ob_records, br_records, knownStation_coor, unknownStation_init_coor):
    """
    start_time : datetime.datetime , 时间段开始时刻
    end_time : datetime.datetime , 时间段结束时刻
    knownStation_ob_records : 基准站多时刻观测数据
    unknownStation_ob_records : 未知站多时刻观测数据
    br_records : 卫星星历数据
    knownStation_coor : 基准站坐标
    unknownStation_init_coor : 未知站坐标
    """
    interval_time = 30
    Tr = start_time
    coordinates = []
    # 第一个历元进行最小二乘解算
    sta2_coor, P, N_float, base_svn, diff_svns, Pse = RTK2.DD_onCarrierPhase_and_Pseudorange_1known(knownStation_ob_records, unknownStation_ob_records, br_records, Tr, knownStation_coor, unknownStation_init_coor, ambi_fix=False)
    # 给初值
    # sta2_coor = unknownStation_init_coor
    # N_float = [i/get_lamb_from_band('L1') for i in Pse]
    # P = np.diag([10, 10, 10]+[1000 for i in range(len(P)-3)])
    coordinates.append(sta2_coor)
    # 进行传递
    if isinstance(sta2_coor, np.ndarray):
        sta2_coor = sta2_coor.tolist()
    if isinstance(N_float, np.ndarray):
        N_float = N_float.tolist()
    amb_trans = ambiguity_transfer('L1',  N_float, base_svn, diff_svns, get_diaglist_from_matrix(P)[3:])
    # amb_trans = ambiguity_transfer('L1', N_float, base_svn, diff_svns, [10 for i in range(len(P) - 3)])
    # amb_trans = ambiguity_transfer('L1', N_float, base_svn, diff_svns, get_diaglist_from_matrix(P[3:, 3:]))
    # DD_GF_collection = DD_records_collection.DD_GF
    DD_GF_collection = ''
    Tr += datetime.timedelta(seconds=interval_time)
    n = 0
    while Tr < end_time:
        # 找到Tr时刻两个观测站观测记录
        knownStation_ob_records_atTr = list(filter(lambda o:o.time == Tr and o.SVN[0] == "G" and o.data != "", knownStation_ob_records))
        unknownStation_ob_records_atTr = list(filter(lambda o:o.time == Tr and o.SVN[0] == "G" and o.data != "", unknownStation_ob_records))
        # 获取符合条件的所有svn
        svn_Tr = []
        knownsta_svn = [x.SVN for x in knownStation_ob_records_atTr]
        unknownsta_svn = [x.SVN for x in unknownStation_ob_records_atTr]
        for svn in knownsta_svn:
            if svn in unknownsta_svn:
                svn_Tr.append(svn)
        # 初始化
        DD_records_collection = DD_records_atTr(Tr, ['L1', 'P2'])
        trans_matrix = DD_records_collection.set_base_satellite_from_ambiguity_transfer(amb_trans, knownStation_ob_records_atTr, unknownStation_ob_records_atTr)
        svn_removed_index, new_in = DD_records_collection.add_satellites(svn_Tr, knownStation_ob_records_atTr, unknownStation_ob_records_atTr, DD_GF_collection)
        print("-2-:", svn_removed_index, new_in)

        # 判断是否发生了参考星的变化
        if np.all(trans_matrix == 0):    # 参考星不变
            # 调整较上个历元被去除的卫星对应的vc阵
            if svn_removed_index:  # 上个历元的卫星有缺失
                P = remove_cross_from_matrix(P, svn_removed_index, svn_removed_index)
        else:      # 参考星变化
            # 调整较上个历元被去除的卫星对应的vc阵
            if svn_removed_index:  # 上个历元的卫星有缺失
                P = remove_cross_from_matrix(P, svn_removed_index, svn_removed_index)
                trans_matrix = np.delete(trans_matrix, svn_removed_index, axis=1)
            P_trans_matrix = RTK.diagonalize_squarematrix(np.eye(3), trans_matrix)
            P = P_trans_matrix @ P @ P_trans_matrix.T
        # 更新由卫星数量增加导致的模糊度参数vc阵的变化
        if new_in != 0:  # 此历元有新增卫星
            eyes = 10000 * np.eye(new_in)    # 前面的常数为新增卫星的方差
            P = RTK.diagonalize_squarematrix(P, eyes)
        # 开始估计
        KM_estimator = Kalman_Filter(DD_records_collection, 'L1', 'P2')
        x, P = KM_estimator.ekf_estimation(P, sta1_coor=knownStation_coor, sta2_coor=sta2_coor, nav_records=br_records)
        print(Tr, x)
        DD_records_collection.DD_ambiguitys['L1'] = x[3:].tolist()
        ambiguitys = DD_records_collection.DD_ambiguitys['L1']
        base_svn = DD_records_collection.base_svn
        satellites = DD_records_collection.qualitified_satellites
        amb_nos1 = [100 for i in range(len(P) - 3)]
        amb_nos2 = [10 for i in range(len(P) - 3)]
        amb_nos3 = [1 for i in range(len(P) - 3)]
        amb_nos4 = [0.000001 for i in range(len(P) - 3)]
        # 模糊度噪声变化的算法
        # if n < 50:
        #     amb_trans = ambiguity_transfer('L1', ambiguitys, base_svn, satellites, get_diaglist_from_matrix(P)[3:])
        #     # amb_trans = ambiguity_transfer('L1', ambiguitys, base_svn, satellites, amb_nos)
        # # elif n < 20:
        # #     amb_trans = ambiguity_transfer('L1', ambiguitys, base_svn, satellites, amb_nos1)
        # # elif n < 30:
        # #     amb_trans = ambiguity_transfer('L1', ambiguitys, base_svn, satellites, amb_nos2)
        # # elif n < 40:
        # #     amb_trans = ambiguity_transfer('L1', ambiguitys, base_svn, satellites, amb_nos2)
        # else:
        #     amb_trans = ambiguity_transfer('L1', ambiguitys, base_svn, satellites, amb_nos4)
            # amb_trans = ambiguity_transfer('L1', ambiguitys, base_svn, satellites, amb_nos)
        # amb_trans = ambiguity_transfer('L1', ambiguitys, base_svn, satellites, get_diaglist_from_matrix(P)[3:])
        amb_trans = ambiguity_transfer('L1', ambiguitys, base_svn, satellites, amb_nos4)
        DD_GF_collection = DD_records_collection.DD_GF

        # 固定模糊度
        # coor = ambiguity_resolution_solution(P, x)
        coor = x[:3].tolist()
        coordinates.append(coor)
        Tr += datetime.timedelta(seconds=interval_time)
        n += 1

    return coordinates


def remove_cross_from_matrix(matrix, row_index_list, col_index_list):
    matrix_row_deleted = np.delete(matrix, row_index_list, axis=0)
    matrix_deleted = np.delete(matrix_row_deleted, col_index_list, axis=1)
    return matrix_deleted



if __name__ == '__main__':

    # station1_observation_file = r"edata\obs\leij3100.20o"    # 已知站点 leij
    # station2_observation_file = r"edata\obs\zim23100.20o"    # 未知站点 zim2
    station2_observation_file = r"edata\obs\zimm3100.20o"    # 未知站点 zimm
    # station1_observation_file = r"edata\obs\wab23100.20o"    # 已知站点 wab2
    station1_observation_file = r"edata\obs\zim23100.20o"  # 已知站点 zim2
    # station1_observation_file = r"edata\obs\zimm3100.20o"  # 已知站点 zimm
    broadcast_file = r"edata\sat_obit\brdc3100.20n"
    # 读入观测文件内容,得到类型对象列表
    knownStation_ob_records = DoFile.read_Rinex2_oFile(station1_observation_file)
    unknownStation_ob_records = DoFile.read_Rinex2_oFile(station2_observation_file)
    br_records = DoFile.read_GPS_nFile(broadcast_file)
    print("数据读取完毕！")
    knownStation_coor = [4331300.1600, 567537.0810, 4633133.5100]  # zim2
    # knownStation_coor = [4327318.2325, 566955.9585, 4636425.9246]  # wab2
    init_coor = [4331297.3480, 567555.6390, 4633133.7280]  # zimm


    strat_time = datetime.datetime(2020, 11, 5, 3, 0, 0)
    end_time = datetime.datetime(2020, 11, 5, 9, 59, 0)
    cal_coors = constaneously_RTK_withfilter(strat_time, end_time, knownStation_ob_records, unknownStation_ob_records, br_records, knownStation_coor, init_coor)
    true_coors = [init_coor for i in range(len(cal_coors))]

    SPP.cal_NEUerrors(true_coors, cal_coors)
    SPP.cal_XYZerrors(true_coors, cal_coors)
    print("neu各方向RMSE:", ResultAnalyse.get_NEU_rmse(true_coors, cal_coors))
    print("坐标RMSE:", ResultAnalyse.get_coor_rmse(true_coors, cal_coors))






























