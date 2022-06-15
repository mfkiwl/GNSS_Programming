# -*- coding: utf-8 -*-
"""

@title:	RTK with filter
@author: iDeal0103
@status:	Active
@type:	Process
@created:	11-March-2022
@post-History:	11-March-2022

comment:
    1. RTK滤波解算法
    2. 附带基线约束的RTK滤波算法

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


# 模糊度固定后的参数估计值更新
def ambiguity_resolution_solution(Q, x):
    Qaa = RTK.get_symmetric_matrix(Q[3:, 3:])
    Qba = Q[:3, 3:]
    N_fixed, sqnorm, Ps, Qzhat, Z, nfixed, mu = LAMBDA.main(x[3:], Qaa)
    # 更新参数估计
    b_hat = x[:3]
    a_hat = x[3:]
    b_check = RTK.MAPmethod(b_hat, a_hat, Qaa, Qba, N_fixed[:, 0])
    return b_check

# 获得对角阵元素
def get_diaglist_from_matrix(matrix):
    """
    获取方阵对角元素并组成列表
    """
    row, col = matrix.shape
    if row != col:
        print("需要一个方阵！")
        raise SystemExit
    else:
        diag_list = []
        for i in range(row):
            diag_list.append(matrix[i, i])
        return diag_list

# 将矩阵中某一行和列去掉
def remove_cross_from_matrix(matrix, row_index_list, col_index_list):
    matrix_row_deleted = np.delete(matrix, row_index_list, axis=0)
    matrix_deleted = np.delete(matrix_row_deleted, col_index_list, axis=1)
    return matrix_deleted

# 判断某个波段是否为载波相位波段
def band_is_carrierphase(band_name):
    flag = False
    if band_name[-1] == 'L':
        flag = True
    return flag


# 双差观测值类
class DD_recordclass():
    def __init__(self, Tr, sat1, sat2, band):
        self.sat1 = sat1
        self.sat2 = sat2
        self.band = band
        self.T = Tr

    def get_DD_obs(self, obs1_sat1, obs1_sat2, obs2_sat1, obs2_sat2):
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

    def get_DD_GF(self, GF_band2, obs1_sat1, obs1_sat2, obs2_sat1, obs2_sat2):
        GF_band1 = self.band
        self.GF_band = [GF_band1, GF_band2]
        # 分别进行双差观测值构造
        obs_sta1sat1_band1 = float(obs1_sat1.data[GF_band1]['observation'])
        obs_sta1sat2_band1 = float(obs1_sat2.data[GF_band1]['observation'])
        obs_sta2sat1_band1 = float(obs2_sat1.data[GF_band1]['observation'])
        obs_sta2sat2_band1 = float(obs2_sat2.data[GF_band1]['observation'])
        DD_obs_band1 = obs_sta2sat2_band1 - obs_sta1sat2_band1 -obs_sta2sat1_band1 + obs_sta1sat1_band1
        obs_sta1sat1_band2 = float(obs1_sat1.data[GF_band2]['observation'])
        obs_sta1sat2_band2 = float(obs1_sat2.data[GF_band2]['observation'])
        obs_sta2sat1_band2 = float(obs2_sat1.data[GF_band2]['observation'])
        obs_sta2sat2_band2 = float(obs2_sat2.data[GF_band2]['observation'])
        DD_obs_band2 = obs_sta2sat2_band2 - obs_sta1sat2_band2 -obs_sta2sat1_band2 +obs_sta1sat1_band2
        self.DD_GF = get_lamb_from_band(GF_band1) * DD_obs_band1 - get_lamb_from_band(GF_band2) * DD_obs_band2



# 记录参考站与若干待差星的GF组合双差观测值数据
class DD_GF_record():
    def __init__(self, base_svn):
        self.base_svn = base_svn
        self.GF_records = {}

    def add_DD_GF(self, svn, value):
        self.GF_records[svn] = value



# 适用于历元间传递双差数据的类
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
            if band_is_carrierphase(band):
                self.DD_ambiguitys[band] = []
                self.DD_ambiguitys_noise[band] = []
        self.init_amb_noise = 10000000     # 模糊度初始值对应的噪声方差
        self.DD_GF = {}

    def set_base_satellite(self, base_SVN, sta1obs_records_Tr, sta2obs_records_Tr):
        self.base_svn = base_SVN
        self.base_sta1obs_record = list(filter(lambda o: o.SVN == base_SVN, sta1obs_records_Tr))[0]
        self.base_sta2obs_record = list(filter(lambda o: o.SVN == base_SVN, sta2obs_records_Tr))[0]

    # 增加新卫星进行双差
    def add_satellite(self, the_svn, sta1obs_records_Tr, sta2obs_records_Tr):
        # 确定GF周跳探测波段
        band1 = 'L1_L'
        band2 = 'L2_L'
        if self.base_svn:
            other_sta1obs_record = list(filter(lambda o: o.SVN == the_svn, sta1obs_records_Tr))[0]
            other_sta2obs_record = list(filter(lambda o: o.SVN == the_svn, sta2obs_records_Tr))[0]
            # 对卫星数据进行检查
            qualified_flag = True
            DD_obses = []
            for band in self.bands:
                DD_record = DD_recordclass(self.Tr, self.base_svn, the_svn, band)
                DD_record.get_DD_obs(self.base_sta1obs_record, other_sta1obs_record, self.base_sta2obs_record, other_sta2obs_record)
                if band_is_carrierphase(band):
                    DD_record.get_DD_GF(band2, self.base_sta1obs_record, other_sta1obs_record, self.base_sta2obs_record, other_sta2obs_record)  # 计算GF双差值
                if DD_record.DD_obs == '':
                    qualified_flag = False
                    break
                DD_obses.append(DD_record)
            # 根据卫星数据合格与否进行是否加入待差卫星的判断, 并增加数据
            if qualified_flag:
                self.qualitified_satellites.append(the_svn)
                for band in self.bands:
                    self.DD_observations_data[band].append(DD_obses[self.bands.index(band)])
                    if band_is_carrierphase(band):
                        lamb = get_lamb_from_band(band)
                        self.DD_ambiguitys[band].append(DD_obses[self.bands.index('L2_C')].DD_obs/lamb)
                        self.DD_ambiguitys_noise[band].append(self.init_amb_noise)
            else:
                self.unqualitified_satellites.append(the_svn)
                print(the_svn + '卫星数据不符合要求！')
        else:
            qualified_flag = False
            print('请设置参考星！')
        return qualified_flag

    def get_DD_GF_record(self, band):
        for data in self.DD_observations_data[band]:
            self.DD_GF[data.sat2] = data.DD_GF

    def add_satellites(self, svns, sta1obs_records_Tr, sta2obs_records_Tr, DD_GF_before={}):
        # 确定GF所使用的波段名称
        band1 = 'L1_L'
        band2 = 'L2_L'
        qualitified_satellites_copy = self.qualitified_satellites.copy()
        if self.qualitified_satellites:
            # 先获得所有数据缺失卫星的原始索引
            svn_satellite_miss = []
            for svn in self.qualitified_satellites:
                if svn not in svns:
                    svn_satellite_miss.append(svn)
            # 将所有数据缺失的上历元合格卫星删除
            self.delete_from_qualitified_svn(svn_satellite_miss)
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
                    if band_is_carrierphase(band):
                        DD_record.get_DD_GF(band2, self.base_sta1obs_record, sta1obs_records_Tr_svn, self.base_sta2obs_record,
                                             sta2obs_records_Tr_svn)    # 计算GF双差值
                    if DD_record.DD_obs == '':
                        qualified_flag = False
                        break
                    DD_obses.append(DD_record)
                # 根据上述卫星数据检查判断是否加入待差卫星
                if qualified_flag:
                    print("观测数据加入卫星", svn)
                    for band in self.bands:
                        self.DD_observations_data[band].append(DD_obses[self.bands.index(band)])
                else:
                    svn_data_incomplete.append(svn)
                    print(svn + '卫星数据不符合要求！')
            # 将数据不完整的对象删除（因为不完整的数据没有加进来，所以不需要删除索引对应的观测数据）
            if svn_data_incomplete:
                self.delete_from_qualitified_svn(svn_data_incomplete)
            # 进行周跳探测, 将认为有周跳的卫星删除
            cycle_slip_svn = []
            self.get_DD_GF_record(band1)
            if DD_GF_before:
                for svn in self.qualitified_satellites:
                    GF_before = DD_GF_before[svn]
                    GF_now = self.DD_GF[svn]
                    if CycleslipsDetection.GF_detector(GF_before, GF_now, threshold=1):
                        cycle_slip_svn.append(svn)
                        print(svn + '卫星发生周跳！deltaGF=', abs(GF_before-GF_now))
                if cycle_slip_svn:
                    self.delete_withrecord_from_qualitified_svn(cycle_slip_svn)
            # 获得所有被移除的上一历元合格卫星的索引
            svn_removed = svn_satellite_miss + svn_data_incomplete + cycle_slip_svn
            svn_removed_index = []
            for svn in svn_removed:
                svn_index = qualitified_satellites_copy.index(svn)
                svn_removed_index.append(svn_index)
            # 将新卫星加入
            new_in = 0
            new_in_svn = []
            for svn in svns:
                if svn not in self.qualitified_satellites and svn not in self.unqualitified_satellites and svn != self.base_svn:
                    new_in_flag = self.add_satellite(svn, sta1obs_records_Tr, sta2obs_records_Tr)
                    if new_in_flag:
                        print("观测数据加入卫星", svn)
                        new_in_svn.append(svn)
                        new_in += 1
            # 获取所有卫星的GF双差观测值
            self.get_DD_GF_record(band1)
        else:
            new_in = 0
            new_in_svn = []
            for svn in svns:
                if svn != self.base_svn:
                    new_in_flag = self.add_satellite(svn, sta1obs_records_Tr, sta2obs_records_Tr)
                    if new_in_flag:
                        new_in_svn.append(svn)
                        new_in += 1
        return svn_removed_index, new_in, new_in_svn

    def delete_from_qualitified_svn(self, svns):
        if svns:
            for svn in svns:
                svn_index = self.qualitified_satellites.index(svn)
                self.qualitified_satellites.remove(svn)
                self.unqualitified_satellites.append(svn)
                for band in self.bands:
                    if band_is_carrierphase(band):
                        # 删除模糊度和模糊度参数噪声
                        self.DD_ambiguitys[band].pop(svn_index)
                        self.DD_ambiguitys_noise[band].pop(svn_index)

    def delete_withrecord_from_qualitified_svn(self, svns):
        if svns:
            for svn in svns:
                svn_index = self.qualitified_satellites.index(svn)
                self.qualitified_satellites.remove(svn)    # 不加入self.unqualified_satellites的原因是后续还可加入该卫星的记录
                for band in self.bands:
                    # 删除双差观测记录（如有）
                    self.DD_observations_data[band].pop(svn_index)
                    if band_is_carrierphase(band):
                        # 删除模糊度和模糊度参数噪声
                        self.DD_ambiguitys[band].pop(svn_index)
                        self.DD_ambiguitys_noise[band].pop(svn_index)

    def set_base_satellite_from_ambiguity_transfer(self, ambiguity_transfer, sta1obs_records_Tr, sta2obs_records_Tr, ele_limit=10):
        base_sta1obs_record = list(filter(lambda o: o.SVN == ambiguity_transfer.base_svn, sta1obs_records_Tr))
        base_sta2obs_record = list(filter(lambda o: o.SVN == ambiguity_transfer.base_svn, sta2obs_records_Tr))
        transfer_matrix = np.array([0])
        # 判断上一历元继承的参考星存在及本历元中数据完整, 否则需要更换参考星
        if not (base_sta1obs_record and base_sta2obs_record):
            # 更换参考星
            print('因为数据不完整更换参考星!原参考星为：' + ambiguity_transfer.base_svn)
            # 选择卫星高度角最大且数据合格的那个
            ele_dict = ambiguity_transfer.satellites_ele.copy()
            for i in range(len(ambiguity_transfer.satellites_ele)):
                the_svn = max(ele_dict, key=ele_dict.get)
                sta1obs_record = list(filter(lambda o: o.SVN == the_svn, sta1obs_records_Tr))
                sta2obs_record = list(filter(lambda o: o.SVN == the_svn, sta2obs_records_Tr))
                if sta1obs_record and sta2obs_record:
                    break
                else:
                    ele_dict.pop(the_svn)
            self.set_base_satellite(the_svn, sta1obs_record, sta2obs_record)
            transfer_matrix = ambiguity_transfer.shift_to_new_basesats(the_svn)
            self.qualitified_satellites = ambiguity_transfer.satellites
            self.DD_ambiguitys[ambiguity_transfer.band] = ambiguity_transfer.ambiguitys
            self.DD_ambiguitys_noise[ambiguity_transfer.band] = ambiguity_transfer.amb_noises
            print('更换参考星!现参考星为：' + self.base_svn)
        # 参考星数据完整，再判断高度角是否满足条件
        else:
            # 判断是否满足高度角截止条件, 不满足则更换参考星
            if ambiguity_transfer.satellites_ele[ambiguity_transfer.base_svn] < ele_limit:
                print('因高度角太小更换参考星!原参考星为：' + ambiguity_transfer.base_svn)
                # 选择卫星高度角最大且数据合格的那个
                ele_dict = ambiguity_transfer.satellites_ele.copy()
                while True:
                    the_svn = max(ele_dict, key=ele_dict.get)
                    sta1obs_record = list(filter(lambda o: o.SVN == the_svn, sta1obs_records_Tr))
                    sta2obs_record = list(filter(lambda o: o.SVN == the_svn, sta2obs_records_Tr))
                    if sta1obs_record and sta2obs_record:
                        break
                    else:
                        ele_dict.pop(the_svn)
                self.set_base_satellite(the_svn, sta1obs_record, sta2obs_record)
                transfer_matrix = ambiguity_transfer.shift_to_new_basesats(the_svn)
                self.qualitified_satellites = ambiguity_transfer.satellites
                self.DD_ambiguitys[ambiguity_transfer.band] = ambiguity_transfer.ambiguitys
                self.DD_ambiguitys_noise[ambiguity_transfer.band] = ambiguity_transfer.amb_noises
                print('更换参考星!现参考星为：' + self.base_svn)
            else:
                self.set_base_satellite(ambiguity_transfer.base_svn, base_sta1obs_record, base_sta2obs_record)
                self.qualitified_satellites = ambiguity_transfer.satellites
                self.DD_ambiguitys[ambiguity_transfer.band] = ambiguity_transfer.ambiguitys
                self.DD_ambiguitys_noise[ambiguity_transfer.band] = ambiguity_transfer.amb_noises
        # set_base_satellite_from_ambiguity_transfer 与 add_satellites 配合使用
        return transfer_matrix

    def set_elebest_base_satellite_from_ambiguity_transfer(self, ambiguity_transfer, sta1obs_records_Tr, sta2obs_records_Tr):
        transfer_matrix = np.array([0])
        # 寻找高度角最高的卫星作为参考星
        ele_dict = ambiguity_transfer.satellites_ele.copy()
        for i in range(len(ambiguity_transfer.satellites_ele)):
            the_svn = max(ele_dict, key=ele_dict.get)
            sta1obs_record = list(filter(lambda o: o.SVN == the_svn, sta1obs_records_Tr))
            sta2obs_record = list(filter(lambda o: o.SVN == the_svn, sta2obs_records_Tr))
            if sta1obs_record and sta2obs_record:
                break
            else:
                ele_dict.pop(the_svn)
        if ambiguity_transfer.base_svn != the_svn:
            self.set_base_satellite(the_svn, sta1obs_record, sta2obs_record)
            transfer_matrix = ambiguity_transfer.shift_to_new_basesats(the_svn)
            self.qualitified_satellites = ambiguity_transfer.satellites
            self.DD_ambiguitys[ambiguity_transfer.band] = ambiguity_transfer.ambiguitys
            self.DD_ambiguitys_noise[ambiguity_transfer.band] = ambiguity_transfer.amb_noises
            print('更换参考星!现参考星为：' + self.base_svn)
        else:
            self.set_base_satellite(ambiguity_transfer.base_svn, sta1obs_record, sta2obs_record)
            self.qualitified_satellites = ambiguity_transfer.satellites
            self.DD_ambiguitys[ambiguity_transfer.band] = ambiguity_transfer.ambiguitys
            self.DD_ambiguitys_noise[ambiguity_transfer.band] = ambiguity_transfer.amb_noises
        return transfer_matrix



class ambiguity_transfer():
    def __init__(self, band, ambiguitys, base_svn, satellites, amb_noises, rover_station=""):
        self.band = band
        self.base_svn = base_svn
        self.ambiguitys = ambiguitys
        self.satellites = satellites
        self.amb_noises = amb_noises
        self.satellites_ele = {}
        self.rover_station = rover_station

    def update_ele(self, rover_station_coor, t, nav_data):
        w, s = TimeSystem.from_datetime_cal_GPSws(t)
        t = TimeSystem.GPSws(w, s)
        # 计算参考星的高度角
        base_sat_coor = SatellitePosition.cal_SatellitePosition_GPS_GPSws(t, self.base_svn, nav_data)
        self.satellites_ele[self.base_svn] = CoorTransform.cal_ele_and_A(rover_station_coor, base_sat_coor)[0] * 180 / math.pi
        # 计算各个待差星位置
        for svn in self.satellites:
            sat_coor = SatellitePosition.cal_SatellitePosition_GPS_GPSws(t, svn, nav_data)
            self.satellites_ele[svn] = CoorTransform.cal_ele_and_A(rover_station_coor, sat_coor)[0] * 180 / math.pi

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

    def delete_amb_svn(self, svn):
        if svn in self.satellites:
            svn_index = self.satellites.index(svn)
            self.satellites.pop(svn_index)
            self.ambiguitys.pop(svn_index)
            self.satellites_ele.pop(svn)
            self.amb_noises.pop(svn_index)
        else:
            print('上一历元不含有该卫星数据！')

    def imitate_ambtransfer(self, ambtransfer):
        new_satellites = []
        new_ambiguitys = []
        new_amb_noises = []
        unabled_svns = []    # 记录无法参与解算的卫星
        for svn in ambtransfer.satellites:
            if svn in self.satellites:
                svn_index = self.satellites.index(svn)
                new_satellites.append(svn)
                new_ambiguitys.append(self.ambiguitys[svn_index])
                new_amb_noises.append(self.amb_noises[svn_index])
            else:
                unabled_svns.append(svn)
                print("缺少卫星"+svn)
        return unabled_svns


# 卡尔曼滤波器
class Kalman_Filter():
    def __init__(self, DD_records_atTr, cp_band, pr_band, mode_int, residual_manager=None, DDambiguity_manager=None, position_manager=None):
        """
        DD_records_atTr : DD_records_atTr class, 双差观测集合
        cp_band : 载波相位
        pr_band : 伪距
        mode_int : 滤波解模式，0为静态，1为动态
        """
        self.DD_records = DD_records_atTr
        self.cp_band = cp_band
        self.pr_band = pr_band
        self.mode = mode_int
        self.residual_manager = residual_manager
        self.DDambiguity_manager = DDambiguity_manager
        self.position_manager = position_manager

    def getF(self):
        n = len(self.DD_records.DD_observations_data[self.cp_band]) + 3
        Fmatrix = np.diag([1 for i in range(n)])
        return Fmatrix

    def getH(self, nav_data, sta1_coor, sta2_coor, x_array):
        H11 = []
        H21 = []
        hx = []
        for DD_data in self.DD_records.DD_observations_data[self.cp_band]:
            # 计算站1星1元素
            ts_sta1sat1, dts_sta1sat1 = SPP.cal_EmitTime_from_datetime(DD_data.T, DD_data.sat1, DD_data.obs1_sat1.data[self.pr_band]['observation'], nav_data, doCRC=True)
            coorX_sta1sat1, coorY_sta1sat1, coorZ_sta1sat1 = SatellitePosition.cal_SatellitePosition_GPS_GPSws(ts_sta1sat1, DD_data.sat1, nav_data)
            dt_sta1sat1 = DD_data.obs1_sat1.data[self.pr_band]['observation'] / c
            Xeci_sta1sat1, Yeci_sta1sat1, Zeci_sta1sat1 = CoorTransform.earth_rotation_correction([coorX_sta1sat1, coorY_sta1sat1, coorZ_sta1sat1], dt_sta1sat1)
            lou_sta1sat1 = CoorTransform.cal_distance(sta1_coor, [Xeci_sta1sat1, Yeci_sta1sat1, Zeci_sta1sat1])

            # 计算站1星2元素
            ts_sta1sat2, dts_sta1sat2 = SPP.cal_EmitTime_from_datetime(DD_data.T, DD_data.sat2, DD_data.obs1_sat2.data[self.pr_band]['observation'], nav_data, doCRC=True)
            coorX_sta1sat2, coorY_sta1sat2, coorZ_sta1sat2 = SatellitePosition.cal_SatellitePosition_GPS_GPSws(ts_sta1sat2, DD_data.sat2, nav_data)
            dt_sta1sat2 = DD_data.obs1_sat2.data[self.pr_band]['observation'] / c
            Xeci_sta1sat2, Yeci_sta1sat2, Zeci_sta1sat2 = CoorTransform.earth_rotation_correction([coorX_sta1sat2, coorY_sta1sat2, coorZ_sta1sat2], dt_sta1sat2)
            lou_sta1sat2 = CoorTransform.cal_distance(sta1_coor, [Xeci_sta1sat2, Yeci_sta1sat2, Zeci_sta1sat2])

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

            H11.append([a_sta2_X, a_sta2_Y, a_sta2_Z])
            H21.append([a_sta2_X, a_sta2_Y, a_sta2_Z])
            hx.append(lou_sta2sat2 - lou_sta1sat2 - lou_sta2sat1 + lou_sta1sat1)

        n = len(self.DD_records.DD_observations_data[self.cp_band])
        H12 = np.eye(n) * get_lamb_from_band(self.cp_band)
        H22 = np.zeros((n, n))
        H = np.block([[np.array(H11), H12], [np.array(H21), H22]])
        hx1 = np.array(2*hx)
        hx2 = get_lamb_from_band(self.cp_band) * np.array(x_array[3:].tolist() + [0 for i in range(len(x_array[3:].tolist()))])
        hx = hx1+hx2
        return H, hx

    def gety(self):
        y1 = []
        y2 = []
        for DD_data in self.DD_records.DD_observations_data[self.cp_band]:
            y1.append(get_lamb_from_band(self.cp_band) * DD_data.DD_obs)
        for DD_data in self.DD_records.DD_observations_data[self.pr_band]:
            y2.append(DD_data.DD_obs)
        y = np.array(y1+y2)
        return y

    def getx(self, sta2_coor):
        # 转换类型
        if isinstance(sta2_coor, np.ndarray):
            sta2_coor = sta2_coor.tolist()
        if isinstance(self.DD_records.DD_ambiguitys[self.cp_band], np.ndarray):
            DD_ambiguitys = self.DD_records.DD_ambiguitys[self.cp_band].tolist()
        else:
            DD_ambiguitys = self.DD_records.DD_ambiguitys[self.cp_band]
        x = np.array(sta2_coor + DD_ambiguitys)
        return x

    def getQ(self):
        Qs = self.DD_records.DD_ambiguitys_noise[self.cp_band]
        # 静态基线模式
        if self.mode == 0:
            Q = np.diag([1e-1, 1e-1, 1e-1] + Qs)
        # 动态基线模式
        elif self.mode == 1:
            Q = np.diag([10, 10, 10] + Qs)
        return Q

    def getR(self, sigma1=0.02, sigma2=3):
        nDD = len(self.DD_records.DD_observations_data[self.cp_band])
        covDD = np.full((nDD, nDD), 1).astype(float)
        for i in range(nDD):
            covDD[i, i] = 2
        covDD1 = 2 * sigma1**2 * covDD
        covDD2 = 2 * sigma2**2 * covDD
        R = RTK.diagonalize_squarematrix(covDD1, covDD2)
        return R

    def ekf_estimation(self, P_before, sta1_coor, sta2_coor, nav_records):
        # 预测
        F = self.getF()
        Q = self.getQ()
        x_pri = self.getx(sta2_coor)
        P_pri = F @ P_before @ F.T + Q
        # 更新
        H, hx = self.getH(nav_records, sta1_coor, sta2_coor, x_pri)
        R = self.getR()
        y = self.gety()
        K = P_pri @ H.T @ np.linalg.inv(H @ P_pri @ H.T + R)
        x_est = x_pri + K @ (y - hx)
        # 输出残差
        H, hx = self.getH(nav_records, sta1_coor, x_pri.tolist()[:3], x_pri)
        if self.residual_manager:
            self.residual_manager.add_epoch_residuals(self.DD_records.Tr, self.DD_records.base_svn, self.DD_records.qualitified_satellites, y-hx)
        if self.DDambiguity_manager:
            self.DDambiguity_manager.add_epoch_ambiguity(self.DD_records.Tr, self.DD_records.base_svn, self.DD_records.qualitified_satellites, x_est[3:])
        if self.position_manager:
            self.position_manager.add_epoch_position(self.DD_records.Tr, x_est[:3])
        P_est = (np.eye(len(x_est)) - K @ H) @ P_pri
        return x_est, P_est


# 附带坐标约束的卡尔曼滤波器
class Kalman_Filter_with_coordinate_constrained():
    def __init__(self, DD_records_atTr, cp_band, pr_band, mode_int, coordinate, residual_manager=None, DDambiguity_manager=None, position_manager=None):
        """
        DD_records_atTr : DD_records_atTr class, 双差观测集合
        cp_band : 载波相位
        pr_band : 伪距
        mode_int : 滤波解模式，0为静态，1为动态
        """
        self.DD_records = DD_records_atTr
        self.cp_band = cp_band
        self.pr_band = pr_band
        self.mode = mode_int
        self.coordinate = coordinate
        self.residual_manager = residual_manager
        self.DDambiguity_manager = DDambiguity_manager
        self.position_manager = position_manager

    def getF(self):
        n = len(self.DD_records.DD_observations_data[self.cp_band]) + 3
        Fmatrix = np.diag([1 for i in range(n)])
        return Fmatrix

    def getH(self, nav_data, sta1_coor, sta2_coor, x_array):
        H11 = []
        H21 = []
        hx = []
        for DD_data in self.DD_records.DD_observations_data[self.cp_band]:
            # 计算站1星1元素
            ts_sta1sat1, dts_sta1sat1 = SPP.cal_EmitTime_from_datetime(DD_data.T, DD_data.sat1, DD_data.obs1_sat1.data[self.pr_band]['observation'], nav_data, doCRC=True)
            coorX_sta1sat1, coorY_sta1sat1, coorZ_sta1sat1 = SatellitePosition.cal_SatellitePosition_GPS_GPSws(ts_sta1sat1, DD_data.sat1, nav_data)
            dt_sta1sat1 = DD_data.obs1_sat1.data[self.pr_band]['observation'] / c
            Xeci_sta1sat1, Yeci_sta1sat1, Zeci_sta1sat1 = CoorTransform.earth_rotation_correction([coorX_sta1sat1, coorY_sta1sat1, coorZ_sta1sat1], dt_sta1sat1)
            lou_sta1sat1 = CoorTransform.cal_distance(sta1_coor, [Xeci_sta1sat1, Yeci_sta1sat1, Zeci_sta1sat1])

            # 计算站1星2元素
            ts_sta1sat2, dts_sta1sat2 = SPP.cal_EmitTime_from_datetime(DD_data.T, DD_data.sat2, DD_data.obs1_sat2.data[self.pr_band]['observation'], nav_data, doCRC=True)
            coorX_sta1sat2, coorY_sta1sat2, coorZ_sta1sat2 = SatellitePosition.cal_SatellitePosition_GPS_GPSws(ts_sta1sat2, DD_data.sat2, nav_data)
            dt_sta1sat2 = DD_data.obs1_sat2.data[self.pr_band]['observation'] / c
            Xeci_sta1sat2, Yeci_sta1sat2, Zeci_sta1sat2 = CoorTransform.earth_rotation_correction([coorX_sta1sat2, coorY_sta1sat2, coorZ_sta1sat2], dt_sta1sat2)
            lou_sta1sat2 = CoorTransform.cal_distance(sta1_coor, [Xeci_sta1sat2, Yeci_sta1sat2, Zeci_sta1sat2])

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

            H11.append([a_sta2_X, a_sta2_Y, a_sta2_Z])
            H21.append([a_sta2_X, a_sta2_Y, a_sta2_Z])
            hx.append(lou_sta2sat2 - lou_sta1sat2 - lou_sta2sat1 + lou_sta1sat1)

        n = len(self.DD_records.DD_observations_data[self.cp_band])
        H12 = np.eye(n) * get_lamb_from_band(self.cp_band)
        H22 = np.zeros((n, n))

        hx1 = np.array(2*hx+sta2_coor)
        hx2 = get_lamb_from_band(self.cp_band) * np.array(x_array[3:].tolist() + [0 for i in range(len(x_array[3:].tolist()))] + [0, 0, 0])
        hx = hx1+hx2
        # 坐标值约束观测方程的添加
        H31 = np.eye(3)
        H32 = np.zeros((3, n))
        # 合并系数阵和观测值常数阵
        H = np.block([[np.array(H11), H12], [np.array(H21), H22], [H31, H32]])
        hx = hx1+hx2
        return H, hx

    def gety(self):
        y1 = []
        y2 = []
        for DD_data in self.DD_records.DD_observations_data[self.cp_band]:
            y1.append(get_lamb_from_band(self.cp_band) * DD_data.DD_obs)
        for DD_data in self.DD_records.DD_observations_data[self.pr_band]:
            y2.append(DD_data.DD_obs)
        y = np.array(y1+y2+self.coordinate)
        return y

    def getx(self, sta2_coor):
        # 转换类型
        if isinstance(sta2_coor, np.ndarray):
            sta2_coor = sta2_coor.tolist()
        if isinstance(self.DD_records.DD_ambiguitys[self.cp_band], np.ndarray):
            DD_ambiguitys = self.DD_records.DD_ambiguitys[self.cp_band].tolist()
        else:
            DD_ambiguitys = self.DD_records.DD_ambiguitys[self.cp_band]
        x = np.array(sta2_coor + DD_ambiguitys)
        return x

    def getQ(self):
        Qs = self.DD_records.DD_ambiguitys_noise[self.cp_band]
        # 静态基线模式
        if self.mode == 0:
            Q = np.diag([1e-1, 1e-1, 1e-1] + Qs)
        # 动态基线模式
        elif self.mode == 1:
            Q = np.diag([1000, 1000, 1000] + Qs)
        return Q

    def getR(self, sigma1=0.02, sigma2=3, sigma3=[0.1, 0.1, 0.1]):
        nDD = len(self.DD_records.DD_observations_data[self.cp_band])
        covDD = np.full((nDD, nDD), 1).astype(float)
        for i in range(nDD):
            covDD[i, i] = 2
        covDD1 = 2 * sigma1**2 * covDD
        covDD2 = 2 * sigma2**2 * covDD
        R = RTK.diagonalize_squarematrix(covDD1, covDD2)
        R_coor = np.diag([sigma3[0]**2, sigma3[1]**2, sigma3[2]**2])
        R = RTK.diagonalize_squarematrix(R, R_coor)
        return R

    def ekf_estimation(self, P_before, sta1_coor, sta2_coor, nav_records):
        # 预测
        F = self.getF()
        Q = self.getQ()
        x_pri = self.getx(sta2_coor)
        P_pri = F @ P_before @ F.T + Q
        # 更新
        H, hx = self.getH(nav_records, sta1_coor, sta2_coor, x_pri)
        R = self.getR()
        y = self.gety()
        K = P_pri @ H.T @ np.linalg.inv(H @ P_pri @ H.T + R)
        x_est = x_pri + K @ (y - hx)
        # 输出残差
        H, hx = self.getH(nav_records, sta1_coor, sta2_coor, x_pri)
        if self.residual_manager:
            self.residual_manager.add_epoch_residuals(self.DD_records.Tr, self.DD_records.base_svn, self.DD_records.qualitified_satellites, y-hx)
        if self.DDambiguity_manager:
            self.DDambiguity_manager.add_epoch_ambiguity(self.DD_records.Tr, self.DD_records.base_svn, self.DD_records.qualitified_satellites, x_est[3:])
        if self.DDambiguity_manager:
            self.position_manager.add_epoch_position(self.DD_records.Tr, x_est[:3])
        P_est = (np.eye(len(x_est)) - K @ H) @ P_pri
        return x_est, P_est



# 附带基线长约束的卡尔曼滤波器
class Kalman_Filter_with_baseline_constrained():
    def __init__(self, DD_records_atTr, cp_band, pr_band, baseline_length, l_sigma, mode_int, residual_manager=None, DDambiguity_manager=None, position_manager=None):
        """
        DD_records_atTr : DD_records_atTr class, 双差观测集合
        band1 : 载波相位
        band2 : 伪距
        baseline_length : 基线长度
        len_sigma : 基线长度测量标准差
        mode_int : 定位模式，0为静态，1为动态
        """
        self.DD_records = DD_records_atTr
        self.cp_band = cp_band
        self.pr_band = pr_band
        self.baseline_length = baseline_length
        self.l_sigma = l_sigma
        self.mode = mode_int
        self.residual_manager = residual_manager
        self.DDambiguity_manager = DDambiguity_manager
        self.position_manager = position_manager
        # print(len(DD_records_atTr.DD_observations_data[band1]), len(DD_records_atTr.DD_ambiguitys_noise[band1]))

    def getF(self):
        n = len(self.DD_records.DD_observations_data[self.cp_band]) + 3
        Fmatrix = np.diag([1 for i in range(n)])
        return Fmatrix

    def getH(self, nav_data, sta1_coor, sta2_coor, x_array):
        H11 = []
        H21 = []
        hx = []
        for DD_data in self.DD_records.DD_observations_data[self.cp_band]:
            # 计算站1星1元素
            ts_sta1sat1, dts_sta1sat1 = SPP.cal_EmitTime_from_datetime(DD_data.T, DD_data.sat1, DD_data.obs1_sat1.data[self.pr_band]['observation'], nav_data, doCRC=True)
            coorX_sta1sat1, coorY_sta1sat1, coorZ_sta1sat1 = SatellitePosition.cal_SatellitePosition_GPS_GPSws(ts_sta1sat1, DD_data.sat1, nav_data)
            dt_sta1sat1 = DD_data.obs1_sat1.data[self.pr_band]['observation'] / c
            Xeci_sta1sat1, Yeci_sta1sat1, Zeci_sta1sat1 = CoorTransform.earth_rotation_correction([coorX_sta1sat1, coorY_sta1sat1, coorZ_sta1sat1], dt_sta1sat1)
            lou_sta1sat1 = CoorTransform.cal_distance(sta1_coor, [Xeci_sta1sat1, Yeci_sta1sat1, Zeci_sta1sat1])

            # 计算站1星2元素
            ts_sta1sat2, dts_sta1sat2 = SPP.cal_EmitTime_from_datetime(DD_data.T, DD_data.sat2, DD_data.obs1_sat2.data[self.pr_band]['observation'], nav_data, doCRC=True)
            coorX_sta1sat2, coorY_sta1sat2, coorZ_sta1sat2 = SatellitePosition.cal_SatellitePosition_GPS_GPSws(ts_sta1sat2, DD_data.sat2, nav_data)
            dt_sta1sat2 = DD_data.obs1_sat2.data[self.pr_band]['observation'] / c
            Xeci_sta1sat2, Yeci_sta1sat2, Zeci_sta1sat2 = CoorTransform.earth_rotation_correction([coorX_sta1sat2, coorY_sta1sat2, coorZ_sta1sat2], dt_sta1sat2)
            lou_sta1sat2 = CoorTransform.cal_distance(sta1_coor, [Xeci_sta1sat2, Yeci_sta1sat2, Zeci_sta1sat2])

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
            H11.append([a_sta2_X, a_sta2_Y, a_sta2_Z])
            H21.append([a_sta2_X, a_sta2_Y, a_sta2_Z])
            hx.append(lou_sta2sat2 - lou_sta1sat2 - lou_sta2sat1 + lou_sta1sat1)

        n = len(self.DD_records.DD_observations_data[self.cp_band])
        H12 = np.eye(n) * get_lamb_from_band(self.cp_band)
        H22 = np.zeros((n, n))
        # 基线长约束观测方程的添加
        l12 = CoorTransform.cal_distance(sta1_coor, sta2_coor)
        H31 = np.array([(sta2_coor[0]-sta1_coor[0])/l12, (sta2_coor[1]-sta1_coor[1])/l12, (sta2_coor[2]-sta1_coor[2])/l12])
        H32 = np.zeros((1, n))
        # 合成整个设计矩阵
        H = np.block([[np.array(H11), H12], [np.array(H21), H22], [H31, H32]])
        hx1 = np.array(2*hx+[l12])
        hx2 = get_lamb_from_band(self.cp_band) * np.array(x_array[3:].tolist() + [0 for i in range(len(x_array[3:].tolist()))] + [0])
        hx = hx1+hx2
        return H, hx

    def gety(self, L):
        y1 = []
        y2 = []
        for DD_data in self.DD_records.DD_observations_data[self.cp_band]:
            y1.append(get_lamb_from_band(self.cp_band) * DD_data.DD_obs)
        for DD_data in self.DD_records.DD_observations_data[self.pr_band]:
            y2.append(DD_data.DD_obs)
        y = np.array(y1+y2+[L])
        return y

    def getx(self, sta2_coor):
        # 转换类型
        if isinstance(sta2_coor, np.ndarray):
            sta2_coor = sta2_coor.tolist()
        if isinstance(self.DD_records.DD_ambiguitys[self.cp_band], np.ndarray):
            DD_ambiguitys = self.DD_records.DD_ambiguitys[self.cp_band].tolist()
        else:
            DD_ambiguitys = self.DD_records.DD_ambiguitys[self.cp_band]
        x = np.array(sta2_coor + DD_ambiguitys)
        return x

    def getQ(self):
        Qs = self.DD_records.DD_ambiguitys_noise[self.cp_band]
        # 静态基线模式
        if self.mode == 0:
            Q = np.diag([0.1, 0.1, 0.1] + Qs)
        # 动态基线模式
        elif self.mode == 1:
            Q = np.diag([1000, 1000, 1000] + Qs)
        return Q

    def getR(self, sigma3, sigma1=0.002, sigma2=3):
        nDD = len(self.DD_records.DD_observations_data[self.cp_band])
        covDD = np.full((nDD, nDD), 1).astype(float)
        for i in range(nDD):
            covDD[i, i] = 2
        covDD1 = 2 * sigma1**2 * covDD
        covDD2 = 2 * sigma2**2 * covDD
        R = RTK.diagonalize_squarematrix(covDD1, covDD2)
        R_l = np.array([[sigma3**2]])
        R = RTK.diagonalize_squarematrix(R, R_l)
        return R

    def ekf_estimation(self, P_before, sta1_coor, sta2_coor, nav_records):
        # 预测
        F = self.getF()
        Q = self.getQ()
        x_pri = self.getx(sta2_coor)
        P_pri = F @ P_before @ F.T + Q
        # 更新
        H, hx = self.getH(nav_records, sta1_coor, sta2_coor, x_pri)
        R = self.getR(self.l_sigma)
        y = self.gety(self.baseline_length)
        K = P_pri @ H.T @ np.linalg.inv(H @ P_pri @ H.T + R)
        x_est = x_pri + K @ (y - hx)
        # 输出残差
        H, hx = self.getH(nav_records, sta1_coor, sta2_coor, x_pri)
        if self.residual_manager:
            self.residual_manager.add_epoch_residuals(self.DD_records.Tr, self.DD_records.base_svn, self.DD_records.qualitified_satellites, y-hx)
        if self.DDambiguity_manager:
            self.DDambiguity_manager.add_epoch_ambiguity(self.DD_records.Tr, self.DD_records.base_svn, self.DD_records.qualitified_satellites, x_est[3:])
        if self.DDambiguity_manager:
            self.position_manager.add_epoch_position(self.DD_records.Tr, x_est[:3])
        P_est = (np.eye(len(x_est)) - K @ H) @ P_pri
        return x_est, P_est


def get_common_svn(knownStation_ob_records_atTr, unknownStation_ob_records_atTr):
    """
    knownStation_ob_records_atTr : 参考站观测记录
    unknownStation_ob_records_atTr : 移动站观测记录
    """
    svn_Tr = []
    knownsta_svn = [x.SVN for x in knownStation_ob_records_atTr]
    unknownsta_svn = [x.SVN for x in unknownStation_ob_records_atTr]
    for svn in knownsta_svn:
        if svn in unknownsta_svn:
            svn_Tr.append(svn)
    return svn_Tr


def constaneously_RTK_withfilter(start_time, end_time, knownStation_ob_records, unknownStation_ob_records, br_records,
                                 knownStation_coor, unknownStation_init_coor, cp_band, pr_band,
                                 interval_time, cycle_slip_detect=False, mode=0, baseline_length=0, l_sigma=0, ele_limit=13):
    """
    start_time : datetime.datetime , 时间段开始时刻
    end_time : datetime.datetime , 时间段结束时刻
    knownStation_ob_records : 基准站多时刻观测数据
    unknownStation_ob_records : 未知站多时刻观测数据
    br_records : 卫星星历数据
    knownStation_coor : 基准站坐标
    unknownStation_init_coor : 未知站坐标
    cp_band : 载波波段
    pr_band : 伪距波段
    interval_time : 间隔时间, s
    mode : 静态基线为0， 动态为1
    baseline_length : 基线长度, m
    l_sigma : 基线长度标准差, m
    """
    Tr = start_time
    # coordinates = []
    fixed_coordinates = []
    float_coordinates = []
    # 残差记录
    DD_residual_manager = PictureResults.plot_DDresidual_bysvn()  # 双差观测值残差
    DD_ambiguity_manager = PictureResults.plot_DDambiguity_bysvn()   # 双差模糊度
    position_manager = PictureResults.plot_position()     # 位置解算结果
    ele_manager = PictureResults.plot_elevation()    # 高度角

    # 第一个历元进行最小二乘解算，给初值
    sta2_coor, P, N_float, base_svn, diff_svns, Pse, residual = RTK2.DD_onCarrierPhase_and_Pseudorange_1known(knownStation_ob_records, unknownStation_ob_records, br_records, Tr, knownStation_coor, unknownStation_init_coor, bands=[cp_band, pr_band], ambi_fix=False)
    DD_residual_manager.add_epoch_residuals(Tr, base_svn, diff_svns, residual)
    DD_ambiguity_manager.add_epoch_ambiguity(Tr, base_svn, diff_svns, N_float)
    position_manager.add_epoch_position(Tr, sta2_coor)
    # sta2_coor = unknownStation_init_coor
    P = np.diag([100**2 for i in range(3)]+[300**2 for i in range(len(N_float))])
    # coordinates.append(sta2_coor)
    fixed_coordinates.append(sta2_coor)
    float_coordinates.append(sta2_coor)
    coor = sta2_coor
    # 进行传递
    if isinstance(sta2_coor, np.ndarray):
        sta2_coor = sta2_coor.tolist()
    if isinstance(N_float, np.ndarray):
        N_float = N_float.tolist()
    amb_trans = ambiguity_transfer(cp_band, N_float, base_svn, diff_svns, get_diaglist_from_matrix(P)[3:])
    amb_trans.update_ele(sta2_coor, Tr, br_records)

    # 进行第二个历元的结算
    DD_GF_collection = ''
    Tr += datetime.timedelta(seconds=interval_time)
    n = 0
    while Tr < end_time:
        # 找到Tr时刻两个观测站观测记录
        knownStation_ob_records_atTr = list(filter(lambda o:o.time == Tr and o.SVN[0] == "G" and o.managed_data_flag[cp_band] and o.managed_data_flag[pr_band] and o.managed_data_flag['L2_L'], knownStation_ob_records))
        unknownStation_ob_records_atTr = list(filter(lambda o:o.time == Tr and o.SVN[0] == "G" and o.managed_data_flag[cp_band] and o.managed_data_flag[pr_band] and o.managed_data_flag['L2_L'], unknownStation_ob_records))
        # 获取符合条件的所有svn
        svn_Tr = get_common_svn(knownStation_ob_records_atTr, unknownStation_ob_records_atTr)

        # 计算各个卫星的高度角
        sats_ele = {}
        svn_ele_qualitified = []
        w, s = TimeSystem.from_datetime_cal_GPSws(Tr)
        t = TimeSystem.GPSws(w, s)
        for svn in svn_Tr:
            sat_coor = SatellitePosition.cal_SatellitePosition_GPS_GPSws(t, svn, br_records)
            sats_ele[svn] = CoorTransform.cal_ele_and_A(coor, sat_coor)[0] * 180 / math.pi
            if sats_ele[svn] > ele_limit:
                svn_ele_qualitified.append(svn)
        ele_manager.add_epoch_elevations(Tr, list(sats_ele.values()), list(sats_ele.keys()))
        # 提取符合高度角要求的观测记录
        knownStation_ob_records_atTr = list(filter(lambda o: o.SVN in svn_ele_qualitified, knownStation_ob_records_atTr))
        unknownStation_ob_records_atTr = list(filter(lambda o: o.SVN in svn_ele_qualitified, unknownStation_ob_records_atTr))

        # 构造双差观测值
        DD_records_collection = DD_records_atTr(Tr, [cp_band, pr_band])
        trans_matrix = DD_records_collection.set_base_satellite_from_ambiguity_transfer(amb_trans, knownStation_ob_records_atTr, unknownStation_ob_records_atTr, 10)
        # trans_matrix = DD_records_collection.set_elebest_base_satellite_from_ambiguity_transfer(amb_trans, knownStation_ob_records_atTr, unknownStation_ob_records_atTr)
        if cycle_slip_detect:
            # 进行周跳探测
            svn_removed_index, new_in, new_in_svn = DD_records_collection.add_satellites(svn_ele_qualitified, knownStation_ob_records_atTr, unknownStation_ob_records_atTr, DD_GF_collection)
        else:
            # 不带周跳探测
            svn_removed_index, new_in, new_in_svn = DD_records_collection.add_satellites(svn_ele_qualitified, knownStation_ob_records_atTr, unknownStation_ob_records_atTr)
        print("-2-:", svn_removed_index, new_in, new_in_svn)

        # 判断是否发生了参考星的变化
        if np.all(trans_matrix == 0):    # 参考星不变
            # 调整较上个历元被去除的卫星对应的vc阵
            if svn_removed_index:  # 上个历元的卫星有缺失
                svn_removed_index2 = [i+3 for i in svn_removed_index]
                P = remove_cross_from_matrix(P, svn_removed_index2, svn_removed_index2)
        else:      # 参考星变化
            # 调整较上个历元被去除的卫星对应的vc阵
            if svn_removed_index:  # 上个历元的卫星要进行删减的
                svn_removed_index2 = [i + 3 for i in svn_removed_index]
                P = remove_cross_from_matrix(P, svn_removed_index2, svn_removed_index2)
                # trans_matrix = np.delete(trans_matrix, svn_removed_index, axis=1)
                trans_matrix = remove_cross_from_matrix(trans_matrix, svn_removed_index, svn_removed_index)
            P_trans_matrix = RTK.diagonalize_squarematrix(np.eye(3), trans_matrix)
            P = P_trans_matrix @ P @ P_trans_matrix.T
        # 更新由卫星数量增加导致的模糊度参数vc阵的变化
        if new_in != 0:  # 此历元有新增卫星
            eyes = 10000 * np.eye(new_in)    # 前面的常数为新增卫星的方差
            P = RTK.diagonalize_squarematrix(P, eyes)
        # 开始估计
        if baseline_length and l_sigma:    # 带基线长度约束的模式
            KM_estimator = Kalman_Filter_with_baseline_constrained(DD_records_collection, cp_band, pr_band, baseline_length, l_sigma, mode_int=mode)
            x, P = KM_estimator.ekf_estimation(P, sta1_coor=knownStation_coor, sta2_coor=coor, nav_records=br_records)
        else:     # 不带基线长度约束的模式
            KM_estimator = Kalman_Filter(DD_records_collection, cp_band, pr_band, mode_int=mode, residual_manager=DD_residual_manager, DDambiguity_manager=DD_ambiguity_manager, position_manager=position_manager)
            # KM_estimator = Kalman_Filter_with_coordinate_constrained(DD_records_collection, cp_band, pr_band, 0, unknownStation_init_coor, residual_manager=DD_residual_manager, DDambiguity_manager=DD_ambiguity_manager, position_manager=position_manager)
            x, P = KM_estimator.ekf_estimation(P, sta1_coor=knownStation_coor, sta2_coor=coor, nav_records=br_records)
        print(Tr, x)
        DD_records_collection.DD_ambiguitys[cp_band] = x[3:].tolist()
        ambiguitys = DD_records_collection.DD_ambiguitys[cp_band]
        base_svn = DD_records_collection.base_svn
        satellites = DD_records_collection.qualitified_satellites
        amb_nos = [0.000001 for i in range(len(P) - 3)]
        amb_trans = ambiguity_transfer(cp_band, ambiguitys, base_svn, satellites, amb_nos)
        DD_GF_collection = DD_records_collection.DD_GF

        # 固定模糊度
        # if amb_fix:
        #     coor = ambiguity_resolution_solution(P, x)
        # else:
        #     coor = x[:3].tolist()

        fixed_coor = ambiguity_resolution_solution(P, x)
        float_coor = x[:3].tolist()
        fixed_coordinates.append(fixed_coor)
        float_coordinates.append(float_coor)

        Tr += datetime.timedelta(seconds=interval_time)
        amb_trans.update_ele(coor, Tr, br_records)
        n += 1

    # 绘制残差图
    DD_residual_manager.plot_residuals(form="scatter")
    DD_ambiguity_manager.plot_DDambiguitys(form="scatter")
    # position_manager.plot_position(form="scatter")
    ele_manager.plot_elevations(form="scatter")

    return float_coordinates, fixed_coordinates, P, amb_trans





if __name__ == '__main__':

    # station1_observation_file = r"edata\obs\leij3100.20o"    # 已知站点 leij
    # station2_observation_file = r"edata\obs\zim23100.20o"    # 未知站点 zim2
    station2_observation_file = r"edata\obs\zimm3100.20o"    # 未知站点 zimm
    # station1_observation_file = r"edata\obs\wab23100.20o"    # 已知站点 wab2
    station1_observation_file = r"edata\obs\zim23100.20o"  # 已知站点 zim2
    # station1_observation_file = r"edata\obs\zimm3100.20o"  # 已知站点 zimm
    broadcast_file = r"edata\sat_obit\brdc3100.20n"

    #
    # station1_observation_file = r"edata\attitude_project\rinex2\cuaa0030.21o"  # 已知站点 cuaa
    # # station1_observation_file = r"edata\attitude_project\rinex2\cucc0030.21o"  # 已知站点 cucc
    # station2_observation_file = r"edata\attitude_project\rinex2\cubb0030.21o"  # 未知站点 cubb
    # # station2_observation_file = r"edata\attitude_project\rinex2\cut00030.21o"  # 未知站点 cut0
    # broadcast_file = r"edata\attitude_project\rinex2\brdc0030.21n"

    # station1_observation_file = r"D:\Desktop\XY503_XY602\A503_40.21o"  # 已知站点 503
    # station2_observation_file = r"D:\Desktop\XY503_XY602\A602_40.21o"  # 未知站点 602
    # broadcast_file = r"D:\Desktop\XY503_XY602\A602_40.21n"
    # knownStation_ob_records = DoFile.read_Rinex3_oFile(station1_observation_file)
    # unknownStation_ob_records = DoFile.read_Rinex3_oFile(station2_observation_file)
    # br_records = DoFile.read_Renix304_nFile(broadcast_file)


    # 读入观测文件内容,得到类型对象列表
    knownStation_ob_records = DoFile.read_Rinex2_oFile(station1_observation_file)
    unknownStation_ob_records = DoFile.read_Rinex2_oFile(station2_observation_file)
    br_records = DoFile.read_GPS_nFile(broadcast_file)
    print("数据读取完毕！")



    # 坐标
    knownStation_coor = [4331300.1600, 567537.0810, 4633133.5100]  # zim2
    # knownStation_coor = [4327318.2325, 566955.9585, 4636425.9246]  # wab2
    init_coor = [4331297.3480, 567555.6390, 4633133.7280]  # zimm

    # knownStation_coor = [-2364336.1554, 4870280.8223, -3360815.9725]  # cuaa
    # # knownStation_coor = [-2364332.0167, 4870284.0337, -3360814.1380]  # cucc
    # init_coor = [-2364334.0912, 4870285.4649, -3360807.7523]  # cubb

    # knownStation_coor = [-2850420.9626, 4651846.2388, 3292931.3812]  # 503
    # init_coor = [-2850405.2169, 4651847.9444, 3292949.6125]  # 602

    start_time = datetime.datetime(2020, 11, 5, 0, 0, 0)
    end_time = datetime.datetime(2020, 11, 5, 23, 59, 0)
    # start_time = datetime.datetime(2021, 1, 3, 16, 0, 0)
    # end_time = datetime.datetime(2021, 1, 3, 20, 59, 0)
    # start_time = datetime.datetime(2021, 2, 9, 0, 0, 0)
    # end_time = datetime.datetime(2021, 2, 9, 23, 59, 0)

    cal_float_coors, cal_fixed_coors, P, ambi_transfer = constaneously_RTK_withfilter(start_time, end_time, knownStation_ob_records, unknownStation_ob_records, br_records, knownStation_coor, init_coor, 'L1_L', 'L2_C', 30, ele_limit=9, cycle_slip_detect=False)
    true_coors = [init_coor for i in range(len(cal_fixed_coors))]

    time_series = PictureResults.get_time_series(start_time, end_time, 30)


    SPP.cal_NEUerrors(true_coors, cal_fixed_coors, T_series=time_series)
    SPP.cal_NEUerrors(true_coors, cal_float_coors, T_series=time_series)
    print("固定解 neu各方向RMSE:", ResultAnalyse.get_NEU_rmse(true_coors, cal_fixed_coors))
    print("浮点解 neu各方向RMSE:", ResultAnalyse.get_NEU_rmse(true_coors, cal_float_coors))



    # RTK2
    # start_time = datetime.datetime(2021, 1, 3, 6, 20, 30)
    # end_time = datetime.datetime(2021, 1, 3, 10, 59, 0)
    # Tr = start_time
    # cal_coors = []
    # residual_manager = PictureResults.plot_DDresidual_bysvn("双差观测值残差序列图")
    # while Tr<end_time:
    #     print(Tr)
    #     CoorXYZ, Q, N_float, the_SVN, final_SVNs, Pse, residual = RTK2.DD_onCarrierPhase_and_Pseudorange_1known(knownStation_ob_records,
    #                                                                                              unknownStation_ob_records,
    #                                                                                              br_records, Tr,
    #                                                                                              knownStation_coor,
    #                                                                                              init_coor, ambi_fix=True)
    #     Xk, Yk, Zk = CoorXYZ
    #     cal_coors.append([Xk, Yk, Zk])
    #     residual_manager.add_epoch_residuals(Tr, the_SVN, final_SVNs, residual)
    #     Tr+=datetime.timedelta(seconds=30)
    # true_coors = [init_coor for i in range(len(cal_coors))]
    # SPP.cal_NEUerrors(true_coors, cal_coors)
    # residual_manager.plot_residuals()


    # SPP.cal_NEUerrors(true_coors, cal_coors)
    # SPP.cal_XYZerrors(true_coors, cal_coors)
    # print("neu各方向RMSE:", ResultAnalyse.get_NEU_rmse(true_coors, cal_coors))
    # print("坐标RMSE:", ResultAnalyse.get_coor_rmse(true_coors, cal_coors))
































