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
import numpy as np


carrier_phase_list = []
pseudorange_list = []
lamb_dictionary = {}

# 双差观测值类
class DD_record():
    def __init__(self, Tr, sat1, sat2, band):
        self.sat1 = sat1
        self.sat2 = sat2
        self.band = band
        self.T = Tr

    def obs_is_null(self, band, obs1_sat1, obs1_sat2, obs2_sat1, obs2_sat2):
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
        if self.obs_is_null(self.band):
            # 构造双差观测值
            obs_sta1sat1 = float(obs1_sat1.data[self.band]['observation'])
            obs_sta1sat2 = float(obs1_sat2.data[self.band]['observation'])
            obs_sta2sat1 = float(obs2_sat1.data[self.band]['observation'])
            obs_sta2sat2 = float(obs2_sat2.data[self.band]['observation'])
            self.DD_obs = obs_sta2sat2 - obs_sta1sat2 - obs_sta2sat1 + obs_sta1sat1
        else:
            self.DD_obs = 0


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
        for band in self.bands:
            self.DD_observations_data[band] = []
            if band in carrier_phase_list:
                self.DD_ambiguitys[band] = []

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
                DD_record = DD_record(self.Tr, self.base_svn, the_svn, band)
                DD_record.get_DD_obs(self.base_sta1obs_record, other_sta1obs_record, self.base_sta2obs_record, other_sta2obs_record)
                if DD_record.DD_obs == 0:
                    qualified_flag = False
                    break
                DD_obses.append(DD_record)
            # 根据卫星数据进行是否加入待差卫星的判断
            if qualified_flag:
                self.qualitified_satellites.append(the_svn)
                for band in self.bands:
                    self.DD_observations_data[band].append(DD_obses[self.bands.index(band)])
                    if band in carrier_phase_list:
                        lamb = lamb_dictionary[band]  # 补一个波长
                        self.DD_ambiguitys.append(DD_obses[self.bands.index(band)].DD_obs/lamb)
            else:
                self.unqualitified_satellites.append(the_svn)
                print('该卫星数据不符合要求！')
        else:
            print('请设置参考星！')

    def add_satellites(self, svns, sta1obs_records_Tr, sta2obs_records_Tr):
        if self.qualitified_satellites:
            # 把不在观测数据中的卫星删除
            for svn in self.qualitified_satellites:
                if svn not in svns:
                    svn_index = self.qualitified_satellites.index(svn)
                    self.qualitified_satellites.remove(svn)
                    self.DD_ambiguitys.pop(svn_index)
                # 对卫星数据进行检查,如有不符合的也从数据中移除
                qualified_flag = True
                DD_obses = []
                for band in self.bands:
                    DD_record = DD_record(self.Tr, self.base_svn, svn, band)
                    DD_record.get_DD_obs(self.base_sta1obs_record, sta1obs_records_Tr, self.base_sta2obs_record,
                                         sta2obs_records_Tr)
                    if DD_record.DD_obs == 0:
                        qualified_flag = False
                        break
                    DD_obses.append(DD_record)
                # 根据卫星数据进行是否加入待差卫星的判断
                if qualified_flag:
                    for band in self.bands:
                        self.DD_observations_data[band].append(DD_obses[self.bands.index(band)])
                else:
                    svn_index = self.qualitified_satellites.index(svn)
                    self.qualitified_satellites.remove(svn)
                    self.DD_ambiguitys.pop(svn_index)
                    print('该卫星数据不符合要求！')
            # 将新卫星加入
            for svn in svns:
                if svn not in self.qualitified_satellites:
                    self.add_satellite(svn, sta1obs_records_Tr, sta2obs_records_Tr)
        else:
            for svn in svns:
                self.add_satellite(svn, sta1obs_records_Tr, sta2obs_records_Tr)





    def set_base_satellite_from_ambiguity_transfer(self, ambiguity_transfer, sta1obs_records_Tr, sta2obs_records_Tr):
        self.base_svn = ambiguity_transfer.new_base_svn
        self.base_sta1obs_record = list(filter(lambda o: o.SVN == self.base_svn, sta1obs_records_Tr))[0]
        self.base_sta2obs_record = list(filter(lambda o: o.SVN == self.base_svn, sta2obs_records_Tr))[0]
        # todo 此处需要再处理参考星数据完整性的问题
        self.qualitified_satellites = ambiguity_transfer.satellites
        self.DD_ambiguitys = ambiguity_transfer.ambiguitys


    # def clone(self):
    #     clone_DDrecords = DD_records_atTr(self.Tr, self.bands, self.sta1, self.sta2)
    #     clone_DDrecords.DD_observations = self.DD_observations_data
    #     clone_DDrecords.qualitified_satellites = self.qualitified_satellites
    #     clone_DDrecords.unqualitified_satellites = self.unqualitified_satellites
    #     return clone_DDrecords

class ambiguity_transfer():
    def __init__(self, band, ambiguitys, satellites):
        self.band = band
        self.ambiguitys = ambiguitys
        self.satellites = satellites

    def transfer_to_new_basesats(self, new_base_svn):
        self.new_base_svn = new_base_svn
        new_base_svn_index = self.satellites.index(new_base_svn)
        transfer_matrix_pre = np.eye(len(self.satellites)-1)
        transfer_matrix = np.insert(transfer_matrix_pre, new_base_svn_index, np.full((len(self.satellites)-1), -1), axis=1)
        self.satellites.remove(new_base_svn)
        self.ambiguitys = (transfer_matrix @ np.array(self.ambiguitys)).tolist()




if __name__ == '__main__':
    # 测试模糊度转移矩阵
    band = 'L1'
    ambiguitys = [2.4, 2.5, 4.5, 5.2]
    sats = ['G02', 'G03', 'G21', 'G12']
    mmm = ambiguity_transfer(band, ambiguitys, sats)
    mmm.transfer_to_new_basesats('G21')






























