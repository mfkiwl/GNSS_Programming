

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




if __name__ == "__main__":
    # station1_observation_file = r"edata\attitude_project\rinex2\cuaa0030.21o"  # 已知站点 cuaa
    # station2_observation_file = r"edata\attitude_project\rinex2\cubb0030.21o"  # 未知站点 cubb
    station1_observation_file = r"edata\obs\zim23100.20o"  # 已知站点 zim2
    station2_observation_file = r"edata\obs\zimm3100.20o"  # 已知站点 zimm

    # broadcast_file = r"edata\attitude_project\rinex2\brdc0030.21n"
    broadcast_file = r"edata\sat_obit\brdc3100.20n"

    # 读入观测文件内容,得到类型对象列表
    knownStation_ob_records = DoFile.read_Rinex2_oFile(station1_observation_file)
    unknownStation_ob_records = DoFile.read_Rinex2_oFile(station2_observation_file)
    br_records = DoFile.read_GPS_nFile(broadcast_file)
    print("数据读取完毕！")

    # station1_coor = [-2364336.1554, 4870280.8223, -3360815.9725]  # cuaa
    # # station2_coor = [-2364334.2691, 4870286.6995, -3360808.7896]  # cubb
    # station2_coor = [-2364334.0912, 4870285.4649, -3360807.7523]  # cubb流动站rtk计算结果
    # # station2_coor = [-2364339.0912, 4870289.4649, -3360809.7523]  # cubb错误

    station1_coor = [4331300.1600, 567537.0810, 4633133.5100]  # zim2
    station2_coor = [4331297.3480, 567555.6390, 4633133.7280]  # zimm

    # RTK2
    # start_time = datetime.datetime(2021, 1, 3, 2, 0, 0)
    # end_time = datetime.datetime(2021, 1, 3, 8, 59, 30)
    start_time = datetime.datetime(2020, 11, 5, 8, 40, 0)
    end_time = datetime.datetime(2020, 11, 5, 10, 0, 0)
    Tr = start_time
    cal_coors = []
    DDambiguity_manager = PictureResults.plot_DDambiguity_bysvn()
    DDobs_residual_manager = PictureResults.plot_DDobs_residual_bysvn()
    prnpoise_manager = PictureResults.plot_obsnoise_bystationsvn()
    while Tr < end_time:
        N_float, residuals, the_SVN, final_SVNs = RTK2.DD_onCPandPR_solve_ambiguity(
            knownStation_ob_records, unknownStation_ob_records, br_records, Tr, station1_coor, station2_coor,
            bands=['L1_L', 'L2_C'],
            DDambiguity_manager=DDambiguity_manager, DDobs_residual_manager=DDobs_residual_manager, prnoise_manager=prnpoise_manager)
        print(Tr, N_float, residuals)
        Tr += datetime.timedelta(seconds=30)
    DDambiguity_manager.plot_DDambiguitys(form="scatter")
    DDobs_residual_manager.plot_residuals(form="scatter")
    # prnpoise_manager.plot_obsnoise()



