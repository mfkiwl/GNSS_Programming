# -*- coding: utf-8 -*-
"""

@title:	pos_plot
@author: iDeal0103
@status: Active
@type:	Process
@created:	30-May-2022
@post-History:	30-May-2022

comment:
    1. plot pos
    2.

"""

# import
import math

import numpy as np

import utils.DoFile as DoFile
import SinglePointPosition as SPP
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import utils.ResultAnalyse as ResultAnalyse
import utils.CoorTransform as CoorTransform
import attitude_determination.TRIAD as TRAID

# 绘制ecef坐标
def plot_pos_xyz(pos_file, timemode=0, posvalmode=1, basecoor=None, xyz2neu=False):
    # 读取文件
    pos_records = DoFile.read_posfile_xyz(pos_file, timemode=timemode, posvalmode=posvalmode, basecoor=basecoor, xyz2neu=xyz2neu)
    # 构造画图数据
    xs = []
    ys = []
    zs = []
    times = []
    for record in pos_records:
        xs.append(record.pos_value[0])
        ys.append(record.pos_value[1])
        zs.append(record.pos_value[2])
        # times.append(record.T.strftime("%H:%M:%S"))
        times.append(record.T)
    # # 绘制方法1调用不同subplot绘图
    # plt.rcParams['font.sans-serif'] = ['SimHei']
    # plt.rcParams['axes.unicode_minus'] = False
    # # 绘制x
    # plt.subplot(311)
    # plt.plot(times, xs, color="r", label="X / m")
    # plt.legend(loc='upper right')
    # plt.grid()
    # # 绘制y
    # plt.subplot(312)
    # plt.plot(times, ys, color="g", label="Y / m")
    # plt.legend(loc='upper right')
    # plt.grid()
    # # 绘制z
    # plt.subplot(313)
    # plt.plot(times, zs, color="b", label="Z / m")
    # plt.legend(loc='upper right')
    # plt.grid()
    # # 显示
    # plt.figure(1)
    # # plt.ylabel('position / m')
    # plt.xlabel('GPST')
    # plt.show()

    # 绘制方法2
    plt.rcParams['font.sans-serif'] = ['SimHei']
    plt.rcParams['axes.unicode_minus'] = False
    fig = plt.figure(tight_layout=True)
    gs = gridspec.GridSpec(3, 1)
    # 绘制x
    ax = fig.add_subplot(gs[0, :])
    ax.plot(times, xs, color="r", label="X/m")
    ax.legend(loc='upper left')
    # ax.set_ylabel('position-X / m')
    # ax.set_xlabel('GPST')
    ax.grid()
    # 绘制y
    ax = fig.add_subplot(gs[1, :])
    ax.plot(times, ys, color="g", label="Y/m")
    ax.legend(loc='upper left')
    # ax.set_ylabel('position-Y / m')
    # ax.set_xlabel('GPST')
    ax.grid()
    # 绘制z
    ax = fig.add_subplot(gs[2, :])
    ax.plot(times, zs, color="b", label="Z/m")
    ax.legend(loc='upper left')
    # ax.set_ylabel('position-Z / m')
    # ax.set_xlabel('GPST')
    ax.grid()
    # 显示结果
    plt.show()


# 绘制neu坐标（也有可能是向量）
def plot_pos_neu(pos_file, timemode=0, posvalmode=1, basecoor=[], xyz2neu=False, fontsize=16):
    # 读取文件
    pos_records = DoFile.read_posfile_xyz(pos_file, timemode=timemode, posvalmode=posvalmode, basecoor=basecoor, xyz2neu=xyz2neu)
    # 构造画图数据
    es = []
    ns = []
    us = []
    times = []
    for record in pos_records:
        es.append(record.pos_value[0])
        ns.append(record.pos_value[1])
        us.append(record.pos_value[2])
        # times.append(record.T.strftime("%H:%M:%S"))
        times.append(record.T)
    # # 绘制方法1调用不同subplot绘图
    # plt.rcParams['font.sans-serif'] = ['SimHei']
    # plt.rcParams['axes.unicode_minus'] = False
    # # 绘制x
    # plt.subplot(311)
    # plt.plot(times, xs, color="r", label="X / m")
    # plt.legend(loc='upper right')
    # plt.grid()
    # # 绘制y
    # plt.subplot(312)
    # plt.plot(times, ys, color="g", label="Y / m")
    # plt.legend(loc='upper right')
    # plt.grid()
    # # 绘制z
    # plt.subplot(313)
    # plt.plot(times, zs, color="b", label="Z / m")
    # plt.legend(loc='upper right')
    # plt.grid()
    # # 显示
    # plt.figure(1)
    # # plt.ylabel('position / m')
    # plt.xlabel('GPST')
    # plt.show()

    # 绘制方法2
    plt.rcParams['font.sans-serif'] = ['SimHei']
    plt.rcParams['axes.unicode_minus'] = False
    fig = plt.figure(tight_layout=False)
    gs = gridspec.GridSpec(3, 1)
    # 绘制n
    ax = fig.add_subplot(gs[0, :])
    ax.plot(times, ns, color="r", label="N/m")
    ax.legend(loc='upper left',fontsize=fontsize)
    # ax.set_ylabel('position-Y / m')
    # ax.set_xlabel('GPST')
    # ax.get_xaxis().set_visible(False)    # 是否显示横坐标
    plt.yticks(fontsize=fontsize)
    plt.xticks(fontsize=1)
    ax.grid()
    # 绘制e
    ax = fig.add_subplot(gs[1, :])
    ax.plot(times, es, color="g", label="E/m")
    ax.legend(loc='upper left', fontsize=fontsize)
    # ax.set_ylabel('position-X / m')
    # ax.set_xlabel('GPST')
    # ax.get_xaxis().set_visible(False)    # 是否显示横坐标
    plt.yticks(fontsize=fontsize)
    plt.xticks(fontsize=1)
    ax.grid()
    # 绘制u
    ax = fig.add_subplot(gs[2, :])
    ax.plot(times, us, color="b", label="U/m")
    ax.legend(loc='upper left', fontsize=fontsize)
    # ax.set_ylabel('position-Z / m')
    # ax.set_xlabel('GPST')
    plt.yticks(fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    ax.grid()
    # 显示结果
    plt.subplots_adjust(hspace=0)
    plt.show()



def plot_pos_moving_baseline(pos_file, timemode=0, l=0 ,fontsize=16):
    # 读取文件
    pos_records = DoFile.read_posfile_xyz(pos_file, timemode=timemode)
    # 构造画图数据
    es = []
    ns = []
    us = []
    ls = []
    times = []
    for record in pos_records:
        es.append(record.pos_value[0])
        ns.append(record.pos_value[1])
        us.append(record.pos_value[2])
        ls.append(math.sqrt(record.pos_value[0]**2 + record.pos_value[1]**2 +record.pos_value[2]**2))
        # times.append(record.T.strftime("%H:%M:%S"))
        times.append(record.T)
    # 绘制方法2
    plt.rcParams['font.sans-serif'] = ['SimHei']
    plt.rcParams['axes.unicode_minus'] = False
    fig = plt.figure(tight_layout=False)     # 设置是否自动调整子图之间横坐标不互相遮盖
    gs = gridspec.GridSpec(4, 1)
    # 绘制x
    ax = fig.add_subplot(gs[0, :])
    ax.plot(times, ns, color="r", label="N/m")
    ax.legend(loc='upper left', fontsize=fontsize)
    plt.xticks(fontsize=1)
    # ax.get_xaxis().set_visible(False)    # 是否显示横坐标
    plt.yticks(fontsize=fontsize)
    # ax.set_ylabel('position-X / m')
    # ax.set_xlabel('GPST')
    ax.grid()
    # 绘制y
    ax = fig.add_subplot(gs[1, :])
    ax.plot(times, es, color="g", label="E/m")
    ax.legend(loc='upper left', fontsize=fontsize)
    plt.xticks(fontsize=1)
    plt.yticks(fontsize=fontsize)
    # ax.get_xaxis().set_visible(False)
    # ax.set_ylabel('position-Y / m')
    # ax.set_xlabel('GPST')
    ax.grid()
    # 绘制z
    ax = fig.add_subplot(gs[2, :])
    ax.plot(times, us, color="b", label="U/m")
    ax.legend(loc='upper left', fontsize=fontsize)
    plt.xticks(fontsize=1)
    plt.yticks(fontsize=fontsize)
    # ax.get_xaxis().set_visible(False)
    # ax.set_ylabel('position-Z / m')
    # ax.set_xlabel('GPST')
    ax.grid(axis='both')
    # 绘制l
    ax = fig.add_subplot(gs[3, :])
    ax.plot(times, ls, color="orange", label="l/m")
    ax.legend(loc='upper left', fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    ax.grid()
    # 显示结果
    plt.subplots_adjust(hspace=0.0)   # 设置子图之间的上下距离
    plt.show()
    print("e-std(m):", ResultAnalyse.get_standard_deviation(es))
    print("n-std(m):", ResultAnalyse.get_standard_deviation(ns))
    print("u-std(m):", ResultAnalyse.get_standard_deviation(us))
    print("l-rmse(m):", ResultAnalyse.get_rmse([l for i in range(len(ls))], ls))


# 绘制定位误差
def plot_pos_xyz_error(pos_file, true_coor, form="scatter", timemode=0, posvalmode=1, basecoor=None):
    # 读取文件
    pos_records = DoFile.read_posfile_xyz(pos_file, timemode=timemode, posvalmode=posvalmode, basecoor=basecoor)
    # 构造画图数据
    cal_coors = []
    true_coors = []
    times = []
    for record in pos_records:
        cal_coors.append(record.pos_value)
        times.append(record.T)
        true_coors.append(true_coor)
    SPP.cal_XYZerrors(cal_coors, true_coors, form=form, T_series="")
    print("neu各方向RMSE:", ResultAnalyse.get_XYZ_rmse(true_coors, cal_coors))
    print("坐标RMSE:", ResultAnalyse.get_coor_rmse(true_coors, cal_coors))

# 绘制定位误差
def plot_pos_neu_error(pos_file, true_coor, form="scatter", timemode=0, posvalmode=1, basecoor=None):
    # 读取文件
    pos_records = DoFile.read_posfile_xyz(pos_file, timemode=timemode, posvalmode=posvalmode, basecoor=basecoor)
    # 构造画图数据
    cal_coors = []
    true_coors = []
    times = []
    for record in pos_records:
        cal_coors.append(record.pos_value)
        times.append(record.T)
        true_coors.append(true_coor)
    SPP.cal_NEUerrors(cal_coors, true_coors, form=form, T_series="")
    print("neu各方向RMSE:", ResultAnalyse.get_NEU_rmse(true_coors, cal_coors))
    print("坐标RMSE:", ResultAnalyse.get_coor_rmse(true_coors, cal_coors))


def plot_courseangle_from_moving_baseline(posfile, timemode=0, fontsize=15):
    pos_records = DoFile.read_posfile_xyz(posfile, timemode=timemode)
    # 计算航向角
    course_angles = []
    times = []
    for record in pos_records:
        unit_baseline = TRAID.get_unit_vector(np.array(record.pos_value))[0]    # [E, N, U]
        course_angle = math.acos(unit_baseline.T @ np.array([0, 1, 0])) / math.pi * 180
        course_angles.append(course_angle)
        times.append(record.T)
    # 绘图
    plt.rcParams['font.sans-serif'] = ['SimHei']
    plt.rcParams['axes.unicode_minus'] = False
    plt.plot(times, course_angles, color="purple", label="course angle/degree")
    plt.legend(loc='upper left', fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.grid()
    plt.show()


def plot_2courseangle_from_moving_baseline(posfile1, posfile2, timemode=0, fontsize=15, label1='', label2=''):
    pos_records1 = DoFile.read_posfile_xyz(posfile1, timemode=timemode)
    pos_records2 = DoFile.read_posfile_xyz(posfile2, timemode=timemode)
    # 计算航向角
    course_angles1 = []
    times1 = []
    course_angles2 = []
    times2 = []
    for record in pos_records1:
        unit_baseline = TRAID.get_unit_vector(np.array(record.pos_value))[0]    # [E, N, U]
        # course_angle = math.acos(unit_baseline.T @ np.array([0, 1, 0])) / math.pi * 180
        course_angle = math.atan2(unit_baseline[0], unit_baseline[1]) / math.pi * 180
        course_angles1.append(course_angle)
        times1.append(record.T)
    for record in pos_records2:
        unit_baseline = TRAID.get_unit_vector(np.array(record.pos_value))[0]    # [E, N, U]
        # course_angle = math.acos(unit_baseline.T @ np.array([0, 1, 0])) / math.pi * 180
        course_angle = math.atan2(unit_baseline[0], unit_baseline[1]) / math.pi * 180
        course_angles2.append(course_angle)
        times2.append(record.T)
    # 绘图
    plt.rcParams['font.sans-serif'] = ['SimHei']
    plt.rcParams['axes.unicode_minus'] = False
    plt.scatter(times1, course_angles1, color="purple", label=label1 + " course angle/degree")
    plt.plot(times1, course_angles1, color="purple")
    plt.scatter(times2, course_angles2, color="lawngreen", label=label2 + " course angle/degree")
    plt.plot(times2, course_angles2, color="lawngreen")
    plt.legend(loc='upper left', ncol=2, fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.grid()
    plt.show()


def plot_2pos_moving_baseline(pos_file1, pos_file2, label1="result1", label2="result2", l=0 ,timemode=0, fontsize=16):
    # 读取文件
    pos_records1 = DoFile.read_posfile_xyz(pos_file1, timemode=timemode)
    pos_records2 = DoFile.read_posfile_xyz(pos_file2, timemode=timemode)
    # 构造画图数据
    es1 = []
    ns1 = []
    us1 = []
    ls1 = []
    es2 = []
    ns2 = []
    us2 = []
    ls2 = []
    times1 = []
    times2 = []
    for record in pos_records1:
        es1.append(record.pos_value[0])
        ns1.append(record.pos_value[1])
        us1.append(record.pos_value[2])
        ls1.append(math.sqrt(record.pos_value[0]**2 + record.pos_value[1]**2 +record.pos_value[2]**2))
        # times.append(record.T.strftime("%H:%M:%S"))
        times1.append(record.T)
    for record in pos_records2:
        es2.append(record.pos_value[0])
        ns2.append(record.pos_value[1])
        us2.append(record.pos_value[2])
        ls2.append(math.sqrt(record.pos_value[0]**2 + record.pos_value[1]**2 +record.pos_value[2]**2))
        times2.append(record.T)
    # 绘制方法2
    plt.rcParams['font.sans-serif'] = ['SimHei']
    plt.rcParams['axes.unicode_minus'] = False
    fig = plt.figure(tight_layout=False)
    gs = gridspec.GridSpec(4, 1)
    # 绘制x
    ax = fig.add_subplot(gs[0, :])
    ax.plot(times1, ns1, color="r", label=label1+" N/m")
    ax.plot(times2, ns2, color="lime", label=label2+" N/m")
    ax.legend(loc='upper left', ncol=2, fontsize=fontsize)
    # ax.set_ylabel('position-X / m')
    # ax.set_xlabel('GPST')
    plt.xticks(fontsize=1)
    # ax.get_xaxis().set_visible(False)
    plt.yticks(fontsize=fontsize)
    ax.grid()
    # 绘制y
    ax = fig.add_subplot(gs[1, :])
    ax.plot(times1, es1, color="g", label=label1+" E/m")
    ax.plot(times2, es2, color="tomato", label=label2+" E/m")
    ax.legend(loc='upper left', ncol=2, fontsize=fontsize)
    # ax.set_ylabel('position-Y / m')
    # ax.set_xlabel('GPST')
    plt.xticks(fontsize=1)
    # ax.get_xaxis().set_visible(False)
    plt.yticks(fontsize=fontsize)
    ax.grid()
    # 绘制z
    ax = fig.add_subplot(gs[2, :])
    ax.plot(times1, us1, color="b", label=label1+" U/m")
    ax.plot(times2, us2, color="gold", label=label2+" U/m")
    ax.legend(loc='upper left', ncol=2, fontsize=fontsize)
    # ax.set_ylabel('position-Z / m')
    # ax.set_xlabel('GPST')
    plt.xticks(fontsize=1)
    # ax.get_xaxis().set_visible(False)
    plt.yticks(fontsize=fontsize)
    ax.grid()
    # 绘制l
    ax = fig.add_subplot(gs[3, :])
    ax.plot(times1, ls1, color="orange", label=label1+" L/m")
    ax.plot(times2, ls2, color="cornflowerblue", label=label2+" L/m")
    ax.legend(loc='upper left', ncol=2, fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    ax.grid()
    # 显示结果
    plt.subplots_adjust(hspace=0.0)   # 设置子图之间的上下距离
    plt.show()
    print("e1-std(m):", ResultAnalyse.get_standard_deviation(es1))
    print("n1-std(m):", ResultAnalyse.get_standard_deviation(ns1))
    print("u1-std(m):", ResultAnalyse.get_standard_deviation(us1))
    print("e2-std(m):", ResultAnalyse.get_standard_deviation(es2))
    print("n2-std(m):", ResultAnalyse.get_standard_deviation(ns2))
    print("u2-std(m):", ResultAnalyse.get_standard_deviation(us2))
    print("l1-rmse(m):", ResultAnalyse.get_rmse([l for i in range(len(ls1))], ls1))
    print("l2-rmse(m):", ResultAnalyse.get_rmse([l for i in range(len(ls2))], ls2))


if __name__ == "__main__":
    # pos_file = r"D:\Desktop\attitude\20210103\cucc0030.pos"
    # pos_file = r"D:\Desktop\attitude\20210103\cucc0030_kinamic_elemask15_float.pos"
    # plot_pos_neu_error(pos_file, [-2364331.9908, 4870284.2312, -3360813.9384], form="plot")


    # pos_file2 = r"D:\Desktop\attitude\20210105\cut00050_movingbase_constraint.pos"
    # pos_file1 = r"D:\Desktop\attitude\20210105\cut00050_movingbase_baseline.pos"
    # pos_file1 = r"D:\Desktop\基线解算\0105\cutc0050_movingbase.pos"
    # pos_file2 = r"D:\Desktop\基线解算\0105\cutc0050_movingbase_constrained.pos"
    # plot_2pos_moving_baseline(pos_file1, pos_file2, "unconstrained", "constrained", timemode=0, l=0.735)
    # plot_pos_xyz_error(pos_file, [4331297.3480, 567555.6390, 4633133.7280], form="plot")

    # pos_file2 = r"D:\Desktop\基线解算\0105\cubb0050_movingbase_constrained.pos"
    # pos_file1 = r"D:\Desktop\基线解算\0105\cubb0050_movingbase.pos"
    # plot_2pos_moving_baseline(pos_file1, pos_file2, "unconstrained", "constrained", timemode=0, l=9.4693)

    pos_file = r"D:\Desktop\基线解算\0105\cubb0050_static.pos"
    plot_pos_neu(pos_file)
    # pos_file = r"D:\Desktop\attitude\20210103\cubb0030_static.pos"
    # plot_pos_neu(pos_file, basecoor=[-2364336.1554, 4870280.8223, -3360815.9725], xyz2neu=True)

    # pos_file = r"D:\Desktop\attitude\20210103\cubb0030_static_baseline.pos"

    # pos_file = r"D:\Desktop\毕业设计\数据\origin\AUX1_movingbase_constraint.pos"
    # pos_file = r"D:\Desktop\毕业设计\数据\origin\AUX1_movingbase_fixed.pos"
    # plot_pos_moving_baseline(pos_file, timemode=0, l=0.735)


    # pos_file = r"D:\Desktop\毕业设计\数据\origin\base-AUX1_fixed.pos"
    # plot_pos_neu(pos_file)

    # pos_file1 = r"D:\Desktop\毕业设计\数据\origin\AUX1_movingbase_constraint.pos"
    # pos_file2 = "D:\Desktop\毕业设计\数据\origin\AUX1_movingbase_fixed.pos"
    # # pos_file1 = r"D:\Desktop\毕业设计\数据\origin\AUX1_baseline_fullytime_constrained.pos"
    # # pos_file2 = r"D:\Desktop\毕业设计\数据\origin\AUX1_baseline_fullytime.pos"
    # plot_2courseangle_from_moving_baseline(pos_file1, pos_file2, fontsize=12, label1="constrained", label2='unconstrained')



















