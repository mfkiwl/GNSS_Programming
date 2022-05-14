# -*- coding: utf-8 -*-
"""

@title:	Style Guide for Python Code
@author: iDeal0103
@status:	Active
@type:	Process
@created:	05-Jan-2022
@post-History:	05-Jan-2022

comment：
    1.画误差序列图
    2.高效画图工具，主要针对残差绘图

"""


# import
import datetime

import utils.DoFile as DoFile
import utils.CoorTransform as CoorTransform
import matplotlib.pyplot as plt
import numpy as np
import random


def get_time_series(start_time, end_time, time_interval):
    time_series = []
    time = start_time
    while time < end_time:
        time_series.append(time)
        time += datetime.timedelta(seconds=time_interval)
    return time_series



def paint_error_sequence_diagram(pefile,igsfile,prnlist,savepath):
    """
    Parameters
    ----------
    pefile : str , 计算所得卫星精密星历数据(带头文件),如 r"E:\大三下\卫星与导航定位\3100_20n.sp3"
    igsfile : str , igs提供的卫星精密星历数据(带头文件),如 r"E:\大三下\卫星与导航定位\igs21304.sp3"
    prnlist : str list, 参与计算的卫星PRN编号列表
    savepath : str , 误差序列图集存储的文件夹路径,如 r"E:\大三下\卫星与导航定位\误差序列图"
    """
    cpe_records=DoFile.read_GPS_sp3File(pefile)
    ipe_records=DoFile.read_GPS_sp3File(igsfile)
    for PRN in prnlist:
        #筛选对象
        cpe_data=list(filter(lambda o:o.PRNs==PRN,cpe_records))
        cpe_data.sort(key=lambda p:p.time)
        ipe_data=list(filter(lambda o:o.PRNs==PRN,ipe_records))
        ipe_data.sort(key=lambda p:p.time)
        #计算数据
        dXs=[]
        dYs=[]
        dZs=[]
        ddTs=[]
        times=[]
        for i in range(len(cpe_data)):
            times.append(cpe_data[i].time)
            dX=cpe_data[i].X-ipe_data[i].X
            dY=cpe_data[i].Y-ipe_data[i].Y
            dZ=cpe_data[i].Z-ipe_data[i].Z
            ddT=cpe_data[i].dT-ipe_data[i].dT
            dXs.append(dX)
            dYs.append(dY)
            dZs.append(dZ)
            ddTs.append(ddT)
        #绘图
        plt.rcParams['font.sans-serif']=['SimHei']
        plt.rcParams['axes.unicode_minus'] = False
        plt.scatter(times,dXs,color="r",label = "delta X / km")
        plt.scatter(times,dYs,color="g",label = "delta Y / km")
        plt.scatter(times,dZs,color="b",label = "delta Z / km")
        plt.scatter(times,ddTs,color="y",label = "delta dT / um")
        plt.legend(loc='upper right')
        plt.title(str(PRN)+" 误差序列图")
        plt.savefig(savepath+"\\"+str(PRN)+"误差序列图")
        plt.show()


def paint_3dimensional_scatterplot_of_positionerrors_inXYZ(cal_coors: list, true_coors, kenamatic_mode: bool = False,
                                                           xaixs_limit: tuple = (-1, 1), yaixs_limit: tuple = (-1, 1),
                                                           zaixs_limit: tuple = (-1, 1)):
    paint_flg = False
    if kenamatic_mode:     # 对象坐标为动态
        if len(cal_coors) != len(true_coors):
            print('坐标计算结果与真实值个数不等！')
            raise SystemExit
        else:
            x_errors = [(cal_coors[i][0] - true_coors[i][0]) for i in range(len(cal_coors))]
            y_errors = [(cal_coors[i][1] - true_coors[i][1]) for i in range(len(cal_coors))]
            z_errors = [(cal_coors[i][2] - true_coors[i][2]) for i in range(len(cal_coors))]
            paint_flg = True
    else:   # 对象坐标为静态
        if len(true_coors) != 1:
            print('静态坐标数量超过1！')
        else:
            x_errors = [(cal_coors[i][0] - true_coors[0][0]) for i in range(len(cal_coors))]
            y_errors = [(cal_coors[i][1] - true_coors[0][1]) for i in range(len(cal_coors))]
            z_errors = [(cal_coors[i][2] - true_coors[0][2]) for i in range(len(cal_coors))]
            paint_flg = True
    # 开始画图
    if paint_flg:
        ax = plt.figure().add_subplot(projection='3d')
        # 绘制误差三维空间分布
        ax.scatter(x_errors, y_errors, z_errors, c="b", label='position errors (unit: m)')
        # 绘制在各平面上误差投影分布
        # ax.scatter(x_errors, y_errors, zs=zaixs_limit[0], zdir='z', c="c", label='position errors in x-y plane (unit: m)')
        # ax.scatter(x_errors, z_errors, zs=yaixs_limit[0], zdir='y', c="m", label='position errors in x-z plane (unit: m)')
        # ax.scatter(y_errors, z_errors, zs=xaixs_limit[0], zdir='x', c="y", label='position errors in y-z plane (unit: m)')
        ax.scatter(x_errors, y_errors, zs=zaixs_limit[0], zdir='z', c="c", marker='+')
        ax.scatter(x_errors, z_errors, zs=yaixs_limit[0], zdir='y', c="m", marker='+')
        ax.scatter(y_errors, z_errors, zs=xaixs_limit[0], zdir='x', c="y", marker='+')
        # 显示参数设置
        ax.legend(loc='upper right')
        if xaixs_limit:
            ax.set_xlim(xaixs_limit[0], xaixs_limit[1])
        if yaixs_limit:
            ax.set_ylim(yaixs_limit[0], yaixs_limit[1])
        if zaixs_limit:
            ax.set_zlim(zaixs_limit[0], zaixs_limit[1])
        ax.set_xlabel('dX (m)')
        ax.set_ylabel('dY (m)')
        ax.set_zlabel('dZ (m)')
    plt.show()


def paint_3dimensional_scatterplot_of_positionerrors_inNEU(cal_coors: list, true_coors, kenamatic_mode: bool = False,
                                                           xaixs_limit: tuple = (-1, 1), yaixs_limit: tuple = (-1, 1),
                                                           zaixs_limit: tuple = (-1, 1)):
    paint_flg = False
    N_errors = []
    E_errors = []
    U_errors = []
    if kenamatic_mode:     # 对象坐标为动态
        if len(cal_coors) != len(true_coors):
            print('坐标计算结果与真实值个数不等！')
            raise SystemExit
        else:
            for i in range(len(cal_coors)):
                n, e, u = CoorTransform.cal_NEU(true_coors[i], cal_coors[i])
                N_errors.append(n)
                E_errors.append(e)
                U_errors.append(u)
            paint_flg = True
    else:   # 对象坐标为静态
        if len(true_coors) != 1:
            print('静态坐标数量超过1！')
        else:
            for i in range(len(cal_coors)):
                n, e, u = CoorTransform.cal_NEU(true_coors[0], cal_coors[i])
                N_errors.append(n)
                E_errors.append(e)
                U_errors.append(u)
            paint_flg = True
    # 开始画图
    if paint_flg:
        ax = plt.figure().add_subplot(projection='3d')
        # 绘制误差三维空间分布
        ax.scatter(N_errors, E_errors, U_errors, c="b", label='position errors in NEU (unit: m)')
        # 绘制在各平面上误差投影分布
        # ax.scatter(N_errors, E_errors, zs=zaixs_limit[0], zdir='z', c="c", marker=',', label='position errors in N-E plane (unit: m)')
        # ax.scatter(N_errors, U_errors, zs=yaixs_limit[0], zdir='y', c="m", marker=',', label='position errors in N-U plane (unit: m)')
        # ax.scatter(E_errors, U_errors, zs=xaixs_limit[0], zdir='x', c="y", marker=',', label='position errors in E-U plane (unit: m)')
        ax.scatter(N_errors, E_errors, zs=zaixs_limit[0], zdir='z', c="c", marker='+')
        ax.scatter(N_errors, U_errors, zs=yaixs_limit[0], zdir='y', c="m", marker='+')
        ax.scatter(E_errors, U_errors, zs=xaixs_limit[0], zdir='x', c="y", marker='+')
        # 显示参数设置
        ax.legend(loc='upper right')
        if xaixs_limit:
            ax.set_xlim(xaixs_limit[0], xaixs_limit[1])
        if yaixs_limit:
            ax.set_ylim(yaixs_limit[0], yaixs_limit[1])
        if zaixs_limit:
            ax.set_zlim(zaixs_limit[0], zaixs_limit[1])
        ax.set_xlabel('dN (m)')
        ax.set_ylabel('dE (m)')
        ax.set_zlabel('dU (m)')
    plt.show()


# 随机生成颜色
def get_cmap(n, name='tab20b'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)

# 单条绘图数值记录类
class plot_record():
    def __init__(self, T, value, label, instruction=""):
        self.T = T
        self.value = value
        self.label = label
        self.instruction = instruction

# 绘图数值记录类的集合管理&画图类
class plot_records_manager():
    def __init__(self, title):
        self.plot_records = []
        self.title = title

    def add_plot_record(self, plot_record):
        self.plot_records.append(plot_record)

    def plot_by_labels(self, labels, form="plot"):
        # 绘图设置
        plt.rcParams['font.sans-serif'] = ['SimHei']
        plt.rcParams['axes.unicode_minus'] = False
        colors = get_cmap(len(labels))
        # 整理记录
        for label in labels:
           # record_list = list(filter(lambda o:o.label == label, self.plot_records)).sort(key=lambda o:o.T)
           record_list = list(filter(lambda o: o.label == label, self.plot_records))
           x = [i.T for i in record_list]
           y = [i.value for i in record_list]
           if form == "plot":
               plt.plot(x, y, color=colors(labels.index(label)), label=label)
           elif form == "scatter":
               plt.scatter(x, y, color=colors(labels.index(label)), label=label)
        plt.title(self.title)
        # plt.legend(loc='best', ncol=5)
        plt.legend(loc='best')
        plt.show()


# 绘制伪距+相位双差残差类
class plot_DDresidual_bysvn():
    def __init__(self, title="双差观测值残差序列图"):
        self.plot_records_manager=plot_records_manager(title)
        self.svn_pairs = []
        self.labels = []

    def add_epoch_residuals(self, T, base_svn, svns, residuals):
        n = len(svns)
        cp_residuals = residuals[:n]
        pr_residuals = residuals[n:]
        for i in range(n):
            svn = svns[i]
            self.add_svnpair_withoutrepeat(base_svn + "-" + svn)
            # 增加载波相位残差对象
            cp_residual = cp_residuals[i]
            cp_record = plot_record(T, cp_residual, base_svn+"-"+svn+"_cp")
            self.plot_records_manager.add_plot_record(cp_record)
            # 增加伪距残差对象
            pr_residual = pr_residuals[i]
            pr_record = plot_record(T, pr_residual, base_svn+"-"+svn+"_pr")
            self.plot_records_manager.add_plot_record(pr_record)

    def add_svnpair_withoutrepeat(self, svn):
        if svn not in self.svn_pairs:
            self.svn_pairs.append(svn)
        else:
            return

    def plot_residuals(self, form="plot"):
        labels_cp = [svn_pair+"_cp" for svn_pair in self.svn_pairs]
        self.plot_records_manager.plot_by_labels(labels_cp, form)
        labels_pr = [svn_pair+"_pr" for svn_pair in self.svn_pairs]
        self.plot_records_manager.plot_by_labels(labels_pr, form)
        # self.plot_records_manager.plot_by_labels(self.svn_pairs)


# 绘制双差观测值（伪距/相位）残差类
class plot_DDobs_residual_bysvn():
    def __init__(self, title="双差观测值残差序列图"):
        self.plot_records_manager=plot_records_manager(title)
        self.svn_pairs = []
        self.labels = []

    def add_epoch_residuals(self, T, base_svn, svns, DDobs_residuals):
        n = len(svns)
        for i in range(n):
            svn = svns[i]
            self.add_svnpair_withoutrepeat(base_svn + "-" + svn)
            # 增加伪距残差对象
            DDobs_residual = DDobs_residuals[i]
            DDobs_record = plot_record(T, DDobs_residual, base_svn+"-"+svn)
            self.plot_records_manager.add_plot_record(DDobs_record)

    def add_svnpair_withoutrepeat(self, svnpais):
        if svnpais not in self.svn_pairs:
            self.svn_pairs.append(svnpais)
        else:
            return

    def plot_residuals(self, form="plot"):
        labels = [svn_pair for svn_pair in self.svn_pairs]
        self.plot_records_manager.plot_by_labels(labels, form)


# 绘制双差模糊度类
class plot_DDambiguity_bysvn():
    def __init__(self, round=False, title="双差模糊度序列"):
        self.plot_records_manager=plot_records_manager(title)
        self.labels = []
        self.round = round

    def add_epoch_ambiguity(self, T, base_svn, dif_svns, ambiguitys):
        n = len(dif_svns)
        for i in range(n):
            svn = dif_svns[i]
            label = base_svn+"-"+svn
            # 增加双差模糊度对象
            DDam_record = plot_record(T, ambiguitys[i], label)
            if label not in self.labels:
                self.labels.append(label)
            self.plot_records_manager.add_plot_record(DDam_record)

    def plot_DDambiguitys(self, form="plot"):
        self.plot_records_manager.plot_by_labels(self.labels, form)


# 绘制坐标序列类
class plot_position():
    def __init__(self, title="坐标序列图"):
        self.plot_records_manager=plot_records_manager(title)
        self.labels = ["X", "Y", "Z"]

    def add_epoch_position(self, T, coor):
        for i in range(3):
            # 增加双差模糊度对象
            DD_position_record = plot_record(T, coor[i], self.labels[i])
            self.plot_records_manager.add_plot_record(DD_position_record)

    def plot_position(self, form="plot"):
        self.plot_records_manager.plot_by_labels(self.labels, form)
        # self.plot_records_manager.plot_by_labels(['X'])
        # self.plot_records_manager.plot_by_labels(['Y'])
        # self.plot_records_manager.plot_by_labels(['Z'])


class plot_elevation():
    def __init__(self, title="高度角序列图"):
        self.plot_records_manager=plot_records_manager(title)
        self.labels = []

    def add_epoch_elevations(self, T, ele_values, svns_label):
        for i in range(len(ele_values)):
            # 增加高度角对象
            ele_record = plot_record(T, ele_values[i], svns_label[i])
            if svns_label[i] not in self.labels:
                self.labels.append(svns_label[i])
            self.plot_records_manager.add_plot_record(ele_record)

    def plot_elevations(self, form="plot"):
        self.plot_records_manager.plot_by_labels(self.labels, form)



# 绘制观测值噪声类
class plot_obsnoise_bystationsvn():
    def __init__(self, round=False, title="观测值噪声序列"):
        self.plot_records_manager=plot_records_manager(title)
        self.labels = []

    def add_epoch_obsnoise(self, T, station, svn, noise):
        label = station+"-"+svn
        # 增加双差模糊度对象
        prnoise_record = plot_record(T, noise, label)
        if label not in self.labels:
            self.labels.append(label)
        self.plot_records_manager.add_plot_record(prnoise_record)

    def plot_obsnoise(self):
        self.plot_records_manager.plot_by_labels(self.labels)




if __name__ == "__main__":
    # 生成坐标
    # cal_coors = []
    # true_coors = []
    # true_coor = [[0.2, 0.2, 0.2]]
    # for i in range(100):
    #     cal_coors.append([random.random(), random.random(), random.random()])
    #     true_coors.append([random.random(), random.random(), random.random()])
    # paint_3dimensional_scatterplot_of_positionerrors_inXYZ(cal_coors, true_coor, kenamatic_mode=False)
    plt.rcParams['font.sans-serif'] = ['SimHei']
    plt.rcParams['axes.unicode_minus'] = False



    x1 = [1, 2, 3, 4]
    x2 = [2, 4, 5, 6]
    y = [2, 5, 8, 9]
    colors = get_cmap(7, "autumn")
    plt.plot(x1, y, color=colors(5), label="x1")
    plt.plot(x2, y, color=colors(10), label="x2")

    # plt.plot(x1, y, 'g', "x1", x2, y, 'r', "x2")
    # plt.plot([1, 2, 3], [1, 2, 3], 'go-', label='line 1')
    # plt.plot([1, 1.5, 3], [1, 4, 9], 'rs', label='line 2')
    plt.legend(loc='upper right')
    plt.title("experiment")
    plt.xlabel("T")
    plt.show()







