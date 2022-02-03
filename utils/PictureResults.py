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
    2.

"""


# import
import utils.DoFile as DoFile
import utils.CoorTransform as CoorTransform
import matplotlib.pyplot as plt
import numpy as np
import random



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


if __name__ == "__main__":
    # 生成坐标
    cal_coors = []
    true_coors = []
    true_coor = [[0.2, 0.2, 0.2]]
    for i in range(100):
        cal_coors.append([random.random(), random.random(), random.random()])
        true_coors.append([random.random(), random.random(), random.random()])
    paint_3dimensional_scatterplot_of_positionerrors_inXYZ(cal_coors, true_coor, kenamatic_mode=False)






