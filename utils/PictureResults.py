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
import matplotlib.pyplot as plt


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
