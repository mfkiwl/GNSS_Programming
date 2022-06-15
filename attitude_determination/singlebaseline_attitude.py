

# import
import datetime

import utils.DoFile as DoFile
import utils.TimeSystem as TimeSystem
import attitude_determination.WOPP as WOPP
import attitude_determination.OPP as OPP
import  attitude_determination.TRIAD as TRAID
from scipy.spatial.transform.rotation import Rotation as rota
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import utils.ResultAnalyse as ResultAnalyse



def plot_attitude(pitch, yaw, roll, times, form="scatter", fontsize=15):
    # 绘图
    plt.rcParams['font.sans-serif'] = ['SimHei']
    plt.rcParams['axes.unicode_minus'] = False
    fig = plt.figure(tight_layout=False)
    gs = gridspec.GridSpec(3, 1)
    # 绘制pitch
    ax = fig.add_subplot(gs[0, :])
    if form=='scatter':
        ax.scatter(times, pitch, color="r", label="pitch/degree")
    elif form=='plot':
        ax.plot(times, pitch, color="r", label="pitch/degree")
    ax.legend(loc='upper left', fontsize=fontsize)
    plt.xticks(fontsize=1)
    plt.yticks(fontsize=fontsize)
    ax.grid()
    # 绘制yaw
    ax = fig.add_subplot(gs[1, :])
    if form == 'scatter':
        ax.scatter(times, yaw, color="g", label="yaw/degree")
    elif form == 'plot':
        ax.plot(times, yaw, color="g", label="yaw/degree")
    ax.legend(loc='upper left', fontsize=fontsize)
    plt.xticks(fontsize=1)
    plt.yticks(fontsize=fontsize)
    ax.grid()
    # 绘制roll
    ax = fig.add_subplot(gs[2, :])
    if form == 'scatter':
        ax.scatter(times, roll, color="b", label="roll/degree")
    elif form=='plot':
        ax.plot(times, roll, color="b", label="roll/degree")
    ax.legend(loc='upper left', fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    ax.grid()

    # 显示其他信息
    # plt.xlabel("epoch", fontdict={'size':12})
    # plt.ylabel('各方向角度', fontdict={'size':12})
    # plt.title("姿态解算结果", fontdict={'size':12})
    plt.subplots_adjust(hspace=0.0)
    plt.show()


if __name__ == "__main__":
    cuaa_coor = [-2364336.1554, 4870280.8223, -3360815.9725]
    # cuaa_coor = [-2364335.4045, 4870281.7841, -3360815.5730]

    # baseline1_records = DoFile.read_posfile_xyz(r"D:\Desktop\attitude\cubb0030.pos", 1, 1, cuaa_coor)
    # baseline2_records = DoFile.read_posfile_xyz(r"D:\Desktop\attitude\cucc0030.pos", 1, 1, cuaa_coor)
    # baseline3_records = DoFile.read_posfile_xyz(r"D:\Desktop\attitude\cut00030.pos", 1, 1, cuaa_coor)
    # baseline1_records = DoFile.read_posfile_xyz(r"D:\Desktop\attitude\20210103\cubb0030_static_combine_baseline.pos", 0, 1, cuaa_coor, xyz2neu=True)
    # baseline2_records = DoFile.read_posfile_xyz(r"D:\Desktop\attitude\20210103\cucc0030_static_combine_baseline.pos", 0, 1, cuaa_coor, xyz2neu=True)
    # baseline3_records = DoFile.read_posfile_xyz(r"D:\Desktop\attitude\20210103\cut00030_static_combine_baseline.pos", 0, 1, cuaa_coor, xyz2neu=True)
    # baseline1_records = DoFile.read_posfile_xyz(r"D:\Desktop\基线解算\0105\cubb0050_static.pos", 0, 0)
    # baseline2_records = DoFile.read_posfile_xyz(r"D:\Desktop\基线解算\0105\cutc0050_static.pos", 0, 0)
    # baseline3_records = DoFile.read_posfile_xyz(r"D:\Desktop\基线解算\0105\cut00050_static.pos", 0, 0)
    baseline1_records = DoFile.read_posfile_xyz(r"D:\Desktop\基线解算\0105\cubb0050_movingbase.pos", 0, 0)
    baseline2_records = DoFile.read_posfile_xyz(r"D:\Desktop\基线解算\0105\cutc0050_movingbase.pos", 0, 0)
    baseline3_records = DoFile.read_posfile_xyz(r"D:\Desktop\基线解算\0105\cut00050_movingbase.pos", 0, 0)



    f1 = np.array([1.8869, 5.8777, 7.1810])
    f2 = np.array([4.1646, 3.4089, 2.0341])
    # f3 = np.array([-2.0202, 4.1602, 7.0336])
    f3 = TRAID.get_unit_vector(np.cross(f1, f2))[0]

    pitch = []
    yaw = []
    roll = []

    epoch = 0
    # time = TimeSystem.GPSws(2139, 68760)
    time = datetime.datetime(2021, 1, 5, 6, 0, 0)
    times = []
    # while time.GpsSecond < 69480:
    while time < datetime.datetime(2021, 1, 5, 9, 0, 0):
        # timemode0
        baseline1_record = list(filter(lambda o: o.T > time-datetime.timedelta(seconds=0.1) and o.T < time+datetime.timedelta(seconds=0.1), baseline1_records))[0]
        baseline2_record = list(filter(lambda o: o.T > time-datetime.timedelta(seconds=0.1) and o.T < time+datetime.timedelta(seconds=0.1), baseline2_records))[0]
        baseline3_record = list(filter(lambda o: o.T > time-datetime.timedelta(seconds=0.1) and o.T < time+datetime.timedelta(seconds=0.1), baseline3_records))[0]
        # timemode1
        # baseline1_record = list(filter(lambda o: round(o.T.GpsSecond) == time.GpsSecond and o.T.GpsWeek == time.GpsWeek, baseline1_records))[0]
        # baseline2_record = list(filter(lambda o: round(o.T.GpsSecond) == time.GpsSecond and o.T.GpsWeek == time.GpsWeek, baseline2_records))[0]
        # baseline3_record = list(filter(lambda o: round(o.T.GpsSecond) == time.GpsSecond and o.T.GpsWeek == time.GpsWeek, baseline3_records))[0]
        # 构造设计基线矩阵F
        F = WOPP.get_matrix_from_vectors([f1.tolist(), f2.tolist(), f3.tolist()])
        # 构造基线矩阵B
        b3 = TRAID.get_unit_vector(np.cross(np.array(baseline1_record.pos_value), np.array(baseline2_record.pos_value)))[0]
        B = WOPP.get_matrix_from_vectors([baseline1_record.pos_value, baseline2_record.pos_value, b3.tolist()])

        # B = WOPP.get_matrix_from_vectors([baseline1_record.pos_value, baseline2_record.pos_value, baseline3_record.pos_value])
        # 构造vc-阵
        QBB=baseline1_record.Q
        QBB = WOPP.diagonalize_squarematrix(QBB, baseline2_record.Q)
        QBB = WOPP.diagonalize_squarematrix(QBB, baseline3_record.Q)
        # QBB = WOPP.diagonalize_squarematrix(QBB, np.eye(3, 3))
        # 计算
        R, Qrr, R_check = WOPP.solve_WOPP_withLagrangianMultipliers(B, F, QBB)
        # R = OPP.solve_OPP_withSVD(F, B, [1, 1, 1])
        # R = TRAID.solve_TRIAD(F[:, 1:], B[:, 1:])
        # print(time.GpsSecond, rota.from_matrix(R).as_euler('zyx', degrees=True))
        print(time, rota.from_matrix(R).as_euler('zyx', degrees=True))
        euler = rota.from_matrix(R).as_euler('zyx', degrees=True)
        pitch.append(euler[0])
        yaw.append(euler[1])
        roll.append(euler[2])

        epoch += 1
        # time.add(30)
        times.append(time)
        time += datetime.timedelta(seconds=30)

    plot_attitude(pitch, yaw, roll, times, form='scatter')
    print("pitch-STD(degree)：", ResultAnalyse.get_standard_deviation(pitch))
    print("yaw-STD(degree)：", ResultAnalyse.get_standard_deviation(yaw))
    print("roll-STD(degree)：", ResultAnalyse.get_standard_deviation(roll))









