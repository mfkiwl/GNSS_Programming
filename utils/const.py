# -*- coding: utf-8 -*-
"""

@title:	Style Guide for Python Code
@author: iDeal0103
@status:	Active
@type:	Process
@created:	14-Apr-2021
@post-History:	14-Apr-2021

comment：
    1.记录地球常数

"""

# 地球重力常数,即GM
miu = 3.986005e14
# 地球自转角速度,单位为rad/s
we = 7.2921151467e-5
# 地球椭球长半轴
a = 6378137.0
# 地球椭球扁率
f = 1/298.257223563
# 光速,单位为m/s
c = 299792458

# GPS信号相关常数
# 频率, 单位MHz
L1_f = 1575.42
L2_f = 1227.60
L5_f = 1176.45
# 波长， 单位m
lamb_L1 = c / (L1_f * 1000000)
lamb_L2 = c / (L2_f * 1000000)
lamb_L5 = c / (L5_f * 1000000)
# 电离层相关
miu1 = 1
miu2 = (L1_f/L2_f)**2
miu5 = (L1_f/L5_f)**2

# 波长
GPS_carrier_phase_list = ['L1C', 'L1S', 'L1L', 'L1X', 'L1P', 'L1W', 'L1Y', 'L1M', 'L1N',
                      'L2C', 'L2D', 'L2S', 'L2L', 'L2X', 'L2P', 'L2W', 'L2Y', 'L2M', 'L2N',
                      'L5I', 'L5Q', 'L5X',
                      'L1', 'L2', 'L5']
GPS_pseudorange_list = ['C1C', 'C1S', 'C1L', 'C1X', 'C1P', 'C1W', 'C1Y', 'C1M',
                      'C2C', 'C2D', 'C2S', 'C2L', 'C2X', 'C2P', 'C2W', 'C2Y', 'C2M',
                      'C5I', 'C5Q', 'C5X',
                      'C1', 'C2', 'C5', 'P1', 'P2']

def get_lamb_from_band(band):
    if '1' in band:
        lamb = lamb_L1
    elif '2' in band:
        lamb = lamb_L2
    elif '5' in band:
        lamb = lamb_L5
    return lamb





