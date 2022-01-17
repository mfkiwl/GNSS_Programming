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






