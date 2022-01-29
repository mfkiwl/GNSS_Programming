# -*- coding: utf-8 -*-
"""

@title:	Style Guide for Python Code
@author: iDeal0103
@status:	Active
@type:	Process
@created:	21-Jan-2022
@post-History:	21-Jan-2022

comment：
    1.Widelane (WL)
    2.

"""

from utils.const import c


def observation_isvalid(obs_record, FreqBandList=['L1', 'L2']):
    '''
    obsrecord_record :  observation_record class , 所使用观测文件读取的纪录类
    FreqBandList : list[str] , 检查数据是否为空的频率波段 (如 'C1','L2','P1' )
    '''
    isnot_null_flags = []
    for band in FreqBandList:
        # 判断波段的数据是否齐全
        if not obs_record.data.__contains__(band):
            isnot_null_flag = False
        elif obs_record.data[band]['observation'] == '':
            isnot_null_flag = False
        else:
            isnot_null_flag = True
        isnot_null_flags.append(isnot_null_flag)
    if False in isnot_null_flags:
        result = False
    else:
        result = True
    return result

def get_widelane_combination(obs_record, band1, band2, f1, f2):
    """


    """
    isvalid_flg = observation_isvalid(obs_record, FreqBandList=[band1, band2])
    if isvalid_flg:
        measurement1 = obs_record.data[band1]['observation']
        measurement2 = obs_record.data[band2]['observation']
        widelane_combination = (f1 * measurement1 - f2 * measurement2) / (f1 - f2)
        lamb_widelane = c / ((f1 -f2) * 1000000)
        return widelane_combination, lamb_widelane


