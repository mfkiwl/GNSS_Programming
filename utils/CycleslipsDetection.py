# -*- coding: utf-8 -*-
"""

@title:	Cycle slips detection
@author: iDeal0103
@status:	Active
@type:	Process
@created:	10-March-2022
@post-History:	10-March-2022

comment:
    1. 周跳探测

"""



def GF_detector(GFbefore, GFnow, threshold=0.05):
    if abs(GFnow - GFbefore) > threshold:
        cycle_slip_flag = True
    else:
        cycle_slip_flag = False
    return cycle_slip_flag





