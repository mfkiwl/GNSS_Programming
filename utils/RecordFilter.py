# -*- coding: utf-8 -*-
"""

@title:	Style Guide for Python Code
@author: iDeal0103
@status:	Active
@type:	Process
@created:	17-Apr-2021
@post-History:	17-Apr-2021

comment：
    1.筛选时间最近的某个PRN号卫星记录
    2.

"""
#import module


#筛选时间最近的某个PRN号卫星记录
def cal_delta_time(time1, time2):
    if time1 > time2:
        delta_time = (time1-time2)
    else:
        delta_time = (time2-time1)
    return delta_time

def find_closest_record(records, time, SVN):  # SVN为卫星编号,如"G**"
    n = 0
    while True:
        if records[n].SVN == SVN:
        # if records[n].serial_no == '{:0>2d}'.format(serial_no):
            closest_record = records[n]
            break
        n += 1
    record_time = records[n].toc
    delta_time_min = cal_delta_time(time, record_time)
    for i in range(n+1, len(records)):
        if records[i].SVN == SVN:
        # if records[i].serial_no == '{:0>2d}'.format(serial_no):
            record_time = records[i].toc
            delta_time = cal_delta_time(time, record_time)
            if delta_time < delta_time_min:
                delta_time_min = delta_time
                closest_record = records[i]
    return closest_record


#将读取n文件获得的GPS_brdc_record列表中非整点数据除去
def GPSBrdcRecord_HourIntegerRecord_Filter(GPS_brdc_records):
    filtered_GPS_brdc_records=[]
    for record in GPS_brdc_records:
        if record.toe%60==0:
            filtered_GPS_brdc_records.append(record)
    return filtered_GPS_brdc_records
    
    
    
    
    
    
    
    
    
    
    




