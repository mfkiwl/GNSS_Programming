# -*- coding: utf-8 -*-
"""

@title:	Style Guide for Python Code
@author: iDeal0103
@status:	Active
@type:	Process
@created:	13-Apr-2021
@post-History:	13-Apr-2021

comment:
    1.UTC-GPS时间转换
    2.GPS有两种时间格式：
        A.datetime型：用datetime.datetime记录,包含(年,月,日,时,分,秒,毫秒),精度只能到1毫秒,
         无法记录小于1毫秒的部分,适合记录星历和观测文件上的时间、或以秒为间隔的时间序列,不适合
         做对时间精确性要求较高的计算(如钟差改正等)
        B.GPSweek+GPSsecond型：用GPSws类型记录,可提供较高精度的GPS秒数据
"""

#import module
import math
import datetime


'''记录GpsWeek和GpsSecond数据'''
class GPSws:
    def __init__(self, GpsWeek="", GpsSecond=""):
        if GpsWeek!="" and GpsSecond!="":
            self.GpsWeek = GpsWeek
            self.GpsSecond = GpsSecond

    def from_GPSdatetime(self, GPSdatetime):
        GPS_start_time = datetime.datetime()
        # todo: 完成GPSdatetime和GPSws之间的转换


    def minus(self, second):
        # 直接在原对象上进行减操作
        self.GpsSecond = self.GpsSecond-second
        while self.GpsSecond < 0:
            self.GpsWeek -= 1
            self.GpsSecond += 86400*7
        return self

    def cal_minus_result(self, second):
        # 只计算减的结果，不对原对象进行操作
        GPSweek_result = self.GpsWeek
        GPSsecond_result = self.GpsSecond-second
        while GPSsecond_result < 0:
            GPSweek_result -= 1
            GPSsecond_result += 86400*7
        return GPSws(GPSweek_result, GPSsecond_result)

    def add(self, second):
        # 直接在原对象上进行减操作
        self.GpsSecond = self.GpsSecond-second
        while self.GpsSecond > 604800:
            self.GpsWeek += 1
            self.GpsSecond -= 86400*7
        return self

    def cal_add_result(self, second):
        # 只计算减的结果，不对原对象进行操作
        GPSweek_result=self.GpsWeek
        GPSsecond_result=self.GpsSecond+second
        while GPSsecond_result > 604800:
            GPSweek_result += 1
            GPSsecond_result -= 86400*7
        return GPSws(GPSweek_result, GPSsecond_result)


'''记录BDSWeek和BDSsecond数据'''
class BDSws:
    def __init__(self, BDSWeek, BDSSecond):
        self.BDSWeek = BDSWeek
        self.BDSSecond = BDSSecond

    def minus(self, second):
        # 直接在原对象上进行减操作
        self.BDSSecond = self.BDSSecond-second
        while self.BDSSecond < 0:
            self.BDSWeek -= 1
            self.BDSSecond += 86400*7
        return self

    def cal_minus_result(self, second):
        # 只计算减的结果，不对原对象进行操作
        BDSweek_result = self.BDSWeek
        BDSsecond_result = self.BDSSecond-second
        while BDSsecond_result < 0:
            BDSweek_result -= 1
            BDSsecond_result += 86400*7
        return BDSws(BDSweek_result, BDSsecond_result)

    def add(self, second):
        # 直接在原对象上进行减操作
        self.BDSSecond = self.BDSSecond-second
        while self.BDSSecond > 604800:
            self.BDSWeek += 1
            self.BDSSecond -= 86400*7
        return self

    def cal_add_result(self, second):
        # 只计算减的结果，不对原对象进行操作
        BDSweek_result=self.BDSWeekQ
        BDSsecond_result=self.BDSSecond+second
        while BDSsecond_result > 604800:
            BDSweek_result += 1
            BDSsecond_result -= 86400*7
        return BDSws(BDSweek_result, BDSsecond_result)

# GPSws和BDSws的转换
def from_GPSws_get_BDSws(GPSws):
    GPSweek = GPSws.GpsWeek
    GPSsecond = GPSws.GpsSecond
    BDSweek = GPSweek - 1356
    BDSsecond = GPSsecond - 14
    return BDSws(BDSweek, BDSsecond)

def from_BDSws_get_GPSws(BDSws):
    BDSweek = BDSws.BDSWeek
    BDSsecond = BDSws.BDSSecond
    GPSweek = BDSweek + 1356
    GPSsecond = BDSsecond +14
    return GPSws(GPSweek, GPSsecond)

# 从GPSws和BDSws得到当周秒 (两者相同)
def get_ephemeris_second(ws):
    if isinstance(ws, GPSws):
        ephemeris_second = ws.GpsSecond
    elif isinstance(ws, BDSws):
        ephemeris_second = ws.BDSSecond
    else:
        return
    return ephemeris_second


# 计算timedelta对象的second值
def cal_deltatime_second(datetime_time):
    day = datetime_time.days
    sec = datetime_time.seconds
    micsec = datetime_time.microseconds
    delta_time = day*86400+sec+micsec/1000000
    return delta_time

# 计算两个GPSws对象的时间差
def cal_GPSws_deltatime_inseconds(GPSws1, GPSws2):
    GpsWeekDelta = GPSws1.GpsWeek - GPSws2.GpsWeek
    GpsSecondDelta = GPSws1.GpsSecond - GPSws2.GpsSecond
    delta_second = GpsWeekDelta * 7 * 86400 + GpsSecondDelta
    return delta_second


# 由儒略日JD计算GPS时间(包括GPS周和GPS秒)
def from_JD_cal_GPStime(JD):
    GpsWeek = int((JD - 2444244.5) / 7)
    GpsSecond = (JD - 2444244.5 - 7 * GpsWeek) * 86400.0
    return GpsWeek, GpsSecond

# 由日期格式时间计算儒略日JD时间
def from_datetime_cal_JD(date_time):
    """
    Parameters
    ----------
        datetime <某个日期格式时间,type->datetime.datetime>
    Returns
    -------
        JD <对应的儒略日时间>
    """
    # 从date_time时间得到年月日时分秒数据
    Y = date_time.year
    M = date_time.month
    D = date_time.day
    hour = date_time.hour
    minute = date_time.minute
    second = date_time.second
    microsecond = date_time.microsecond
    # 将minute和second化入hour
    H = hour+minute/60+(second+microsecond/1000000)/3600
    # H = hour + minute / 60
    if M <= 2:
        y = Y-1
        m = M+12
    else:
        y = Y
        m = M
    JD = int(365.25*y)+int(30.6001*(m+1))+D+H/24+1720981.5
    return JD


# 由儒略日JD时间计算GPS_date_time时间的年月日
def from_JD_cal_datetime(JD):
    a = int(JD+0.5)
    b = a+1537
    c = int((b-122.1)/365.25)
    d = int(365.25*c)
    e = int((b-d)/30.6001)
    D = b-d-int(30.60001*e)+math.modf(JD+0.5)[0]
    M = e-1-12*int(e/14)
    Y = c-4715-int((7+M)/10)
    # 将D化为日(day)时(hour)分(minute)秒(second)
    day = int(math.modf(D)[1])
    seconds = math.modf(D)[0]*86400   # 日的小数部分化为秒
    date_time = datetime.datetime(Y, M, day)+datetime.timedelta(seconds=seconds)
    return date_time


#由GPS的date_time时间计算GPS周和GPS秒
def from_datetime_cal_GPSws(datetime):
    JD = from_datetime_cal_JD(datetime)
    GPS_week, GPS_second = from_JD_cal_GPStime(JD)
    return GPS_week, GPS_second


#由UTC时间计算BD周和BD秒
def from_UTCtime_cal_BDtime_2(UTCtime):
    JD=from_datetime_cal_JD(UTCtime)
    GPS_week,GPS_second = from_JD_cal_GPStime(JD)
    BD_week=GPS_week-1356
    BD_second=GPS_second-14
    return BD_week,BD_second


def from_GPStime_cal_JD(GPS_week, GPS_second):
    JD = GPS_second/86400.0 + 7*GPS_week + 2444244.5
    return JD

#由GPS周和GPS秒计算GPS的datetime时间
def from_GPSws_cal_datetime_2(GPSws):
    GPS_week = GPSws.GpsWeek
    GPS_second = GPSws.GpsSecond
    JD = from_GPStime_cal_JD(GPS_week, GPS_second)
    date = from_JD_cal_datetime(JD)
    return date

#由UTC时间计算GPS系统时间
def from_UTCtime_cal_GPStime_1(UTCtime):
    if UTCtime>datetime.datetime(2017,1,1):
        GPStime=UTCtime+datetime.timedelta(seconds=18)
    elif UTCtime>datetime.datetime(2015,7,1):
        GPStime=UTCtime+datetime.timedelta(seconds=17)
    elif UTCtime>datetime.datetime(2012,7,1):
        GPStime=UTCtime+datetime.timedelta(seconds=16)
    elif UTCtime>datetime.datetime(2009,1,1):
        GPStime=UTCtime+datetime.timedelta(seconds=15)
    elif UTCtime>datetime.datetime(2006,1,1):
        GPStime=UTCtime+datetime.timedelta(seconds=14)
    elif UTCtime>datetime.datetime(1999,1,1):
        GPStime=UTCtime+datetime.timedelta(seconds=13)
    elif UTCtime>datetime.datetime(1997,7,1):
        GPStime=UTCtime+datetime.timedelta(seconds=12)
    elif UTCtime>datetime.datetime(1996,1,1):
        GPStime=UTCtime+datetime.timedelta(seconds=11)
    elif UTCtime>datetime.datetime(1994,7,1):
        GPStime=UTCtime+datetime.timedelta(seconds=10)
    elif UTCtime>datetime.datetime(1993,7,1):
        GPStime=UTCtime+datetime.timedelta(seconds=9)
    elif UTCtime>datetime.datetime(1992,7,1):
        GPStime=UTCtime+datetime.timedelta(seconds=8)
    elif UTCtime>datetime.datetime(1991,1,1):
        GPStime=UTCtime+datetime.timedelta(seconds=7)
    elif UTCtime>datetime.datetime(1990,1,1):
        GPStime=UTCtime+datetime.timedelta(seconds=6)
    elif UTCtime>datetime.datetime(1988,1,1):
        GPStime=UTCtime+datetime.timedelta(seconds=5)
    elif UTCtime>datetime.datetime(1985,7,1):
        GPStime=UTCtime+datetime.timedelta(seconds=4)
    elif UTCtime>datetime.datetime(1983,7,1):
        GPStime=UTCtime+datetime.timedelta(seconds=3)
    elif UTCtime>datetime.datetime(1982,7,1):
        GPStime=UTCtime+datetime.timedelta(seconds=2)
    elif UTCtime>datetime.datetime(1981,7,1):
        GPStime=UTCtime+datetime.timedelta(seconds=1)
    elif UTCtime>datetime.datetime(1980,1,1):
        GPStime=UTCtime+datetime.timedelta(seconds=0)
    return GPStime

#由GPS时间计算UTC系统时间    //未完成
# def from_GPStime_cal_UTCtime_1(GPStime):
#     if GPStime>datetime.datetime(2017,):
#         UTCtime=GPStime+datetime.timedelta(seconds=18)
#     elif GPStime>datetime.datetime(2015,7,1):
#         UTCtime=UTCtime+datetime.timedelta(seconds=17)
#     elif UTCtime>datetime.datetime(2012,7,1):
#         GPStime=UTCtime+datetime.timedelta(seconds=16)
#     elif UTCtime>datetime.datetime(2009,1,1):
#         GPStime=UTCtime+datetime.timedelta(seconds=15)
#     elif UTCtime>datetime.datetime(2006,1,1):
#         GPStime=UTCtime+datetime.timedelta(seconds=14)
#     elif UTCtime>datetime.datetime(1999,1,1):
#         GPStime=UTCtime+datetime.timedelta(seconds=13)
#     elif UTCtime>datetime.datetime(1997,7,1):
#         GPStime=UTCtime+datetime.timedelta(seconds=12)
#     elif UTCtime>datetime.datetime(1996,1,1):
#         GPStime=UTCtime+datetime.timedelta(seconds=11)
#     elif UTCtime>datetime.datetime(1994,7,1):
#         GPStime=UTCtime+datetime.timedelta(seconds=10)
#     elif UTCtime>datetime.datetime(1993,7,1):
#         GPStime=UTCtime+datetime.timedelta(seconds=9)
#     elif UTCtime>datetime.datetime(1992,7,1):
#         GPStime=UTCtime+datetime.timedelta(seconds=8)
#     elif UTCtime>datetime.datetime(1991,1,1):
#         GPStime=UTCtime+datetime.timedelta(seconds=7)
#     elif UTCtime>datetime.datetime(1990,1,1):
#         GPStime=UTCtime+datetime.timedelta(seconds=6)
#     elif UTCtime>datetime.datetime(1988,1,1):
#         GPStime=UTCtime+datetime.timedelta(seconds=5)
#     elif UTCtime>datetime.datetime(1985,7,1):
#         GPStime=UTCtime+datetime.timedelta(seconds=4)
#     elif UTCtime>datetime.datetime(1983,7,1):
#         GPStime=UTCtime+datetime.timedelta(seconds=3)
#     elif UTCtime>datetime.datetime(1982,7,1):
#         GPStime=UTCtime+datetime.timedelta(seconds=2)
#     elif UTCtime>datetime.datetime(1981,7,1):
#         GPStime=UTCtime+datetime.timedelta(seconds=1)
#     elif UTCtime>datetime.datetime(1980,1,1):
#         GPStime=UTCtime+datetime.timedelta(seconds=0)
#     return GPStime
    
#由UTC时间计算BD系统时间   
def from_UTCtime_cal_BDtime_1(UTCtime):
    GPStime=from_UTCtime_cal_GPStime_1(UTCtime)
    BDtime=GPStime-datetime.timedelta(seconds=14)
    return BDtime



