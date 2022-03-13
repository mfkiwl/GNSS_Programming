import utils.TimeSystem as TimeSystem
import utils.RecordFilter as RecordFilter
import math
import utils.const as const

# 以某条记录(GPSweek和GPSsecond)为基础,计算卫星钟差(包含相对论效应)
def cal_ClockError_BDS_GPSws_withRelativisticEffect(time, SVN, brs):
    """
    Parameters
    ----------
        time : GPSws,所求时刻的GPSws类的GPS时间
        SVN : str,卫星的SVN号
        brs : list[GPS_brdc_record class],所依据的广播星历记录
    Returns
    -------
        clockerror : 卫星钟差,单位s
    """
    # 筛选出最接近的一条记录
    # 由GPSws类数据得到datetime.datetime类型时间,便于筛选数据
    datetime_time = TimeSystem.from_GPSws_cal_datetime_2(time)
    BDSws = TimeSystem.from_GPSws_get_BDSws(time)
    br = RecordFilter.find_closest_record(brs, datetime_time, SVN)
    # (1)Time from ephemeris epoch toe
    tk = BDSws.BDSSecond-br.toe
    # (2)Calculate the semimajor axis
    a = br.sqrt_a**2
    # (3)Compute mean motion – rad/s
    n0 = math.sqrt(const.miu/a**3)
    # (4)Correct mean motion
    n = n0+br.delta_n
    # (5)Mean anomaly Mk at epoch t
    Mk = br.M0+n*tk
    # (6)Eccentricity anomaly Ek at epoch t
    Ek = cal_Ek(Mk, br.e)
    # 计算钟差偏移和漂移部分
    clockerror_biasdrift = br.a0+br.a1*tk+br.a2*tk**2
    # 计算钟差相对论效应部分
    clockerror_re = -4.443e-10*br.e*math.sqrt(a)*math.sin(Ek)
    # 合并钟差
    clockerror = clockerror_biasdrift+clockerror_re
    return clockerror