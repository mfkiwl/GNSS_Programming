# -*- coding: utf-8 -*-
"""

@title:	Style Guide for Python Code
@author: iDeal0103
@status:	Active
@type:	Process
@created:	13-Apr-2021
@post-History:	13-Apr-2021

comment:
    1.读取brdc广播文件（GPS）
    2.读取精密星历文件（igs文件)
    3.读取观测文件(O文件)
    4.读取
    
todo:
    1.读取BDS等其他系统卫星的广播星历
    2.读取多源广播星历
    3.
    
    
"""

#将utils加入sys.path
# import sys
# import os
# work_path="../"
# utils_path=os.path.join(work_path, "utils")
# if utils_path not in sys.path:
#     sys.path.append(utils_path)

# import
import datetime
import math
# import TimeSystem


#将文件数据中的科学计数法底数"D",变为"e" 
#如 3.234D+10 -> 3.234e+10
def parse_Dstring(D_string):
    e_string=D_string.replace("D","e")
    e_data=float(e_string)
    return e_data

#
def timeString2UTC_1(timeString):
    """
    Parameters
    ----------
        timeString : GPS_n文件中的toc字符部分,如" 20 11  7  0  0  0.0"
    Returns
    -------
        UTCtime : datetime.datetime,对应toc的UTC时刻
    """
    tslist=timeString.split()
    tslist[0]="20"+tslist[0]
    UTCtime=datetime.datetime(int(tslist[0]),int(tslist[1]),
        int(tslist[2]),int(tslist[3]),int(tslist[4]))
    +datetime.timedelta(seconds=float(tslist[5]))
    return UTCtime

def timeString2UTC_2(timeString):
    """
    Parameters
    ----------
        timeString : renix文件中的toc字符部分,如" 2020 11  7  0  0  0.0"
    Returns
    -------
        UTCtime : datetime.datetime,对应toc的UTC时刻
    """
    tslist=timeString.split()
    tslist[0]=tslist[0]
    UTCtime=datetime.datetime(int(tslist[0]),int(tslist[1]),
        int(tslist[2]),int(tslist[3]),int(tslist[4]))
    +datetime.timedelta(seconds=float(tslist[5]))
    return UTCtime
    
def get_standard_SVN(system, no):
    if isinstance(no, str):
        if no != "":
            SVN = str(system) + '{:0>2d}'.format(int(no))
    elif isinstance(no, int):
        SVN = str(system) + '{:0>2d}'.format(no)
    return SVN


#  n文件
"""
将对应行数的数据作为参数读入并进行初始化解析

已满足的数据有：
    GPS
    
"""
#GPS广播星历记录解析类
class GPS_brdc_record():
    def __init__(self, data):
        #解析第一行
        #self.serial_no=int(data[0][0:2])
        self.SVN=get_standard_SVN("G", data[0][0:2])
        self.toc=timeString2UTC_1(data[0][2:22])
        self.a0=parse_Dstring(data[0][22:41])
        self.a1=parse_Dstring(data[0][41:60])
        self.a2=parse_Dstring(data[0][60:79])
        #解析第二行
        self.IODE=parse_Dstring(data[1][3:22])
        self.Crs=parse_Dstring(data[1][22:41])
        self.delta_n=parse_Dstring(data[1][41:60])
        self.M0=parse_Dstring(data[1][60:79])
        #解析第三行
        self.Cuc=parse_Dstring(data[2][3:22])
        self.e=parse_Dstring(data[2][22:41])
        self.Cus=parse_Dstring(data[2][41:60])
        self.sqrt_a=parse_Dstring(data[2][60:79])
        #解析第四行
        self.toe=parse_Dstring(data[3][3:22])
        self.Cic=parse_Dstring(data[3][22:41])
        self.omega0=parse_Dstring(data[3][41:60])
        self.Cis=parse_Dstring(data[3][60:79])
        #解析第五行
        self.i0=parse_Dstring(data[4][3:22])
        self.Crc=parse_Dstring(data[4][22:41])
        self.w=parse_Dstring(data[4][41:60])
        self.omega_dot=parse_Dstring(data[4][60:79])
        #解析第六行
        self.i_dot=parse_Dstring(data[5][3:22])
        self.Code_L2=parse_Dstring(data[5][22:41])
        self.GPS_week=parse_Dstring(data[5][41:60])
        self.L2_flag=parse_Dstring(data[5][60:79])
        #解析第七行
        self.SV_accuracy=parse_Dstring(data[6][3:22])
        self.SV_health=parse_Dstring(data[6][22:41])
        self.TGD=parse_Dstring(data[6][41:60])
        self.IODC=parse_Dstring(data[6][60:79])




# 读取GPS的n文件,返回GPS_brdc_record类的记录对象
def read_GPS_nFile(GPS_nFile, health_control=True):
    """
    Parameters
    ----------
        GPS_nFile : 文件路径(需包含 r"",如r"E:\大三下\卫星与导航定位\data\brdc3120.20n")   
    Returns
    -------
        list[GPS_brdc_record class]
    """
    nf=open(GPS_nFile,"r")
    GPS_brdc_records=[]     #存储GPS_brdc_record的列表
    #读取文件头
    linedata=nf.readline()
    while linedata.strip()!="END OF HEADER":
        linedata=nf.readline()
    #计入数据文件部分
    linedata=nf.readline()
    while linedata.strip()!="":
        data=[]
        for j in range(8):
            data.append(linedata)
            linedata=nf.readline()
        record=GPS_brdc_record(data)
        if health_control:
            if record.SV_health == 0:
                GPS_brdc_records.append(record)
            else:
                print(str(record.toe)+"时刻卫星%s不够健康"%record.SVN)
        else:
            GPS_brdc_records.append(record)
    nf.close()
    return GPS_brdc_records

    
# Renix3.04广播星历记录解析类
# 定义存储各系统导航星历单条记录的行数
navigation_record_lines_number = {'G': 8, 'C': 8, 'E': 8, 'R': 4, 'J': 8, 'S': 4, 'I': 8}

class Renix304_navigation_record_GPS():
    def __init__(self, data):
        self.system = 'G'
        # 解析第一行
        self.SVN=data[0][0:3]
        self.toc=timeString2UTC_2(data[0][3:23])
        self.a0=float(data[0][23:42])
        self.a1=float(data[0][42:61])
        self.a2=float(data[0][61:80])
        # 解析第二行
        self.IODE=float(data[1][4:23])
        self.Crs=float(data[1][23:42])
        self.delta_n=float(data[1][42:61])
        self.M0=float(data[1][61:80])
        # 解析第三行
        self.Cuc=float(data[2][4:23])
        self.e=float(data[2][23:42])
        self.Cus=float(data[2][42:61])
        self.sqrt_a=float(data[2][61:80])
        # 解析第四行
        self.toe=float(data[3][4:23])
        self.Cic=float(data[3][23:42])
        self.omega0=float(data[3][42:61])
        self.Cis=float(data[3][61:80])
        # 解析第五行
        self.i0=float(data[4][4:23])
        self.Crc=float(data[4][23:42])
        self.w=float(data[4][42:61])
        self.omega_dot=float(data[4][61:80])
        # 解析第六行
        self.i_dot=float(data[5][4:23])
        self.Code_L2=float(data[5][23:42])
        self.GPS_week=float(data[5][42:61])
        self.L2_flag=float(data[5][61:80])
        #解析第七行
        self.SV_accuracy=float(data[6][4:23])
        self.SV_health=float(data[6][23:42])
        self.TGD=float(data[6][42:61])
        self.IODC=float(data[6][61:80])


class Renix304_navigation_record_BDS():
    def __init__(self, data):
        self.system = 'C'
        # 解析第一行
        self.SVN = data[0][0:3]
        self.toc = timeString2UTC_2(data[0][3:23])
        self.a0 = float(data[0][23:42])
        self.a1 = float(data[0][42:61])
        self.a2 = float(data[0][61:80])
        # 解析第二行
        self.AODE = float(data[1][4:23])
        self.Crs = float(data[1][23:42])
        self.delta_n = float(data[1][42:61])
        self.M0 = float(data[1][61:80])
        # 解析第三行
        self.Cuc = float(data[2][4:23])
        self.e = float(data[2][23:42])
        self.Cus = float(data[2][42:61])
        self.sqrt_a = float(data[2][61:80])
        # 解析第四行
        self.toe = float(data[3][4:23])
        self.Cic = float(data[3][23:42])
        self.omega0 = float(data[3][42:61])
        self.Cis = float(data[3][61:80])
        # 解析第五行
        self.i0 = float(data[4][4:23])
        self.Crc = float(data[4][23:42])
        self.w = float(data[4][42:61])
        self.omega_dot = float(data[4][61:80])
        # 解析第六行
        self.i_dot = float(data[5][4:23])
        self.spare1 = float(data[5][23:42])
        self.BDS_week = float(data[5][42:61])
        self.spare2 = float(data[5][61:80])


# 读取Renix3.04的n文件,返回Renix304_brdc_record类的记录对象
def read_Renix304_nFile(Renix304_nFile, health_control=True):
    """
    Parameters
    ----------
        Renix304_nFile : 文件路径(需包含 r"",如r"E:\大三下\卫星与导航定位\data\brdc3120.20n")   
    Returns
    -------
        list[Renix304_brdc_record class]
    """
    nf = open(Renix304_nFile, "r")
    Renix304_brdc_records = []      # 存储GPS_brdc_record的列表
    # 读取文件头
    linedata=nf.readline()
    while linedata.strip() != "END OF HEADER":
        linedata=nf.readline()
    # 计入数据文件部分
    linedata=nf.readline()
    while linedata.strip() != "":
        data = []
        system = linedata[0]
        lines_number = navigation_record_lines_number[system]
        for j in range(lines_number):
            data.append(linedata)
            linedata=nf.readline()
        # GPS
        if system == 'G':
            record = Renix304_navigation_record_GPS(data)
            if health_control:
                if record.SV_health == 0:
                    Renix304_brdc_records.append(record)
                else:
                    print(str(record.toe)+"时刻卫星%s不够健康"%record.SVN)
            else:
                Renix304_brdc_records.append(record)
        # BDS
        elif system == 'C':
            record = Renix304_navigation_record_BDS(data)
            Renix304_brdc_records.append(record)
        # GALILEO
        elif system == 'E':
            pass
        # GLONASS
        elif system == 'R':
            pass
        # QZSS
        elif system == 'J':
            pass
        # SBAS
        elif system == 'S':
            pass
        # IRNSS
        elif system == 'I':
            pass
    nf.close()
    return Renix304_brdc_records



# sp3文件
"""
将对应行数的数据作为参数读入并进行初始化解析

已满足的数据有：
    GPS
    
"""
#GPS精密星历数据记录类
class satelliteobit_XYZdT:
    def __init__(self,time,system,PRN,X,Y,Z,dT):
        self.time=time
        self.system=system
        self.PRN=PRN
        self.X=X
        self.Y=Y
        self.Z=Z
        self.dT=dT

# 读取GPS的sp3文件,返回satelliteobit_XYZdT类的记录对象
def read_GPS_sp3File(GPS_sp3File):
    pef=open(GPS_sp3File, "r")
    for i in range(3):
        line=pef.readline()
    sat_num=int(line[4:6])
    for i in range(3, 22):
        line=pef.readline()
    all_records=pef.readlines()
    pef.close()
    #读取并解析数据
    pe_records=[]
    record_cursor=0
    while record_cursor<len(all_records)-1:
        time_record=all_records[record_cursor]
        year=int(time_record[3:7])
        month=int(time_record[8:10])
        day=int(time_record[11:13])
        hour=int(time_record[14:16])
        minute=int(time_record[17:19])
        second=int(eval(time_record[20:]))
        time=datetime.datetime(year,month,day,hour,minute,second)
        for i in range(1,sat_num+1):
            pp=record_cursor+i
            pe_system=all_records[pp].split()[0][1]  #卫星系统名称
            pe_RPN=all_records[pp].split()[0][2:]  #卫星PRN号
            pe_X=eval(all_records[pp].split()[1])
            pe_Y=eval(all_records[pp].split()[2])
            pe_Z=eval(all_records[pp].split()[3])
            pe_dT=eval(all_records[pp].split()[4])
            pe_sa=satelliteobit_XYZdT(time,pe_system,pe_RPN,pe_X,pe_Y,pe_Z,pe_dT)
            pe_records.append(pe_sa)
        record_cursor+=sat_num+1
    return pe_records





#  o文件
"""
将对应行数的数据作为参数读入并进行初始化解析

可支持混合数据
    
"""
# 观测记录解析类
class observation_record:
    def __init__(self, SVN, time, data):
        self.SVN=SVN
        self.time=time
        self.data=data

# 读取Renix2的观测文件,返回observation_record类的记录对象
def read_Rinex2_oFile(Rinex2_oFile):
    """
    Parameters
    ----------
        Rinex2_oFile : 文件路径(需包含 r"",如r"edata\obs\wab23100.20o")
    Returns
    -------
        list[observation_record class]
        
    """
    of=open(Rinex2_oFile, "r")
    observation_records = []     # 存储observation_record的列表
    # 读取文件头中的观测记录对象描述
    header = []
    line = of.readline()
    header.append(line.strip())
    # 反复读入行直到观测类型记录行
    while line.strip()[-15:] != "TYPES OF OBSERV":
        line = of.readline()
        header.append(line.strip())
    # o文件中观测记录的观测对象个数
    observation_num=int(line[:6])
    # 文件头第12行处开始的观测类型行数
    header_ob_linenum=math.ceil(observation_num/9)    # 文件头中观测类型行数
    lastline_header_ob_num=observation_num%9    # 文件头中观测类型最后一行的个数
    # o文件中每条观测记录的行数
    gap=math.ceil(observation_num/5)
    # lastline_main_ob_num=observation_num%5
    # 读取文件中的观测数类型
    observation_types=[]
    # 只有一行观测类型的情况
    if header_ob_linenum == 1:
        cursor=6
        for i in range(lastline_header_ob_num):
            type_name=line[cursor:cursor+6].strip()
            observation_types.append(type_name)
            cursor+=6
    # 有多行观测类型的情况
    else:
        for j in range(header_ob_linenum-1):
            cursor=6
            for i in range(9):
                type_name=line[cursor:cursor+6].strip()
                observation_types.append(type_name)
                cursor+=6
            line=of.readline()
            header.append(line.strip())
            cursor=6
        for k in range(lastline_header_ob_num):
            type_name=line[cursor:cursor+6].strip()
            observation_types.append(type_name)
            cursor+=6
    # 读取到文件头的尾部
    line=of.readline()
    header.append(line.strip())
    while line.strip()!="END OF HEADER":
        line=of.readline()
        header.append(line.strip())
    # 开始读取主体数据
    main_data=of.readlines()     # 全部主体数据,以列表存储
    of.close()
    firstline_cursor=0
    while firstline_cursor < len(main_data):
        if main_data[firstline_cursor].split()[-1] == "COMMENT" or len(main_data[firstline_cursor].strip())==4:   # 第二个条件为splice之前的“4  1”这样的内容行
            while main_data[firstline_cursor].split()[-1] == "COMMENT" or len(main_data[firstline_cursor].strip())==4:
                firstline_cursor += 1
        # print(main_data[firstline_cursor])
        year=int("20"+main_data[firstline_cursor][1:3])
        month=int(main_data[firstline_cursor][4:6])
        day=int(main_data[firstline_cursor][7:9])
        hour=int(main_data[firstline_cursor][10:12])
        minute=int(main_data[firstline_cursor][13:15])
        second=eval(main_data[firstline_cursor][15:26])
        time=datetime.datetime(year, month, day, hour, minute)+datetime.timedelta(seconds=second)  # 只能精确到1us,观测文件上最小位为0.1us
        # 此处可以加一个对卫星状态的判断
        sat_num=int(main_data[firstline_cursor][29:32])  # 观测卫星个数,即PRN个数
        SVN_row_num=math.ceil(sat_num / 12)   # 得到卫星SVN记录行数
        last_row_col=sat_num%12   # 最后一行卫星SVN的个数
        if last_row_col==0:
            last_row_col=12
        cursor= firstline_cursor + SVN_row_num    # 确定卫星观测数据所起始行
        n_of_line=1     # 在每行所读SVN的个数
        n=0       # 所读SVN所在的行数
        while n<SVN_row_num-1:  # 未读到SVN记录的最后一行
           the_svn = main_data[firstline_cursor + n][29 + 3 * n_of_line:32 + 3 * n_of_line]
           the_data={}
           # 读入每一个观测类型的数据
           for i in range(len(observation_types)):
               l=(i+1)//5   # 目前读取到的卫星观测值在记录块中的行数
               k=(i+1)%5   # 目前读取到的卫星观测数据所在行的位置数
               if k==0:
                   l=l-1
                   k=5
               data_detail={}
               if main_data[cursor+l][16*(k-1):16*k-2].strip() != "":
                   data_detail['observation'] = float(eval(main_data[cursor+l][16*(k-1):16*k-2]))
               else:
                   data_detail['observation'] = ""
               data_detail['LLI'] = main_data[cursor+l][16*k-2:16*k-1].strip()
               data_detail['Signal strength'] = main_data[cursor+l][16*k-1:16*k].strip()
               the_data["%s"%observation_types[i]] = data_detail
           obrecord = observation_record(the_svn, time, the_data)
           observation_records.append(obrecord)
           n_of_line+=1
           cursor+=gap
           if n_of_line==13:
               n+=1
               n_of_line=1
        # 读取最后一行SVN对应的数据
        while n_of_line<=last_row_col:
           the_svn= main_data[firstline_cursor + n][29 + 3 * n_of_line:32 + 3 * n_of_line]
           the_data={}
           for i in range(len(observation_types)):
               l=(i+1) // 5   # 目前读取到的卫星观测值在记录块中的行数
               k=(i+1) % 5   # 目前读取到的卫星观测数据所在行的位置数
               if k == 0:
                   l=l-1
                   k=5
               data_detail={}
               if main_data[cursor+l][16*(k-1):16*k-2].strip() != "":
                   data_detail['observation']=float(eval(main_data[cursor+l][16*(k-1):16*k-2]))
               else:
                   data_detail['observation']=""
               data_detail['LLI']=main_data[cursor+l][16*k-2:16*k-1].strip()
               data_detail['Signal strength']=main_data[cursor+l][16*k-1:16*k].strip()
               the_data["%s"%observation_types[i]]=data_detail
           obrecord=observation_record(the_svn, time, the_data)
           observation_records.append(obrecord)
           n_of_line+=1
           cursor+=gap
        firstline_cursor += (SVN_row_num + gap * sat_num)
    return observation_records


# 读取Rinex3的观测文件,返回observation_record类的记录对象
def read_Rinex3_oFile(Rinex3_oFile):
    """
    Parameters
    ----------
        Rinex3_oFile : 文件路径(需包含 r"",如r"edata\obs\wab23100.20o")
    Returns
    -------
        list[observation_record class]

    """
    of = open(Rinex3_oFile, "r")
    observation_records = []  # 存储observation_record的列表
    # 读取文件头中的观测记录对象描述
    header = []
    line = of.readline()
    header.append(line.strip())
    # 反复读入行直到观测类型记录行
    while line.strip()[-9:] != "OBS TYPES":
        line = of.readline()
        header.append(line.strip())
    # 初始化观测值类型
    observation_types = {}
    C_observation_types = []
    E_observation_types = []
    G_observation_types = []
    R_observation_types = []
    S_observation_types = []
    observation_types['C'] = C_observation_types
    observation_types['E'] = E_observation_types
    observation_types['G'] = G_observation_types
    observation_types['R'] = R_observation_types
    observation_types['S'] = S_observation_types
    while line.strip()[-9:] == "OBS TYPES":
        system = line[0]   # 获取当前系统名
        observation_num = int(line[3:6])
        left_obs_type_lines = observation_num // 13
        # 读取首行观测类型
        if left_obs_type_lines == 0:   # 观测值类型小于13个
            for i in range(observation_num):
                observation_types[system].append(line[6+i*4:10+i*4].strip())
        else:
            # 首行
            for i in range(13):
                observation_types[system].append(line[6+i*4:10+i*4].strip())
            # 非首尾行
            for j in range(left_obs_type_lines-1):
                line = of.readline()
                for i in range(13):
                    observation_types[system].append(line[6+i*4:10+i*4].strip())
            # 尾行
            line = of.readline()
            for i in range(observation_num%13):
                observation_types[system].append(line[6+i*4:10+i*4].strip())
        line = of.readline()
    # 读取到文件头的尾部
    header.append(line.strip())
    while line.strip() != "END OF HEADER":
        line = of.readline()
        header.append(line.strip())
    # 开始读取主体数据
    line = of.readline()
    while line[0] == ">":
        # 获取当前历元时间
        year = int(line[1:6])
        month = int(line[6:9])
        day = int(line[9:12])
        hour = int(line[12:15])
        minute = int(line[15:18])
        second = eval(line[18:29])
        time = datetime.datetime(year, month, day, hour, minute) + datetime.timedelta(
            seconds=second)  # 只能精确到1us,观测文件上最小位为0.1us
        epoch_flag = int(line[29:32])
        obs_satellite_num = int(line[32:35])
        # receiver_clock_offset = eval(line[41:56])
        for j in range(obs_satellite_num):
            line = of.readline()
            SVN = line[0:3]
            the_system = SVN[0]
            the_data = {}
            for i in range(len(observation_types[the_system])):
                data_detail = {}
                if line[3+i*16:17+i*16].strip() != "":
                    data_detail['observation'] = float(eval(line[3+i*16:17+i*16]))
                else:
                    data_detail['observation'] = ""
                data_detail['LLI'] = line[17+i*16:18+i*16].strip()
                data_detail['Signal strength'] = line[18+i*16:19+i*16].strip()
                the_data[observation_types[the_system][i]] = data_detail
            obrecord = observation_record(SVN, time, the_data)
            observation_records.append(obrecord)
        line = of.readline()
        if line == "":
            break
    return observation_records


#  clk文件
"""
将对应行数的数据作为参数读入并进行初始化解析

已满足的数据有：
    GPS
    
"""
#clk观测记录解析类
class clk_satellite_record:
    def __init__(self,system,PRN,time,data):
        self.system_PRN=system
        self.PRN=PRN
        self.time=time
        self.data=data
        
class clk_receiver_record:
    def __init__(self,name,time,data):
        self.receiver=name
        self.time=time
        self.data=data

# 读取GPS的clk文件,返回GPS_observation_record类的记录对象
def read_GPS_clkFile(GPS_clkfile):    
    clkf=open(GPS_clkfile,"r")
    clk_satellite_records=[]     #存储clk_satellite_record的列表
    clk_receiver_records=[]     #存储clk_receiver_record的列表
    line=clkf.readline()
    #跳过头文件,如果后续有需要可以精细读取头文件
    while line.strip()!="END OF HEADER":
        line=clkf.readline()
    #读取各条记录
    line=clkf.readline()
    while line.strip()!="":
        #接收站记录
        if line[:2]=="AR":
            r_name=line[3:7].strip()
            year=int(line[8:12])
            month=int(line[12:15])
            day=int(line[15:18])
            hour=int(line[18:21])
            minute=int(line[21:24])
            second=float(line[24:34])
            time=datetime.datetime(year, month, day,hour,minute)+datetime.timedelta(seconds=second)
            data_num=int(line[34:37])
            data=[]
            for i in range(data_num):
                data.append(float(line[(40+20*i):(59+20*i)]))
            clk_receiver_records.append(clk_receiver_record(r_name,time,data))
            line=clkf.readline()
        #卫星记录
        elif line[:2]=="AS":
            system=line[3:4]
            PRN=line[4:7].strip()
            year=int(line[8:12])
            month=int(line[12:15])
            day=int(line[15:18])
            hour=int(line[18:21])
            minute=int(line[21:24])
            second=float(line[24:34])
            time=datetime.datetime(year, month, day,hour,minute)+datetime.timedelta(seconds=second)
            data_num=int(line[34:37])
            data=[]
            for i in range(data_num):
                data.append(float(line[(40+20*i):(59+20*i)]))
            clk_satellite_records.append(clk_satellite_record(system,PRN,time,data))
            line=clkf.readline()
    clkf.close()
    return clk_satellite_records, clk_receiver_records

if __name__ == "__main__" :
    obs = read_Rinex3_oFile(r"D:\Tongji_study\my_GNSS\GNSS_Programming\edata\obs\renix3\LEIJ00DEU_R_20213120000_01D_30S_MO.rnx")















