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
    
    

#  n文件
"""
将对应行数的数据作为参数读入并进行初始化解析

已满足的数据有：
    GPS
    
"""
#GPS广播星历记录解析类
class GPS_brdc_record():
    def __init__(self,data):
        #解析第一行
        #self.serial_no=int(data[0][0:2])
        self.serial_no=data[0][0:2]
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


# 读取GPS的n文件,返回GPS_brdc_record类的记录对象
def read_GPS_nFile(GPS_nFile):
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
        GPS_brdc_records.append(record)
    nf.close()
    return GPS_brdc_records

    
#Renix3.04广播星历记录解析类
class Renix304_brdc_record():
    def __init__(self,data):
        #解析第一行
        #self.serial_no=int(data[0][0:2])
        self.serial_no=data[0][0:2]
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


# 读取Renix3.04的n文件,返回Renix304_brdc_record类的记录对象
def read_Renix304_nFile(Renix304_nFile):
    """
    Parameters
    ----------
        Renix304_nFile : 文件路径(需包含 r"",如r"E:\大三下\卫星与导航定位\data\brdc3120.20n")   
    Returns
    -------
        list[Renix304_brdc_record class]
    """
    nf=open(Renix304_nFile,"r")
    Renix304_brdc_records=[]     #存储GPS_brdc_record的列表
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
        record=Renix304_brdc_record(data)
        Renix304_brdc_records.append(record)
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
    pef=open(GPS_sp3File,"r")
    for i in range(3):
        line=pef.readline()
    sat_num=int(line[4:6])
    for i in range(3,22):
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

已满足的数据有：
    GPS
    
"""
#GPS观测记录解析类
class GPS_observation_record:
    def __init__(self,system,PRN,time,data):
        self.system=system
        self.PRN=PRN
        self.time=time
        self.data=data

# 读取GPS的o文件,返回GPS_observation_record类的记录对象
def read_GPS_oFile(GPS_oFile):
    """
    Parameters
    ----------
        GPS_nFile : 文件路径(需包含 r"",如r"E:\大三下\卫星与导航定位\data\brdc3120.20n")   
    Returns
    -------
        list[GPS_brdc_record class]
        
    """
    of=open(GPS_oFile, "r")
    GPS_observation_records = []     # 存储GPS_observation_record的列表
    #读取文件头中的观测记录对象描述
    header = []
    line=of.readline()
    header.append(line.strip())
    # 反复读入行直到观测类型记录行
    while line.strip()[-15:]!="TYPES OF OBSERV":
        line=of.readline()
        header.append(line.strip())
    # o文件中观测记录的观测对象个数
    observation_num=int(line[:6])
    # 文件头第12行处开始的观测类型行数
    header_ob_linenum=math.ceil(observation_num/9)    #文件头中观测类型行数
    lastline_header_ob_num=observation_num%9    #文件头中观测类型最后一行的个数
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
    #读取到文件头的尾部
    line=of.readline()
    header.append(line.strip())
    while line.strip()!="END OF HEADER":
        line=of.readline()
        header.append(line.strip())
    #开始读取主体数据
    main_data=of.readlines()     #全部主体数据,以列表存储
    of.close()
    firstline_cursor=0
    while firstline_cursor<len(main_data):
        if main_data[firstline_cursor+1]=="RINEX FILE SPLICE                                           COMMENT\n":
            firstline_cursor+=2
        year=int("20"+main_data[firstline_cursor][1:3])
        month=int(main_data[firstline_cursor][4:6])
        day=int(main_data[firstline_cursor][7:9])
        hour=int(main_data[firstline_cursor][10:12])
        minute=int(main_data[firstline_cursor][13:15])
        second=eval(main_data[firstline_cursor][15:26])
        time=datetime.datetime(year,month,day,hour,minute)+datetime.timedelta(seconds=second)
        #此处可以加一个对卫星状态的判断
        sat_num=int(main_data[firstline_cursor][29:32])  #观测卫星个数,即PRN个数
        PRN_row_num=math.ceil(sat_num/12)   #得到卫星PRN记录行数
        last_row_col=sat_num%12   #最后一行卫星PRN对应的个数
        cursor=firstline_cursor+PRN_row_num    #确定卫星观测数据所起始行
        n_of_line=1     #在每行所读prn的个数
        n=0       #所读prn所在的行数
        while n<PRN_row_num-1:  #未读到PRN记录的最后一行
           the_sprn=main_data[firstline_cursor+n][29+3*n_of_line:32+3*n_of_line]
           the_data={}
           #读入每一个观测类型的数据
           for i in range(len(observation_types)):
               l=(i+1)//5   #目前读取到的卫星观测值在记录块中的行数
               k=(i+1)%5   #目前读取到的卫星观测数据所在行的位置数
               if k==0:
                   l=l-1
                   k=5
               data_detail={}
               if main_data[cursor+l][16*(k-1):16*k-2].strip()!="":  
                   data_detail['observation']=float(eval(main_data[cursor+l][16*(k-1):16*k-2]))
               else:
                   data_detail['observation']=""
               data_detail['LLI']=main_data[cursor+l][16*k-2:16*k-1]
               data_detail['Signal strength']=main_data[cursor+l][16*k-1:16*k]
               the_data["%s"%observation_types[i]]=data_detail
           obrecord=GPS_observation_record(the_sprn[0],the_sprn[1:],time,the_data)
           GPS_observation_records.append(obrecord)
           n_of_line+=1
           cursor+=gap
           if n_of_line==13:
               n+=1
               n_of_line=1
        while n_of_line<=last_row_col:
           the_sprn=main_data[firstline_cursor+n][29+3*n_of_line:32+3*n_of_line]
           the_data={}
           for i in range(len(observation_types)):
               l=(i+1)//5   #目前读取到的卫星观测值在记录块中的行数
               k=(i+1)%5   #目前读取到的卫星观测数据所在行的位置数
               if k==0:
                   l=l-1
                   k=5
               data_detail={}
               if main_data[cursor+l][16*(k-1):16*k-2].strip()!="":  
                   data_detail['observation']=float(eval(main_data[cursor+l][16*(k-1):16*k-2]))
               else:
                   data_detail['observation']=""
               data_detail['LLI']=main_data[cursor+l][16*k-2:16*k-1].strip()
               data_detail['Signal strength']=main_data[cursor+l][16*k-1:16*k].strip()
               the_data["%s"%observation_types[i]]=data_detail
           obrecord=GPS_observation_record(the_sprn[0],the_sprn[1:],time,the_data)
           GPS_observation_records.append(obrecord)
           n_of_line+=1
           cursor+=gap
        firstline_cursor+=(PRN_row_num+gap*sat_num)
    return GPS_observation_records
    

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
    return clk_satellite_records,clk_receiver_records
















