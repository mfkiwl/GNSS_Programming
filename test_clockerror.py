
import utils.RecordFilter as rf
import utils.SatellitePosition as sp
import utils.DoFile as df
import utils.TimeSystem as ts
import datetime
import matplotlib.pyplot as plt

broadcast_file = r"edata\sat_obit\brdc3100.20n"
br_records = df.read_GPS_nFile(broadcast_file)
serial_no = '08'
dts = []
Tr = datetime.datetime(2020, 11, 5, 0, 00, 0)

while Tr < datetime.datetime(2020, 11, 5, 23, 20, 16):
    w, s = ts.from_datetime_cal_GPSws(Tr)
    Tr_GPSws = ts.GPSws(w, s)
    br_records = df.read_GPS_nFile(broadcast_file)
    # br_records = rf.GPSBrdcRecord_HourIntegerRecord_Filter(br_records)
    dt = sp.cal_ClockError_GPS_GPSws(Tr_GPSws, serial_no, br_records)
    dts.append(dt)
    Tr += datetime.timedelta(seconds=30)

print(dts)
plt.plot(dts)
plt.show()

