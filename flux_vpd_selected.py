import os
import numpy as np
import pandas as pd
import xarray as xr
import netCDF4 as nc
import calendar


#判断闰年/平年。
def year_days(year):
    try:
        if calendar.isleap(year):
            year_n = 366
        else:
            year_n = 365
    except:
        print("判断是否闰年错误！")
    return year_n


fluxdir = 'D:/data/PLUMBER2/flux_no1/'    #PLUMBER2 Flux 文件
metdir = 'D:/data/PLUMBER2/met_no1/'      #PLUMBER2 Met 文件
files = os.listdir(fluxdir)
outsite_year = {}    #用于存储符合标准的站点和年份


for file in files:
    # 判断站点观测分辨率是半小时还是一小时
    data_time = nc.Dataset(fluxdir+file,'r')
    ts   = data_time['time'][:] 
    if ts[1] == 1800:
        day_n = 48
    elif ts[1] == 3600:
        day_n = 24
    else:
        print('Error!')


    qle_year = []   #潜热可用年份
    qh_year = []    #感热可用年份
    data = xr.open_dataset(fluxdir+file)
    # 对qc信息以年为单位进行分组
    h_group = data['Qh_qc'].resample(time = 'Y')
    lh_group = data['Qle_qc'].resample(time = 'Y')
    # 挑选感热qc>1比例不超过10%的年份
    for name,h in h_group:
        name = pd.to_datetime(name).to_pydatetime()
        year = name.year    # 年份
        days = year_days(year)    # 天数
        h = h.values.flatten()
        h[h==1] = 0
        gap_n = len(np.where(h>0)[0])    # qc>1的数量
        gappercent = gap_n/(days*day_n)    # qc>1的比例
        # 把qc>1的比例设为不超过10%。
        if gappercent <= 0.1:
            qle_year.append(year)
    # 挑选潜热qc>1比例不超过10%的年份
    for name,lh in lh_group:
        name = pd.to_datetime(name).to_pydatetime()
        year = name.year
        days = year_days(year)
        lh = lh.values.flatten()
        lh[lh==1] = 0
        gap_n = len(np.where(lh>0)[0])
        gappercent = gap_n/(days*day_n)
        if gappercent <= 0.1:
            qh_year.append(year)

    # 通过潜热、感热筛选后可用的年份
    available_year = list(set(qle_year) & set(qh_year))    #求取交集
    available_year.sort()


    # 对于OzFlux站点，不需要对VPD进行筛选。
    if file[0:2] == 'AU':
        if len(available_year) != 0:
            outsite_year[file[0:6]] = available_year    #存储可用站点与年份
        continue    #挑出本次循环，处理下一个站点。
    
    # 挑选VPD qc>0 比例不超过10%的年份
    file = file.replace('Flux.nc', 'Met.nc') 
    data = xr.open_dataset(metdir+file)    # 读取气象数据
    vpd_year = []    #用于存储VPD gap-filled大于10%的年份
    vpd_group = data['VPD_qc'].resample(time = 'Y')
    for name,vpd in vpd_group:
        name = pd.to_datetime(name).to_pydatetime()
        year = name.year
        days = year_days(year)
        vpd = vpd.values.flatten()
        gap_n = len(np.where(vpd>0)[0])
        gappercent = gap_n/(days*day_n)
        if gappercent <= 0.1:
            vpd_year.append(year)

    
    all_avai = list(set(available_year) & set(vpd_year))    #求取交集
    all_avai.sort()
    if len(all_avai) != 0:
        outsite_year[file[0:6]] = all_avai



print(outsite_year)    








