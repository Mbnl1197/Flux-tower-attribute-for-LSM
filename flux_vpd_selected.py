import os
import numpy as np
import pandas as pd
import xarray as xr
import netCDF4 as nc
import calendar


# determine if a year is a leap year or a common year
def year_days(year):
    try:
        if calendar.isleap(year):
            year_n = 366
        else:
            year_n = 365
    except:
        print("Error in leap year judgement")
    return year_n


fluxdir      = '/stu01/shijh21/data/forcingPLUMBER2/flux_no1/'    # PLUMBER2 Flux file
metdir       = '/stu01/shijh21/data/forcingPLUMBER2/met_no1/'     # PLUMBER2 Met file
files        = os.listdir(fluxdir)
outsite_year = {}    # Used for storing compliant sites and years


for file in files:

    #determine if the observation resolution of the site is
    # half an hour or one hour
    data_time = nc.Dataset(fluxdir+file,'r')
    ts   = data_time['time'][:]
    if ts[1] == 1800:
        day_n = 48
    elif ts[1] == 3600:
        day_n = 24
    else:
        print('Error!')


    qle_year = []    # selected years for latent heat flux
    qh_year  = []    # selected years for sensible heat flux
    data     = xr.open_dataset(fluxdir+file)

    # Group qc information by year
    h_group  = data['Qh_qc'].resample(time  = 'Y')
    lh_group = data['Qle_qc'].resample(time = 'Y')

    # Select years with sensible heat qc>1 not exceeding 10%.
    for name,h in h_group:
        name = pd.to_datetime(name).to_pydatetime()
        year = name.year          # year
        days = year_days(year)    # day
        h = h.values.flatten()
        h[h==1] = 0
        gap_n = len(np.where(h>0)[0])    # number of qc>1
        gappercent = gap_n/(days*day_n)  # ratio of qc>1

        # set qc>1 ratio <= 10%
        if gappercent <= 0.1:
            qle_year.append(year)

    # Select years with latent heat qc>1 not exceeding 10%.
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

    # avalable years after screening by sensible and latent qc filter
    available_year = list(set(qle_year) & set(qh_year))    # overlapped year
    available_year.sort()


    # for OzFlux sites, no need VPD screening
    if file[0:2] == 'AU':
        if len(available_year) != 0:
            outsite_year[file[0:6]] = available_year       # store the available years
        continue     # exit cycle and process the next site

    # select year that VPD qc>0 ratio <= 10%
    file = file.replace('Flux.nc', 'Met.nc')
    data = xr.open_dataset(metdir+file)    # read met data
    vpd_year = []                          # store VPD gap-filled > 10% data
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


    all_avai = list(set(available_year) & set(vpd_year))   # overlapped year
    all_avai.sort()
    if len(all_avai) != 0:
        outsite_year[file[0:6]] = all_avai

outsite_year = sorted(outsite_year.items())
print(outsite_year)
