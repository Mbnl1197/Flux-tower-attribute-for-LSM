import os
import re
import numpy as np
import pandas as pd
import netCDF4 as nc

############################################################
#
# Add soil depth, the reference measurement height of 
# air temperature and humidity into the NetCDF file.
#
# Author: Jiahao Shi, 06/2024
#
############################################################

################ read feature file ###############
data = pd.read_csv('./creat_nc_read/soil_depth-reference_height.csv',index_col=0)

sites = data.index.unique()

# PLUMBER2 data path
forname_files = os.listdir('D:/data/PLUMBER2/met/')


############################################################
#                 loop for each site                       
print('#############################################################################')
print('Loop for each site, add ssoil depth, the reference measurement height of air temperature and hunidity into the NetCDF file.')
print('#############################################################################')
############################################################
for site in sites:

    print(f'Processing site {site}! Adding feature.')

    # read site attribute information from file created by Creat_nc.py
    metname = re.compile(site + r'.*' + '.nc')
    metfile = [metfile for metfile in forname_files if re.match(metname,metfile)][0]
    if metfile[-11:-7] == 'Flux':
        file = nc.Dataset(f'./siteinfo_out_nc/{site}_OzFlux_Veg_Soil_Topography_ReferenceHeight.nc', 'r+')
    elif metfile[-11:-7] == '2015':
        file = nc.Dataset(f'./siteinfo_out_nc/{site}_FLUXNET2015_Veg_Soil_Topography_ReferenceHeight.nc', 'r+')
    elif metfile[-11:-7] == 'uile':
        file = nc.Dataset(f'./siteinfo_out_nc/{site}_LaThuile_Veg_Soil_Topography_ReferenceHeight.nc', 'r+')
    else:
        print('error!!!!!!!!!!!!')

    # read site topography data
    sitedata = data.loc[site]

    ######################################################################
    ############### add soil depth ##############################

    # create var and set value
    s_d = file.createVariable('Soil_depth', 'f4')
    s_d[()] = float(sitedata.loc['soil_depth'])

    # add attribute info
    s_d.long_name = 'Soil depth'
    s_d.unit = 'cm'
    if 'Flux' in sitedata.loc['soil_depth_source'] or 'FLUX' in sitedata.loc['soil_depth_source']:    
        s_d.source = sitedata['soil_depth_source'] + ' (' + sitedata['soil_depth_web'] + ')'
    else:
        s_d.source = sitedata['soil_depth_source'] + ', ' + sitedata['soil_depth_web']

    ######################################################################
    ################ add slope attribute data if it exists ###############
    # create var and set value
    measurement_height_t = file.createVariable('Reference_height_t', 'f4')
    measurement_height_t[()] = float(sitedata.loc['TA_Height'])
    measurement_height_q = file.createVariable('Reference_height_q', 'f4')
    measurement_height_q[()] = float(sitedata.loc['RH_Height'])

    # add attribute info
        # temperature
    measurement_height_t.long_name = 'Reference measurement height of air temperature or flux'
    measurement_height_t.units     = 'm'
    if 'Flux' in sitedata.loc['Height_source'] or 'FLUX' in sitedata.loc['Height_source']:
        measurement_height_t.source = sitedata.loc['Height_source'] +' (' + sitedata.loc['Height_web'] + ')'
    else:
        measurement_height_t.source = sitedata.loc['Height_source'] +', ' + sitedata.loc['Height_web']
    measurement_height_t.measured_variable = sitedata.loc['measurement_variable_1']
        # humidity
    measurement_height_q.long_name = 'Reference measurement height of air humidity or flux'
    measurement_height_q.units     = 'm'
    if 'Flux' in sitedata.loc['Height_source'] or 'FLUX' in sitedata.loc['Height_source']:
        measurement_height_q.source = sitedata.loc['Height_source'] +' (' + sitedata.loc['Height_web'] + ')'
    else:
        measurement_height_q.source = sitedata.loc['Height_source'] +', ' + sitedata.loc['Height_web']
    measurement_height_q.measured_variable = sitedata.loc['measurement_variable_2']



    file.close()


# END
print('#############################################################################')
print('Processing completed! Please check the files in the siteinfo_out_nc directory.')
print('#############################################################################')





























