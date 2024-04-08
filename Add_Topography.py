import os
import re
import numpy as np
import pandas as pd
import netCDF4 as nc

############################################################
#
# Add site elevation, slope, and aspect.
# attribute values into the NetCDF file.
#
# Author: Jiahao Shi, 04/2024
#
############################################################

################ read topography feature file ###############
data = pd.read_csv('./creat_nc_read/topography.csv',index_col=0)

sites = data.index.unique()

# PLUMBER2 data path
forname_files = os.listdir('/stu01/shijh21/data/forcingPLUMBER2/met/')


############################################################
#                 loop for each site                       
print('#############################################################################')
print('Loop for each site, add site elevation, slope, and aspect attribute values into the NetCDF file.')
print('#############################################################################')
############################################################
for site in sites:

    print(f'Processing site {site}! Adding topography feature.')

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
    ############### add elevation attribute ##############################
    # create var and set value
    elev = file.createVariable('Elevation', 'f4')
    elev[()] = float(sitedata.loc['elevation'])

    # add attribute info
    elev.long_name = 'Site elevation'
    elev.unit = 'm'
    if 'Flux' in sitedata.loc['elevation_source'] or 'FLUX' in sitedata.loc['elevation_source']:    
        elev.source = sitedata['elevation_source'] + ' (' + sitedata['elevation_web'] + ')'
    else:
        elev.source = sitedata['elevation_source'] + ', ' + sitedata['elevation_web']
    
    ######################################################################
    ################ add slope attribute data if it exists ###############
    if not pd.isna(sitedata['slope']):
        # create var and set value
        slope = file.createVariable('Slope', 'S100')
        slope[()] = sitedata.loc['slope']
        
        # add attribute info
        slope.long_name = 'Site slope'
        slope.content_reads = sitedata.loc['slope']
        slope.description = '''Slope is given as a percentage, or `Flat` in the case of flat land.'''
        if 'Flux' in sitedata.loc['slope_source'] or 'FLUX' in sitedata.loc['slope_source']:    
            slope.source = sitedata['slope_source'] + ' (' + sitedata['slope_web'] + ')'
        else:
            slope.source = sitedata['slope_source'] + ', ' + sitedata['slope_web']

    ######################################################################
    ################ add aspect attribute data if it exists ###############   
    if not pd.isna(sitedata['aspect']):
        # create var and set value
        aspect = file.createVariable('Aspect', 'S100')
        aspect[()] = sitedata.loc['aspect']
        
        # add attribute info
        aspect.long_name = 'Site aspect'
        aspect.content_reads = sitedata.loc['aspect']
        aspect.description = 'Flat: The site is flat and the footprint/exposure is not in any specific direction; N: North; NNE: North-northeast; NE: Northeast; ENE: East-northeast; E: East; ESE: East-southeast; SE: Southeast; SSE: South-southeast; S: South; SSW: South-southwest; SW: Southwest; WSW: West-southwest; W: West; WNW: West-northwest; NW: Northwest; NNW: North-northwest.'
        if 'Flux' in sitedata.loc['aspect_source'] or 'FLUX' in sitedata.loc['aspect_source']:    
            aspect.source = sitedata['aspect_source'] + ' (' + sitedata['aspect_web'] + ')'
        else:
            aspect.source = sitedata['aspect_source'] + ', ' + sitedata['aspect_web']


    file.close()


# END
print('#############################################################################')
print('Processing completed! Please check the files in the siteinfo_out_nc directory.')
print('#############################################################################')






















