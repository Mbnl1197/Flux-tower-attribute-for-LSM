import os
import re
import numpy as np
import pandas as pd
import netCDF4 as nc

############################################################
#
# Add site soil Bulk Density (BD) and Soil Organic Carbon (SOC)
# attribute values into the NetCDF file.
#
# Author: Jiahao Shi, 10/2023
#
############################################################


################ read BD，SOC attribute file ###############
data = pd.read_csv('./creat_nc_read/bd_soc.csv',index_col=0)

sites = data.index.unique()

# PLUMBER2 data path
forname_files = os.listdir('/stu01/shijh21/data/forcingPLUMBER2/met/')


############################################################
#                 loop for each site
print('#############################################################################')
print('Loop for each site, add site soil bulk density (BD) and soil organic carbon (SOC) attribute values into the NetCDF file.')
print('#############################################################################')
############################################################
for site in sites:

    print(f'Processing site {site}! Adding BD and SOC.')

    # read site attribute information from file created by Creat_nc.py
    metname = re.compile(site + r'.*' + '.nc')
    metfile = [metfile for metfile in forname_files if re.match(metname,metfile)][0]
    if metfile[-11:-7] == 'Flux':
        file = nc.Dataset(f'./siteinfo_out_nc/{site}_OzFlux_Veg_Soil_ReferenceHeight.nc', 'r+')
    elif metfile[-11:-7] == '2015':
        file = nc.Dataset(f'./siteinfo_out_nc/{site}_FLUXNET2015_Veg_Soil_ReferenceHeight.nc', 'r+')
    elif metfile[-11:-7] == 'uile':
        file = nc.Dataset(f'./siteinfo_out_nc/{site}_LaThuile_Veg_Soil_ReferenceHeight.nc', 'r+')
    else:
        print('error!!!!!!!!!!!!')

    # read site BD and SOC data
    sitedata = data.loc[site]

    # add SOC attribute data if it exists
    if not pd.isna(sitedata['soc'][0]):

        # create var and set value
        soc = file.createVariable('Soil_OC', 'f4', ('soil_layer'))
        soc[0] = sitedata['soc'][0]

        # add attribute info
        soc.long_name = 'Soil organic carbon concentration'
        soc.units = '%'
        if 'Flux' in sitedata['source'][0]:
            soc.source = sitedata['source'][0] + ' (' + sitedata['web'][0] + ')'
        else:
            soc.source = sitedata['source'][0] + ', ' + sitedata['web'][0]
        soc.layer_1_depth = sitedata['soc'][1]

        # process SOC data for multiple-layer cases
        if not pd.isna(sitedata['soc1'][0]):
            soc[1] = sitedata['soc1'][0]
            soc.layer_2_depth = sitedata['soc1'][1]
        else:
            soc[1] = np.nan
            soc.layer_2_depth = 'Na'
        if not pd.isna(sitedata['soc2'][0]):
            soc[2] = sitedata['soc2'][0]
            soc.layer_3_depth = sitedata['soc2'][1]
        else:
            soc[2] = np.nan
            soc.layer_3_depth = 'Na'
        if not pd.isna(sitedata['soc3'][0]):
            soc[3] = sitedata['soc3'][0]
            soc.layer_4_depth = sitedata['soc3'][1]
        else:
            soc[3] = np.nan
            soc.layer_4_depth = 'Na'

    # add BD attribute data if it exists
    if not pd.isna(sitedata['bd'][0]):

        # create var and set value
        bd = file.createVariable('Soil_BD', 'f4', ('soil_layer'))
        bd[0] = sitedata['bd'][0]

        # add attribute info
        bd.long_name = 'Soil bulk density'
        bd.units = 'g cm-3'
        if 'Flux' in sitedata['source'][0]:
            bd.source = sitedata['source'][0] + ' (' + sitedata['web'][0] + ')'
        else:
            bd.source = sitedata['source'][0] + ', ' + sitedata['web'][0]
        bd.layer_1_depth = sitedata['bd'][1]

        # process BD data for multiple-layer cases
        if not pd.isna(sitedata['bd1'][0]):
            bd[1] = sitedata['bd1'][0]
            bd.layer_2_depth = sitedata['bd1'][1]
        else:
            bd[1] = np.nan
            bd.layer_2_depth = 'Na'
        if not pd.isna(sitedata['bd2'][0]):
            bd[2] = sitedata['bd2'][0]
            bd.layer_3_depth = sitedata['bd2'][1]
        else:
            bd[2] = np.nan
            bd.layer_3_depth = 'Na'
        if not pd.isna(sitedata['bd3'][0]):
            bd[3] = sitedata['bd3'][0]
            bd.layer_4_depth = sitedata['bd3'][1]
        else:
            bd[3] = np.nan
            bd.layer_4_depth = 'Na'
    file.close()

    # set SOC data reference source for site CN-Du2
    if site == 'CN-Du2':
        file = nc.Dataset(f'./siteinfo_out_nc/{site}_FLUXNET2015_Veg_Soil_ReferenceHeight.nc', 'r+')
        Soil_OC = file.variables['Soil_OC']
        Soil_OC.setncattr('source','ChinaFlux (http://www.chinaflux.org/)')
        file.close()

# END
print('#############################################################################')
print('Processing completed! Please check the files in the siteinfo_out_nc directory.')
print('#############################################################################')