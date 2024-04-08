import os
import re
import numpy as np
import pandas as pd
import netCDF4 as nc
import datetime as dt

############################################################
#
# Read the filtered site years and collected site attribute
# data for writing into an NC file
#
# Author: Jiahao Shi, 10/06/2023
#
############################################################

data_year = pd.read_csv('./creat_nc_read/select_year.csv',index_col = 0)
data_attr = pd.read_csv('./creat_nc_read/siteinfo.csv',index_col = 0)
data_qc   = pd.read_csv('./creat_nc_read/lai_qc.csv',index_col = 0)
igbpdata  = pd.read_csv('./creat_nc_read/igbp.csv',index_col = 0)
sites     = data_attr.index.unique()

# Read PLUMBER2 site data to obtain the dataset described by the site.
# (OzFlux, FLUXNET2015, LaThuile)
forname_files = os.listdir('/stu01/shijh21/data/forcingPLUMBER2/met/')

# Creating Output Folder
os.system('mkdir -p ./siteinfo_out_nc')

############################################################
# process each site (loop)
print('#############################################################################')
print('Pocess each site (loop)!!!')
print('#############################################################################')
############################################################

######## Read site-related attribute information ###########
for site in sites:

    print(f'Processing site {site}!')

    # read IGBP type
    igbp_short = igbpdata.loc[site].IGBP_short
    igbp_long  = igbpdata.loc[site].IGBP_long
    igbp_index = igbpdata.loc[site].IGBP_index

    # read site PFT, soil texture, maximum LAI, canopy height, reference
    # height, lai, lon
    sitedata = data_attr.loc[site]
    sitedata = sitedata.set_index('variable')

    # read the seleted years
    siteyear  = []
    siteyears = data_year.loc[site][0]
    siteyears = siteyears.split(',')
    for year in siteyears:
        siteyear.append(int(year))
    siteyear.sort()

################create NC file ############################

    # set filename
    metname = re.compile(site + r'.*' + '.nc')
    metfile = [metfile for metfile in forname_files if re.match(metname,metfile)][0]
    if metfile[-11:-7] == 'Flux':
        newfile = nc.Dataset(f'./siteinfo_out_nc/{site}_OzFlux_Veg_Soil_Topography_ReferenceHeight.nc', 'w', format='NETCDF4')
    elif metfile[-11:-7] == '2015':
        newfile = nc.Dataset(f'./siteinfo_out_nc/{site}_FLUXNET2015_Veg_Soil_Topography_ReferenceHeight.nc', 'w', format='NETCDF4')
    elif metfile[-11:-7] == 'uile':
        newfile = nc.Dataset(f'./siteinfo_out_nc/{site}_LaThuile_Veg_Soil_Topography_ReferenceHeight.nc', 'w', format='NETCDF4')
    else:
        print('error!!!!!!!!!!!!')

    # create dimensions
    long       = newfile.createDimension('longitude', size=1)     # longitute
    lati       = newfile.createDimension('latitude', size=1)      # latitude
    PFT        = newfile.createDimension('pft', size=16)          # PFT
    particle   = newfile.createDimension('particle_size', size=3) # soil particle size
    soil_layer = newfile.createDimension('soil_layer', size=4)    # soil layer
    year       = newfile.createDimension('year', 21)              # data year range 1997-2018，total 21

    # create variables
    lon                = newfile.createVariable('longitude', 'f4', dimensions = 'longitude')
    lat                = newfile.createVariable('latitude', 'f4', dimensions = 'latitude')
    PCT_PFT            = newfile.createVariable('PCT_PFT','f4',('pft'))
    soil_tex           = newfile.createVariable('Soil_TEX','f4',('particle_size','soil_layer'))
    LAI_Max            = newfile.createVariable('LAI_Max','f4')
    canopy_height      = newfile.createVariable('Canopy_height','f4',)
    measurement_height = newfile.createVariable('Reference_height','f4')
    pft                = newfile.createVariable('pft','i4',('pft'))
    particle_size      = newfile.createVariable('particle_size','i4',('particle_size'))
    soil_layer         = newfile.createVariable('soil_layer','i4',('soil_layer'))
    year_var           = newfile.createVariable('year', np.int32, ('year',))
    IGBP_index         = newfile.createVariable('IGBP','i4')
    selected_year      = newfile.createVariable('year_qc', np.int32, ('year',))


############################# set values ##########################################
    pft[:] = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]    # 15 PFT types
    particle_size[:] = [1,2,3]                          # 3 particle types, sand/silt/clay
    soil_layer[:]    = [1, 2, 3, 4]                     # 4 soil layer depth
    year_var[:]      = np.arange(1997, 2018)            # data year range
    IGBP_index[()]   = igbp_index                       # IGBP index
    year_range       = np.arange(1997,2018)             # year that selected
    siteyear         = np.array(siteyear)
    b                = np.in1d(year_range,siteyear)
    b                = b.astype(int)
    selected_year[:] = b

    # set for maximum LAI
    LAI_Max[()] = float(sitedata.loc['LAI'][0])         # max LAI value
    if not pd.isna(sitedata.loc['LAI'][3]):             # get the year of max LAI value and its value
        year1    = sitedata.loc['LAI_year_range'][3]
        lai1     = sitedata.loc['LAI'][3]
        laiyear1 = str(year1) + ' (' + str(lai1) + ')'
    else:
        laiyear1 = 'Na'
    if not pd.isna(sitedata.loc['LAI'][4]):
        year2    = sitedata.loc['LAI_year_range'][4]
        lai2     = sitedata.loc['LAI'][4]
        laiyear2 = str(year2) + ' (' + str(lai2) + ')'
    else:
        laiyear2 = 'Na'
    if not pd.isna(sitedata.loc['LAI'][5]):
        year3    = sitedata.loc['LAI_year_range'][5]
        lai3     = sitedata.loc['LAI'][5]
        laiyear3 = str(year3) + ' (' + str(lai3) + ')'
    else:
        laiyear3 = 'Na'

    if pd.isna(sitedata.loc['LAI'][3]):                 # set the year of max LAI value and its value
        LAI_Max.LAI_Max_year = 'Na'
    if laiyear1 != 'Na' and laiyear2 == 'Na':
        LAI_Max.LAI_Max_year = laiyear1
    if (laiyear1 != 'Na' and laiyear2 != 'Na') and laiyear3 == 'Na':
        LAI_Max.LAI_Max_year = laiyear1 + ', ' + laiyear2
    if laiyear1 != 'Na' and laiyear2 != 'Na' and laiyear3 != 'Na':
        LAI_Max.LAI_Max_year = laiyear1 + ', ' + laiyear2 + ', ' + laiyear3


    # set soil texture for each layers
    soil_tex[:,0] = sitedata.value[16:19].values
    soil_tex.layer_1_depth = sitedata.loc['tex_depth'][0]
    if not pd.isna(sitedata.loc['sand'][3]):
        soil_tex[:,1] = sitedata.sup_info1[16:19].values
        soil_tex.layer_2_depth = sitedata.loc['tex_depth'][3]
    else:
        soil_tex[:,1] = np.nan
        soil_tex.layer_2_depth = 'Na'
    if not pd.isna(sitedata.loc['sand'][4]):
        soil_tex[:,2] = sitedata.sup_info2[16:19].values
        soil_tex.layer_3_depth = sitedata.loc['tex_depth'][4]
    else:
        soil_tex[:,2] = np.nan
        soil_tex.layer_3_depth = 'Na'
    if not pd.isna(sitedata.loc['sand'][5]):
        soil_tex[:,3] = sitedata.sup_info3[16:19].values
        soil_tex.layer_4_depth = sitedata.loc['tex_depth'][5]
    else:
        soil_tex[:,3] = np.nan
        soil_tex.layer_4_depth = 'Na'


    # set lat, lon
    lat[:] = float(sitedata.loc['lat'][0])
    lon[:] = float(sitedata.loc['lon'][0])

    # set PCT_PFT
    PCT_PFT[:] = sitedata.value[0:16].values

    # set canopy height
    canopy_height[()] = float(sitedata.loc['canopy height'][0])

    # set reference height
    measurement_height[()] = float(sitedata.loc['reference height'][0])


    # add global attributes
    newfile.site_name = site
    newfile.QC_flag_descriptions = '0: Data collected from site-related materials; 1: Set for LAI_Max only. Data collected from site-related materials, but the reference source only gives a site LAI value without attribute descriptions such as observation period; 2: Filled using relevant global data products.'
    newfile.title = 'Flux tower site attribute data sets for land surface modeling'
    newfile.contact = 'Jiahao Shi (shijh26@mail2.sysu.edu.cn), Hua Yuan (yuanh25@mail.sysu.edu.cn)'
    newfile.institution = 'Land-Atmosphere Interaction Research Group at Sun Yat-sen University (http://globalchange.bnu.edu.cn)'
    # newfile.Creation_date = dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    newfile.Creation_date = '2023-12-25'


    # add vars' attributes
    lon.long_name = 'Longitude'
    lon.units     = 'degrees east'
    lat.long_name = 'Latitude'
    lat.units     = 'degrees north'

    pft.long_name = 'Plant Functional Types'
    pft.pft_index = '0:bare soil; 1:Needleleaf evergreen tree,temperature; 2:Needleleaf evergreen tree,boreal; 3:Needleleaf deciduous tree; 4:Broadleaf evergreen tree,tropical; 5:Broadleaf evergreen tree,temperate; 6:Broadleaf deciduous tree,tropical; 7:Broadleaf deciduous tree,temperate; 8:Broadleaf deciduous tree,boreal; 9:Broadleaf evergreen shrub,temperate; 10:Broadleaf deciduous shrub,temperate; 11:Broadleaf deciduous shrub,boreal; 12:C3 grass,arctic; 13:C3 grass; 14:C4 grass; 15:Crop'
    particle_size.long_name      = 'Size of soil particles (sand/silt/clay)'
    particle_size.particle_index = 'particle_size_1:sand; particle_size_2:silt; particle_size_3:clay'

    soil_layer.long_name = 'Soil layer'
    year_var.long_name   = 'Year range of all sites'
    year_var.units       = 'years'

    # for PCT_PFT
    PCT_PFT.long_name = 'Percent PFT (Plant Functional Type) cover'
    PCT_PFT.units     = '%'
    if 'Flux' in sitedata.loc['pft0'][1] or 'FLUX' in sitedata.loc['pft0'][1]:
        PCT_PFT.source = sitedata.loc['pft0'][1] +' (' + sitedata.loc['pft0'][2] + ')'
    else:
        PCT_PFT.source = sitedata.loc['pft0'][1] +', ' + sitedata.loc['pft0'][2]

    if PCT_PFT.source[0:9] == 'Satellite':
        PCT_PFT.qc = 2
    else:
        PCT_PFT.qc = 0

    # for soil texture
    soil_tex.long_name = 'Soil texture (sand/silt/clay)'
    soil_tex.units     = '%'
    if 'Flux' in sitedata.loc['sand'][1] or 'FLUX' in sitedata.loc['sand'][1]:
        soil_tex.source = sitedata.loc['sand'][1] +' (' + sitedata.loc['sand'][2] + ')'
    else:
        soil_tex.source = sitedata.loc['sand'][1] +', ' + sitedata.loc['sand'][2]

    if soil_tex.source[1:10] == 'Shangguan':
        soil_tex.qc = 2
    else:
        soil_tex.qc = 0

    # for maximum LAI
    LAI_Max.long_name = 'Maximum leaf area index'
    LAI_Max.units     = 'm^2/m^2'
    if 'Flux' in sitedata.loc['LAI'][1] or 'FLUX' in sitedata.loc['LAI'][1]:
        LAI_Max.source = sitedata.loc['LAI'][1] +' (' + sitedata.loc['LAI'][2] + ')'
    else:
        LAI_Max.source = sitedata.loc['LAI'][1] +', ' + sitedata.loc['LAI'][2]
    LAI_Max.year_range = sitedata.loc['LAI_year_range'][0]
    LAI_Max.qc = int(data_qc.loc[site].qc)

    # for canopy height
    canopy_height.long_name = 'Canopy height'
    canopy_height.units     = 'm'
    if 'Flux' in sitedata.loc['canopy height'][1] or 'FLUX' in sitedata.loc['canopy height'][1]:
        canopy_height.source = sitedata.loc['canopy height'][1] +' (' + sitedata.loc['canopy height'][2] +')'
    else:
        canopy_height.source = sitedata.loc['canopy height'][1] +', ' + sitedata.loc['canopy height'][2]

    # for measurement height
    measurement_height.long_name = 'Reference measurement height of wind speed or flux'
    measurement_height.units     = 'm'
    if 'Flux' in sitedata.loc['reference height'][1] or 'FLUX' in sitedata.loc['reference height'][1]:
        measurement_height.source = sitedata.loc['reference height'][1] +' (' + sitedata.loc['reference height'][2] + ')'
    else:
        measurement_height.source = sitedata.loc['reference height'][1] +', ' + sitedata.loc['reference height'][2]
    measurement_height.measured_variable = sitedata.loc['reference height'][3]

    # for selected year
    selected_year.long_name   = 'Selected years of high quality data'
    selected_year.description = 'The selected high-quality year is represented by the value 1'

    # for IGBP type
    IGBP_index.long_name      = 'IGBP index number'
    IGBP_index.IGBP_veg_long  = igbp_long
    IGBP_index.IGBP_veg_short = igbp_short

    ## close file
    newfile.close()

################# set attributes for several specific sites ##################
    if site in ['BE-Vie', 'US-MMS', 'US-SRG', 'US-SRM', 'US-Var', 'US-FPe']:
        if metfile[-11:-7] == 'Flux':
            file = nc.Dataset(f'./siteinfo_out_nc/{site}_OzFlux_Veg_Soil_Topography_ReferenceHeight.nc', 'r+')
        elif metfile[-11:-7] == '2015':
            file = nc.Dataset(f'./siteinfo_out_nc/{site}_FLUXNET2015_Veg_Soil_Topography_ReferenceHeight.nc', 'r+')
        elif metfile[-11:-7] == 'uile':
            file = nc.Dataset(f'./siteinfo_out_nc/{site}_LaThuile_Veg_Soil_Topography_ReferenceHeight.nc', 'r+')
        else:
            print('error!!!!!!!!!!!!')

        PCT_PFT = file.variables['PCT_PFT']
        if site in [ 'US-MMS', 'US-SRG']:
            PCT_PFT.setncattr('source','[1] ' + sitedata.loc['pft0'][1] +', ' + sitedata.loc['pft0'][2] + '; [2] ' + sitedata.loc['pft0'][3] +' (' + sitedata.loc['pft0'][4] + ')')
        elif site in ['BE-Vie', 'US-SRM', 'US-Var']:
            PCT_PFT.setncattr('source','[1] ' + sitedata.loc['pft0'][1] +', ' + sitedata.loc['pft0'][2] + '; [2] ' + sitedata.loc['pft0'][3] +', ' + sitedata.loc['pft0'][4])
        elif site == 'US-FPe':
            PCT_PFT.setncattr('source','[1] ' + sitedata.loc['pft0'][1] +', ' + sitedata.loc['pft0'][2] + '; [2] ' + sitedata.loc['pft0'][3] +', ' + sitedata.loc['pft0'][4])
        else:
            print('error!!!!!!!!!!!')

        file.close()

# END
print('#############################################################################')
print('Processing completed! Please check the files in the siteinfo_out_nc directory.')
print('#############################################################################')
