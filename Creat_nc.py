import os
import re
import numpy as np
import pandas as pd
import netCDF4 as nc
import datetime as dt

############################################################
#读取筛选到的站点年份、收集到的站点属性数据，用于写入NC文件
############################################################

data_year = pd.read_csv('./select_year.csv',index_col=0)
data_attr = pd.read_excel('./siteinfo.xlsx',index_col=0)
data_qc = pd.read_excel('./lai_qc.xlsx',index_col=0)
igbpdata = pd.read_csv('./igbp.csv',index_col=0)
sites = data_attr.index.unique()

#读取PLUMBER2站点数据，用于获取站点所述数据集（OzFlux, FLUXNET2015, LaThuile）。
forname_files = os.listdir('D:/data/PLUMBER2/select_site_met/')


############################################################
###############################循环处理每个站点##############
############################################################
for site in sites:
##############################读取站点相关属性信息#############################
    #读取IGBP信息
    igbp_short = igbpdata.loc[site].IGBP_short
    igbp_long = igbpdata.loc[site].IGBP_long
    igbp_index = igbpdata.loc[site].IGBP_index
    #读取站点属性信息（PFT, soil texture, LAI, , canopy height, reference height, lai, lon）
    sitedata = data_attr.loc[site]
    sitedata = sitedata.set_index('variable')
    #读取站点筛选的年份
    siteyear = []
    siteyears = data_year.loc[site][0]
    siteyears = siteyears.split(',')
    for year in siteyears:
        siteyear.append(int(year))
    siteyear.sort()
##############################创建NC文件#############################
    # 根据站点所属数据集设置文件名
    metname = re.compile(site + r'.*' + '.nc')
    metfile = [metfile for metfile in forname_files if re.match(metname,metfile)][0]
    if metfile[-11:-7] == 'Flux':
        newfile = nc.Dataset(f'./siteinfo_out1/{site}_OzFlux_Veg_Soil_ReferenceHeight.nc', 'w', format='NETCDF4')
    elif metfile[-11:-7] == '2015':
        newfile = nc.Dataset(f'./siteinfo_out1/{site}_FLUXNET2015_Veg_Soil_ReferenceHeight.nc', 'w', format='NETCDF4')
    elif metfile[-11:-7] == 'uile':
        newfile = nc.Dataset(f'./siteinfo_out1/{site}_LaThuile_Veg_Soil_ReferenceHeight.nc', 'w', format='NETCDF4')
    else:
        print('error!!!!!!!!!!!!')
    #创建维度
    long = newfile.createDimension('longitude', size=1)           #经度
    lati = newfile.createDimension('latitude', size=1)            #维度
    PFT = newfile.createDimension('pft', size=16)                 #PFT类型 
    particle = newfile.createDimension('particle_size', size=3)   #土壤颗粒大小
    soil_layer = newfile.createDimension('soil_layer', size=4)    #土壤层
    year = newfile.createDimension('year', 21)                    #筛选后数据在1997-2018，共21年
    #创建变量
    lon = newfile.createVariable('longitude', 'f4', dimensions='longitude')
    lat = newfile.createVariable('latitude', 'f4', dimensions='latitude')
    PCT_PFT = newfile.createVariable('PCT_PFT','f4',('pft'))
    soil_tex = newfile.createVariable('Soil_TEX','f4',('particle_size','soil_layer'))
    LAI_Max = newfile.createVariable('LAI_Max','f4')
    canopy_height = newfile.createVariable('Canopy_height','f4',)
    measurement_height = newfile.createVariable('Reference_height','f4')
    pft = newfile.createVariable('pft','i4',('pft'))
    particle_size = newfile.createVariable('particle_size','i4',('particle_size'))
    soil_layer = newfile.createVariable('soil_layer','i4',('soil_layer'))
    year_var = newfile.createVariable('year', np.int32, ('year',))
    IGBP_index = newfile.createVariable('IGBP','i4')
    selected_year = newfile.createVariable('year_qc', np.int32, ('year',))


#############################为变量赋值############################################
    pft[:] = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]    #15种PFT类型
    particle_size[:] = [1,2,3]                          #3种土壤颗粒，sand/silt/clay
    soil_layer[:] = [1, 2, 3, 4]                        #4个土壤层深度
    year_var[:] = np.arange(1997, 2018)                 #数据集所覆盖的21个年份
    IGBP_index[()] = igbp_index                         #为IGBP索引赋值
    year_range = np.arange(1997,2018)                   #为挑选的年份selected_year赋值
    siteyear = np.array(siteyear)
    b = np.in1d(year_range,siteyear)
    b = b.astype(int)
    selected_year[:] = b
    ###########为LAI最大值赋值
    LAI_Max[()] = float(sitedata.loc['LAI'][0])         #LAI最大值赋值
    if not pd.isna(sitedata.loc['LAI'][3]):             #取出特点年份的LAI最大值
        year1 = sitedata.loc['LAI_year_range'][3]
        lai1 = sitedata.loc['LAI'][3]
        laiyear1 = str(lai1) + '_' + str(year1)
    else:
        laiyear1 = 'Na'
    if not pd.isna(sitedata.loc['LAI'][4]):
        year2 = sitedata.loc['LAI_year_range'][4]
        lai2 = sitedata.loc['LAI'][4]
        laiyear2 = str(lai2) + '_' + str(year2)
    else:
        laiyear2 = 'Na'
    if not pd.isna(sitedata.loc['LAI'][5]):
        year3 = sitedata.loc['LAI_year_range'][5]
        lai3 = sitedata.loc['LAI'][5]
        laiyear3 = str(lai3) + '_' + str(year3)
    else:
        laiyear3 = 'Na'

    if pd.isna(sitedata.loc['LAI'][3]):                 #设置特定年份LAI最大值
        LAI_Max.LAI_Max_year = 'Na'
    if laiyear1 != 'Na' and laiyear2 == 'Na':
        LAI_Max.LAI_Max_year = laiyear1 
    if (laiyear1 != 'Na' and laiyear2 != 'Na') and laiyear3 == 'Na':
        LAI_Max.LAI_Max_year = laiyear1 + '; ' + laiyear2
    if laiyear1 != 'Na' and laiyear2 != 'Na' and laiyear3 != 'Na':
        LAI_Max.LAI_Max_year = laiyear1 + '; ' + laiyear2 + '; ' + laiyear3


    ###########为不同土壤层的土壤质地赋值
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


    #为经纬度赋值
    lat[:] = float(sitedata.loc['lat'][0])
    lon[:] = float(sitedata.loc['lon'][0])
    #为PFT比例赋值
    PCT_PFT[:] = sitedata.value[0:16].values
    #为冠层高度赋值
    canopy_height[()] = float(sitedata.loc['canopy height'][0])
    #为测量高度赋值
    measurement_height[()] = float(sitedata.loc['reference height'][0])


    #添加全局属性信息
    newfile.site_name = site
    newfile.QC_flag_descriptions = '0: Data collected from site-related materials; 1: Set for LAI only. Data collected from site-related materials, but the reference source only gives a site LAI value without attribute descriptions such as observation period; 2: Filled using relevant global data products.'
    newfile.title = 'Flux tower attribute data sets for land surface and climate modelling'
    newfile.contact = 'Jiahao,Shi (shijh26@mail2.sysu.edu.cn)'
    newfile.institution = 'Land-Atmosphere Interaction Research Group at Sun Yat-sen University (http://globalchange.bnu.edu.cn)'
    newfile.Creation_date = dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


    #添加变量属性信息
    lon.long_name = 'Longitude'
    lon.units = 'degrees east'
    lat.long_name = 'Latitude'
    lat.units = 'degrees north'

    pft.long_name = 'Plant functional types'
    pft.pft_index = '0:bare soil; 1:Needleleaf evergreen tree,temperature; 2:Needleleaf evergreen tree,boreal; 3:Needleleaf deciduous tree; 4:Broadleaf evergreen tree,tropical; 5:Broadleaf evergreen tree,temperate; 6:Broadleaf deciduous tree,tropical; 7:Broadleaf deciduous tree,temperate; 8:Broadleaf deciduous tree,boreal; 9:Broadleaf evergreen shrub,temperate; 10:Broadleaf deciduous shrub,temperate; 11:Broadleaf deciduous shrub,boreal; 12:C3 grass,arctic; 13:C3 grass; 14:C4 grass; 15:Crop'
    particle_size.long_name = 'Size of soil particles (sand/silt/clay)'
    particle_size.particle_index = 'particle_size_1:sand; particle_size_2:silt; particle_size_3:clay'
    
    soil_layer.long_name = 'Soil layer'
    year_var.long_name = 'Year range of all sites'
    year_var.units = 'years'
    
    PCT_PFT.long_name = 'Percent PFT(Plant Functional Type) cover'
    PCT_PFT.units = '%'
    if 'Flux' in sitedata.loc['pft0'][1] or 'FLUX' in sitedata.loc['pft0'][1]:
        PCT_PFT.source = sitedata.loc['pft0'][1] +' (' + sitedata.loc['pft0'][2] + ')'
    else:
        PCT_PFT.source = sitedata.loc['pft0'][1] +', ' + sitedata.loc['pft0'][2]

    if PCT_PFT.source[0:9] == 'Satellite':
        PCT_PFT.qc = 2
    else:
        PCT_PFT.qc = 0

    soil_tex.long_name = 'Soil texture (sand/silt/clay)'
    soil_tex.units = '%'
    if 'Flux' in sitedata.loc['sand'][1] or 'FLUX' in sitedata.loc['sand'][1]:
        soil_tex.source = sitedata.loc['sand'][1] +' (' + sitedata.loc['sand'][2] + ')'
    else:
        soil_tex.source = sitedata.loc['sand'][1] +', ' + sitedata.loc['sand'][2]

    if soil_tex.source[1:10] == 'Shangguan':
        soil_tex.qc = 2
    else:
        soil_tex.qc = 0
    
    LAI_Max.long_name = 'Maximum leaf area index '
    LAI_Max.units = 'm2/m2'
    if 'Flux' in sitedata.loc['LAI'][1] or 'FLUX' in sitedata.loc['LAI'][1]:
        LAI_Max.source = sitedata.loc['LAI'][1] +' (' + sitedata.loc['LAI'][2] + ')'
    else:
        LAI_Max.source = sitedata.loc['LAI'][1] +', ' + sitedata.loc['LAI'][2]
    LAI_Max.year_range = sitedata.loc['LAI_year_range'][0]
    LAI_Max.qc = int(data_qc.loc[site].qc)

    canopy_height.long_name = 'Canopy height'
    canopy_height.units = 'm'
    if 'Flux' in sitedata.loc['canopy height'][1] or 'FLUX' in sitedata.loc['canopy height'][1]:
        canopy_height.source = sitedata.loc['canopy height'][1] +' (' + sitedata.loc['canopy height'][2] +')'
    else:
        canopy_height.source = sitedata.loc['canopy height'][1] +', ' + sitedata.loc['canopy height'][2]

    measurement_height.long_name = 'Reference measurement height of wind speed or flux'
    measurement_height.units = 'm'
    if 'Flux' in sitedata.loc['reference height'][1] or 'FLUX' in sitedata.loc['reference height'][1]:
        measurement_height.source = sitedata.loc['reference height'][1] +' (' + sitedata.loc['reference height'][2] + ')'
    else:
        measurement_height.source = sitedata.loc['reference height'][1] +', ' + sitedata.loc['reference height'][2]
    measurement_height.measured_variable = sitedata.loc['reference height'][3]

    selected_year.long_name = 'Selected years of high quality data'
    selected_year.description = 'The selected high-quality year is represented by the value 1'

    IGBP_index.long_name = 'IGBP index number'
    IGBP_index.IGBP_veg_long = igbp_long
    IGBP_index.IGBP_veg_short = igbp_short

    ## close file
    newfile.close()

#################对有两处引用的站点进行规范##########################################################
    if site in ['BE-Vie', 'US-MMS', 'US-SRG', 'US-SRM', 'US-Var']:
        if metfile[-11:-7] == 'Flux':
            file = nc.Dataset(f'./siteinfo_out1/{site}_OzFlux_Veg_Soil_ReferenceHeight.nc', 'r+')
        elif metfile[-11:-7] == '2015':
            file = nc.Dataset(f'./siteinfo_out1/{site}_FLUXNET2015_Veg_Soil_ReferenceHeight.nc', 'r+')
        elif metfile[-11:-7] == 'uile':
            file = nc.Dataset(f'./siteinfo_out1/{site}_LaThuile_Veg_Soil_ReferenceHeight.nc', 'r+')
        else:
            print('error!!!!!!!!!!!!')

        PCT_PFT = file.variables['PCT_PFT']
        if site in [ 'US-MMS', 'US-SRG']:
            PCT_PFT.setncattr('source','[1] ' + sitedata.loc['pft0'][1] +', ' + sitedata.loc['pft0'][2] + '; [2] ' + sitedata.loc['pft0'][3] +' (' + sitedata.loc['pft0'][4] + ')')
        elif site in ['BE-Vie', 'US-SRM', 'US-Var']:
            PCT_PFT.setncattr('source','[1] ' + sitedata.loc['pft0'][1] +', ' + sitedata.loc['pft0'][2] + '; [2] ' + sitedata.loc['pft0'][3] +', ' + sitedata.loc['pft0'][4])
        else:
            print('error!!!!!!!!!!!')

        file.close()






