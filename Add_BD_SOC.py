#为含有站点土壤容重（BD）和有机碳浓度（SOC）的站点添加相应变量到NC文件
import os
import re
import numpy as np
import pandas as pd
import netCDF4 as nc


################读取BD，SOC属性文件################################
data = pd.read_csv('./bd_soc.csv',index_col=0)
sites = data.index.unique()

#PLUMBER2数据存储路径
forname_files = os.listdir('D:/data/PLUMBER2/select_site_met/')


############################################################
###############循环处理每个站点##############################
############################################################
for site in sites:

    #读取站点属性信息数据（由Creat_nc生成的文件）
    metname = re.compile(site + r'.*' + '.nc')
    metfile = [metfile for metfile in forname_files if re.match(metname,metfile)][0]
    if metfile[-11:-7] == 'Flux':
        file = nc.Dataset(f'./siteinfo_adbdsoc/{site}_OzFlux_Veg_Soil_ReferenceHeight.nc', 'r+')
    elif metfile[-11:-7] == '2015':
        file = nc.Dataset(f'./siteinfo_adbdsoc/{site}_FLUXNET2015_Veg_Soil_ReferenceHeight.nc', 'r+')
    elif metfile[-11:-7] == 'uile':
        file = nc.Dataset(f'./siteinfo_adbdsoc/{site}_LaThuile_Veg_Soil_ReferenceHeight.nc', 'r+')
    else:
        print('error!!!!!!!!!!!!')
    
    #读取站点BD和SOC数据
    sitedata = data.loc[site]

    #如果站点存在SOC属性数据，则添加变量
    if not pd.isna(sitedata['soc'][0]):
        #创建变量并赋值
        soc = file.createVariable('Soil_OC', 'f4', ('soil_layer'))
        soc[0] = sitedata['soc'][0]
        #添加属性信息
        soc.long_name = 'Soil organic carbon concentration'
        soc.units = '%'
        if 'Flux' in sitedata['source'][0]:
            soc.source = sitedata['source'][0] + ' (' + sitedata['web'][0] + ')'
        else:
            soc.source = sitedata['source'][0] + ', ' + sitedata['web'][0]
        soc.layer_1_depth = sitedata['soc'][1]
        #对多层SOC数据进行处理
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

    #如果站点存在BD属性数据，则添加变量
    if not pd.isna(sitedata['bd'][0]):
        #创建变量并赋值
        bd = file.createVariable('Soil_BD', 'f4', ('soil_layer'))
        bd[0] = sitedata['bd'][0]
        #添加属性信息
        bd.long_name = 'Soil bulk density'
        bd.units = 'g cm-3'
        if 'Flux' in sitedata['source'][0]:
            bd.source = sitedata['source'][0] + ' (' + sitedata['web'][0] + ')'
        else:
            bd.source = sitedata['source'][0] + ', ' + sitedata['web'][0]
        bd.layer_1_depth = sitedata['bd'][1]
        #对多层BD数据进行处理
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

    #设置CN-Du2站点SOC参考来源
    if site == 'CN-Du2':
        file = nc.Dataset(f'./siteinfo_adbdsoc/{site}_FLUXNET2015_Veg_Soil_ReferenceHeight.nc', 'r+')
        Soil_OC = file.variables['Soil_OC']
        Soil_OC.setncattr('source','ChinaFlux (http://www.chinaflux.org/)')
        file.close()




