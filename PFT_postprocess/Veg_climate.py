##########################################
#用于读取koppen数据，并根据Poulter et al., (2011)对植被进行气候分类
###############################################33
from osgeo import gdal, osr
import numpy as np
import pandas as pd


def get_tif_info(tif_path):
    if tif_path.endswith('.tif') or tif_path.endswith('.TIF'):
        dataset = gdal.Open(tif_path)
        pcs = osr.SpatialReference()
        pcs.ImportFromWkt(dataset.GetProjection())
        gcs = pcs.CloneGeogCS()
        extend = dataset.GetGeoTransform()
        # im_width = dataset.RasterXSize #栅格矩阵的列数
        # im_height = dataset.RasterYSize #栅格矩阵的行数
        shape = (dataset.RasterYSize, dataset.RasterXSize)
    else:
        raise "Unsupported file format"

    img = dataset.GetRasterBand(1).ReadAsArray()  # (height, width)
    # img(ndarray), gdal数据集、地理空间坐标系、投影坐标系、栅格影像大小
    return img, dataset, gcs, pcs, extend, shape
def longlat_to_xy(gcs, pcs, lon, lat):
    ct = osr.CoordinateTransformation(gcs, pcs)
    coordinates = ct.TransformPoint(lon, lat)
    return coordinates[0], coordinates[1], coordinates[2]

def xy_to_rowcol(extend, x, y):
    a = np.array([[extend[1], extend[2]], [extend[4], extend[5]]])
    b = np.array([x - extend[0], y - extend[3]])
    row_col = np.linalg.solve(a, b)
    row = int(np.floor(row_col[1]))
    col = int(np.floor(row_col[0]))
    return row, col

##################读取影像指定坐标（rowcol/lonlat/xy）的值
def get_value_by_coordinates(tif_pah, coordinates, coordinate_type='rowcol'):
    img, dataset, gcs, pcs, extend, shape = get_tif_info(tif_pah)

    if coordinate_type == 'rowcol':
        value = img[coordinates[0], coordinates[1]]
    elif coordinate_type == 'lonlat':
        x, y, _ = longlat_to_xy(gcs, pcs, coordinates[0], coordinates[1])
        row, col = xy_to_rowcol(extend, x, y)
        value = img[row, col]
    elif coordinate_type == 'xy':
        row, col = xy_to_rowcol(extend, coordinates[0], coordinates[1])
        value = img[row, col]
    else:
        raise 'coordinated_type error'
    return value


file = './Beck_KG_V1/Beck_KG_V1_present_0p0083.tif'  # koppen气候分区文件
data = pd.read_csv('./site_latlon.csv',index_col=0)  #读取站点经纬度
sites = data.index
climate = pd.DataFrame(columns=['climate','veg_cli'])#定义输出的dataframe，用于存储站点气候类型


##############################循环处理每个站点
for site in sites:
    lat = np.round(data.loc[site].lat,4)
    lon = np.round(data.loc[site].lon,4)
    value = get_value_by_coordinates(file, [lon, lat], coordinate_type='lonlat')

    #根据Poulter et al., (2011)，使用koppen气候数据对植被进行气候分类
    if value <= 4 or value == 6:
        veg_cli = 'tropical'
    elif value <= 16:
        veg_cli = 'temperature'
    elif value >= 17:
        veg_cli = 'boreal'

    #存储站点koppen气候类型序号以及对应植被气候类型
    climate.loc[site] = [value,veg_cli]








