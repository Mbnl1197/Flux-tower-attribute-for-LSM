######################################################################
# Read Köppen data and categorize the climate type of vegetation 
# according to Poulter et al., (2011). 
######################################################################
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
        # im_width = dataset.RasterXSize 
        # im_height = dataset.RasterYSize 
        shape = (dataset.RasterYSize, dataset.RasterXSize)
    else:
        raise "Unsupported file format"

    img = dataset.GetRasterBand(1).ReadAsArray()  # (height, width)
    # img(ndarray), gdal dataset, geospatial coordinate system, projected coordinate system, raster image size
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

####################################################################
# Retrieves the value of the specified coordinates (rowcol/lonlat/xy) of the image
#######################################################################
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

# Köppen_climate file
file = './Beck_KG_V1/Beck_KG_V1_present_0p0083.tif'  

# site lat and lon
data = pd.read_csv('./site_latlon.csv',index_col=0) 
sites = data.index

# Define output dataframe for storing site climate types
climate = pd.DataFrame(columns=['Köppen_climate','Veg_cli'])

############################################################
# process each site (loop)
print('#############################################################################')
print('Pocess each site (loop)!!!')
print('#############################################################################')
############################################################
for site in sites:

    print(f'Processing site {site}!')

    lat = np.round(data.loc[site].lat,4)
    lon = np.round(data.loc[site].lon,4)
    value = get_value_by_coordinates(file, [lon, lat], coordinate_type='lonlat')

    # Climate classification of vegetation using Köppen climate data according to Poulter et al., (2011)
    if value <= 4 or value == 6:
        veg_cli = 'tropical'
    elif value <= 16:
        veg_cli = 'temperature'
    elif value >= 17:
        veg_cli = 'boreal'

    # Store the site Köppen climate type number and the corresponding vegetation climate type.
    climate.loc[site] = [value,veg_cli]


# OUT
climate.to_csv('./Veg_cli.csv')

# END
print('#############################################################################')
print('Processing completed, please check the Veg_cli.csv in the current folder')
print('#############################################################################')






