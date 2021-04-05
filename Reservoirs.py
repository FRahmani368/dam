import os
import pandas as pd
import numpy as np
import rasterio
import rasterio.mask
from rasterio.merge import merge
import fiona
import glob
import geopandas as gpd
import multiprocessing
import timeit
import gdal
import shutil
from matplotlib import pyplot
import matplotlib.pyplot as plt
import haversine as hs  # this is for distance calculation between two points (pip install haversine)
import pysheds
# import xarray as xr
from pysheds.grid import Grid  # use conda install -c conda-forge pysheds
from src.createOutput import (
    outputDir,
    find_watershed_outlet,
    general_purpose_watershed,
    create_FlowAccu_tif_file
)

## defining paths for shapefiles
## shapefiles are in different folders (different categories)

path_dict = {}
path_dict['output_dir'] = r"/data/fzr5082/"
path_dict['path_shp'] = r"/data/wxt42/hydroDL_reservoir_data/total_shp/CONUS_basin_shp/"  # shapefiles directory
path_dict['path_dams_shp'] = r"/data/fzr5082/NID/maj_dams_shp/us-dams.shp"  # dams point shapefile
path_dict['path_dams_info'] = r"/data/fzr5082/NID/NID2019_U.feather"  # dams excel file information
path_dict['gages_path'] = r"/data/wxt42/hydroDL_reservoir_data/gagesII_9322_point_shapefile/gagesII_9322_sept30_2011.shp"  # streamflow gauges point shapefile
path_dict['topography_path'] = r"/data/shared_data/NED_10m/"

# creating directories
path_dict = outputDir(path_dict)

## Reading shapefiles directory
# dirs_shp = os.listdir(path_shp)




# It is a shapefile that illustrates each dam as a point
dams_shp = gpd.read_file(path_dict['path_dams_shp'])
dams_shp_prj = dams_shp.to_crs('epsg:4269')  # Do projection for the shapefile
## some information of the dams in dams_shp is in the following file (dams_info)
dams_info = pd.read_feather(path_dict['path_dams_info'])
print("dams information file was read")

# Reading flow gages shapefile which is used for finding the outlet of a watershed
gages = gpd.read_file(path_dict['gages_path'])
gages_prj = gages.to_crs('epsg:4269')  # Do projection for the shapefile


##################### MAIN PART #########################
def reservoirs(wtshd):
    data = []  # to store all variables we need and append it to pandas dataframe file at the end
    watershed = gpd.read_file(wtshd)
    watershed_prj = watershed.to_crs('epsg:4269')
    data.append(watershed_prj['GAGE_ID'][0])
    data.append(watershed_prj['AREA'][0])

    flowAccu_temp_path, AccuProcess = create_FlowAccu_tif_file(wtshd, path_dict)
    ### if AccuProcess = True --> it means there is at least a tif file that has overlap with the watershed
    ### if AccuProcess = False --> it means there is not any tif file that has overlap with the watershed
    if AccuProcess == True:
        watershed_outlet = find_watershed_outlet(watershed_prj, flowAccu_temp_path, gages_prj)
    else:
        watershed_outlet = ["_", "_", "_", "_"]

    data.extend(watershed_outlet)

    ## finding all dams inside a watershed
    selected_dams = gpd.clip(dams_shp_prj, watershed_prj)
    No_dams = len(selected_dams)  # number of dams inside a watershed
    data.append(No_dams)

    if (No_dams) > 0:  # at least one dam in the watershed

        # now calculate the Minimum and Mean Distance between the dams and watershed_outlet
        # and to get "General_Purpose" of watershed
        distance_outlet_dams = []
        dam_ID = []
        for i, point in enumerate(selected_dams['geometry']):
            x = point.xy[0][0]
            y = point.xy[1][0]
            dam_ID.append(selected_dams['NID_ID_Cod'].values[i])

            # calculating distance between outlet and dams
            if (watershed_outlet[2] != "_") & (watershed_outlet[3] != "_"):
                dam_point = (x, y)
                outlet_point = (watershed_outlet[2], watershed_outlet[3])
                distance = hs.haversine(dam_point, outlet_point, unit=hs.Unit.METERS)
                distance_outlet_dams.append(distance)
            else:
                distance_outlet_dams.append(np.nan)

        data.append(np.nanmin(distance_outlet_dams))  # distance of nearest dam from the outlet
        data.append(np.nanmean(distance_outlet_dams))  # average distance of dams from outlet

        # To calculate General Purpose of each watershed based on the dams inside it
        dam = dams_info.loc[dams_info["NIDID"].isin(dam_ID)]
        if len(dam) > 0:  # sometimes some dams in shape file are not in the excel file
            general_purpose = general_purpose_watershed(dam)
        else:
            general_purpose = '_'

    else:  # No dam in watershed
        data.extend([np.nan, np.nan])  # this is for dams distance and watershed outlet (means no dam in watershed)
        general_purpose = '_'

    data.append(general_purpose)
    print("Watershed:       ", os.path.split(wtshd)[-1])
    return data



shp_lst = glob.glob(os.path.join(path_dict['path_shp'], '*.shp'))
result = []
for wtshd in shp_lst[18:]:
    data = reservoirs(wtshd)
    result.append(data)


###########################################################
#### multiprocessing part is not working yet ##############

# num_cores = multiprocessing.cpu_count()
# # num_cores = 1
# print('Number cores:', num_cores, '\n')
# p = multiprocessing.Pool(num_cores)
# start = timeit.default_timer()
# result = p.map(reservoirs, shp_lst)
# p.close()
# p.join()


# appending data that is collected for this watershed to Reservoirs dataframe
###### difining a pandas file to store averything inside it.
###### This is my results
columns = [
    'GAGE_ID', 'AREA', "STAID", "flow_Accu", "LAT_GAGE", "LNG_GAGE", 'MAJ_NDAMS', "NEAREST_DIS", "AVE_DIS",
    'general_purpose'
]
Reservoirs = pd.DataFrame(result, columns=columns)
Reservoirs.to_csv(os.path.join(path_dict['output_dir'], 'Reservoirs.csv'))
print('END')
