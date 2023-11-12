import os
import pandas as pd
import numpy as np
# import pyarrow
import glob
import geopandas as gpd
# import multiprocessing
# import timeit
# import ogr
# import gdal
# import shutil
# from matplotlib import pyplot
# import matplotlib.pyplot as plt
import haversine as hs  # this is for distance calculation between two points (pip install haversine)
# import pysheds
# # import xarray as xr
# from pysheds.grid import Grid  # use conda install -c conda-forge pysheds
from src.createOutput import (
    outputDir,
    find_watershed_outlet,
    general_purpose_watershed,
    create_FlowAccu_tif_file,
    finding_NDAMS
)

## defining paths for shapefiles
## shapefiles are in different folders (different categories)

path_dict = {}
# path_dict['output_dir'] = r"/data/fzr5082/"
# path_dict['output_dir'] = r"/scratch/fzr5082"   # r"/data/fzr5082/"
# path_dict['path_shp'] = r"/data/wxt42/hydroDL_reservoir_data/total_shp/CONUS_basin_shp/"  # shapefiles directory
# # path_dict['path_shp'] = r"/data/wxt42/raw_data/shape_file/SRB_HUC10_all/"  # shapefiles directory
# path_dict['path_dams_shp'] = r"/data/fzr5082/NID/maj_dams_shp/us-dams.shp"  # dams point shapefile
# path_dict['path_dams_info'] = r"/data/fzr5082/NID/NID2019_U.feather"  # dams excel file information
# path_dict['gages_path'] = r"/data/wxt42/hydroDL_reservoir_data/gagesII_9322_point_shapefile/gagesII_9322_sept30_2011.shp"  # streamflow gauges point shapefile
# path_dict['topography_path'] = r"/data/shared_data/NED_10m/"
path_dict['output_dir'] = r"/data/fzr5082/"
path_dict['path_shp'] = r"/data/wxt42/hydroDL_reservoir_data/total_shp/CONUS_basin_shp/"  # shapefiles directory
# path_dict['path_shp'] = r"/data/wxt42/raw_data/shape_file/SRB_HUC10_all/"  # shapefiles directory
path_dict['path_dams_shp'] = r"/data/fzr5082/NID/maj_dams_shp/us-dams.shp"  # dams point shapefile
path_dict['path_dams_info'] = r"/data/fzr5082/NID/NID2019_U.feather"  # dams excel file information
path_dict['gages_path'] = r"/data/wxt42/hydroDL_reservoir_data/gagesII_9322_point_shapefile/gagesII_9322_sept30_2011.shp"  # streamflow gauges point shapefile
path_dict['topography_path'] = r"/data/shared_data/NED_10m/"
# path_dict['topography_path'] = r"/data/shared_data/DEM/SRTMGL1v003_30m/"
# Do we need Flow accumulation (DEM part) part? If yes --> needDEM = True, if No --> needDEM = False
needDEM = False
# creating directories
path_dict = outputDir(path_dict)

## Reading shapefiles directory
# dirs_shp = os.listdir(path_shp)




# It is a shapefile that illustrates each dam as a point
dams_shp = gpd.read_file(path_dict['path_dams_shp'])
dams_shp_prj = dams_shp.to_crs('epsg:4269')  # Do projection for the shapefile
## some information of the dams in dams_shp is in the following file (dams_info)
dams_info = pd.read_feather(path_dict['path_dams_info'])
dams_info_m = dams_info.loc[dams_info['OTHERSTRUCTUREID']!=dams_info['OTHERSTRUCTUREID']].reset_index(drop=True)
dams_info_gdf = gpd.GeoDataFrame(
    dams_info_m, geometry=gpd.points_from_xy(dams_info_m.LONGITUDE, dams_info_m.LATITUDE))

dams_info_gdf.set_crs(epsg=4269, inplace=True)
print("dams information file was read")
### remove maj dams that are not in dams_info_gdf
tempsites = dams_shp_prj['NID_ID_Cod'].tolist()
x = dams_info_gdf.loc[dams_info_gdf['NIDID'].isin(tempsites),'NIDID'].tolist()
dams_shp_prj = dams_shp_prj.loc[dams_shp_prj['NID_ID_Cod'].isin(x)]


# Reading flow gages shapefile which is used for finding the outlet of a watershed
gages = gpd.read_file(path_dict['gages_path'])
gages_prj = gages.to_crs('epsg:4269')  # Do projection for the shapefile


##################### MAIN PART #########################
def reservoirs(wtshd):
    data = []  # to store all variables we need and append it to pandas dataframe file at the end
    watershed = gpd.read_file(wtshd)
    geom_buffer = watershed.buffer(0)    # it was 400 before, Farshid changed it to 0
    watershed1 = gpd.read_file(wtshd)
    watershed1['geometry'] = geom_buffer
    watershed_prj = watershed1.to_crs('epsg:4269')
    watershed_prj_path = os.path.join(path_dict['tempFolder_path'], 'watershed_prj.shp')
    watershed_prj.to_file(watershed_prj_path)
    watershed_prj_noBuffer = watershed.to_crs('epsg:4269')

    # to get basin's ID:
    # the column_name is different in shapefile directories:
    if "GAGE_ID" in watershed.columns:    # this is for Conus shapefiles,
        data.append(watershed_prj_noBuffer['GAGE_ID'][0])
    elif "HUC10" in watershed.columns:   # this is for running SRB-HUC10 shapefiles
        data.append(watershed_prj_noBuffer['HUC10'][0])

    # to get basin's area:
    # The column_name is different in different shapefiles
    if "AREA" in watershed.columns:
       data.append(watershed_prj_noBuffer['AREA'][0])
    elif "AreaSqKm" in watershed.columns:
        data.append(watershed_prj_noBuffer['AreaSqKm'][0] * 1e6)    # converting sqkm to sqm

    if (needDEM == True) & (data[1] < 24e9):  # 25e9 is the maximum area that we can musaic,
        # larger than this, we get memory issue. needs to get fixed in the future
        flowAccu_temp_path, AccuProcess = create_FlowAccu_tif_file(watershed_prj_path, path_dict)
        ### if AccuProcess = True --> it means there is at least a tif file that has overlap with the watershed
        ### if AccuProcess = False --> it means there is not any tif file that has overlap with the watershed
    else:
        AccuProcess = False

    if AccuProcess == True:
        watershed_outlet = find_watershed_outlet(watershed_prj, flowAccu_temp_path, gages_prj)
    else:
        watershed_outlet = ["_", "_", "_", "_"]

    data.extend(watershed_outlet)

    ## finding all major dams inside a watershed
    selected_maj_dams = gpd.clip(dams_shp_prj, watershed_prj_noBuffer)    # it was dams_shp_prj
    # selected_maj_dams = selected_dams.loc[(selected_dams['DAM_HEIGHT'] >= 50) |
    #                                       (selected_dams['NORMAL_STORAGE'] >= 5000)]
    MAJ_NDAMS = len(selected_maj_dams)  # number of dams inside a watershed
    data.append(MAJ_NDAMS)

    if (MAJ_NDAMS) > 0:  # at least one dam in the watershed

        # now calculate the Minimum and Mean Distance between the dams and watershed_outlet
        # and to get "General_Purpose" of watershed
        distance_outlet_dams = []
        dam_ID = []
        for i, point in enumerate(selected_maj_dams['geometry']):
            x = point.xy[0][0]
            y = point.xy[1][0]
            dam_ID.append(selected_maj_dams['NID_ID_Cod'].values[i])    #NID_ID_Cod

            # calculating distance between outlet and dams
            if (watershed_outlet[2] != "_") & (watershed_outlet[3] != "_"):
                dam_point = (x, y)
                outlet_point = (watershed_outlet[2], watershed_outlet[3])
                distance = hs.haversine(dam_point, outlet_point, unit=hs.Unit.METERS)
                distance_outlet_dams.append(distance)
            else:
                distance_outlet_dams.append(-999)   # it was np.nan before

        data.append(np.nanmin(distance_outlet_dams))  # distance of nearest dam from the outlet
        data.append(np.nanmean(distance_outlet_dams))  # average distance of dams from outlet

        # To calculate General Purpose of each watershed based on the major dams inside it
        dam = dams_info_gdf.loc[dams_info_gdf["NIDID"].isin(dam_ID)]
        if len(dam) > 0:  # sometimes some major dams in shape file are not in the excel file
            general_purpose, max_norm_stor, std_norm_stor = general_purpose_watershed(dam, dams_shp_prj)
        else:   # sometimes the dataset has some problem
            general_purpose = -1
            max_norm_stor = 0
            std_norm_stor = 0

    else:  # No dam in watershed
        data.extend([-999, -999])  # this is for dams distance and watershed outlet (means no dam in watershed)
        general_purpose = -1
        max_norm_stor = 0
        std_norm_stor = 0

    data.append(general_purpose)
    data.append(max_norm_stor)
    data.append(std_norm_stor)
    ### finding NDAMS_2009 and normal storage for all dams:
    NDAMS, STOR_NOR_2009 = finding_NDAMS(watershed_prj_noBuffer, dams_info_gdf, path_dict, data)

    data.append(NDAMS)
    data.append(STOR_NOR_2009)

    print("Watershed:       ", os.path.split(wtshd)[-1])
    return data


result = []
shp_lst = glob.glob(os.path.join(path_dict['path_shp'], '*.shp'))
#########################################################
## single cpu processing
count = 1
for wtshd in shp_lst:   # [1878:]
    data = reservoirs(wtshd)
    result.append(data)
    print(count)
    count = count + 1
    columns = [
        'GAGE_ID', 'AREA', "STAID", "flow_Accu", "LAT_GAGE", "LNG_GAGE", 'MAJ_NDAMS', "NEAREST_DIS", "AVE_DIS",
        'general_purpose', "max_normal_storage", 'std_norm_stor', 'NDAMS_2009', 'STOR_NOR_2009'
    ]
    Reservoirs = pd.DataFrame(result, columns=columns)
    Reservoirs.to_csv(os.path.join(path_dict['output_dir'], 'Reservoirs.csv'))

###########################################################
#### multiprocessing part  ##############
# comment single cpu processing part above
# def parallel_extraction(item):
#     wtshd = shp_lst[item]
#     data = reservoirs(wtshd)
#     #result.append(data)
#     result = np.expand_dims(np.array(data),axis = 0)
#     # print(count)
#     # count = count + 1
#     columns = [
#         'GAGE_ID', 'AREA', "STAID", "flow_Accu", "LAT_GAGE", "LNG_GAGE", 'MAJ_NDAMS', "NEAREST_DIS", "AVE_DIS",
#         'general_purpose', "max-normal_storage", 'std_norm_stor', 'NDAMS_2009', 'STOR_NOR_2009'
#     ]
#
#     Reservoirs = pd.DataFrame(np.array(result), columns=columns)
#     Reservoirs.to_csv(os.path.join(path_dict['output_dir'], os.path.split(wtshd)[1].split(".shp")[0]+'.csv'))
#     shutil.rmtree(path_dict['tempFolder_ID_path'], ignore_errors=False, onerror=None)
#
# items = [x for x in np.arange(len(shp_lst))]
# p = multiprocessing.Pool(5)
#
# # b = p.map(rs.getNLDAS2, items)
# b = p.map(parallel_extraction, items)
# p.close()
# p.join()

print('END')
