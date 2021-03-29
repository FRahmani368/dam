import os
import pandas as pd
import numpy as np
import glob
import geopandas as gpd
import gdal
import rasterio
import haversine as hs # this is for distance calculation between two points
import pysheds
import xarray as xr
from pysheds.grid import Grid   # use conda install -c conda-forge pysheds

## defining paths for shapefiles
## shapefiles are in different folders (different categories)

path_shp = os.getcwd() + r'\tot'
dirs_shp = os.listdir(path_shp)

###### difining a pandas file to store averything inside it.
###### This is my results
Reservoirs = pd.DataFrame(columns=[
    'GAGE_ID', 'AREA', "STAID", "flow_Accu", "LAT_GAGE", "LNG_GAGE", 'MAJ_NDAMS', "NEAREST_DIS", "AVE_DIS", 'generral_purpose'
    ])

# It is a shapefile that illustrates each dam as a point
path_dams_shp = r'NID\us-dams.shp'
dams_shp = gpd.read_file(path_dams_shp)

## some information of the dams in dams_shp is in the following file (dams_info)
dams_info = pd.read_feather(r'NID\NID2019_U.feather')
print("dams information file was read")



##################### MAIN PART #########################

# Do projection for the shapefile
dams_shp_prj = dams_shp.to_crs({'init':'epsg:4326'})

# Reading the Raster file of flow accumulation (D8 method)
flowAccu_path = r"Topography\FlowAccuD8.tif"
flowAccu_prj_path = os.path.join(os.path.sep, os.getcwd(), os.path.split(flowAccu_path)[0], os.path.split(flowAccu_path)[1].split('.tif')[0] + '_prj.tif')
gdal.Warp(flowAccu_prj_path, flowAccu_path, dstSRS='EPSG:4326')
# flowdir_raster = xr.open_rasterio(flowdir_path)

# Reading flow gages shapefile which is used for finding the outlet of a watershed
gages_path = r"gagesII_9322_point_shapefile\gagesII_9322_sept30_2011.shp"
gages = gpd.read_file(gages_path)
gages_prj = gages.to_crs({'init':'epsg:4326'})


def general_purpose_watershed(dam):
    purpose_list = pd.DataFrame([["I", "H", "C", "N", "S", "R", "P", "F", "D", "T", "O", "G"],
                                 np.zeros(12), np.zeros(12)]).transpose()
    purpose_list.columns = ['purpose', 'priority', 'normal_stor']

    if len(dam) == 1:
        general_purpose = dam['PURPOSES'].values[0][0]
    else:
        pur = dam['PURPOSES'].tolist()
        norm_stor = dam['NORMAL_STORAGE'].tolist()
        for i in range(len(dam)):
            purpose_list.loc[purpose_list['purpose'] == pur[i], 'normal_stor'] += norm_stor[i]
            purpose_list.loc[purpose_list['purpose'] == pur[i], 'priority'] += 1

        # sort by normal storage
        purpose_list = purpose_list.sort_values(by='normal_stor')
        datatemp = purpose_list.loc[purpose_list['normal_stor'] == purpose_list['normal_stor'][0]]
        # sort by priority
        datatemp = datatemp.sort_values(by='priority')
        general_purpose = (datatemp.iloc[0]['purpose'])
    return general_purpose


def find_watershed_outlet(watershed_prj, flowAccu_prj_path, gages_prj):
    gages_prj_clip = gpd.clip(gages_prj, watershed_prj)
    if len(gages_prj_clip) > 0:    # it means at least one gage is inside the waterhed
        flowAccu_prj = rasterio.open(flowAccu_prj_path)
        gage_accu_val = []
        for i, point in enumerate(gages_prj_clip['geometry']):
            x = point.xy[0][0]
            y = point.xy[1][0]
            row, col = flowAccu_prj.index(x, y)
            flow_accu_val = flowAccu_prj.read(1)[row, col]

            gage_accu_val.append([gages_prj_clip['STAID'].values[i], flow_accu_val, x, y])

        gage_accu_val_pd = pd.DataFrame(gage_accu_val, columns=["STAID", "flow_Accu", "LAT_GAGE", "LNG_GAGE"])
        gage_accu_val_pd = gage_accu_val_pd.sort_values(by="flow_Accu", ascending=False)
        watershed_outlet = gage_accu_val_pd.iloc[0].tolist()
        return watershed_outlet
    else:
        return ["_", "_", "_", "_"]


for dirs in dirs_shp:
    path = os.path.join(os.path.sep, path_shp, dirs)
    shp_list = glob.glob(path + '\*.shp')
    for wtshd in shp_list[0:10]:
        data = []   # to store all variables we need and append it to pandas dataframe file at the end
        watershed = gpd.read_file(wtshd)
        # reproject watershed to epsg:4326 using gdal
        watershed_prj = watershed.to_crs({'init':'epsg:4326'})
        data.append(watershed_prj['GAGE_ID'][0])
        data.append(watershed_prj['AREA'][0])


        #### clipping flow Accumulation raster file by watershed_prj using gdal
        watershed_prj.to_file("watershed_prj.shp")    # first save the shape file
        tempfile_name = os.path.join(os.getcwd(), "flowAccu_prj_clip.tif")      # output raster file path
        watershed_prj_path = os.path.join(os.getcwd(), "watershed_prj.shp")      # the shape file path
        gdal.Warp(tempfile_name, flowAccu_prj_path, cutlineDSName=watershed_prj_path, cropToCutline=True)
        watershed_outlet = find_watershed_outlet(watershed_prj, flowAccu_prj_path, gages_prj)
        data.extend(watershed_outlet)

        ## finding all dams inside a watershed
        selected_dams = gpd.clip(dams_shp_prj, watershed_prj)
        No_dams = len(selected_dams)   # number of dams inside a watershed
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

            data.append(np.nanmin(distance_outlet_dams))   #  distance of nearest dam from the outlet
            data.append(np.nanmean(distance_outlet_dams))      # average distance of dams from outlet

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

        # appending data that is collected for this watershed to Reservoirs dataframe
        data = pd.Series(data, index=Reservoirs.columns)
        Reservoirs = Reservoirs.append(data, ignore_index=True)
        # print(watershed)

        # removing the projected shapefile which was stored in temporary_folder
        # arcpy.management.Delete(watershed_prj)

print('END')