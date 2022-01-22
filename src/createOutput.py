import os
import shutil
import pandas as pd
import numpy as np
import geopandas as gpd
import glob
import rasterio
import rasterio.mask
import rasterio.shutil
from rasterio.plot import show
from rasterio.merge import merge
import fiona
import haversine as hs # this is for distance calculation between two points
import pysheds
from pysheds.grid import Grid   # use conda install -c conda-forge pysheds
import shutil
from matplotlib import pyplot
import matplotlib.pyplot as plt
import ogr
from shapely.geometry import box
import pycrs


def outputDir(path_dict):
    # creating output directory:
    if not os.path.isdir(path_dict['output_dir']):
        os.mkdir(path_dict['output_dir'])
        os.mkdir(os.path.join(path_dict['output_dir'], 'tempFolder'))
        path_dict['tempFolder_path'] = os.path.join(path_dict['output_dir'], 'tempFolder')
    else:
        if not os.path.isdir(os.path.join(path_dict['output_dir'], 'tempFolder')):
            os.mkdir(os.path.join(path_dict['output_dir'], 'tempFolder'))
        path_dict['tempFolder_path'] = os.path.join(path_dict['output_dir'], 'tempFolder')
    return path_dict

def find_rasterBounds(path_dict):
    filename = 'Bounds_' + path_dict['topography_path'].split(os.sep)[-2] + '.feather'
    if not os.path.isfile(os.path.join(path_dict['rasterBounds_path'], filename)):
        rs_bounds = pd.DataFrame(columns=['filename', 'left', 'bottom', 'right', 'top', 'crs'])
        topo_tif_lst = glob.glob(path_dict['topography_path'] + '*.tif')
        # shpdir_lst = glob.glob(path_dict['path_shp'] + '*')
        # for shpdir in shpdir_lst:
        #     shpfile_lst = glob.glob(shpdir + '*.shp')
        #     for shp in shpfile_lst:
        for tif_file in topo_tif_lst:

            rs = rasterio.open(tif_file)
            bounds = [os.path.split(tif_file)[1], rs.bounds[0], rs.bounds[1], rs.bounds[2], rs.bounds[3], rs.crs.data]   # left, bottom, right, top
            datatemp = pd.Series(bounds, index=rs_bounds.columns)
            rs_bounds = rs_bounds.append(datatemp, ignore_index=True)
        rs_bounds.to_feather(os.path.join(path_dict['rasterBounds_path'], filename))
    else:
        rs_bounds = pd.read_feather(os.path.join(path_dict['rasterBounds_path'], filename))

    return rs_bounds

def general_purpose_watershed(dam, dams_shp_prj):
    # purpose_list = pd.DataFrame([["I", "H", "C", "N", "S", "R", "P", "F", "D", "T", "O", "G"],
    #                              np.zeros(12), np.zeros(12)]).transpose()
    purpose_list = pd.DataFrame([[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
                                 ["_", "I", "H", "C", "N", "S", "R", "P", "F", "D", "T", "O", "G"],
                                 np.zeros(13), np.zeros(13)]).transpose()
    purpose_list.columns = ['purpose_code', 'purpose', 'priority', 'normal_stor']

    if len(dam) == 1:
        if dam['PURPOSES'].values[0] is None:
            pur = (dams_shp_prj.loc[dams_shp_prj['NID_ID_Cod'] == dam['NIDID'].values[0], 'Purposes'].values[0])
            # general_purpose = purpose_list.loc[purpose_list['purpose']==(dams_shp_prj.loc[dams_shp_prj['NID_ID_Cod'] == dam['NIDID'].values[0], 'Purposes'].values[0]),
            #                                    'purpose_code'].values[0]
            if pur is None:
                general_purpose = purpose_list.loc[purpose_list['purpose']=="_", 'purpose_code'].values[0]
            else:     # if there were more than one purpose for a dam
                general_purpose = purpose_list.loc[purpose_list['purpose'] == pur[0], 'purpose_code'].values[0]
        else:
            general_purpose = purpose_list.loc[purpose_list['purpose']==dam['PURPOSES'].values[0][0], 'purpose_code'].values[0]

        ## finding the max(normal_storage)
        max_norm_stor = dam['NORMAL_STORAGE'].values[0]
        if max_norm_stor is None:
            max_norm_stor = (dams_shp_prj.loc[dams_shp_prj['NID_ID_Cod'] == dam['NIDID'].values[0], 'Normal_Sto'].values[0])
            if max_norm_stor is None:
                max_norm_stor = 0
        # finding the standard deviation of normal storage of dams. here, as there is one dam --> std = 0
        std_norm_stor = 0
    else:
        pur = dam['PURPOSES'].tolist()
        for ii, i in enumerate(pur):
            if i is None:
                pur[ii] = dams_shp_prj.loc[dams_shp_prj['NID_ID_Cod'] == dam['NIDID'].tolist()[ii], 'Purposes'].values[0]
                if pur[ii] is None:
                  continue
                elif len(pur[ii]) > 1:
                    pur[ii] = pur[ii][0]
        norm_stor = dam['NORMAL_STORAGE'].tolist()

        for ii, i in enumerate(norm_stor):
            if i is None:
                norm_stor[ii] = dams_shp_prj.loc[dams_shp_prj['NID_ID_Cod'] == dam['NIDID'].tolist()[ii], 'Normal_Sto'].values[0]
                if len(norm_stor[ii]) > 1:
                    norm_stor[ii] = norm_stor[ii][0]
        ### now removing None in pur and norm_stor
        for ii, i in enumerate(pur):
            if i is None:
                del pur[ii]
                del norm_stor[ii]
        for i in range(len(pur)):
            for j in range(len(pur[i])):
                purpose_list.loc[purpose_list['purpose'] == pur[i][j], 'normal_stor'] += norm_stor[i]
                purpose_list.loc[purpose_list['purpose'] == pur[i][j], 'priority'] += j

        # sort by normal storage
        purpose_list = purpose_list.sort_values(by='normal_stor', ascending=False).reset_index(drop=True)
        datatemp = purpose_list.loc[purpose_list['normal_stor'] == purpose_list['normal_stor'][0]]
        # sort by priority
        datatemp = datatemp.sort_values(by='priority').reset_index(drop=True)
        general_purpose = (datatemp.iloc[0]['purpose_code'])

        ## finding the max(normal_storage)
        max_norm_stor = np.nanmax(norm_stor)
        # finding the standard deviation of normal storage of dams
        std_norm_stor = np.nanstd(norm_stor)

    return general_purpose, max_norm_stor, std_norm_stor

def find_watershed_outlet(watershed_prj, flowAccu_prj_path, gages_prj):
    gages_prj_clip = gpd.clip(gages_prj, watershed_prj)
    if len(gages_prj_clip) > 0:    # it means at least one gage is inside the waterhed
        flowAccu_prj = rasterio.open(flowAccu_prj_path)
        gage_accu_val = []
        p1, p2 = np.where(flowAccu_prj.read(1) == flowAccu_prj.read(1).max())
        lat_lon_best = flowAccu_prj.transform * (p1, p2)
        for i, point in enumerate(gages_prj_clip['geometry']):
            x = point.xy[0][0]
            y = point.xy[1][0]
            distance = hs.haversine(lat_lon_best, (x, y), unit=hs.Unit.METERS)
            row, col = flowAccu_prj.index(x, y)
            flow_accu_val = flowAccu_prj.read(1)[row, col]
            ### based flow accumulation value of the gages that we have observed data
            # gage_accu_val.append([gages_prj_clip['STAID'].values[i], flow_accu_val, x, y])

        # gage_accu_val_pd = pd.DataFrame(gage_accu_val, columns=["STAID", "flow_Accu", "LAT_GAGE", "LNG_GAGE"])
        # gage_accu_val_pd = gage_accu_val_pd.sort_values(by="flow_Accu", ascending=False)
        # watershed_outlet = gage_accu_val_pd.iloc[0].tolist()
            gage_accu_val.append([gages_prj_clip['STAID'].values[i], distance, x, y])
        gage_accu_val_pd = pd.DataFrame(gage_accu_val, columns=["STAID", "distance", "LAT_GAGE", "LNG_GAGE"])
        gage_accu_val_pd = gage_accu_val_pd.sort_values(by="distance", ascending=True)
        watershed_outlet = gage_accu_val_pd.iloc[0].tolist()
        return watershed_outlet
    else:
        return ["_", "_", "_", "_"]


def getFeatures(gdf):
    """Function to parse features from GeoDataFrame in such a manner that rasterio wants them"""
    import json
    return [json.loads(gdf.to_json())['features'][0]['geometry']]

def create_FlowAccu_tif_file(shapefile, path_dict):
    topo_tif_lst = glob.glob(path_dict['topography_path'] + '*.tif')
    # shpdir_lst = glob.glob(path_dict['path_shp'] + '*')
    # for shpdir in shpdir_lst[2:3]:
    #     shpfile_lst = glob.glob(shpdir + '/' + '*.shp')
    #     for shapefile in shpfile_lst:
    shape = gpd.read_file(shapefile)
    shape = shape.to_crs('epsg:4269')
    shape_path_temp = os.path.join(path_dict['tempFolder_path'], 'watershed.shp')
    shape.to_file(shape_path_temp)
    with fiona.open(shape_path_temp) as shp:
        shape_geo = [feature['geometry'] for feature in shp]
        rs_lst = []
        for tif_file in topo_tif_lst:
            with rasterio.open(tif_file) as tile:
                try:
                    out_image, out_transform = rasterio.mask.mask(tile, shape_geo, crop=True)
                except ValueError:
                    continue
                else:      # it means there was not any except
                    rs_lst.append(tif_file)

    temp_mosaic_path = os.path.join(path_dict['tempFolder_path'], 'temp_mosaic.tif')
    ### to  release the memory
    tile = None
    if len(rs_lst) > 1:   # need to mosaic
        ### firstly, doing clip, and then mosaic the clipped files
        #read shapefile
        with fiona.open(shape_path_temp, "r") as shpfile:
            shapes = [feature["geometry"] for feature in shpfile]

            # read raster files one by one
            temp_path = []
            for i, rs in enumerate(rs_lst):
                with rasterio.open(rs) as src:
                    bbox = box(np.maximum(src.bounds[0], shpfile.bounds[0]),  # left
                               np.maximum(src.bounds[1], shpfile.bounds[1]),  # bottom
                               np.minimum(src.bounds[2], shpfile.bounds[2]),  # right
                               np.minimum(src.bounds[3], shpfile.bounds[3])    # top
                               )
                    geo = gpd.GeoDataFrame({'geometry': bbox}, index=[0], crs="EPSG:4269")
                    geo = geo.to_crs(crs=src.crs.data['init'])
                    coords = getFeatures(geo)
                    out_img, out_transform = rasterio.mask.mask(dataset=src, shapes=coords, crop=True)
                    out_meta = src.meta.copy()
                    # we need to parse the epsg value from the CRS so that we can create a proj4 -string using pyCRS
                    # library (to ensure that the projection information is saved correctly)
                    epsg_code = int(src.crs.data['init'][5:])
                    # print(epsg_code)
                    out_meta.update({"driver": "GTiff",
                                     "height": out_img.shape[1],
                                     "width": out_img.shape[2],
                                     "transform": out_transform,
                                     "crs": pycrs.parse.from_epsg_code(epsg_code).to_proj4()}
                                    )
                    out_tif = os.path.join(path_dict['tempFolder_path'], 'topo_temp_clip_' + str(i) + '.tif')
                    temp_path.append(out_tif)
                    with rasterio.open(out_tif, "w", **out_meta) as dest:
                        dest.write(out_img)
                    # clipped = rasterio.open(out_tif)
                    # show(clipped, cmap='terrain')
                    dest.close()
                    ### to release the memory
                    dest = None
                    out_img = None







                    # out_image, out_transform = rasterio.mask.mask(src, shapes, crop=True)
                    # out_meta = src.meta
                    # # make a new file
                    # out_meta.update({"driver": "GTiff",
                    #                  "height": out_image.shape[1],
                    #                  "width": out_image.shape[2],
                    #                  "transform": out_transform,
                    #                  "crs": "epsg:4269"})
                    # path = os.path.join(path_dict['tempFolder_path'], 'topo_temp_clip_' + str(i) + '.tif')
                    # temp_path.append(path)
                    # with rasterio.open(path, "w", **out_meta) as dest:
                    #     dest.write(out_image)
                    # dest.close()
                    # ### to release the memory
                    # dest = None
                    # out_image = None

            #### now we do mosaic:

        src_files_to_mosaic = []
        for fp in temp_path:
            src = rasterio.open(fp)
            src_files_to_mosaic.append(src)
        out_image_mos, out_trans = rasterio.merge.merge(src_files_to_mosaic)
        # make a new file
        out_meta = src.meta.copy()
        out_meta.update(
            {
                "driver": "GTiff",
                "height": out_image_mos.shape[1],
                "width": out_image_mos.shape[2],
                "transform": out_trans,
                "crs": "epsg:4269"
            }
        )
        with rasterio.open(temp_mosaic_path, "w", **out_meta) as dest:
            dest.write(out_image_mos)
        dest.close()
        ### to release memory
        dest = None
        out_image_mos = None
        AccuProcess = True
        ###########################################



       #  src_files_to_mosaic = []
       #  for fp in rs_lst:
       #      src = rasterio.open(fp)
       #      src_files_to_mosaic.append(src)
       #  out_image_mos, out_trans = rasterio.merge.merge(src_files_to_mosaic)
       # # make a new file
       #  out_meta = src.meta.copy()
       #  out_meta.update(
       #      {
       #          "driver": "GTiff",
       #          "height": out_image_mos.shape[1],
       #          "width": out_image_mos.shape[2],
       #          "transform": out_trans,
       #          "crs": "epsg:4269"
       #      }
       #  )
       #  with rasterio.open(temp_mosaic_path, "w", **out_meta) as dest:
       #      dest.write(out_image_mos)
       #  dest.close()
       #  ### to release memory
       #  dest = None
       #  out_image_mos = None
       #  AccuProcess = True

    elif len(rs_lst) == 1:
        # copy the tif file to tempfolder to work on
        # shutil.copyfile(rs_lst[0], temp_mosaic_path)
        rasterio.shutil.copy(rs_lst[0], temp_mosaic_path)
        # make a new file
        AccuProcess = True
    else:
        print("There was no Topography .TIFF file for watershed:   " + os.path.split(shapefile)[1], '\n')
        AccuProcess = False
        return temp_mosaic_path, AccuProcess   ##  it means the code stops here. because there is not any .tif file

    ## clip shapefile and mosaic
    try:
        with fiona.open(shape_path_temp, "r") as shpfile:
            shapes = [feature["geometry"] for feature in shpfile]
            with rasterio.open(temp_mosaic_path) as src:
                out_image, out_transform = rasterio.mask.mask(src, shapes, crop=True)
                out_meta = src.meta
                # make a new file
                out_meta.update({"driver": "GTiff",
                                 "height": out_image.shape[1],
                                 "width": out_image.shape[2],
                                 "transform": out_transform,
                                 "crs": "epsg:4269"})
                topo_clip_temp_path = os.path.join(path_dict['tempFolder_path'], 'topo_temp_clip.tif')
                with rasterio.open(topo_clip_temp_path, "w", **out_meta) as dest:
                    dest.write(out_image)
                dest.close()
                ### to release the memory
                dest = None
                out_image = None
    except:
        AccuProcess = False
        return temp_mosaic_path, AccuProcess
    #  calculating flow direction and flow accumulation (using pysheds package)
    try:
        grid = Grid.from_raster(topo_clip_temp_path, data_name='dem')
        grid.fill_depressions(data='dem', out_name='flooded_dem')
        grid.dem = None
        grid.resolve_flats(data='flooded_dem', out_name='inflated_dem')
        grid.flooded_dem = None
        # dirmap = (1, 2, 3, 4, 5, 6, 7, 8)
        dirmap = (64, 128, 1, 2, 4, 8, 16, 32)  # ESRI default
        grid.flowdir(data='inflated_dem', out_name='dir', dirmap=dirmap)
        grid.inflated_dem = None
        grid.accumulation(data='dir', out_name='acc')
        grid.dir = None
        flowAccu_temp_path = os.path.join(path_dict['tempFolder_path'], 'flowAccu_temp.tif')
        # grid.to_raster('dir', flowAccu_temp_path)
        grid.to_raster('acc', flowAccu_temp_path, dtype=np.int32)
        grid = None
        return flowAccu_temp_path, AccuProcess
    except:
        AccuProcess = False
        return flowAccu_temp_path, AccuProcess


def finding_NDAMS(watershed_prj, dams_info_gdf_prj, path_dict,  data):
    dams_clip = gpd.clip(dams_info_gdf_prj, watershed_prj)
    NDAMS = len(dams_clip)
    # STOR_NOR_2009 = 1.233 * (np.nansum(dams_clip['NORMAL_STORAGE'].tolist()))/(watershed_prj['AREA'][0] * 1e-6) # 1.233 Acre-feet to Megaliter
    STOR_NOR_2009 = 1.233 * (np.nansum(dams_clip['NORMAL_STORAGE'].tolist())) / (
                data[1] * 1e-6)  # 1.233 Acre-feet to Megaliter
    return NDAMS, STOR_NOR_2009



