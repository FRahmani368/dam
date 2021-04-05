import os
import shutil
import pandas as pd
import numpy as np
import geopandas as gpd
import glob
import rasterio
import rasterio.mask
from rasterio.merge import merge
import fiona
import pysheds
from pysheds.grid import Grid   # use conda install -c conda-forge pysheds
import shutil
from matplotlib import pyplot
import matplotlib.pyplot as plt
import ogr


def outputDir(path_dict):
    # creating output directory:
    if not os.path.isdir(path_dict['output_dir']):
        os.mkdir(path_dict['output_dir'])
        os.mkdir(os.path.join(path_dict['output_dir'], 'flowAccu'))
        path_dict['flowAccu_path'] = os.path.join(path_dict['output_dir'], 'flowAccu')
        os.mkdir(os.path.join(path_dict['output_dir'], 'tempFolder'))
        path_dict['tempFolder_path'] = os.path.join(path_dict['output_dir'], 'tempFolder')
        os.mkdir(os.path.join(path_dict['output_dir'], 'rasterBounds'))
        path_dict['rasterBounds_path'] = os.path.join(path_dict['output_dir'], 'rasterBounds')
    else:
        if not os.path.isdir(os.path.join(path_dict['output_dir'], 'flowAccu')):
            os.mkdir(os.path.join(path_dict['output_dir'], 'flowAccu'))
        path_dict['flowAccu_path'] = os.path.join(path_dict['output_dir'], 'flowAccu')

        if not os.path.isdir(os.path.join(path_dict['output_dir'], 'tempFolder')):
            os.mkdir(os.path.join(path_dict['output_dir'], 'tempFolder'))
        path_dict['tempFolder_path'] = os.path.join(path_dict['output_dir'], 'tempFolder')

        if not os.path.isdir(os.path.join(path_dict['output_dir'], 'rasterBounds')):
            os.mkdir(os.path.join(path_dict['output_dir'], 'rasterBounds'))
        path_dict['rasterBounds_path'] = os.path.join(path_dict['output_dir'], 'rasterBounds')

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
    if len(rs_lst) > 1:   # need to mosaic
        src_files_to_mosaic = []
        for fp in rs_lst:
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
        AccuProcess = True
    elif len(rs_lst) == 1:
        # copy the tif file to tempfolder to work on
        shutil.copyfile(rs_lst[0], temp_mosaic_path)
        AccuProcess = True
    else:
        print("There was no Topography .TIFF file for watershed:   " + os.path.split(shapefile)[1], '\n')
        AccuProcess = False
        return temp_mosaic_path, AccuProcess   ##  it means the code stops here. because there is not any .tif file

    ## clip shapefile and mosaic
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

    #  calculating flow direction and flow accumulation (using pysheds package)
    grid = Grid.from_raster(topo_clip_temp_path, data_name='dem')
    grid.fill_depressions(data='dem', out_name='flooded_dem')
    grid.resolve_flats(data='flooded_dem', out_name='inflated_dem')

    # dirmap = (1, 2, 3, 4, 5, 6, 7, 8)
    dirmap = (64, 128, 1, 2, 4, 8, 16, 32)  # ESRI default
    grid.flowdir(data='inflated_dem', out_name='dir', dirmap=dirmap)
    grid.accumulation(data='dir', out_name='acc')
    flowAccu_temp_path = os.path.join(path_dict['tempFolder_path'], 'flowAccu_temp.tif')
    grid.to_raster('dir', flowAccu_temp_path)
    return flowAccu_temp_path, AccuProcess