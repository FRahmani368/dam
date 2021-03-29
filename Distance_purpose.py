# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 21:45:54 2021

@author: fzr5082
"""
# IMPORTANT NOTES:
    # 1) Please change line 28 to set the workspace
    # 2) Please install haversine Package using "pip install haversine"
    # 3) I know that it is possible (and maybe easier) to calculate distance between two point by arcpy but in our project
    # I wanted to test haversine for my own project. 
    # 4) In Line 151, You can choose the number of watersheds that you want to run the code. Now it is not set.
    # and you can change it
###############

import arcpy
import os
import glob
import pandas as pd
import numpy as np
import haversine as hs # this is for distance calculation between two points 

# This is for resetting the variables  automatically.
from IPython import get_ipython
get_ipython().magic('reset -sf')


arcpy.env.workspace = os.getcwd()
arcpy.env.overwriteOutput = True

## defining paths for shapefiles
## shapefiles are in different folders (different categories)
path_shp = arcpy.env.workspace + r'\tot'  
dirs_shp = os.listdir(path_shp)

# first: deleting all rasters in the workspace, if any.
# This code deletes everythin at the end, tis part is just to be in the safe side 
raster_list = arcpy.ListRasters()
for ras in raster_list:
    arcpy.Delete_management(ras)

###### difining a pandas file to store averything inside it.
###### This is my results
Reservoirs = pd.DataFrame(columns=[
    'GAGE_ID', 'AREA', "STAID", "flow_Accu", "LAT_GAGE", "LNG_GAGE", 'MAJ_NDAMS', "NEAREST_DIS", "AVE_DIS", 'generral_purpose'
    ])

# It is a shapefile that illustrates each dam as a point
dams_shp = r'NID\us-dams.shp'  

## some information of the dams in dams_shp is in the following file (dams_info)
dams_info = pd.read_feather(r'NID\NID2019_U.feather')
print("dams information file was read")



def general_purpose_watershed(dam):
    purpose_list = pd.DataFrame([["I", "H", "C", "N", "S", "R", "P", "F", "D", "T", "O", "G"], 
                                np.zeros(12), np.zeros(12)]).transpose()
    purpose_list.columns = ['purpose', 'priority', 'normal_stor']
    
    if len(dam) == 1:
        general_purpose=dam['PURPOSES'].values[0][0]
    else:
        pur = dam['PURPOSES'].tolist()
        norm_stor = dam['NORMAL_STORAGE'].tolist()
        for i in range(len(dam)):
            purpose_list.loc[purpose_list['purpose']==pur[i], 'normal_stor'] += norm_stor[i]
            purpose_list.loc[purpose_list['purpose']==pur[i], 'priority'] += 1
        
        # sort by normal storage
        purpose_list = purpose_list.sort_values(by='normal_stor')
        datatemp = purpose_list.loc[purpose_list['normal_stor']==purpose_list['normal_stor'][0]]
        # sort by priority
        datatemp = datatemp.sort_values(by='priority')
        general_purpose = (datatemp.iloc[0]['purpose'])
    return general_purpose

def Project(shapefile, coord_system = 4326):
    sr = arcpy.SpatialReference(coord_system)
    new_shapefile = "new_" + shapefile.split("\\")[-1]
    arcpy.Project_management(shapefile, new_shapefile, sr)
    #new_shapefile = arcpy.DefineProjection_management(shapefile, sr)
    return new_shapefile

def Project_raster(rasterfile, newfile_name, coord_system=4326):
    sr = arcpy.SpatialReference(coord_system)
    new_rasterfile = rasterfile.name
    #if os.path.isfile(new_rasterfile):
     #   arcpy.Delete_management(new_rasterfile)
    arcpy.ProjectRaster_management(rasterfile, new_rasterfile, sr)
    new_rasterfile1 = arcpy.Raster(new_rasterfile)
    return new_rasterfile1

def find_watershed_outlet(watershed_prj, flowdir_prj, gages_prj):
    raster_cliped = arcpy.ia.Clip(flowdir_prj, watershed_prj)
    flow_accu = arcpy.sa.FlowAccumulation(raster_cliped)
    # finding all gages inside the watershed
    selected_gages = arcpy.SelectLayerByLocation_management(gages_prj, "CONTAINED_BY", watershed_prj)
    No_selected_gages = arcpy.GetCount_management(selected_gages)
    # finding the gage with the highest flow accumulation number
    if int(No_selected_gages[0])>0:
        gage_accu_val = []
        with arcpy.da.SearchCursor(selected_gages,["SHAPE@", "STAID", "LAT_GAGE", "LNG_GAGE"]) as cursor:
            for row in cursor:
                flow_accu_val = 0
                arcpy.sa.ExtractValuesToPoints(row[0], flow_accu, "outVal")
                with arcpy.da.SearchCursor("outVal.shp", ["RASTERVALU"]) as cursorPoint:
                    for rowPoint in cursorPoint:
                        flow_accu_val = rowPoint[0]
                del cursorPoint
                gage_accu_val.append([row[1], flow_accu_val, row[2], row[3]])
        del cursor
        gage_accu_val_pd = pd.DataFrame(gage_accu_val, columns=["STAID", "flow_Accu", "LAT_GAGE", "LNG_GAGE"])
        gage_accu_val_pd = gage_accu_val_pd.sort_values(by="flow_Accu", ascending=False)
        watershed_outlet = gage_accu_val_pd.iloc[0].tolist()
        return watershed_outlet
    else:
        return ["_", "_", "_", "_"]



##################### MAIN PART #########################

# Do projection for the shapefile
dams_shp_prj = Project(dams_shp)


# Reading the Raster file of elevations in the US
# topo_path = r"Topography\elev_0.01.tif"
# topo_raster = arcpy.Raster(topo_path)
# Do projection for the raster file
# topo_raster_prj = Project_raster(topo_raster, newfile_name="topoprj.tif")


# Reading Flow direction raster file which has been created by Arcgis Pro tool using elevation file
# FlowDirection() tool. This tool is available only in ARCGIS Enterprise
# So I created it first, and saved it as a .TIFF file
flowdir_path = r"Topography\flowdir\FlowDir_elev1.tif"
flowdir = arcpy.Raster(flowdir_path)
# Do projection
flowdir_prj = Project_raster(flowdir, newfile_name="dirprj.tif")

# Reading flow gages shapefile which is used for finding the outlet of a watershed
gages_path = r"gagesII_9322_point_shapefile\gagesII_9322_sept30_2011.shp"
gages_prj = Project(gages_path)

for dirs in dirs_shp:
    path = os.path.join(os.path.sep, path_shp, dirs)
    shp_list = glob.glob(path + '\*.shp')
    for watershed in shp_list:
        watershed_prj = Project(watershed)
        
        # to store results for each watershed in it and append it to "Reservoirs" pandas dataframe
        data = [] 
        
        with arcpy.da.SearchCursor(watershed_prj, ['GAGE_ID', 'AREA']) as cursor:
            for row in cursor:
                data.append(row[0])
                data.append(row[1])
        del cursor
        
        # Now going tofind watershed outlet by claculating the flow accumulation for the watershed
        watershed_outlet = find_watershed_outlet(watershed_prj, flowdir_prj, gages_prj)
        data.extend(watershed_outlet)
        
        ## finding all dams inside a watershed
        selected_dams = arcpy.SelectLayerByLocation_management(dams_shp_prj, "CONTAINED_BY", watershed_prj)
        ## Counting number of dams inside the watershed
        No_dams = arcpy.GetCount_management(selected_dams)
        data.append(int(No_dams[0]))
        
        if (int(No_dams[0]))>0:  # at least one dam in the watershed
        
            # now calculate the Minimum and Mean Distance between the dams and watershed_outlet
            #and to get "Geeneral_Purpose" of watershed
            distance_outlet_dams = [] 
            dam_ID = []
            with arcpy.da.SearchCursor(selected_dams, ['NID_ID_Cod', "Latitude", "Longitude"]) as cursor:
                for row in cursor:
                    dam_ID.append(row[0])
                    # calculating distance between outlet and dams
                    if (watershed_outlet[2] != "_") or (watershed_outlet[3] != "_"): 
                        dam_point = (row[1], row[2])                  
                        outlet_point = (watershed_outlet[2], watershed_outlet[3])
                        distance = hs.haversine(dam_point, outlet_point, unit=hs.Unit.METERS) 
                        distance_outlet_dams.append(distance)
                    else:
                        distance_outlet_dams.append(np.nan)
            del cursor
            data.append(np.nanmin(distance_outlet_dams))
            data.append(np.nanmean(distance_outlet_dams))
            
            #To calculate General Purpose of each watershed based on the dams inside it
            dam = dams_info.loc[dams_info["NIDID"].isin(dam_ID)]
            if len(dam) >0: # sometimes some dams in shape file are not in the excel file
                general_purpose = general_purpose_watershed(dam)
            else:
                general_purpose = '_'
                
        else: # No dam in watershed
            data.extend([np.nan, np.nan])   # this is for dams distance and watershed outlet (means no dam in watershed)
            general_purpose = '_'
        
        data.append(general_purpose)
        
        # appending data that is collected for this watershed to Reservoirs dataframe
        data = pd.Series(data, index=Reservoirs.columns)
        Reservoirs = Reservoirs.append(data, ignore_index=True)
        print(watershed)
        # removing the projected shapefile which was stored in temporary_folder
        arcpy.management.Delete(watershed_prj)

# Deleting shapefiles and rasters
arcpy.management.Delete(dams_shp_prj)
arcpy.management.Delete(gages_prj)
arcpy.management.Delete("outVal.shp")
raster_list = arcpy.ListRasters()
for ras in raster_list:
    arcpy.Delete_management(ras)
    

# Saving the result as a csv file in workspace
Reservoirs.to_csv(r'Reservoirs.csv')
print("END")





