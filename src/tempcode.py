#
#
# shp_lst = [#path_dict['path_shp'] + '07191500.shp',
#            # path_dict['path_shp'] + '07050700.shp',
#            path_dict['path_shp'] + '05330000.shp',
#            path_dict['path_shp'] + '10333000.shp',
#            path_dict['path_shp'] + '07052500.shp',
#            path_dict['path_shp'] + '13081500.shp',
#            path_dict['path_shp'] + '06657000.shp',
#            path_dict['path_shp'] + '02469761.shp',
#            path_dict['path_shp'] + '02359170.shp']

#
#
# ################# test the code
# shapefile = shp_lst[1]
# topo_hgt_lst = glob.glob(path_dict['topography_path'] + '*[44-48]*.hgt')
# shape = gpd.read_file(shapefile)
# intersection_lst = []
# for topo in (topo_hgt_lst):
#     raster = gdal.Open(topo)
#     vector = ogr.Open(shapefile)    # vector and raster should be in the same coordinate system.
#
#     # Get raster geometry
#     transform = raster.GetGeoTransform()
#     pixelWidth = transform[1]
#     pixelHeight = transform[5]
#     cols = raster.RasterXSize
#     rows = raster.RasterYSize
#
#     xLeft = transform[0]
#     yTop = transform[3]
#     xRight = xLeft + cols * pixelWidth
#     yBottom = yTop - rows * pixelHeight
#
#     ring = ogr.Geometry(ogr.wkbLinearRing)
#     ring.AddPoint(xLeft, yTop)
#     ring.AddPoint(xLeft, yBottom)
#     ring.AddPoint(xRight, yTop)
#     ring.AddPoint(xRight, yBottom)
#     ring.AddPoint(xLeft, yTop)
#     rasterGeometry = ogr.Geometry(ogr.wkbPolygon)
#     rasterGeometry.AddGeometry(ring)
#
#     # Get vector geometry
#     layer = vector.GetLayer()
#     feature = layer.GetFeature(0)
#     vectorGeometry = feature.GetGeometryRef()
#     if rasterGeometry.Intersect(vectorGeometry):
#         # if True, it means they have overlap
#         intersection_lst.append(topo)
#
#
# shape = shape.to_crs('epsg:4269')
# shape_path_temp = os.path.join(path_dict['tempFolder_path'], 'watershed.shp')
# shape.to_file(shape_path_temp)