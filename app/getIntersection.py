from qgis.core import *
import pyproj
import urllib.parse
import matplotlib
import matplotlib.pyplot as plt
import branca.colormap as cm
import folium
matplotlib.use('agg')
import time
import numpy as np
import pandas as pd
import geopandas as gpd
import movingpandas as mpd
from shapely.geometry import Point
from scipy.spatial import cKDTree
import math
import os

#convert the movingpandas trajectorycollection to a geopandas dataframe to work with
def convertToGeoPandasFrame(data):
    #read the testing dataset
    traj = data

    #convert the trajectory collection to a dataframe
    frame = traj.to_point_gdf()

    #put the crs in the same epsg the strava dataset will be in
    if frame.crs.to_epsg() != 3857:
        frame = frame.to_crs('EPSG:3857')
    else:
        frame = frame
    #this ends the program if the dataset is too large, eventually I will subset the data based on the track id or on area
    # in place because later on the blocking of the data fails when the raster is too large
    if ((int(frame.total_bounds[2]) - int(frame.total_bounds[0])) > 500000) | ((int(frame.total_bounds[3]) - int(frame.total_bounds[1])) > 500000):
        print("Study size is too large, please subset your data to an area smaller than 500,000 metres by 500,000 metres and try again. Try using a smaller number of individuals from the dataset. ")
        exit()

    return frame

#initialize the qgis environment
def qgsAppInit():
    ##set the path to the qgis bin
    try:
        if (os.environ.get('LOCAL_DEV', 'off') == 'on'):
            # local development
            QgsApplication.setPrefixPath("./Library/bin/", True)
        else:
            # running in production (aka on MoveApps)
            # kudos: https://gis.stackexchange.com/a/263853
            os.environ["QT_QPA_PLATFORM"] = "offscreen"
            QgsApplication.setPrefixPath("/usr", False)
    except Exception as e:
        print("Error setting application path. \n" + str(e))

    ##initialize the qgis session
    try:
        qgs = QgsApplication([], False)
        qgs.initQgis()
    except Exception as e:
        print("Error initializing QGIS for analysis. \n" + str(e))

    return qgs

#exit the qgis environment
def qgsAppExit(qgs):
    time.sleep(1)
    qgs.exitQgis()

#format the strava url based on the user input
def getStravaLayer(urlKeysEntry):

    #basic checks for missing elements of the cookie information
    if urlKeysEntry == None:
        print("No information passed. Please input correct parameter and try again. ")
        exit()
    elif (urlKeysEntry.find("Key-Pair-Id") == -1):
        print("Key-Pair-Id incorrect.")
        exit()
    elif (urlKeysEntry.find("Policy") == -1):
        print("Policy incorrect.")
        exit()
    elif (urlKeysEntry.find("Signature") == -1):
        print("Signature incorrect.")
        exit()

    #after the word red, the url needs to be encoded up to the end of the signature key in the cookies
    urlStart = "type=xyz&url=https://heatmap-external-a.strava.com/tiles-auth/all/blue/"
    urlZoomSet = "{z}/{x}/{y}.png"
    urlQ = "?"
    urlKeys = urlKeysEntry
    urlZoomParams = "&zmax=12&zmin=0"

    #configure the url and initialize the conection to Strava
    try:
        #format the url
        url = urlStart + urllib.parse.quote(urlZoomSet) + urlQ + urllib.parse.quote(urlKeys) + urlZoomParams
                
        #create the xyz tile connection
        strava_raster = QgsRasterLayer(url, "Red", "wms")

    except Exception as e:
        print("Error initializing the connection to Strava. \nThere is a problem with the url provided. Please ensure url is formatted correctly and try again.\n" + str(e))
        exit()

    return strava_raster

#get the strava data and the bands of the raster
def extractRasterBands(strava_raster, geoframe):

    #this should always be true
    if (strava_raster.isValid() == True):
        try:
            #set the pipe and provider based on the strava raster layer
            pipe = QgsRasterPipe()
            provider = strava_raster.dataProvider()
            pipe.set(provider.clone())
        except Exception as e:
            print("Error writing file. " + str(e))

        #create a qgs rectangle based on the geoframe
        try:
            input_qgs_rect = QgsRectangle(geoframe.total_bounds[0], geoframe.total_bounds[1], geoframe.total_bounds[2], geoframe.total_bounds[3])
        except Exception as e:
            print("Error obtaining the shapefile. Is it in a shapefile? " + str(e))

        #computes the resolution, should be a whole number
        try:
            pixel_size = 100  #meters (approximately)
            xres = int((input_qgs_rect.xMaximum() - input_qgs_rect.xMinimum()) / pixel_size)
            yres = int((input_qgs_rect.yMaximum() - input_qgs_rect.yMinimum()) / pixel_size)
        except Exception as e:
            print("Error formatting resolution. " + str(e))
            exit()
    else:
        print("Raster layer was invalid, please ensure the url parameters you entered were correct. ")
        exit()

    #block the raster to get the data based on the qgis rectangle
    # this fails if the study areas is too large - threshold for failure is unknown
    try:
        raster_bands = provider.block(1, input_qgs_rect, xres, yres)
    except Exception as e:
        print("Error converting to byte array. " + str(e))
        exit()
    except QgsException as e:
        print("Timer exception. ")
        exit()
    try:
        byte_array = raster_bands.data()
    except Exception as e:
        print("Error converting to byte array. " + str(e))
        exit()

    #convert to a numpy array
    numpy_array = np.frombuffer(byte_array.data(), dtype=np.uint8)
    #reshape the array based on the number of bands
    reshaped_array = np.reshape(numpy_array, (yres, xres, 4))

    return reshaped_array, xres, yres, input_qgs_rect

#perform operations on the array to format it for the TrajectoryCollection
def arrayOperations(geoframe, reshaped_array, xres, yres, input_qgs_rect):           

    #get the crs of the original geoframe to return it to the user as it was given
    geo_crs = geoframe.crs

    #compute the appropriate pixel size
    pixel_x_size = (input_qgs_rect.xMaximum() - input_qgs_rect.xMinimum()) / (reshaped_array.shape[1])
    pixel_y_size = (input_qgs_rect.yMaximum() - input_qgs_rect.yMinimum()) / (reshaped_array.shape[0])
    xmin = input_qgs_rect.xMinimum()
    ymin = input_qgs_rect.yMinimum()
    xmax = input_qgs_rect.xMaximum()
    ymax = input_qgs_rect.yMaximum()
    data = []

    #make a grid of x and y coordinates
    x_coords = np.arange(0, xres) * pixel_x_size + xmin
    y_coords = np.arange(yres - 1, -1, -1) * pixel_y_size + ymin

    #create a meshgrid
    x_mesh, y_mesh = np.meshgrid(x_coords, y_coords)

    #flatten into 1D arrays
    x_flattened = x_mesh.flatten()
    y_flattened = y_mesh.flatten()

    #calculate midpoints
    x_midpoints = x_flattened + (pixel_x_size / 2)
    y_midpoints = y_flattened + (pixel_y_size / 2)

    #create shapely points
    geometries = [Point(x, y) for x, y in zip(x_midpoints, y_midpoints)]

    #create a dictionary:key set for the bands
    band_data = {
            'band_1': reshaped_array[:, :, 2].flatten(),
            'band_2': reshaped_array[:, :, 1].flatten(),
            'band_3': reshaped_array[:, :, 0].flatten(),
            'band_4': reshaped_array[:, :, 3].flatten()
        }

    #append the geometry and band data to a geodataframe
    data = {'geometry': geometries, **band_data}
    df = gpd.GeoDataFrame(data, crs=geo_crs)

    # remove rows with missing geometry
    df = df[~df['geometry'].isnull()]

    #reproject just in case
    if df.crs != 'EPSG:3857':
        df = df.to_crs('EPSG:3857')

    return df, pixel_x_size, pixel_y_size

#combine the dataframes together
def combineFrames(bandGdf, geoframe, pixel_x_size, pixel_y_size):

    maxBandValue = 255

    #extract coords
    geoframe_coords = [(geom.x, geom.y) for geom in geoframe.geometry]

    #extract coords
    bandGdf_coords = [(geom.x, geom.y) for geom in bandGdf.geometry]

    #create kdtree for searching
    tree = cKDTree(bandGdf_coords)

    #get nearest neighbors and distances
    distances, indices = tree.query(geoframe_coords)

    #assign the nearest neighbor distances to a new column
    geoframe['nearest_distance'] = distances

    #extract band values
    band1_values = bandGdf.loc[indices, 'band_1'].tolist()
    band2_values = bandGdf.loc[indices, 'band_2'].tolist()
    band3_values = bandGdf.loc[indices, 'band_3'].tolist()
    band4_values = bandGdf.loc[indices, 'band_4'].tolist()

    #append the band values to the geoframe
    geoframe['band_1'] = band1_values
    geoframe['band_2'] = band2_values
    geoframe['band_3'] = band3_values
    geoframe['band_4'] = band4_values
    
    #set radius for the distance to include/exclude
    hypotenuse = math.sqrt((pixel_x_size*pixel_x_size) + (pixel_y_size*pixel_y_size))

    #calculate intensity based on the band values and distances
    selected_columns = ["band_1", "band_2", "band_3"]
    geoframe['intensity'] = (geoframe[selected_columns].sum(axis=1) / (maxBandValue*3))

    #if the distance is too far, set the values to na
    criterion = geoframe['nearest_distance'] > (hypotenuse / 2)
    geoframe.loc[criterion, ['band_1', 'band_2', 'band_3', 'intensity']] = np.nan

    return geoframe

#produce graphs and figures
def produceOutputs(intersection_gdf, moveapps_io):
    try:
        #filter based on only above 0 band values
        criteria = (intersection_gdf['band_1'] > 0) | (intersection_gdf['band_2'] > 0) | (intersection_gdf['band_3'] > 0)

        #create a filtered version of only points that intersected the heatmap
        filtered_gdf = intersection_gdf[criteria]

    ############### bar graphs
        has_intersection = len(filtered_gdf)
        no_intersection = len(intersection_gdf) - len(filtered_gdf)

        fig, ax = plt.subplots()

        intersect_labs = ['# intersecting pixels', '# not intersecting pixels']
        intersect_values = [has_intersection, no_intersection]
        bar_labels = ['# intersecting', '# not intersecting']
        bar_container = ax.bar(intersect_labs, intersect_values)
        bar_colours = ['red', 'blue']

        ax.bar(intersect_labs, intersect_values, label=bar_labels, color=bar_colours)
        ax.bar_label(bar_container, label_type='center')
        plt.savefig(moveapps_io.create_artifacts_file('intersections.png'), bbox_inches='tight')

    ########################## create histograms for band values
        #check that it is doing the transformation correctly
        columns = ['band_1', 'band_2', 'band_3', 'intensity']
        titles = ['Band 1 RGB Value', 'Band 2 RGB Value', 'Band 3 RGB Value', 'Human Activity Intensity Distribution']
        filenames = ['band1.png', 'band2.png', 'band3.png', 'intensity.png']

        for column, title, filename in zip(columns, titles, filenames):
            fig, ax = plt.subplots()
            ax.hist(intersection_gdf[column], bins=20, log=True)
            ax.set_xlabel(title)
            ax.set_ylabel("Log of Count")
            plt.savefig(moveapps_io.create_artifacts_file(filename), bbox_inches='tight')
            plt.close(fig)
    ########################### making map
        #convert band valeus to hex rgb
        def rgb_to_hex(rgb):
            return "#{:02X}{:02X}{:02X}".format(int(rgb['band_1']), int(rgb['band_2']), int(rgb['band_3']))

        #get the unique number of tracks
        unique_track_ids = intersection_gdf['trackId'].unique()
        intersection_gdf = intersection_gdf.to_crs("EPSG:4326")
        #create the colouring column for points on map
        intersection_gdf['colouring'] = intersection_gdf.apply(rgb_to_hex, axis=1)
        #remove black points
        intersection_gdf = intersection_gdf[intersection_gdf['colouring'] != '#000000'].copy()

        cent = intersection_gdf.dissolve().centroid

        white = '#FFFFFF'
        blue = '#0000FF'
        black = '#000000'
        colours = [black, blue, white]
        colourmap = cm.LinearColormap(colours, vmin=0, vmax=1, max_labels = 2)
        colourmap.caption = "Human Activity Intensity Low (0) to High (1)"

        m = folium.Map(
            location=[cent.y[0], cent.x[0]],
            zoom_start=11,
            tiles="OpenStreetMap",
            attr="(c) OpenStreetMap",
        )
        
        k = 0
        for individual_id in unique_track_ids:
            k += 1
            if k<11:
                feature_group = folium.FeatureGroup(individual_id)
                working_gdf = intersection_gdf[intersection_gdf["trackId"] == individual_id]
                i = 0
                for item, row in working_gdf.iterrows():
                    #Place the markers with the popup labels and data
                    feature_group.add_child(
                        folium.Circle(
                            location=(row.geometry.y, row.geometry.x),
                            popup="\nIntensity: " + str(row["intensity"]),
                            color = row["colouring"], 
                            fill = True,
                            fillOpacity = 1
                        )
                    )
                    i = i + 1
                feature_group.add_to(m)
                i = 0

        folium.LayerControl(position="topright").add_to(m)
        m.add_child(colourmap)
        
        m.save(moveapps_io.create_artifacts_file('map.html'))

    except Exception as e:
        print("Couldn't generate figures. " + str(e))
        
#revert back to movingpandas TrajectoryCollection
def gpdToMpd(gpd):
    gpd = gpd.to_crs('EPSG:4326')
    return mpd.TrajectoryCollection(gpd, "trackId", t="timestamps", crs=gpd.crs)
    