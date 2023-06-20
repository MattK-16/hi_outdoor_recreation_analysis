import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('agg')
import numpy as np
import geopandas as gpd
import movingpandas as mpd
from shapely.geometry import Point
from scipy.spatial import cKDTree
import math

def convertToPointGpd(data):
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
    # based on some testing, 652 million metres squared is still a serviceable area,
    if (((int(frame.total_bounds[2]) - int(frame.total_bounds[0])) * (int(frame.total_bounds[3]) - int(frame.total_bounds[1]))) > 652000000000):
        print("Study size is too large, please subset your data to an area smaller than 650 million metres squared and try again. Try using a smaller number of individuals from the dataset. ")
        exit()

    return frame

#perform operations on the array to format it for the TrajectoryCollection
def cellPointArray(geoframe, reshaped_array, xres, yres, input_qgs_rect):           

    #get the crs of the original geoframe to return it to the user as it was given
    geo_crs = geoframe.crs

    #compute the appropriate pixel size
    pixel_x_size = (input_qgs_rect.xMaximum() - input_qgs_rect.xMinimum()) / (reshaped_array.shape[1])
    pixel_y_size = (input_qgs_rect.yMaximum() - input_qgs_rect.yMinimum()) / (reshaped_array.shape[0])
    xmin = input_qgs_rect.xMinimum()
    ymin = input_qgs_rect.yMinimum()
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
def computePointNearest(bandGdf, geoframe, pixel_x_size, pixel_y_size):

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
def producePointOutputs(intersection_gdf, moveapps_io):
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
        
    ########################### 
    except Exception as e:
        print("Couldn't generate figures. " + str(e))

#revert back to movingpandas TrajectoryCollection
def gpdToMpd(gpd):
    return mpd.TrajectoryCollection(gpd, "trackId", t="timestamps", crs=gpd.crs)