import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('agg')
import movingpandas as mpd
import geopandas as gpd
import numpy as np
from rtree import index
from shapely.geometry import Point, box
from rtree import index
import math

#convert to a line gdf
def convertToLineGDF(data):
    traj = data
    frame = traj.to_line_gdf()
    crs_frame = traj.to_point_gdf()

    crsf = crs_frame.crs
    frame = frame.set_crs(crsf)
    
    if frame.crs.to_epsg() != 3857:
        frame = frame.to_crs('EPSG:3857')
    else:
        frame = frame
    #this creates multiline strings by each trackId
    #frame = frame.dissolve(by='trackId').reset_index()
    if (((int(frame.total_bounds[2]) - int(frame.total_bounds[0])) * (int(frame.total_bounds[3]) - int(frame.total_bounds[1]))) > 652000000000):
        print("Study size is too large, please subset your data to an area smaller than 650 billion metres squared and try again. Try using a smaller number of individuals from the dataset. ")
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

    selected_columns = ["band_1", "band_2", "band_3"]
    df['intensity'] = (df[selected_columns].sum(axis=1) / (255*3))

    #reproject just in case
    if df.crs != 'EPSG:3857':
        df = df.to_crs('EPSG:3857')

    return df, pixel_x_size, pixel_y_size

#combine the dataframes together
def computePointNearest(bandGdf, geoframe, pixel_x_size, pixel_y_size):
    extended_bbox = []
    filtered_points_count = []
    mean_intensity = []
    hypotenuse = math.sqrt((pixel_x_size*pixel_x_size) + (pixel_y_size*pixel_y_size))
    #geoframe = geoframe.dissolve(by="trackId").reset_index()

    #use an rtree
    idx = index.Index()
    for i, geom in enumerate(bandGdf['geometry']):
        idx.insert(i, geom.bounds)

    #perform calculations based on a buffer around each line segment
    for geom in geoframe['geometry']:
        bbox = box(*geom.bounds)
        buffer_geom = bbox.buffer(hypotenuse)
        #range query
        candidate_indices = list(idx.intersection(buffer_geom.bounds))
        #fliter the points
        filtered_points = bandGdf.iloc[candidate_indices][bandGdf.iloc[candidate_indices].geometry.within(buffer_geom)]
        #get the values and append them to the lists
        extended_bbox.append(buffer_geom)
        filtered_points_count.append(len(filtered_points))
        mean_intensity.append(filtered_points['intensity'].mean())

    #assign the final values for visualisations
    geoframe['extended_bbox'] = extended_bbox
    geoframe['filtered_points_count'] = filtered_points_count
    geoframe['mean_intensity'] = mean_intensity
    geoframe['length'] = [geom.length for geom in geoframe['geometry']]

    return geoframe


#produce graphs and figures
def produceOutputs(intersection_gdf, moveapps_io):
    try:
        print(intersection_gdf)
        #filter based on only above 0 band values
        criteria = (intersection_gdf['mean_intensity'] > 0)

        #create a filtered version of only points that intersected the heatmap
        filtered_gdf = intersection_gdf[criteria]

    ############### bar graphs
        has_intersection = len(filtered_gdf)
        no_intersection = len(intersection_gdf) - len(filtered_gdf)

        fig, ax = plt.subplots()

        intersect_labs = ['# mean intensity > 0', '# mean intensity <= 0']
        intersect_values = [has_intersection, no_intersection]
        bar_labels = ['# mean intensity > 0', '# mean intensity <= 0']
        bar_container = ax.bar(intersect_labs, intersect_values)
        bar_colours = ['red', 'blue']

        ax.bar(intersect_labs, intersect_values, label=bar_labels, color=bar_colours)
        ax.bar_label(bar_container, label_type='center')
        plt.savefig(moveapps_io.create_artifacts_file('segmentIntensity.png'), bbox_inches='tight')

    ########################## create histograms for band values
        #check that it is doing the transformation correctly
        columns = ['length', 'mean_intensity', 'mean_intensity']
        titles = ['Segment Length Distribution (m)', 'Mean Intensity Distribution']
        filenames = ['segmentLength.png', 'meanIntensity.png']

        for column, title, filename in zip(columns, titles, filenames):
            fig, ax = plt.subplots()
            ax.hist(intersection_gdf[column], bins=20, log=False)
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