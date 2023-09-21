from sdk.moveapps_spec import hook_impl
from movingpandas import TrajectoryCollection
import time
import os
from . import getIntersection as mv


class App(object):

    def __init__(self, moveapps_io):
        self.moveapps_io = moveapps_io

    #app run
    @hook_impl
    def execute(self, data: TrajectoryCollection, config: dict) -> TrajectoryCollection:
        print("Executing...")
        start = time.time()
        geoframe = mv.convertToGeoPandasFrame(data)
        qgs = mv.qgsAppInit()
        strava_raster = mv.getStravaLayer(config.get('keys'))
        numpy_array, xres, yres, input_qgs_rect = mv.extractRasterBands(strava_raster, geoframe)
        banded_array, pixel_x, pixel_y = mv.arrayOperations(geoframe, numpy_array, xres, yres, input_qgs_rect)
        intersected_gpd = mv.combineFrames(banded_array, geoframe, pixel_x, pixel_y)
        mv.produceOutputs(intersected_gpd, self.moveapps_io)
        trajAltered = mv.gpdToMpd(intersected_gpd)
        mv.qgsAppExit(qgs)
        end = time.time()
        print("Completed in " + str(end - start) + " seconds")

        return trajAltered

    
