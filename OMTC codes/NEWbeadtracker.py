import sys
import os
import csv
import subprocess

from ij.plugin.frame import RoiManager
from ij import IJ, ImagePlus, ImageStack, WindowManager
from ij.gui import WaitForUserDialog

from fiji.plugin.trackmate import Model
from fiji.plugin.trackmate import Settings
from fiji.plugin.trackmate import TrackMate
from fiji.plugin.trackmate import SelectionModel
from fiji.plugin.trackmate import Logger
from fiji.plugin.trackmate.tracking import LAPUtils
from fiji.plugin.trackmate.gui.displaysettings import DisplaySettingsIO
import fiji.plugin.trackmate.visualization.hyperstack.HyperStackDisplayer as HyperStackDisplayer
import fiji.plugin.trackmate.features.FeatureFilter as FeatureFilter
from fiji.plugin.trackmate import Model, Settings, TrackMate, SelectionModel, Logger
from fiji.plugin.trackmate.detection import LogDetectorFactory, DogDetectorFactory
from fiji.plugin.trackmate.tracking.sparselap import SparseLAPTrackerFactory, SimpleSparseLAPTrackerFactory
from fiji.plugin.trackmate.features import FeatureFilter
from fiji.plugin.trackmate.io import TmXmlReader, TmXmlWriter
from fiji.plugin.trackmate.util import TMUtils
from fiji.plugin.trackmate.detection import DetectorKeys
from fiji.plugin.trackmate.features import FeatureAnalyzer, ModelFeatureUpdater, SpotFeatureCalculator
from fiji.plugin.trackmate.features.spot import SpotContrastAndSNRAnalyzerFactory, SpotContrastAndSNRAnalyzer
from fiji.plugin.trackmate.features.track import TrackDurationAnalyzer, TrackSpeedStatisticsAnalyzer

########################## MAKE SURE ROI 1 CLICK --> 35x35 (only 'Add to ROI Manager' checked) and bead's oscillation is ~centered ################################################################################################################################################################

# Must do the folowing before everything to avoid errors  
reload(sys)
sys.setdefaultencoding('utf-8')

######################################### Path (USER INPUT) #################################################################################################################################################################################

userDirectory = "D:/RAW DATA/"

experiment = "ET_OMTCTFM_12kPa_SC_PC3CSS_PC31nM_48hour_mOrange_heat_13072022"   

condition = "PC3_1nMR1881_12kPa_SC"

#############################################################################################################################################################################################################################################

directory_experiment = userDirectory + experiment
directory_condition = directory_experiment + "/"+condition

# Rename func
def rename_roi(bead_directory, condition, series_num, bead_num):
	
	bead_name = condition +"_series_" + series_num + "_bead_"+ bead_num
	rm.select(i-1)
	RM.runCommand("Rename", bead_name )

# Iterating over condition folder w/ presaved roi zips (already have selected beads)			
folder = os.listdir(directory_condition)
for f in folder: # Goes through condition folder
	if len(f) == 8:
		series = f[7]
	else:
		series = f[7:9]
	print(series)
	
	directory_series = directory_condition + "/series_" + str(series)
	npath = userDirectory + experiment + "//"+ condition + "//" + "series_" + str(series)
	
	imp = IJ.openImage(directory_series + "/OMTC_" + str(series) + '.tif')
	imp.show()
	IJ.openImage(directory_series +"/RoiSet.zip")
	
	# ROI
	RM = RoiManager.getInstance()         
	rm = RM.getRoiManager()  
	
	# Iterating over ROI      
	roiCount = RM.getCount()
	print(roiCount)
	
	for i in range(1, roiCount+1): # Goes through the ROI manager
		
		# Prepare ROI 
		if not os.path.isdir(npath+"//bead_"+str(i)):
			bead_directory=os.makedirs(npath+"//bead_"+str(i))
		bead_directory=os.path.join(npath+"//bead_"+str(i))
		path = userDirectory + experiment + "//"+ condition + "//" + "series_" + str(series) + "//" + "bead_" + str(i) + "//"
		filename = condition + "_" + "series_" + str(series) + "_bead_" + str(i) + "_SPOTS"
		completename = os.path.join(path,  filename + ".csv")
		IJ.selectWindow("OMTC_"+str(series)+".tif")
		rm.select(i-1)
		IJ.run("Duplicate...", "duplicate")
		IJ.run("Split Channels")
	  	IJ.selectWindow("C1-OMTC_"+ str(series) + "-1.tif")
		C1_tiff = "C1_" + str(condition) + "_series_" + str(series) + "_bead_" + str(i) + ".tif"
		IJ.saveAs("Tiff", bead_directory + "/" + C1_tiff )
		IJ.getImage().close()	
		IJ.selectWindow("C2-OMTC_"+ str(series) + "-1.tif")
		IJ.run("Invert", "stack")
		IJ.run("Smooth", "stack")
		IJ.run("Enhance Contrast", "saturated=0.35")
		#IJ.run("Gaussian Blur...", "sigma=1 stack")
		
		# Create model object 
		model = Model()
		
		# Logger
		model.setLogger(Logger.IJ_LOGGER)
		
		# Create settigns object
		imp = IJ.openImage("C2-OMTC_"+ str(series) + "-1.tif")
		imp = WindowManager.getCurrentImage()
		settings = Settings(imp)
		print(settings)
		
		# Configure Detector
		settings.detectorFactory = LogDetectorFactory()
		settings.detectorSettings = {
		
			'DO_SUBPIXEL_LOCALIZATION' : True,
			'RADIUS' : 2.25, # this is radius, trackmate GUI asks for diameter
			'TARGET_CHANNEL' : 1,
			'THRESHOLD' : 0.0,
			'DO_MEDIAN_FILTERING' : True,
		}
		
		# Configure Spot filters
		filter1 = FeatureFilter('QUALITY', 3.5, True) # 3.5
		settings.addSpotFilter(filter1)
		filter2 = FeatureFilter('CONTRAST_CH1', .135, False) # .135 or .2
		settings.addSpotFilter(filter2)
		filter8 = FeatureFilter('POSITION_X', 5.0, True) # 5.0
		settings.addSpotFilter(filter8)
		filter9 = FeatureFilter('POSITION_Y', 15.0, False) # 15.0
		settings.addSpotFilter(filter9)
		
		# Configure Track filters
		filter3 = FeatureFilter('TRACK_DURATION', 5.0, True) # 5.0
		settings.addTrackFilter(filter3)
		filter4 = FeatureFilter('TRACK_MEAN_QUALITY', 4.0, True) # 4.0
		settings.addTrackFilter(filter4)
		filter5 = FeatureFilter('TRACK_MEAN_X', 14.0, False) # 14.0
		settings.addTrackFilter(filter5)
		filter6 = FeatureFilter('TRACK_MEAN_Y', 5.0, True) # 5.0
		settings.addTrackFilter(filter6)
		#filter7 = FeatureFilter('MAX_DISTANCE_TRAVELED', .8, True) 
		#settings.addTrackFilter(filter7)
	
		 
		# Configure Tracker
		settings.trackerFactory = SimpleSparseLAPTrackerFactory()
		settings.trackerSettings = LAPUtils.getDefaultLAPSettingsMap()
		settings.trackerSettings['LINKING_MAX_DISTANCE'] = 10.0
		settings.trackerSettings['GAP_CLOSING_MAX_DISTANCE'] = 15.0
		settings.trackerSettings['MAX_FRAME_GAP'] = 2
		settings.trackerSettings['ALLOW_TRACK_SPLITTING'] = False
		settings.trackerSettings['ALLOW_TRACK_MERGING'] = False
		settings.trackerSettings['ALLOW_GAP_CLOSING'] = True
	
		# All feature analyzers
		settings.addAllAnalyzers()
		
		# Instantiate plugin
		IJ.selectWindow("C2-OMTC_"+ str(series) + "-1.tif")
		trackmate = TrackMate(model, settings)
		
		# Check in
		ok = trackmate.checkInput()
		if not ok:
			sys.exit(str(trackmate.getErrorMessage()))
			 
		ok = trackmate.process()
		if not ok:
			sys.exit(str(trackmate.getErrorMessage()))
			    
		# Get results
		model.getLogger().log('Found ' + str(model.getTrackModel().nTracks(True)) + ' tracks.')
		sm = SelectionModel(model)
		ds = DisplaySettingsIO.readUserDefault()
		displayer =  HyperStackDisplayer(model, sm, imp, ds)
		displayer.render()
		displayer.refresh()
		fm = model.getFeatureModel()
		model.getLogger().log(str(model))
		
		# Export csv
		fieldnames = ['index','Time','x','y']
		with open(completename, 'w') as f:
			writer = csv.DictWriter(f, fieldnames = fieldnames)
			writer.writeheader()
		
			for id in model.getTrackModel().trackIDs(True):
			# Get all the spots of the current track.
				track = model.getTrackModel().trackSpots(id)
				for spot in track:
					sid = spot.ID()
					# Fetch spot features directly from spot.
					x = spot.getFeature('POSITION_X')
					y = spot.getFeature('POSITION_Y')
					t = spot.getFeature('POSITION_T')
					index = fm.getTrackFeature(id, 'TRACK_INDEX')
					writer.writerow({'index': index, 'Time': t, 'x' : x, 'y': y})
	                	
		rename_roi(bead_directory, condition, str(series) , str(i))
		IJ.selectWindow("C2-OMTC_"+ str(series) + "-1.tif")
		IJ.saveAs("Tiff", bead_directory + "/C2_" + condition + "_series_" + str(series) + "_bead_" + str(i) + ".tif")
		IJ.getImage().close() 
		WindowManager.getWindow("OMTC_" + str(series) + ".tif")
		
	WindowManager.getWindow("OMTC_" + str(series) + ".tif")
	IJ.getImage().close()
	
	rm.runCommand("Deselect"); # deselect ROIs to save them all
	rm.runCommand("Save", directory_series + "/RoiSet.zip")
	IJ.selectWindow('ROI Manager')
	IJ.run("Close")
	
IJ.run("Close All", "")

print('DONE')
        
        
	
        
        


