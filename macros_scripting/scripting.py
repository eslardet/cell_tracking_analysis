import sys
import os
import time
import inspect

from ij import IJ
from ij import ImagePlus
from ij import WindowManager
import ij.plugin.ContrastEnhancer as ContrastEnhancer
 
from fiji.plugin.trackmate import Model
from fiji.plugin.trackmate import Settings
from fiji.plugin.trackmate import TrackMate
from fiji.plugin.trackmate import SelectionModel
from fiji.plugin.trackmate import Logger
from fiji.plugin.trackmate.detection import HessianDetectorFactory
from fiji.plugin.trackmate.tracking.jaqaman import SparseLAPTrackerFactory
from fiji.plugin.trackmate.gui.displaysettings import DisplaySettingsIO
import fiji.plugin.trackmate.visualization.hyperstack.HyperStackDisplayer as HyperStackDisplayer
import fiji.plugin.trackmate.features.FeatureFilter as FeatureFilter
import fiji.plugin.trackmate.action.ExportTracksToXML as ExportTracksToXML
##from fiji.plugin.trackmate.action import CloseGapsByLinearInterpolationAction
import java.io.File as File
import fiji.plugin.trackmate.action.AbstractTMAction
import fiji.plugin.trackmate.action.CaptureOverlayAction as CaptureOverlayAction
##import fiji.plugin.trackmate.action.CloseGapsByLinearInterpolationAction as CloseGapsByLinearInterpolationAction


##print(inspect.getmembers(fiji.plugin.trackmate.action))


def tracking(plate, phase, well, pic, inputFolder, outputFolder,
			radius, quality_threshold, max_link, gap_closing, max_frame_gap, max_gap, min_displacement, min_track_duration, max_track_start,
			display=False, save_xml=True):
	# We have to do the following to avoid errors with UTF8 chars generated in 
	# TrackMate that will mess with our Fiji Jython.
	reload(sys)
	sys.setdefaultencoding('utf-8')
	 
	# Get image
	filename = "VID" + plate + '_' + phase + '_' + well + '_' + pic
	image = inputFolder + plate + "/" + filename + '.tif'
	
	if os.path.exists(image) == False:
		print("Image does not exist!")
		print(image)
		return


	imp = IJ.openImage(image)
	title = imp.getTitle()
	print(title)
	#imp.show()
		
	folder = outputFolder + plate + "/"
	if not os.path.exists(folder):
		os.makedirs(folder)
	outFile = folder + filename + '.xml'
	##if os.path.exists(outFile):
	##	print("Already tracked!")
	##	return
	
	#----------------
	# Configure image
	#----------------
	# Swap Z and T dimensions if T=1
	dims = imp.getDimensions() # default order: XYCZT
	if dims[4] == 1:
		imp.setDimensions(dims[2], dims[4], dims[3])
	print(dims)
	## Resize stack to delete first slice and any number of the last slices
	# imp2 = imp.getStack()
	# imp.deleteSlice(1)
	# for i in range(imp.getSize()):
	# 	if i > end_frame-1:
	# 		imp.deleteLastSlice()
	# print(imp.getSize())

	# imp = ImagePlus(title, imp2)
	# print(imp.getTitle())
	# dims = imp.getDimensions()
	# print(dims)
	
	#----------------------------
	# Create the model object now
	#----------------------------
	 
	# Empty model
	model = Model()
	
	# Send all messages to ImageJ log window.
	model.setLogger(Logger.IJ_LOGGER)
	
	#------------------------
	# Prepare settings object
	#------------------------
	
	settings = Settings(imp)
	
	## Configure detector - We use the Strings for the keys
	settings.detectorFactory = HessianDetectorFactory()
	settings.detectorSettings = {
	    'DO_SUBPIXEL_LOCALIZATION' : True,
	    'NORMALIZE' : True,
	    'RADIUS' : radius,
	    'RADIUS_Z': 0.0,
	    'TARGET_CHANNEL' : 1,
	    'THRESHOLD' : quality_threshold
	}  
	
	## Configure spot filters - Classical filter on quality
	#filter1 = FeatureFilter('QUALITY', 30, True)
	#settings.addSpotFilter(filter1)
	 
	# Configure tracker - We want to allow merges and fusions
	settings.trackerFactory = SparseLAPTrackerFactory()
	settings.trackerSettings = settings.trackerFactory.getDefaultSettings()
	settings.trackerSettings['ALLOW_TRACK_SPLITTING'] = False
	settings.trackerSettings['ALLOW_TRACK_MERGING'] = False
	settings.trackerSettings['LINKING_MAX_DISTANCE'] = max_link
	settings.trackerSettings['ALLOW_GAP_CLOSING'] = gap_closing
	settings.trackerSettings['MAX_FRAME_GAP'] = max_frame_gap
	settings.trackerSettings['GAP_CLOSING_MAX_DISTANCE'] = max_gap
	
	# Add ALL the feature analyzers known to TrackMate. They will 
	# yield numerical features for the results, such as speed, mean intensity etc.
	settings.addAllAnalyzers()
	
	# Configure track filters
	 
	filter1 = FeatureFilter('TRACK_DISPLACEMENT', min_displacement, True)
	filter2 = FeatureFilter('TRACK_DURATION', min_track_duration, True)
	filter3 = FeatureFilter('TRACK_START', max_track_start, False)
	settings.addTrackFilter(filter1)
	settings.addTrackFilter(filter2)
	settings.addTrackFilter(filter3)
	
	# #-------------------
	# # Instantiate plugin
	# #-------------------
	 
	trackmate = TrackMate(model, settings)
	 
	# #--------
	# # Process
	# #--------
	 
	ok = trackmate.checkInput()
	if not ok:
		sys.exit(str(trackmate.getErrorMessage()))
	    
	ok = trackmate.process()
	if not ok:
	    sys.exit(str(trackmate.getErrorMessage()))
	
	#-----------------
	# Save to xml file
	#-----------------
	if save_xml == True:
		folder = outputFolder + plate + "/"
		if not os.path.exists(folder):
			os.makedirs(folder)
		outFile = File(folder + filename + '.xml')
		# selectionModel = SelectionModel( model )
		# ds = DisplaySettingsIO.readUserDefault()
		# Parent = imp.getParent()
		# CloseGaps.execute(model, selectionModel, ds, Parent)
		ExportTracksToXML.export(model, settings, outFile)
		
	if display == True:
		# A selection.
		selectionModel = SelectionModel( model )
		 
		# Read the default display settings.
		ds = DisplaySettingsIO.readUserDefault()
		 
		displayer =  HyperStackDisplayer( model, selectionModel, imp, ds )
		displayer.render()
		displayer.refresh()
		model.getLogger().log( str( model ) )

inputFolder = "/Volumes/T7 Shield/Incucyte_data/processed_data/"
outputFolder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/xml/"

radius = 3.0
quality_threshold = 0.53
max_link = 50.0
gap_closing = False
max_frame_gap = 0
max_gap = 0.0
min_displacement = 0.5
min_track_duration = 12.0
max_track_start = 200.0

display=False
save_xml=True

all_cells = ['B', 'C', 'D', 'E', 'F', 'G']

plate = '1738'
phase = 'green'
##well = 'B3'
pic = '1'

t0 = time.time()
for cell in ['G']:
	for i in range(2,12):
		t1 = time.time()
		well = cell + str(i)
		tracking(plate, phase, well, pic, inputFolder, outputFolder,
	   			radius, quality_threshold, max_link, gap_closing, max_frame_gap, max_gap, min_displacement, min_track_duration, max_track_start,
	     		display, save_xml)
		print(time.time() - t1)
print(time.time()-t0)
