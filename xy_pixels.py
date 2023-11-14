from matplotlib import pyplot as plt
import numpy as np
import csv
from analysis_functions_xml import *
from well_plate_dictionary import *
from analysis_functions_image import *

plate = "1737"
well = "F4"
phase = "green"
frames = [0,216]




t0 = time.time()
# intensity = get_xy_pixels_all(plate, well, phase)
# intensity = get_xy_pixels_array(plate, well, phase)
im = get_xy_pixels(plate, well, phase, frame=72)

print(time.time()-t0)

display_image(im)