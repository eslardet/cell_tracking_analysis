import os
import numpy as np
import time
import math
from matplotlib import pyplot as plt
from analysis_functions_xml import *
from analysis_functions_image import *
from well_plate_dictionary import *


plate = 1737
well = "D6"
phase = "green"
frame = 123

im = get_intensity(plate, well, phase)[frame]
im_mean = im - np.mean(im)
im_scaled = im_mean/np.max(im_mean)
im_thresh = im_scaled > 0.75


display_image(im)
display_image(im_thresh)
display_image(im>0.25)