import os
import numpy as np
import time
import math
from matplotlib import pyplot as plt
from analysis_functions_tracks import *
from analysis_functions_image import *
from well_plate_dictionary import *


plate = 1737
well = "D6"
phase = "green"
frame = 200

im = get_intensity(plate, well, phase)[frame]
im_mean = im - np.mean(im)
im_scaled = im_mean/np.max(im_mean)
im_thresh = im_scaled > 0.75

# im_thresh = im > 0.25

# print(sum(im_thresh_flat))
# display_image(im)
display_image(im_thresh)
# display_image(im>0.25)

# display_image(np.multiply(im, im_thresh))


# red_thresh = 0.15
# green_thresh = 0.75
# slice = 108
# im_red = get_image(plate, well, "red")[slice]
# im_green = get_image(plate, well, "green")[slice]

# ## First segment with thresholding
# im_red_thresh = im_red/np.max(im_red) > red_thresh
# # im_green_thresh = im_green/np.max(im_green) > green_thresh
# im_mean = im_green - np.mean(im_green)
# im_scaled = im_mean/np.max(im_mean)
# im_green_thresh = im_scaled > green_thresh


# coeff = measure.manders_coloc_coeff(im_green/np.max(im_green), im_red_thresh)
# print(coeff)