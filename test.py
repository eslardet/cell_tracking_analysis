import numpy as np
from analysis_functions_tracks import *
from analysis_functions_image import *
from skimage import io
from matplotlib import pyplot as plt
from matplotlib import cm

plate = 1728
phase = "green"
well = "F4"
frame = 100

# folder = "/Volumes/T7 Shield/Incucyte_data/raw_data/exp_109/"
# filename = folder + "VID" + str(plate) + '_' + str(phase) + "_" + well + '_1.tif'
# im_all = io.imread(filename)
# im = im_all[frame]
# # io.imshow(im)
# io.imshow(im/np.max(im), cmap=cm.gray)
# plt.show()

# im = get_image(plate, well, phase)[frame]
# display_image(im)

print([[i,i+12] for i in range(0,72,12)])