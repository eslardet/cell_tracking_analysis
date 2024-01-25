import os, sys, inspect
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
sys.path.insert(0, parentdir)

# from PIL import Image, ImageEnhance
from skimage import io
from matplotlib import pyplot as plt
import numpy as np
from analysis_functions_tracks import *
from well_plate_dictionary import *
from analysis_functions_image import *
import os
import time
from scipy.signal import correlate, correlate2d, fftconvolve
from scipy.fftpack import fft2, ifft2
from scipy.spatial import distance_matrix
from scipy.spatial.distance import cdist
from scipy.optimize import curve_fit


cell_type = 'Microglia & T-cells'
group_list = ["Uninfected healthy control", "AC lPVL", "AC hPVL", "HAM"]
stim_list = ["With stimulation"]
t_cell_list = ["Specific CD4"]
plate_list = ["1737"]

phase = "green"
clean_green = True
frame = 108
frame_list = np.arange(108,215,9)
r_min = 5
r_max = 10
r_bin_num = 10
ylim = 40


# plot_corr_exp_scatter(cell_type, plate_list, group_list, stim_list, t_cell_list, phase, frame, r_min, r_max, r_bin_num, save_plot=True, show_plot=False)

# for cell_type in ["Microglia & T-cells"]:
#     t0 = time.time()
#     for stim in all_stim:
#         for t_cell in all_t_cell:
#             stim_list = [stim]
#             t_cell_list = [t_cell]
#             # plot_corr_exp_scatter(cell_type, plate_list, group_list, stim_list, t_cell_list, phase, frame, r_min, r_max, r_bin_num, ylim=ylim, 
#             #                       save_plot=True, show_plot=False, clean_green=True)
#             plot_corr_exp_scatter_av(cell_type, plate_list, group_list, stim_list, t_cell_list, phase, frame_list, r_min, r_max, r_bin_num, ylim=ylim, 
#                                      save_plot=True, show_plot=False, clean_green=True)
#     print(time.time()-t0)


plate = 1737
phase = "green"

well = "B2"
threshold = True
# frames = [108,215]
# frames = np.arange(5,216,10)

# im = get_image(plate, well, phase)[-1]
# im_thresh = im/np.max(im) > 0.99
# display_image(im_thresh)

# plot_corr_log(plate, well, phase, frames=frames, r_min=10, r_max=50, r_bin_num=25, threshold=threshold, corr_av=True, save_plot=True, show_plot=False, x_lim=True)
plot_exponents_time(plate, well, phase, frames=np.arange(5, 216, 10), r_min=10, r_max=20, threshold=threshold, corr_av=True, save_plot=True, show_plot=False)

# for well in all_wells:
#     path = get_tif_file(plate, well, phase)
#     if os.path.exists(path):
#         path_plot = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/plots/correlation_functions/" + str(plate) + "/" + str(plate) + '_' + well + '_log.png'
#         # if not os.path.exists(path_plot):
#         #     print(well)
#         # else:
#         if phase == "red":
#             plot_exponents_time(plate, well, phase, frames=np.arange(5, 216, 10), r_min=10, r_max=50, threshold=False, corr_av=True, save_plot=True, show_plot=False)
#         elif phase == "green":
#             plot_exponents_time(plate, well, phase, frames=np.arange(5, 216, 10), r_min=5, r_max=10, threshold=True, corr_av=True, save_plot=True, show_plot=False)

