# from PIL import Image, ImageEnhance
from skimage import io
from matplotlib import pyplot as plt
import numpy as np
from analysis_functions_xml import *
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
plate_list = ["1723"]

phase = "red"
frame = 215
r_min = 10
r_max = 50
r_bin_num = 25

# plot_corr_exp_scatter(cell_type, plate_list, group_list, stim_list, t_cell_list, phase, frame, r_min, r_max, r_bin_num, save_plot=True, show_plot=False)

# for cell_type in ["Microglia & T-cells"]:
#     t0 = time.time()
#     for stim in all_stim:
#         for t_cell in all_t_cell:
#             stim_list = [stim]
#             t_cell_list = [t_cell]
#             plot_corr_exp_scatter(cell_type, plate_list, group_list, stim_list, t_cell_list, phase, frame, r_min, r_max, r_bin_num, save_plot=True, show_plot=False)
#     print(time.time()-t0)


plate = 1737
well = "E5"
phase = "red"
threshold = False
frames = [1,72,144,216]

# im = get_image(plate, well, phase)[-1]
# im_thresh = im/np.max(im) > 0.99
# display_image(im_thresh)

# plot_corr_log(plate, well, phase, frames=frames, r_min=10, r_max=20, r_bin_num=10, threshold=threshold, corr_av=True, save_plot=True, show_plot=False, x_lim=True)
plot_exponents_time(plate, well, phase, frames=np.arange(1, 216, 5), r_min=10, r_max=50, threshold=threshold, corr_av=True, save_plot=True, show_plot=False)

# for well in all_wells:
#     path = get_tif_file(plate, well, phase)
#     if os.path.exists(path):
#         path_plot = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/plots/correlation_functions/" + str(plate) + "/" + str(plate) + '_' + well + '_log.png'
#         if not os.path.exists(path_plot):
#             print(well)
#             # plot_corr_log(plate, well, phase, frames=[1,216], r_min=10, r_max=50, r_bin_num=25, corr_av=True, save_plot=True, show_plot=False)
#             # plot_corr_log(plate, well, phase, frames=[1,216], r_min=0, r_max=5, r_bin_num=5, corr_av=True, save_plot=True, show_plot=False, x_lim=True)

    # plot_exponents_time(plate, well, phase, frames=np.arange(1, 212, 10), r_min=10, r_max=20, corr_av=True, save_plot=True, show_plot=False)


## Green
# for well in all_wells:
#     path = get_tif_file(plate, well, phase)
#     if os.path.exists(path):
#         path_plot = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/plots/correlation_functions/" + str(plate) + "/" + str(plate) + '_' + well + '_log.png'
#         if not os.path.exists(path_plot):
#             print(well)
#         else:
#             plot_exponents_time(plate, well, phase, frames=np.arange(1, 212, 10), r_min=10, r_max=20, threshold=threshold, corr_av=True, save_plot=True, show_plot=False)

