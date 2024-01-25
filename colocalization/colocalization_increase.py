import os, sys, inspect
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
sys.path.insert(0, parentdir)

import numpy as np
# from analysis_functions_xml import *
from analysis_functions_image import *



cell_type = 'Microglia & T-cells'
group_list = ["Uninfected healthy control", "AC lPVL", "AC hPVL", "HAM"]
# group_list = ["Uninfected healthy control"]
stim = "With stimulation"
t_cell = "Specific CD4"

plate_list = [1728, 1737]
stim_list = ["No stimulation", "With stimulation"]
t_cell_list = ["Non-specific CD4", "Non-specific CD8", "Specific CD8"]

# slice_compare = [3, 108]
# slice_compare = [3, 216]
# slice_compare = [3, 36]

# rescale_list = ["1737 C4", "1737 G2"]


# t0 = time.time()
# plot_manders_increase(cell_type, plate_list, group_list, stim, t_cell, slice_compare, rescale_list,
#                       time_av=False, show_plot=False, save_plot=False, save_data=True)
# # # plot_manders_double_time(cell_type, group_list, stim_list, t_cell_list, plate_list, slice_range)
# # plot_manders_exponent(cell_type, plate_list, group_list, stim_list, t_cell_list, slice_range)
# print(time.time()-t0)

for cell_type in ["Microglia & T-cells"]:
    for stim in stim_list:
        for t_cell in t_cell_list:
            for slice_compare in [[3, 36], [3, 72], [3, 108], [3, 144], [3, 180], [3, 216]]:
                print(cell_type, stim, t_cell, slice_compare)
                t0 = time.time()
                plot_manders_increase(cell_type, plate_list, group_list, stim, t_cell, slice_compare,
                                     green_thresh_default=0.75, time_av=False, show_plot=False, save_plot=False, save_data=True)
                # plot_manders_exponent(cell_type, plate_list, group_list, stim_list, t_cell_list, slice_range)
                print(time.time()-t0)

# plate = 1728
# well = "B2"
# slice = 1
# im_red = get_image(plate, well, "red")[slice]
# im_green = get_image(plate, well, "green")[slice]
# red_thresh = 0.15
# # ## First segment with thresholding
# im_red_thresh = im_red/np.max(im_red) > red_thresh

# io.imshow(im_red_thresh)
# plt.show()