import numpy as np
# from analysis_functions_xml import *
from analysis_functions_image import *


plate = 1737
well = "B2"

slice_range = np.arange(1, 108, 5)

# t0 = time.time()
# plot_manders_vs_time(plate, well, slice_range)
# print(time.time()-t0)

for well in ["E5"]:
# for well in all_wells:
    t0 = time.time()
    plot_manders_vs_time(plate, well, slice_range)
    # print(get_manders_coeff(plate, well, -1))
    print(time.time()-t0)


cell_type = 'Microglia & T-cells'
group_list = ["Uninfected healthy control", "AC lPVL", "AC hPVL", "HAM"]
# group_list = ["Uninfected healthy control"]
stim_list = ["With stimulation"]
t_cell_list = ["Specific CD4"]

slice_compare = [1, 72]
# slice_compare = [1, 215]
# slice_compare = [1, 36]

# t0 = time.time()
# plot_manders_increase(cell_type, group_list, stim_list, t_cell_list, slice_compare, green_thresh=0.25, red_thresh=0.15, time_av=False, show_plot=False, save_plot=True)
# print(time.time()-t0)

# for cell_type in ["Astrocytes & T-cells", "Microglia & T-cells"]:
#     t0 = time.time()
#     for stim in all_stim:
#         for t_cell in all_t_cell:
#             stim_list = [stim]
#             t_cell_list = [t_cell]
#             for slice_compare in [[1, 36]]:
#                 t0 = time.time()
#                 plot_manders_increase(cell_type, group_list, stim_list, t_cell_list, slice_compare, green_thresh=0.25, red_thresh=0.15, time_av=False, show_plot=False, save_plot=True)
#                 print(time.time()-t0)

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