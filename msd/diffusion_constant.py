import os, sys, inspect
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
sys.path.insert(0, parentdir)
import numpy as np
import time
import math
from matplotlib import pyplot as plt
from analysis_functions_tracks import *
from well_plate_dictionary import *
import warnings


cell_type = 'Microglia & T-cells'
group = "AC lPVL"
stim = "With stimulation"
t_cell = "Specific CD4"

# warnings.filterwarnings("ignore")
# get_d_average(cell_type, group, stim, t_cell)

cell_type = 'Microglia & T-cells'
group_list = ["Uninfected healthy control", "AC lPVL", "AC hPVL", "HAM"]
stim_list = ["With stimulation"]
t_cell_list = ["Specific CD4"]

# plot_d_violin(cell_type, group_list, stim_list, t_cell_list, save_plot=True, show_plot=False, y_upper=True)
# plot_d_box(cell_type, group_list, stim_list, t_cell_list, save_plot=True, show_plot=False, y_upper=True)

# print(get_d_average(cell_type, group, stim, t_cell, average_type="mean_track"))

# print(get_d(plate="1723", well="C2"))
# print(get_d(plate="1737", well="C2"))

# print(get_wells(group="AC lPVL", stim="With stimulation", t_cell="Specific CD4"))

# for cell_type in ["Microglia & T-cells"]:
#     t0 = time.time()
#     for stim in all_stim:
#         for t_cell in all_t_cell:
#             stim_list = [stim]
#             t_cell_list = [t_cell]
#             plot_folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/plots/diffusion_const"
#             # plot_name = cell_type + "_" + stim + "_" + t_cell + ".png"
#             # if os.path.exists(os.path.join(plot_folder, plot_name)) == 0:
#             # plot_d_violin(cell_type, group_list, stim_list, t_cell_list, save_plot=True, show_plot=False, y_upper=True)
#             plot_d_box(cell_type, group_list, stim_list, t_cell_list, save_plot=True, show_plot=False, y_upper=True)
#             # else:
#             #     print("Already plotted " + plot_name)
#     print(time.time()-t0)




