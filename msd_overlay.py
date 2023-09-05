import os
import numpy as np
import time
from matplotlib import pyplot as plt
from analysis_functions_xml import *
from well_plate_dictionary import *

cell_type = 'Astrocytes & T-cells'
group_list = ["Uninfected healthy control"]
stim_list = ["With stimulation"]
t_cell_list = ["Non-specific CD4"]

phase = 'green'
log = False
remove_outliers = False

def plot_msd_overlay(cell_type, group_list, stim_list, t_cell_list, phase='green', log=False, remove_outliers=False):

    fig, ax = plt.subplots()

    plate_list = get_plates(cell_type)
    for plate in plate_list:
        for group in group_list:
            for stim in stim_list:
                for t_cell in t_cell_list:
                    well_list = get_wells(group, stim, t_cell)
                    for well in well_list:
                        file_path = get_xml_file(plate, well, phase)
                        if os.path.exists(file_path):
                            msd = get_msd(file_path, max_frames=200, plot=False, remove_outliers=remove_outliers)
                            plot_label = get_well_info(well)[0] + ", " + get_well_info(well)[1] + ", " + get_well_info(well)[2]
                            plot_label = plate + ", " + well
                            if log == True:
                                ax.loglog(np.arange(len(msd))/3, msd, label=plot_label)
                            else:
                                ax.plot(np.arange(len(msd))/3, msd, label=plot_label)

    ax.set_xlabel(r"$\tau$ (hours)")
    ax.set_ylabel("MSD (pixels)")
    plt.legend()
    plt.show()

def plot_msd_all_tracks(plate, well):
    file_path = get_xml_file(plate, well)
    msd_all = get_msd_individual_tracks(file_path)
    fig, ax = plt.subplots()
    for msd in msd_all:
        n_frames = len(msd)
        max_t = n_frames // 3
        # ax.plot(np.arange(len(msd))/3, msd)
        ax.plot(np.arange(max_t)/3, msd[:max_t])

    ax.set_xlabel(r"$\tau$ (hours)")
    ax.set_ylabel("MSD (pixels)")
    plt.show()

# plot_msd_overlay(cell_type, group_list, stim_list, t_cell_list, phase, log, remove_outliers=False)
# plot_msd_overlay(cell_type, group_list, stim_list, t_cell_list, phase, log, remove_outliers=True)


plot_msd_all_tracks(plate='1736', well='F4')