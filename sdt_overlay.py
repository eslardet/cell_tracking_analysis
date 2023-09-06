import os
import numpy as np
import time
from matplotlib import pyplot as plt
from analysis_functions_xml import *
from well_plate_dictionary import *



def plot_sdt_overlay(cell_type, group_list, stim_list, t_cell_list, phase='green', log=False, max_frames=200):

    fig, ax = plt.subplots()

    plate_list = get_plates(cell_type)
    for group in group_list:
        for stim in stim_list:
            for t_cell in t_cell_list:
                for plate in plate_list:
                    well_list = get_wells(group, stim, t_cell)
                    for well in well_list:
                        file_path = get_xml_file(plate, well, phase)
                        if os.path.exists(file_path):
                            sdt = get_sdt(file_path, max_frames=max_frames, plot=False)
                            plot_label = get_well_info(well)[0] + ", " + get_well_info(well)[1] + ", " + get_well_info(well)[2]
                            # plot_label = plate + ", " + well
                            if log == True:
                                ax.loglog(np.arange(len(sdt))/3, sdt, label=plot_label)
                            else:
                                ax.plot(np.arange(len(sdt))/3, sdt, label=plot_label)

    ax.set_xlabel(r"$\tau$ (hours)")
    ax.set_ylabel("SDT (pixels^2)")
    plt.legend()
    plt.show()

def plot_sdt_all_tracks(plate, well):
    file_path = get_xml_file(plate, well)
    sdt_all = get_sdt_individual_tracks(file_path)
    fig, ax = plt.subplots()
    for sdt in sdt_all:
        n_frames = len(sdt)
        max_t = n_frames // 3
        # ax.plot(np.arange(max_t)/3, sdt[:max_t])

        ax.plot(np.arange(len(sdt))/3, sdt/3)

    ax.set_xlabel(r"$\tau$ (hours)")
    ax.set_ylabel("SDT (pixels^2)")
    plt.show()

cell_type = 'Microglia & T-cells'
group_list = ["Uninfected healthy control", "HAM"]
stim_list = ["With stimulation"]
t_cell_list = ["Non-specific CD4"]

phase = 'green'
log = False
remove_outliers = False

plot_sdt_overlay(cell_type, group_list, stim_list, t_cell_list, phase, log, max_frames=72)


# plot_sdt_all_tracks(plate='1736', well='F2')