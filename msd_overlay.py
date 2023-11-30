import os
import numpy as np
import time
from matplotlib import pyplot as plt
from analysis_functions_tracks import *
from well_plate_dictionary import *



def plot_msd_overlay(cell_type, group_list, stim_list, t_cell_list, phase='green', log=False, remove_outliers=False):

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
                            # msd = get_msd(file_path, max_frames=200, plot=False, remove_outliers=remove_outliers)
                            msd = get_msd_filter(file_path)
                            plot_label = get_well_info(well)[0] + ", " + get_well_info(well)[1] + ", " + get_well_info(well)[2]
                            plot_label = plate + ", " + well
                            if log == True:
                                ax.loglog(np.arange(len(msd))/3, msd, label=plot_label)
                            else:
                                ax.plot(np.arange(len(msd))/3, msd, label=plot_label)
                                # ax.plot(np.arange(len(msd2))/3, msd2, label="filtered")
                            msd_1 = msd[3:15]
                            fit_1 = np.polyfit(np.log(np.arange(len(msd))[3:15]/3), np.log(msd_1), 1)
                            msd_2 = msd[15:90]
                            fit_2 = np.polyfit(np.log(np.arange(len(msd))[15:90]/3), np.log(msd_2), 1)
                            print(fit_1[0], fit_2[0])
    ax.plot(np.arange(217)/3, np.arange(217)/3 *100)
    # ax.plot(np.arange(217)/3, np.arange(217)**2/3 *100)
    ax.set_xlabel(r"$\tau$ (hours)")
    ax.set_ylabel("MSD (pixels^2)")
    plt.legend()
    plt.show()

def plot_msd_all_tracks(plate, well, show_mean=False):
    file_path = get_xml_file(plate, well)
    msd_all = get_msd_individual_tracks(file_path)
    fig, ax = plt.subplots()
    for msd in msd_all:
        n_frames = len(msd)
        max_t = n_frames // 3
        # ax.plot(np.arange(max_t)/3, msd[:max_t])
        ax.plot(np.arange(len(msd))/3, msd)
    if show_mean == True:
        msd, sd = get_msd(file_path, max_frames=200, plot=False, remove_outliers=remove_outliers, return_sd=True)
        msd = np.array(msd)
        ax.plot(np.arange(len(msd))/3, msd, color='k') # Mean MSD track
        ax.fill_between(np.arange(len(msd))/3, msd-sd, msd+sd, color='k', alpha=0.2) # STD as conifidence interval
    # ax.plot(np.arange(217)/3, np.arange(217)/3 *100, 'k')
    ax.set_xlabel(r"$\tau$ (hours)")
    ax.set_ylabel("MSD (pixels^2)")
    # ax.set_ylim(0,50000)
    plt.show()

def plot_msd_all_tracks_subset(plate, well, show_mean=False, d_threshold=400):
    file_path = get_xml_file(plate, well)
    msd_all = get_msd_individual_tracks(file_path)
    fig, ax = plt.subplots()
    for msd in msd_all:
        n_frames = len(msd)
        max_t = n_frames // 3
        # ax.plot(np.arange(max_t)/3, msd[:max_t])
        if n_frames > 20:
            min_frame = 15
            if n_frames < 72:
                max_frame = n_frames
            else:
                max_frame = 72
            d = np.polyfit(np.arange(max_frame-min_frame)/3, msd[min_frame:max_frame], 1)[0]/4
            if d > d_threshold:
                ax.plot(np.arange(len(msd))/3, msd, color='r')
            else:
                ax.plot(np.arange(len(msd))/3, msd, color='b')
    if show_mean == True:
        msd, sd = get_msd(file_path, max_frames=200, plot=False, remove_outliers=remove_outliers, return_sd=True)
        msd = np.array(msd)
        ax.plot(np.arange(len(msd))/3, msd, color='k') # Mean MSD track
        ax.fill_between(np.arange(len(msd))/3, msd-sd, msd+sd, color='k', alpha=0.2) # STD as conifidence interval
    # ax.plot(np.arange(217)/3, np.arange(217)/3 *100, 'k')
    ax.set_xlabel(r"$\tau$ (hours)")
    ax.set_ylabel("MSD (pixels^2)")
    # ax.set_ylim(0,50000)
    plt.show()

cell_type = 'Microglia & T-cells'
group_list = ["Uninfected healthy control"]
# group_list = all_groups
stim_list = ["With stimulation"]
t_cell_list = ["Specific CD4"]

phase = 'green'
log = True
remove_outliers = False

# plot_msd_overlay(cell_type, group_list, stim_list, t_cell_list, phase, log, remove_outliers)


# plot_msd_all_tracks(plate='1727', well='B2', show_mean=True)
plot_msd_all_tracks_subset(plate='1723', well='B3', show_mean=False, d_threshold=400)

# plot_msd_all_tracks_subset(plate='1737', well='C2', show_mean=False, d_threshold=400)