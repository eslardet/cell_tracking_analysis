import os
import numpy as np
import time
from matplotlib import pyplot as plt
from analysis_functions_xml import *
from well_plate_dictionary import *
import warnings

cell_type = 'Astrocytes & T-cells'
group_list = ["Uninfected healthy control"]
stim_list = ["With stimulation"]
t_cell_list = ["Non-specific CD4"]

phase = 'green'

def get_d(plate, well):
    file_path = get_xml_file(plate, well)
    msd = get_msd(file_path)
    n_frames = len(msd)
    max_t = n_frames // 3
    d = np.polyfit(np.arange(max_t)/3, msd[:max_t], 1)[0]/4  
    return d

def get_d_all_tracks(plate, well):
    file_path = get_xml_file(plate, well)
    msd_all = get_msd_individual_tracks(file_path)
    # fig, ax = plt.subplots()
    all_d = []
    for msd in msd_all:
        n_frames = len(msd)
        if n_frames > 12:
            max_t = n_frames // 3
            d = np.polyfit(np.arange(max_t)/3, msd[:max_t], 1)[0]/4
            all_d.append(d)
    return all_d

def get_d_average(cell_type, group, stim, t_cell, average_type="both"):
    """
    average_type: "all_track" or "mean_track" or "both"
    """
    plate_list = get_plates(cell_type)
    well_list = get_wells(group, stim, t_cell)
    all_d = []
    d_av = []
    for plate in plate_list:
        for well in well_list:
            file_path = get_xml_file(plate, well)
            if os.path.exists(file_path):
                all_d += get_d_all_tracks(plate, well)
                d_av.append(get_d(plate, well))
    if average_type == "all_track":
        return np.mean(all_d)
    elif average_type == "mean_track":
        return np.mean(d_av)
    else:
        return np.mean(all_d), np.mean(d_av)


def set_axis_style(ax, labels):
    ax.get_xaxis().set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_xlim(0.25, len(labels) + 0.75)

def plot_d_violin(cell_type, group_list, stim_list, t_cell_list, save_plot=False, show_plot=True):
    fig, ax = plt.subplots(figsize=(12,6))

    plot_d = []
    x_ticks = []
    d_mean = []
    plate_list = get_plates(cell_type)
    for group in group_list:
        for stim in stim_list:
            for t_cell in t_cell_list:
                well_list = get_wells(group, stim, t_cell)
                all_d = []
                d_av = []
                for plate in plate_list:
                    for well in well_list:
                        file_path = get_xml_file(plate, well)
                        if os.path.exists(file_path):
                            d_tracks = get_d_all_tracks(plate, well)
                            all_d += d_tracks
                            d_av.append(get_d(plate, well))
                plot_d.append(all_d)
                x_ticks.append(group + "\n " + stim + "\n " + t_cell)

                d_mean.append(np.mean(d_av))
    ax.violinplot(plot_d, showmeans=False, showextrema=True, showmedians=False, points=500)
    ax.scatter(np.arange(1, len(d_mean)+1), d_mean)
    ax.set_title(cell_type)
    set_axis_style(ax,x_ticks)
    ax.set_ylim(0,200)
    ax.set_ylabel("Diffusion constant (pixels^2/hour)")

    if save_plot == True:  
        plot_folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/plots/diffusion_const"
        plot_name = cell_type + "_" + stim + "_" + t_cell + ".png"
        plt.savefig(os.path.join(plot_folder, plot_name))

    if show_plot == True:
        plt.show()


cell_type = 'Astrocytes & T-cells'
group = "Uninfected healthy control"
stim = "With stimulation"
t_cell = "Non-specific CD4"

# warnings.filterwarnings("ignore")
# get_d_average(cell_type, group, stim, t_cell)

cell_type = 'Microglia & T-cells'
group_list = ["Uninfected healthy control", "AC lPVL", "AC hPVL", "HAM"]
stim_list = ["With stimulation"]
t_cell_list = ["Specific CD4"]

for stim in all_stim:
    for t_cell in all_stim:
        stim_list = [stim]
        t_cell_list = [t_cell]
        plot_d_violin(cell_type, group_list, stim_list, t_cell_list, save_plot=True, show_plot=False)
