import os, sys, inspect
parentdir = os.path.diranme(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
sys.path.insert(0, parentdir)
import numpy as np
import time
import math
from matplotlib import pyplot as plt
from analysis_functions_tracks import *
from well_plate_dictionary import *
import warnings
# from diffusion_constant import get_d

def get_d(plate, well, d_threshold=400, d_negative=True):
    file_path = get_xml_file(plate, well)
    msd = get_msd_filter(file_path, d_threshold=d_threshold, d_negative=d_negative)
    n_frames = len(msd)
    # max_t = n_frames // 3
    min_frame = 15
    max_frame = 72
    d = np.polyfit(np.arange(max_frame-min_frame)/3, msd[min_frame:max_frame], 1)[0]/4  
    return d

def get_d_all_tracks(plate, well, d_threshold=400, d_negative=True):
    file_path = get_xml_file(plate, well)
    msd_all = get_msd_individual_tracks(file_path)
    # fig, ax = plt.subplots()
    all_d = []
    for msd in msd_all:
        n_frames = len(msd)
        if n_frames > 20:
            # max_t = n_frames // 3
            min_frame = 15
            if n_frames < 72:
                max_frame = n_frames
            else:
                max_frame = 72
            d = np.polyfit(np.arange(max_frame-min_frame)/3, msd[min_frame:max_frame], 1)[0]/4
            if d < d_threshold:
                if d_negative == False:
                    if d > 0:
                        all_d.append(d)
                else:
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
                print(plate, well, get_d(plate, well))
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

def plot_d_violin(cell_type, group_list, stim_list, t_cell_list, save_plot=False, show_plot=True, y_upper=None):
    fig, ax = plt.subplots(figsize=(12,6))

    plot_d = []
    x_ticks = []
    d_mean = []
    d_sd = []
    plot_d_mean = []
    plot_d_sd = []
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
                print(len(all_d), group, stim, t_cell)
                plot_d.append(all_d)
                plot_d_mean.append(np.mean(all_d))
                plot_d_sd.append(np.std(all_d))
                x_ticks.append(group + "\n " + stim + "\n " + t_cell)
                d_mean.append(np.mean(d_av))
                d_sd.append(np.std(d_av))
    ax.violinplot(plot_d, showmeans=False, showextrema=True, showmedians=True, points=500)
    # ax.errorbar(np.arange(1, len(plot_d_mean)+1), plot_d_mean, yerr=plot_d_sd, fmt='o', capsize=3)
    # ax.scatter(np.arange(1, len(d_mean)+1), d_mean)
    ax.errorbar(np.arange(1, len(d_mean)+1), d_mean, yerr=d_sd, fmt='o', capsize=3, color='k')
    ax.set_title(cell_type)
    set_axis_style(ax,x_ticks)
    if y_upper == True:
        max_d = np.max(d_mean)
        y_upper = int(math.ceil(max_d/100)) * 100
        ax.set_ylim(0,y_upper)
    ax.set_ylabel("Diffusion constant (pixels^2/hour)")

    if save_plot == True:  
        plot_folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/plots/diffusion_const_violin"
        plot_name = cell_type + "_" + stim + "_" + t_cell + ".png"
        plt.savefig(os.path.join(plot_folder, plot_name))

    if show_plot == True:
        plt.show()


def plot_d_box(cell_type, group_list, stim_list, t_cell_list, d_threshold=400, d_negative=True, save_plot=False, show_plot=True, y_upper=False):
    fig, ax = plt.subplots(figsize=(12,6))

    plot_d = []
    x_ticks = []
    d_mean = []
    d_sd = []
    plot_d_mean = []
    plot_d_sd = []
    plate_list = get_plates(cell_type)
    pos = 1
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
                            # Individual tracks
                            d_tracks = get_d_all_tracks(plate, well, d_threshold, d_negative)
                            all_d += d_tracks

                            # Mean tracks
                            d_av.append(get_d(plate, well, d_threshold, d_negative=True))
                            d_c = get_d(plate, well, d_threshold, d_negative=True)
                            ax.scatter(pos, d_c, color='k')        
                print(len(all_d), group, stim, t_cell)

                # Individual tracks
                plot_d.append(all_d)
                # plot_d_mean.append(np.mean(all_d))
                # plot_d_sd.append(np.std(all_d))

                # Mean tracks
                d_mean.append(np.mean(d_av))
                d_sd.append(np.std(d_av))

                x_ticks.append(group + "\n " + stim + "\n " + t_cell)
                pos += 1
    ax.boxplot(plot_d, notch=True)
    # ax.errorbar(np.arange(1, len(plot_d_mean)+1), plot_d_mean, yerr=plot_d_sd, fmt='o', capsize=3)
    # ax.scatter(np.arange(1, len(d_mean)+1), d_mean)
    # ax.errorbar(np.arange(1, len(d_mean)+1), d_mean, yerr=d_sd, fmt='o', capsize=3, color='k')
    ax.set_title(cell_type)
    set_axis_style(ax,x_ticks)
    if y_upper == True:
        # max_d = np.max(d_mean)
        max_d = d_threshold / 2
        y_upper = int(math.ceil(max_d/100)) * 100
        ax.set_ylim(0,y_upper)
        # ax.set_ylim(0,d_threshold)
    ax.set_ylabel("Diffusion constant (pixels^2/hour)")

    if save_plot == True:  
        plot_folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/plots/diffusion_const_box/filter_" + str(d_threshold)
        if not os.path.exists(plot_folder):
            os.makedirs(plot_folder)
        plot_name = cell_type + "_" + stim + "_" + t_cell + ".png"
        plt.savefig(os.path.join(plot_folder, plot_name))

    if show_plot == True:
        plt.show()

cell_type = 'Microglia & T-cells'
group = "AC lPVL"
stim = "With stimulation"
t_cell = "Specific CD4"

# warnings.filterwarnings("ignore")
# get_d_average(cell_type, group, stim, t_cell)

cell_type = 'Microglia & T-cells'
group_list = ["Uninfected healthy control", "AC lPVL", "AC hPVL", "HAM"]
# group_list = ["Uninfected healthy control"]
stim_list = ["With stimulation"]
t_cell_list = ["Specific CD4"]

# plot_d_violin(cell_type, group_list, stim_list, t_cell_list, save_plot=True, show_plot=False, y_upper=False)
# plot_d_box(cell_type, group_list, stim_list, t_cell_list, d_threshold=400, d_negative=False, save_plot=True, show_plot=False, y_upper=True)

# print(get_d_average(cell_type, group, stim, t_cell, average_type="mean_track"))

# d_threshold = 100
# print(get_d(plate="1723", well="C2", d_threshold=d_threshold))
# print(get_d(plate="1737", well="C2", d_threshold=d_threshold))

# print(get_wells(group="AC lPVL", stim="With stimulation", t_cell="Specific CD4"))

for cell_type in ["Astrocytes & T-cells"]:
    t0 = time.time()
    for stim in all_stim:
        for t_cell in all_t_cell:
            stim_list = [stim]
            t_cell_list = [t_cell]
            # plot_folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/plots/diffusion_const/middle_of_track_5-24hr"
            # plot_name = cell_type + "_" + stim + "_" + t_cell + ".png"
            # if os.path.exists(os.path.join(plot_folder, plot_name)) == 0:
            # plot_d_violin(cell_type, group_list, stim_list, t_cell_list, save_plot=True, show_plot=False, y_upper=True)
            plot_d_box(cell_type, group_list, stim_list, t_cell_list, d_threshold=400, d_negative=False, save_plot=True, show_plot=False, y_upper=True)
            # else:
            #     print("Already plotted " + plot_name)
    print(time.time()-t0)


plate = "1737"
well = "C4"

# for well in all_wells:
#     d_all = get_d_all_tracks(plate, well, d_threshold=100000)

#     # print(len(d_all))
#     # # print(np.median([i for i in d_all if i>100]))
#     # # print(np.max(d_all))
#     print(well, sum(i<0 for i in d_all)/len(d_all))
