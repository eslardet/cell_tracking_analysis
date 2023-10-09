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



def plot_corr_log(plate, well, phase, frames, r_min, r_max, r_bin_num=100, corr_av=True, save_plot=False, show_plot=True, x_lim=False):
    im = get_image(plate, well, phase)
    
    fig, ax = plt.subplots()

    for i in frames:
        radius, corr = get_corr_binned(im[i], corr_av=corr_av, corr_r_min=0, corr_r_max=r_max, r_bin_num=r_bin_num)

        ax.plot(radius, corr, label="Frame=" + str(i))

        # xi, coeff = fit_log_line(radius, corr, 0, r_min)
        # print(xi, coeff)
        # r_plot = np.linspace(0, r_min, 100)
        # ax.plot(r_plot, func(r_plot, xi, coeff), '--', label=r"$\xi_0=$" + str(round(xi,2)))

        xi, coeff = fit_log_line(radius, corr, r_min, r_max)
        print(xi, coeff)
        r_plot = np.linspace(r_min, r_max, 100)
        ax.plot(r_plot, exp_func(r_plot, xi, coeff), '--', label=r"$\xi=$" + str(round(xi,2)))

        ax.set_yscale('log')

    ax.set_xlabel("r (pixels)")
    ax.set_ylabel("C(r)")
    # ax.set_ylim([0,1])
    if x_lim == True:
        ax.set_xlim(0,r_max)
    else:
        ax.set_xlim(0,50)
    ax.legend()
    # ax.set_ylim(10**-2, 1)
    if save_plot == True:
        plot_folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/plots/correlation_functions/" + str(plate) + "_" + phase + "/"
        if not os.path.exists(plot_folder):
            os.makedirs(plot_folder)
        plt.savefig(plot_folder + str(plate) + '_' + well + '_log.png', bbox_inches='tight')

    if show_plot == True:
        plt.show()

    plt.close()


def plot_exponents_time(plate, well, phase, frames, r_min, r_max, corr_av=True, save_plot=False, show_plot=True):
    im_all = get_image(plate, well, phase)
    
    exp = []
    for i in frames:
        im = im_all[i]
        exp.append(get_exponent(im, r_min, r_max, corr_av))

    fig, ax = plt.subplots()
    ax.plot(frames/3, exp, '-o')
    ax.set_xlabel("Time (hours)")
    ax.set_ylabel(r"Correlation length $\xi$")

    if save_plot == True:
        plot_folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/plots/correlation_exponents_time/" + str(plate) + "_" + phase + "/"
        if not os.path.exists(plot_folder):
            os.makedirs(plot_folder)
        plt.savefig(plot_folder + str(plate) + '_' + well + '.png', bbox_inches='tight')

    if show_plot == True:
        plt.show()

    plt.close()


def plot_corr_exp_scatter(cell_type, plate_list, group_list, stim_list, t_cell_list, phase, frame, r_min, r_max, r_bin_num=100, save_plot=False, show_plot=True):
    fig, ax = plt.subplots(figsize=(12,6))

    x_ticks = []
    exp_mean = []
    exp_sd = []

    pos = 1
    for group in group_list:
        for stim in stim_list:
            for t_cell in t_cell_list:
                well_list = get_wells(group, stim, t_cell)
                all_exp = []
                for plate in plate_list:
                    for well in well_list:
                        file_path = get_tif_file(plate, well, phase)
                        if os.path.exists(file_path):
                            im = get_image(plate, well, phase)
                            radius, corr = get_corr_binned(im[frame], corr_r_min=r_min, corr_r_max=r_max, r_bin_num=r_bin_num)
                            xi, coeff = fit_log_line(radius, corr, r_min, r_max)
                            all_exp.append(xi)
                            ax.scatter(pos, xi, marker='^', color='tab:blue')
                            ax.annotate(str(plate) + ", " + well, (pos, xi)) 
                print(len(all_exp), group, stim, t_cell)

                # Mean tracks
                exp_mean.append(np.mean(all_exp))
                exp_sd.append(np.std(all_exp))

                x_ticks.append(group + "\n " + stim + "\n " + t_cell)
                pos += 1
    ax.errorbar(np.arange(1, len(exp_mean)+1), exp_mean, yerr=exp_sd, fmt='o', capsize=3, color='k')
    ax.set_title(cell_type)
    set_axis_style(ax,x_ticks)
    ax.set_ylabel(r"Correlation length, $\xi$ (pixels))")
    ax.set_ylim(0,200)

    ax.set_title(cell_type + ", " + stim + ", " + t_cell + ", frame=" + str(frame))

    if save_plot == True:  
        plot_folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/plots/correlation_exponents"
        if not os.path.exists(plot_folder):
            os.makedirs(plot_folder)
        plot_name = cell_type + "_" + stim + "_" + t_cell + ".png"
        plt.savefig(os.path.join(plot_folder, plot_name))

    if show_plot == True:
        plt.show()

    plt.close()

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

for cell_type in ["Microglia & T-cells"]:
    t0 = time.time()
    for stim in all_stim:
        for t_cell in all_t_cell:
            stim_list = [stim]
            t_cell_list = [t_cell]
            plot_corr_exp_scatter(cell_type, plate_list, group_list, stim_list, t_cell_list, phase, frame, r_min, r_max, r_bin_num, save_plot=True, show_plot=False)
    print(time.time()-t0)


plate = 1723
well = "F4"
phase = "red"
# plot_corr_log(plate, well, phase, frames=[1,216], r_min=0, r_max=50, r_bin_num=50, corr_av=True, save_plot=True, show_plot=False, x_lim=True)
# plot_exponents_time(plate, well, phase, frames=np.arange(1, 212, 10), r_min=10, r_max=30, corr_av=True, save_plot=True, show_plot=True)

# for well in all_wells:
#     path = get_tif_file(plate, well, phase)
#     if os.path.exists(path):
#         path_plot = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/plots/correlation_functions/" + str(plate) + "/" + str(plate) + '_' + well + '_log.png'
#         if not os.path.exists(path_plot):
#             print(well)
#             # plot_corr_log(plate, well, phase, frames=[1,216], r_min=10, r_max=50, r_bin_num=25, corr_av=True, save_plot=True, show_plot=False)
#             # plot_corr_log(plate, well, phase, frames=[1,216], r_min=0, r_max=5, r_bin_num=5, corr_av=True, save_plot=True, show_plot=False, x_lim=True)

#             plot_exponents_time(plate, well, phase, frames=np.arange(1, 212, 10), r_min=10, r_max=30, corr_av=True, save_plot=True, show_plot=False)


