# from PIL import Image, ImageEnhance
from skimage import io
# from skimage.transform import resize
from matplotlib import pyplot as plt
from matplotlib import cm
import numpy as np
import time
from well_plate_dictionary import *
import os
from scipy.optimize import curve_fit
from scipy.signal import fftconvolve
from skimage import measure, segmentation, filters
from skimage.measure import manders_coloc_coeff, intersection_coeff, pearson_corr_coeff

def set_axis_style(ax, labels):
    ax.get_xaxis().set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_xlim(0.25, len(labels) + 0.75)

def get_image(plate, well, phase):
    folder = "/Volumes/T7 Shield/Incucyte_data/processed_data/" + str(plate) + "/"
    filename = folder + "VID" + str(plate) + '_' + str(phase) + "_" + well + '_1.tif'
    im = io.imread(filename)
    return im

def get_tif_file(plate, well, phase):
    folder = "/Volumes/T7 Shield/Incucyte_data/processed_data/" + str(plate) + "/"
    filename = folder + "VID" + str(plate) + '_' + str(phase) + "_" + well + '_1.tif'
    return filename

def get_xy_pixels(plate, well, phase, frame):
    folder = "/Volumes/T7 Shield/Incucyte_data/t-cell-intensity/" + str(plate) + "/"
    filename = folder + "VID" + str(plate) + '_' + str(phase) + "_" + well + '_1.txt'
    intensity = np.zeros((972, 1296))
    with open(filename) as f:
        for line in f:
            l = line.split("\t")
            if l[2] == str(frame):
                intensity[int(l[1]), int(l[0])] = float(l[3])
    return intensity

def get_xy_pixels_all(plate, well, phase):
    folder = "/Volumes/T7 Shield/Incucyte_data/t-cell-intensity/" + str(plate) + "/"
    filename = folder + "VID" + str(plate) + '_' + str(phase) + "_" + well + '_1.txt'
    intensity = np.zeros((217, 972, 1296))
    frame = 0
    with open(filename) as f:
        for line in f:
            l = line.split("\t")
            if l[2] == str(frame):
                intensity[frame, int(l[1]), int(l[0])] = float(l[3])
            else:
                frame += 1
                if l[2] == str(frame):
                    intensity[frame, int(l[1]), int(l[0])] = float(l[3])
    return intensity

def display_image(im):
    io.imshow(im/np.max(im), cmap=cm.gray)
    plt.show()

def get_intensity(plate, well, phase):
    im = get_image(plate, well, phase)
    im = im/np.max(im)
    return im

def scaled_threshold(im, threshold=0.75):
    im_mean = im - np.mean(im)
    im_scaled = im_mean/np.max(im_mean)
    im_thresh = im_scaled > threshold

    return im_thresh

def get_coverage(plate, well, phase, threshold=10000):
    im = get_image(plate, well, phase)

    frames = np.arange(im.shape[0])

    coverage = []
    for i in frames:
        a = im[i] > threshold
        coverage.append(a.sum()/a.size * 100)

    times = frames / 3
    coverage = np.array(coverage)

    return times, coverage

def resize_image(im, scale_factor):
    h,w = im.shape
    h_scaled = h//scale_factor
    w_scaled = w//scale_factor
    print(h_scaled,w_scaled)

    im_scaled = np.zeros((h_scaled, w_scaled))

    for i in range(h_scaled):
        for j in range(w_scaled):
            im_scaled[i,j] = np.mean(im[i*scale_factor:(i+1)*scale_factor, j*scale_factor:(j+1)*scale_factor])
    return im_scaled


def get_distance_matrix(im):
    """
    Output matrix is distance shift matrix in terms of x, y pixels
    """
    h,w = im.shape
    n = max(h,w)
    x = np.tile(np.arange(0,n), (n,1))
    y = x.transpose()
    dist = np.sqrt(x**2+y**2)
    return dist[:h,:w]



def get_corr_fft(im, fluc=True, take_av=True):
    h,w = im.shape
    if fluc == True:
        im = im - np.mean(im)
    corr = np.round(fftconvolve(im,im[::-1,::-1], mode='full')[h-1:,w-1:],0)
    avs = np.outer(np.arange(h,0,-1), np.arange(w,0,-1))
    if take_av == True:
        corr = corr/avs
    corr = corr/corr[0,0]
    return corr

# def get_corr_r(plate, well, phase, frame):
#     im = get_image(plate, well, phase)[frame]
#     corr = get_corr_fft(im).flatten()
#     dist = get_distance_matrix(im).flatten()

#     return corr, dist

def get_corr_im(im, take_av=True):
    # im = get_image(plate, well, phase)[frame]
    corr = get_corr_fft(im, take_av).flatten()
    dist = get_distance_matrix(im).flatten()

    return corr, dist

def get_corr_binned(im, corr_r_min=0, corr_r_max=250, r_bin_num=200, corr_av=True):
    corr_all, dist = get_corr_im(im, corr_av)

    bin_size = (corr_r_max-corr_r_min) / r_bin_num
    r_plot = np.linspace(corr_r_min, corr_r_max+1, num=r_bin_num, endpoint=False) + bin_size/2
    corr_plot = []
    for i in range(r_bin_num):
        lower = r_plot[i]
        try:
            upper = r_plot[i+1]
        except:
            upper = corr_r_max+1
        idx = np.where((dist>lower)&(dist<upper))[0]
        if len(idx) != 0:
            corr = np.mean(corr_all[idx])
            corr_plot.append(corr)

    return r_plot, corr_plot

def exp_func(r, xi, coeff):
    return np.exp(coeff)*np.exp(-r/xi)

def fit_exp_curve(r_plot, corr):
    xi, coeff = curve_fit(exp_func, r_plot, corr)[0]
    return xi, coeff

def fit_log_line(r_plot, corr, r_min, r_max):
    r_plot = np.array(r_plot)
    corr = np.array(corr)
    idx1 = np.where(r_plot<r_max)[0]
    idx2 = np.where(r_plot>r_min)[0]
    idx = list(set(idx1) & set(idx2))
    corr = np.abs(corr[idx])
    r_plot = r_plot[idx]

    alpha, coeff = np.polyfit(r_plot, np.log(corr), deg=1)
    xi = -1/alpha
    return xi, coeff

def get_exponent(im, r_min, r_max, corr_av=True):
    radius, corr = get_corr_binned(im, corr_r_min=r_min, corr_r_max=r_max, r_bin_num=int(r_max-r_min), corr_av=corr_av)
    xi, coeff = fit_log_line(radius, corr, r_min, r_max)
    return xi


## Correlation Plotting ##

def plot_corr_log(plate, well, phase, frames, r_min, r_max, r_bin_num=100, threshold=False, corr_av=True, save_plot=False, show_plot=True, x_lim=False):
    im = get_image(plate, well, phase)
    # if threshold == True:
        # im = im/np.max(im) > 0.99

    fig, ax = plt.subplots()

    for i in frames:
        im_frame = im[i]
        if threshold == True:
            im_mean = im_frame - np.mean(im_frame)
            im_scaled = im_mean/np.max(im_mean)
            im_frame = im_scaled > 0.75
        radius, corr = get_corr_binned(im_frame, corr_av=corr_av, corr_r_min=0, corr_r_max=r_max, r_bin_num=r_bin_num)

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


def plot_exponents_time(plate, well, phase, frames, r_min, r_max, threshold=False, corr_av=True, save_plot=False, show_plot=True):
    im_all = get_image(plate, well, phase)
    
    exp = []
    for i in frames:
        im = im_all[i]
        if threshold == True:
            # im = im/np.max(im) > 0.99
            im_mean = im - np.mean(im)
            im_scaled = im_mean/np.max(im_mean)
            im = im_scaled > 0.75
        exp.append(get_exponent(im, r_min, r_max, corr_av))

    fig, ax = plt.subplots()
    ax.plot(frames/3, exp, '-o', label=str(plate), color="tab:" + phase)
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

## Colocalization ##

def get_manders_coeff(plate, well, slice, green_thresh=0.25, red_thresh=0.15):
    im_red = get_image(plate, well, "red")[slice]
    im_green = get_image(plate, well, "green")[slice]

    ## First segment with thresholding
    im_red_thresh = im_red/np.max(im_red) > red_thresh
    # im_green_thresh = im_green/np.max(im_green) > green_thresh
    im_mean = im_green - np.mean(im_green)
    im_scaled = im_mean/np.max(im_mean)
    im_green_thresh = im_scaled > green_thresh
    

    coeff = measure.manders_coloc_coeff(im_green_thresh, im_red_thresh)
    return coeff

def get_manders_increase(plate, well, compare, green_thresh=0.25, red_thresh=0.15):
    coeff0 = get_manders_coeff(plate, well, compare[0], green_thresh, red_thresh)
    coeff1 = get_manders_coeff(plate, well, compare[1], green_thresh, red_thresh)
    
    percentage_increase = (coeff1-coeff0)/coeff0*100

    return percentage_increase

## Colocalization Plotting ##

def plot_manders_vs_time(plate, well, slice_range, green_thresh=0.25, red_thresh=0.15, time_av=False, show_plot=False, save_plot=True):
    green_path = get_tif_file(plate, well, "green")
    red_path = get_tif_file(plate, well, "red")
    if not os.path.exists(green_path) or not os.path.exists(red_path):
        print("File does not exist: " + str(plate) + ", " + str(well))
        return
    coeff = []
    for i in slice_range:
        if time_av == True:
            av = 0
            for j in range(10):
                av += get_manders_coeff(plate, well, i+j, green_thresh, red_thresh)
            coeff.append(av/10)
        else:
            coeff.append(get_manders_coeff(plate, well, i, green_thresh, red_thresh))

    fig,ax = plt.subplots()

    ax.plot(slice_range/3, coeff, 'o-')
    ax.set_xlabel("Time (hours)")
    ax.set_ylabel("Manders' colocalization coefficient")
    ax.set_ylim([0,1])

    if save_plot == True:
        plot_folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/plots/colocalization/" + str(plate) + "/"
        if not os.path.exists(plot_folder):
            os.makedirs(plot_folder)
        plt.savefig(plot_folder + str(plate) + '_' + well + '.png', bbox_inches='tight')
    if show_plot == True:
        plt.show()
    plt.close()

def plot_manders_increase(cell_type, group_list, stim_list, t_cell_list, slice_compare, green_thresh=0.25, red_thresh=0.15, time_av=False, show_plot=False, save_plot=True):
    plot_folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/plots/colocalization_increase/"
    if not os.path.exists(plot_folder):
        os.makedirs(plot_folder)
    plot_name = cell_type + "_" + stim_list[-1] + "_" + t_cell_list[-1] + "_" + str(slice_compare[1]) + ".png"
    if os.path.exists(os.path.join(plot_folder, plot_name)):
        print("Already plotted!")
        print(os.path.join(plot_folder, plot_name))
        return
    
    fig, ax = plt.subplots(figsize=(12,6))

    x_ticks = []

    plate_list = get_plates(cell_type)
    pos = 1
    for group in group_list:
        for stim in stim_list:
            for t_cell in t_cell_list:
                well_list = get_wells(group, stim, t_cell)
                for plate in plate_list:
                    for well in well_list:
                        green_path = get_tif_file(plate, well, "green")
                        red_path = get_tif_file(plate, well, "red")
                        if os.path.exists(green_path) and os.path.exists(red_path):
                    
                            percentage_increase = get_manders_increase(plate, well, slice_compare, green_thresh, red_thresh)
                            if percentage_increase > 0:
                                ax.scatter(pos, percentage_increase, color='k')
                                ax.annotate(str(plate) + ", " + well, (pos, percentage_increase))
                x_ticks.append(group + "\n " + stim + "\n " + t_cell)
                pos += 1
    ax.set_ylabel("Percentage increase in Manders' Coefficient (%)")
    ax.set_ylim(0,800)
    set_axis_style(ax,x_ticks)

    ax.set_title(cell_type + ", compare slices " + str(slice_compare[0]) + " and " + str(slice_compare[1]))

    if save_plot == True:
        plot_folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/plots/colocalization_increase/"
        if not os.path.exists(plot_folder):
            os.makedirs(plot_folder)
        plot_name = cell_type + "_" + stim + "_" + t_cell + "_" + str(slice_compare[1]) + ".png"
        plt.savefig(os.path.join(plot_folder, plot_name))

    if show_plot == True:
        plt.show()

    plt.close()