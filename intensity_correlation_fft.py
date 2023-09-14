# from PIL import Image, ImageEnhance
from skimage import io
from matplotlib import pyplot as plt
import numpy as np
from analysis_functions_xml import *
from well_plate_dictionary import *
from image_analysis_functions import *
import os
import time
from scipy.signal import correlate, correlate2d, fftconvolve
from scipy.fftpack import fft2, ifft2
from scipy.spatial import distance_matrix
from scipy.spatial.distance import cdist
from scipy.optimize import curve_fit


def get_corr_slow(im):
    # im = get_image(plate, well, phase)[frame]
    h,w = im.shape
    corr = correlate2d(im,im, mode='full')[h-1:,w-1:]
    avs = np.outer(np.arange(h,0,-1), np.arange(w,0,-1))
    corr = corr/avs
    corr = corr/corr[0,0]
    corr = corr.flatten()
    dist = get_distance_matrix(im).flatten()

    return corr, dist

def plot_corr_slow(im):
    # corr_all, dist = get_corr_r(plate, well, phase, frame)
    corr_all, dist = get_corr_slow(im)
    corr_r_max = np.max(dist)
    corr_r_min = 0
    r_bin_num = 10

    bin_size = (corr_r_max-corr_r_min) / r_bin_num
    r_plot = np.linspace(corr_r_min, corr_r_max, num=r_bin_num, endpoint=False) + bin_size/2

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

    fig, ax = plt.subplots()
    ax.plot(r_plot, corr_plot)
    plt.show()



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
    r_plot = np.linspace(corr_r_min, corr_r_max, num=r_bin_num, endpoint=False) + bin_size/2
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

def func(r, xi, coeff):
    return np.exp(coeff)*np.exp(-r/xi)

def fit_curve(r_plot, corr):
    xi, coeff = curve_fit(func, r_plot, corr)[0]
    return xi, coeff

def fit_log_line(r_plot, corr, r_min, r_max):
    r_plot = np.array(r_plot)
    corr = np.array(corr)
    idx1 = np.where(r_plot<r_max)[0]
    idx2 = np.where(r_plot>r_min)[0]
    idx = list(set(idx1) & set(idx2))
    corr = corr[idx]
    r_plot = r_plot[idx]

    alpha, coeff = np.polyfit(r_plot, np.log(corr), deg=1)
    xi = -1/alpha
    return xi, coeff


def plot_corr_log(im, frames, r_min, r_max, corr_av=True, save_plot=False, show_plot=True, x_lim=False):
    fig, ax = plt.subplots()

    for i in frames:
        radius, corr = get_corr_binned(im[i], corr_av=corr_av, corr_r_min=r_min, corr_r_max=r_max*2, r_bin_num=100)

        ax.plot(radius, corr, label="Frame=" + str(i))
        xi, coeff = fit_log_line(radius, corr, r_min, r_max)
        print(xi, coeff)
        r_plot = np.linspace(r_min, r_max, 100)
        # # r_plot = np.logspace(np.log(10), np.log(50), 100, base=np.exp(1))
        ax.plot(r_plot, func(r_plot, xi, coeff), '--', label=r"$\xi=$" + str(round(xi,2)))

        ax.set_yscale('log')

    ax.set_xlabel("r (pixels)")
    ax.set_ylabel("C(r)")
    # ax.set_ylim([0,1])
    if x_lim == True:
        ax.set_xlim(0,r_max*2)
    ax.legend()

    if save_plot == True:
        plate = "1737"
        well = "F4"
        plot_folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/plots/correlation_functions/"
        plt.savefig(plot_folder + str(plate) + '_' + well + '_log.png', bbox_inches='tight')

    if show_plot == True:
        plt.show()

def get_exponent(im, r_min, r_max, corr_av=True):
    radius, corr = get_corr_binned(im, corr_r_min=r_min, corr_r_max=r_max, r_bin_num=int(r_max-r_min), corr_av=corr_av)
    xi, coeff = fit_log_line(radius, corr, r_min, r_max)
    return xi


def plot_exponents_time(im_all, frames, r_min, r_max, corr_av=True, save_plot=False, show_plot=True):
    exp = []
    for i in frames:
        im = im_all[i]
        exp.append(get_exponent(im, r_min, r_max, corr_av))

    fig, ax = plt.subplots()
    ax.plot(frames/3, exp, '-o')
    ax.set_xlabel("Time (hours)")
    ax.set_ylabel(r"Correlation length $\xi$")
    if save_plot == True:
        plot_folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/plots/correlation_exponents/" + str(plate) + "/"
        if not os.path.exists(plot_folder):
            os.makedirs(plot_folder)
        plt.savefig(plot_folder + str(plate) + '_' + well + '.png', bbox_inches='tight')

    if show_plot == True:
        plt.show()

    plt.close()

plate = 1737
well = "F4"
phase = "red"

im = get_image(plate, well, phase)

# plot_corr_log(im, frames=[0,216], r_min=10, r_max=50, corr_av=True, save_plot=False, show_plot=True)

# plot_exponents_time(im_all=im, frames=np.arange(1, 212, 10), r_min=10, r_max=50, corr_av=True, save_plot=True, show_plot=True)

# for well in all_wells:
#     if not os.path.exists(get_tif_file(plate, well, phase)):
#         print("File does not exist: " + str(plate) + ", " + str(well))
#     else:
#         im_all = get_image(plate, well, phase)
#         frames = np.arange(1, 212, 10)
#         plot_exponents_time(im_all=im_all, frames=np.arange(1, 212, 10), r_min=10, r_max=50, corr_av=True, save_plot=True, show_plot=False)