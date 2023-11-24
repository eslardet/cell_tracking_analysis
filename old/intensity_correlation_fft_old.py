# from PIL import Image, ImageEnhance
from skimage import io
from matplotlib import pyplot as plt
import numpy as np
from analysis_functions_tracks import *
from well_plate_dictionary import *
from analysis_functions_image import *
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

def get_corr_binned(im, corr_r_max=500, r_bin_num=100, corr_av=True):
    corr_all, dist = get_corr_im(im, corr_av)
    corr_r_min = 0

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

def plot_corr(im_all, frame, corr_r_max=500, r_bin_num=100, corr_av=True, show_fit=True, show_plot=True, save_plot=False):
    im = im_all[frame]
    r_plot, corr_plot = get_corr_binned(im, corr_r_max=corr_r_max, r_bin_num=r_bin_num, corr_av=corr_av)

    fig, ax = plt.subplots()
    ax.plot(r_plot, corr_plot)
    ax.set_xlabel("Distance (pixels)")
    ax.set_ylabel("Correlation")

    if show_fit == True:
        r_plot, corr_plot = get_corr_binned(im, corr_r_max=corr_r_max, r_bin_num=r_bin_num, corr_av=corr_av)
        xi, coeff, C_inf = fit_curve(r_plot, corr_plot)
        ax.plot(r_plot, func(r_plot, xi, coeff, C_inf), '--', label=r"$\xi=$" + str(np.round(1/xi,2)))
        ax.legend()

    if save_plot == True:
        plot_folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/plots/correlation_functions/"
        plt.savefig(plot_folder + str(plate) + '_' + well + '_' + str(frame) + '.png', bbox_inches='tight')

    if show_plot == True:
        plt.show()


def func(r, xi, coeff, C_inf):
    return np.exp(coeff)*np.exp(-r/xi) + C_inf

def fit_curve(r_plot, corr):
    xi, coeff, C_inf = curve_fit(func, r_plot, corr)[0]
    return xi, coeff, C_inf

def get_exponent(im, corr_r_max=500, r_bin_num=100, corr_av=True):
    r_plot, corr_plot = get_corr_binned(im, corr_r_max=corr_r_max, r_bin_num=r_bin_num, corr_av=corr_av)
    xi, coeff, C_inf = fit_curve(r_plot, corr_plot)
    return 1/xi

def plot_exponents_time(im_all, frames, corr_r_max=500, r_bin_num=100, corr_av=True, save_plot=False, show_plot=True):
    exp = []
    for i in frames:
        im = im_all[i]
        exp.append(get_exponent(im, corr_r_max=corr_r_max, r_bin_num=r_bin_num, corr_av=corr_av))

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

plot_corr(im, frame=216, corr_r_max=250, r_bin_num=100, corr_av=True, save_plot=True, show_plot=False)

# for well in all_wells:
#     if not os.path.exists(get_tif_file(plate, well, phase)):
#         print("File does not exist: " + str(plate) + ", " + str(well))
#     else:
#         im_all = get_image(plate, well, phase)
#         frames = np.arange(1, 202, 10)
#         plot_exponents_time(im_all, frames, save_plot=True, show_plot=False)