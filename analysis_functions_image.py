# from PIL import Image, ImageEnhance
from skimage import io
# from skimage.transform import resize
from matplotlib import pyplot as plt
import numpy as np
import time
from well_plate_dictionary import *
import os
from scipy.optimize import curve_fit
from scipy.signal import fftconvolve

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

def display_image(im):
    io.imshow(im/np.max(im))
    plt.show()

def get_intensity(plate, well, phase):
    im = get_image(plate, well, phase)
    im = im/np.max(im)
    return im

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