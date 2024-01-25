import os, sys, inspect
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
sys.path.insert(0, parentdir)

# from PIL import Image, ImageEnhance
from skimage import io
from matplotlib import pyplot as plt
import numpy as np
from well_plate_dictionary import *
from analysis_functions_image import *
import os
import time
from scipy.signal import correlate, correlate2d, fftconvolve, convolve
# from scipy.fftpack import fft2, ifft2
# from scipy.spatial import distance_matrix
# from scipy.spatial.distance import cdist
from intensity_correlation import get_corr_fft, get_corr_slow, get_corr_im

def get_corr_manual(im, corr_r_min, corr_r_max, r_bin_num):
    h,w = im.shape
    corr_list = []
    r_list = []
    for i in np.arange(h):
        for j in  np.arange(w):
            for k in np.arange(i,h):
                for l in np.arange(j,w):
                    corr = im[i,j]*im[k,l]
                    r = np.sqrt((i-k)**2+(j-l)**2)
                    corr_list.append(corr)
                    r_list.append(r)
    c0 = 0
    for i in np.arange(h):
        for j in np.arange(w):
            c0 += im[i,j]**2
    c0 /= h*w

    corr_all = np.array(corr_list)
    rij_all = np.array(r_list)
    # corr_r_max = np.max(rij_all)
    # corr_r_min = 0
    # r_bin_num = 50

    bin_size = (corr_r_max-corr_r_min) / r_bin_num
    r_plot = np.linspace(corr_r_min, corr_r_max, num=r_bin_num, endpoint=False) + bin_size/2

    corr_plot = []
    for i in range(r_bin_num):
        lower = r_plot[i]
        try:
            upper = r_plot[i+1]
        except:
            upper = corr_r_max+1
        idx = np.where((rij_all>lower)&(rij_all<upper))[0]
        if len(idx) != 0:
            corr = np.mean(corr_all[idx])
            corr_plot.append(np.abs(corr/c0))

    return r_plot, corr_plot

def get_corr_slow(im, corr_r_min, corr_r_max, r_bin_num, take_av=True):
    # im = get_image(plate, well, phase)[frame]
    h,w = im.shape
    corr = correlate2d(im,im, mode='full')[h-1:,w-1:]
    avs = np.outer(np.arange(h,0,-1), np.arange(w,0,-1))
    if take_av == True:
        corr = corr/avs
    corr = corr/corr[0,0]
    corr = corr.flatten()
    rij_all = get_distance_matrix(im).flatten()

    # corr_r_max = np.max(rij_all)
    # corr_r_min = 0
    # r_bin_num = 50

    bin_size = (corr_r_max-corr_r_min) / r_bin_num
    r_plot = np.linspace(corr_r_min, corr_r_max, num=r_bin_num, endpoint=False) + bin_size/2

    corr_plot = []
    for i in range(r_bin_num):
        lower = r_plot[i]
        try:
            upper = r_plot[i+1]
        except:
            upper = corr_r_max+1
        idx = np.where((rij_all>lower)&(rij_all<upper))[0]
        if len(idx) != 0:
            c = np.mean(corr[idx])
            corr_plot.append(np.abs(c))

    return r_plot, corr_plot

def get_corr_fft(im, corr_r_min, corr_r_max, r_bin_num, take_av=True):
    h,w = im.shape
    corr = np.round(fftconvolve(im,im[::-1,::-1], mode='full')[h-1:,w-1:],0)
    avs = np.outer(np.arange(h,0,-1), np.arange(w,0,-1))
    if take_av == True:
        corr = corr/avs
    corr = corr/corr[0,0]
    corr = corr.flatten()

    rij_all = get_distance_matrix(im).flatten()

    # corr_r_max = np.max(rij_all)
    # # corr_r_max = 100
    # corr_r_min = 0
    # r_bin_num = 50

    bin_size = (corr_r_max-corr_r_min) / r_bin_num
    r_plot = np.linspace(corr_r_min, corr_r_max, num=r_bin_num, endpoint=False) + bin_size/2

    corr_plot = []
    for i in range(r_bin_num):
        lower = r_plot[i]
        try:
            upper = r_plot[i+1]
        except:
            upper = corr_r_max+1
        idx = np.where((rij_all>lower)&(rij_all<upper))[0]
        if len(idx) != 0:
            c = np.mean(corr[idx])
            corr_plot.append(np.abs(c))

    return r_plot, corr_plot


n = 100
im = get_image(1737, "F4", "red")[-1][:n,:n]
im = im - np.mean(im)

corr_r_min = 0
corr_r_max = 20
r_bin_num = 20


fig, ax = plt.subplots()

t0 = time.time()
r_plot, corr_plot = get_corr_manual(im, corr_r_min, corr_r_max, r_bin_num)
t1 = time.time() - t0
ax.plot(r_plot, corr_plot, label="manual, t=" + str(np.round(t1,2)) + "s")

t0 = time.time()
r_plot, corr_plot = get_corr_slow(im, corr_r_min, corr_r_max, r_bin_num, take_av=True)
t1 = time.time() - t0
ax.plot(r_plot, corr_plot, label="correlate2d, t=" + str(np.round(t1,2)) + "s")

t0 = time.time()
r_plot, corr_plot = get_corr_fft(im, corr_r_min, corr_r_max, r_bin_num, take_av=True)
t1 = time.time() - t0
ax.plot(r_plot, corr_plot, label="fftconvolve, t=" + str(np.round(t1,2)) + "s")

# t0 = time.time()
# r_plot, corr_plot = get_corr_fft(im, corr_r_min, corr_r_max, r_bin_num, take_av=False)
# t1 = time.time() - t0
# ax.plot(r_plot, corr_plot, label="fftconvolve, no av")

# t0 = time.time()
# r_plot, corr_plot = get_corr_fft(im, corr_r_min, corr_r_max, r_bin_num, take_av=True)
# t1 = time.time() - t0
# ax.plot(r_plot, corr_plot, label="fftconvolve, with av")


ax.set_yscale("log")
# ax.set_xbound(0,100)
ax.legend()
plt.show()




# n = 5
# m = 5
# # im = np.arange(0,n*m).reshape((n,m))
# # im = np.tile(np.arange(1,6), (5,1))
# # im = np.arange(0,n*m).reshape((n,m))/10
# # im = np.random.randint(1,10, (n,m))
# # im = np.random.rand(n,m)
# im = get_image(1737, "F4", "red")[-1][:100,:100]
# im = im - np.mean(im)

# # im = np.array([[2049, 1195, 1195],
# #  [1423, 1480, 1081],
# #  [854, 1935, 1935]])

# n,m = im.shape
# print(n,m)

# t0 = time.time()

# corr = np.round(correlate2d(im,im, mode='full')[n-1:,m-1:],0)
# print(time.time()-t0)

# t0 = time.time()
# corr2 = np.round(fftconvolve(im,im[::-1,::-1], mode='full')[n-1:,m-1:],0)
# # avs = np.outer(np.arange(n,0,-1), np.arange(m,0,-1))
# print(time.time()-t0)

# print(np.all(corr==corr2))



# # def corr_fft(im):
# #     h,w = im.shape
# #     corr = np.round(fftconvolve(im,im[::-1,::-1], mode='full')[h-1:,w-1:],0)
# #     avs = np.outer(np.arange(h,0,-1), np.arange(w,0,-1))
# #     corr = corr/avs
# #     corr = corr/corr[0,0]
# #     return corr

# # def corr_corr(im):
# #     h,w = im.shape
# #     corr = correlate2d(im,im, mode='full')[h-1:,w-1:]
# #     return corr


# # corr3 = corr_corr(im)[0]
# # corr4 = corr_fft(im)[0]

# # print(np.all(corr3==corr4))