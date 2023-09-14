# from PIL import Image, ImageEnhance
from skimage import io
from matplotlib import pyplot as plt
import numpy as np
from well_plate_dictionary import *
from image_analysis_functions import *
import os
import time
from scipy.signal import correlate, correlate2d, fftconvolve
from scipy.fftpack import fft2, ifft2
from scipy.spatial import distance_matrix
from scipy.spatial.distance import cdist

plate = 1728
well = "G4"
phase = "red"

# im = get_image(plate, well, phase)[-1]

im = np.random.rand(100,100)

n=10
im = np.arange(0,n**2).reshape((n,n))
print(im)
# corr = correlate2d(im,im, mode='full')[n-1:,n-1:]
corr = np.round(fftconvolve(im,im[::-1,::-1], mode='full')[n-1:,n-1:],0)
print(corr)

avs = np.outer(np.arange(n,0,-1), np.arange(n,0,-1))
print(corr/avs)

dist = get_distance_matrix(im)

fig, ax = plt.subplots()
c0 = corr[0,0]/avs[0,0]
# ax.scatter(dist, corr/corr[0,0])
# ax.scatter(dist, corr/avs/c0)

corr_all = np.array(corr/avs).flatten()
rij_all = np.array(dist).flatten()
corr_r_max = np.max(rij_all)
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
    idx = np.where((rij_all>lower)&(rij_all<upper))[0]
    if len(idx) != 0:
        corr = np.mean(corr_all[idx])
        corr_plot.append(corr/c0)

ax.scatter(r_plot, corr_plot)



corr_list = []
r_list = []
for i in np.arange(n):
    for j in  np.arange(n):
        for k in np.arange(i,n):
            for l in np.arange(j,n):
                corr = im[i,j]*im[k,l]
                r = np.sqrt((i-k)**2+(j-l)**2)
                corr_list.append(corr)
                r_list.append(r)

corr_all = np.array(corr_list)
rij_all = np.array(r_list)
corr_r_max = np.max(rij_all)
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
    idx = np.where((rij_all>lower)&(rij_all<upper))[0]
    if len(idx) != 0:
        corr = np.mean(corr_all[idx])
        corr_plot.append(corr/c0)

ax.scatter(r_plot, corr_plot)
plt.show()

# im = np.arange(0,n)
# print(im)
# print(np.flip(np.correlate(im,im, mode='full')[:n]))


# x = np.random.rand(3,3)
# x = np.pad(x, ([0,2],[0,2]), mode='constant')
# x = fft2(x)
# # x[:,1] = np.conj(x[:,1])
# corr = ifft2(np.prod(x,x))

# print(np.abs(corr))
