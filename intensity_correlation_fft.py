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

im = get_image(plate, well, phase)[-1]

# def corr_fft(im):
# x = im
# x = np.fft2(np.pad(x,([0,0],[0,0],[0,34],[0,34]),mode='constant'))

# xp = (x - np.average(x))/np.std(x)
# f= fft2(xp)
# p = np.abs(f)**2
# pi = ifft2(p)
# corr = np.real(pi)[:int(len(xp)/2)]/int(len(xp))


# plt.plot(np.arange(len(corr)), corr)
# plt.show()

x = im

x = np.random.rand(3,3)
x.flatten()

y = np.random.rand(3,3)
y.flatten()

# fftconvolve(x,x, mode='same')
# t0 = time.time()
# a = correlate(x,y)
# print(time.time()-t0)

# t0 = time.time()
# b = fftconvolve(x,x)
# print(time.time()-t0)

# t0 = time.time()
# count = 0
# h,w = x.shape
# for i in np.arange(h):
#     for j in  np.arange(w):
#         for k in np.arange(h):
#             for l in np.arange(w):
#                     corr = x[i,j]*x[k,l]
#                     count += 1
# print(count)
# print(time.time()-t0)

# print(a)
# print(b)

# ## To get pairwise distances
# h,w = im.shape
# n = max(h,w)
# x = np.tile(np.arange(0,n), (n,1))
# y = x.transpose()
# # Look up distance by computing difference between coordinates
# dist = np.sqrt(x**2+y**2)

x = np.random.rand(3,3)
x = np.pad(x, ([0,2],[0,2]), mode='constant')
x = fft2(x)
# x[:,1] = np.conj(x[:,1])
corr = ifft2(np.prod(x,x))

# print(np.abs(corr))
