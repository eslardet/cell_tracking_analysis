# from PIL import Image, ImageEnhance
from skimage import io
# from skimage.transform import resize
from matplotlib import pyplot as plt
import numpy as np
import time
from well_plate_dictionary import *
import os

plate = 1723
exp = 107
phase = "red"
labels = ["Only Microglia", "Uninfected Healthy", "AC lPVL", "AC hPVL", "HAM"]

def get_image(plate, well, phase):
    folder = "/Volumes/T7 Shield/Incucyte_data/processed_data/" + str(plate) + "/"
    filename = folder + "VID" + str(plate) + '_' + str(phase) + "_" + well + '_1.tif'
    im = io.imread(filename)
    return im

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


# im = resize_image(get_image("1723", "G2", "red")[-1], 4)
# display_image(im)