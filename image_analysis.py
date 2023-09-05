# from PIL import Image, ImageEnhance
from skimage import io
from matplotlib import pyplot as plt
import numpy as np
from well_plate_dictionary import *
import os

plate = 1723
exp = 107
phase = "red"
labels = ["Only Microglia", "Uninfected Healthy", "AC lPVL", "AC hPVL", "HAM"]


def get_intensity(exp, plate, well):
    folder = "/Volumes/T7 Shield/Incucyte_data/processed_data/"
    filename = folder + "exp_" + str(exp) + "/" + "VID" + str(plate) + '_' + str(phase) + "_" + well + '_1.tif'
    im = io.imread(filename)
    im = im/np.max(im)
    return im

def get_coverage(exp, plate, well, threshold=10000):
    folder = "/Volumes/T7 Shield/Incucyte_data/processed_data/"
    filename = folder + "exp_" + str(exp) + "/" + "VID" + str(plate) + '_' + str(phase) + "_" + well + '_1.tif'

    im = io.imread(filename)

    frames = np.arange(im.shape[0])

    coverage = []
    for i in frames:
        a = im[i] > threshold
        coverage.append(a.sum()/a.size * 100)

    times = frames / 3
    coverage = np.array(coverage)

    return times, coverage
    
## Plot coverage over time, averaged over each patient group
fig, ax = plt.subplots()
patient_groups = patient_group_dict.keys()
for group in patient_groups:
    coverage_all = []
    for well in patient_group_dict[group]:
        folder = "/Volumes/T7 Shield/Incucyte_data/processed_data/"
        filename = folder + "exp_" + str(exp) + "/" + "VID" + str(plate) + '_' + str(phase) + "_" + well + '_1.tif'
        if os.path.exists(filename):
            times, coverage = get_coverage(exp, plate, well)
            coverage_all.append(coverage)
    mean_coverage = np.mean(coverage_all, axis=0)
    ax.plot(times[2:], mean_coverage[2:], label=group)
ax.set_xlabel("Time (hrs)")
ax.set_ylabel("Percentage coverage (%)")
ax.set_ylim(bottom=0)
ax.legend()
plt.show()