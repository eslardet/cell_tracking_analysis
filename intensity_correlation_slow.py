# from PIL import Image, ImageEnhance
from skimage import io
from matplotlib import pyplot as plt
import numpy as np
from well_plate_dictionary import *
from analysis_functions_image import *
import os
import time
from numba import jit

t0 = time.time()

plate = 1723
exp = 107
phase = "red"
well = "G4"
labels = ["Only Microglia", "Uninfected Healthy", "AC lPVL", "AC hPVL", "HAM"]


# for well in ["F6", "G4", "C2", "C4", "C8"]:
#     folder = "/Volumes/T7 Shield/Incucyte_data/processed_data/"
#     filename = folder + "exp_" + str(exp) + "/" + "VID" + str(plate) + '_' + str(phase) + "_" + well + '_1.tif'
#     im = io.imread(filename)
#     # io.imshow(im[-1])
#     # print(np.max(im[-1]))
#     # print(np.mean(im[-1]/np.max(im[-1])))
#     # plt.show()

#     # Select image and normalize intensities to 0-1
#     # im_select = im[-1]/np.max(im[-1])
#     im_norm = im[-1]/np.max(im[-1])
#     # im_fluc = im_norm - np.mean(im_norm)
#     im_fluc = im_norm

#     c0 = np.mean(im_fluc**2)

#     corr_list = []
#     r_list = []

#     h,w = im_fluc.shape

#     for i in np.random.choice(range(h),200):
#         for j in  np.random.choice(range(w),200):
#             for k in np.random.choice(range(i-200,i+200),20):
#                 if 0 <= k < h:
#                     for l in np.random.choice(range(j-200, j+200),20):
#                         if 0 <= l < w:
#                             corr = im_fluc[i,j]*im_fluc[k,l]
#                             r = np.sqrt((i-k)**2+(j-l)**2)
#                             corr_list.append(corr)
#                             r_list.append(r)

#     corr_all = np.array(corr_list)
#     rij_all = np.array(r_list)
#     corr_r_max = np.max(rij_all)
#     corr_r_min = 0
#     r_bin_num = 20

#     bin_size = (corr_r_max-corr_r_min) / r_bin_num
#     r_plot = np.linspace(corr_r_min, corr_r_max, num=r_bin_num, endpoint=False) + bin_size/2

#     corr_plot = []
#     for i in range(r_bin_num):
#         lower = r_plot[i]
#         try:
#             upper = r_plot[i+1]
#         except:
#             upper = corr_r_max+1
#         idx = np.where((rij_all>lower)&(rij_all<upper))[0]
#         if len(idx) != 0:
#             corr = np.mean(corr_all[idx])
#             corr_plot.append(corr / c0)

#     print(time.time() - t0)

#     fig, ax = plt.subplots()
#     info = get_well_info(well)
#     ax.plot(r_plot, np.abs(corr_plot), '-', label=info)
#     ax.set_xlabel("r (pixels)")
#     ax.set_ylabel("C(r)")
#     ax.set_ylim([0,1])
#     ax.legend()

#     plot_name = "VID" + str(plate) + '_' + str(phase) + "_" + well + ".png"
#     plot_folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/plots/intensity_correlation"
#     plt.savefig(os.path.join(plot_folder, plot_name))

#     plt.close()

#     io.imshow(im[-1])
#     plot_name = "VID" + str(plate) + '_' + str(phase) + "_" + well + "_image.png"
#     plt.savefig(os.path.join(plot_folder, plot_name))
#     # plt.show()


@jit
def corr_calc(im):
    h,w = im.shape
    for i in np.arange(h):
        for j in  np.arange(w):
            for k in np.arange(10):
                if 0 <= k < h:
                    for l in np.arange(10):
                        if 0 <= l < w:
                            corr = im[i,j]*im[k,l]
                            r = np.sqrt((i-k)**2+(j-l)**2)
                            # corr_list.append(corr)
                            # r_list.append(r)

t0 = time.time()
im = get_image(plate, well, phase)[-1]
im = im/np.max(im)
# corr_calc(im)
# print(time.time() - t0)

# x = np.tile(np.arange(0,n), (n,1))
# x = np.random.rand(n,n)
# print(x)

