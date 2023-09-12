import numpy as np
from image_analysis_functions import *
from scipy import ndimage as ndi
from skimage import measure, segmentation, filters
from skimage.measure import manders_coloc_coeff, intersection_coeff


def get_manders_coeff(plate, well, slice, green_thresh=0.25, red_thresh=0.15):
    im_red = get_image(plate, well, "red")[slice]
    im_green = get_image(plate, well, "green")[slice]

    # ## First segment with thresholding
    im_red_thresh = im_red/np.max(im_red) > red_thresh
    im_green_thresh = im_green/np.max(im_green) > green_thresh

    coeff = measure.manders_coloc_coeff(im_green_thresh, im_red_thresh)

    return coeff

def plot_manders_vs_time(plate, well, slice_range, green_thresh=0.25, red_thresh=0.15):
    coeff = []
    for i in slice_range:
        coeff.append(get_manders_coeff(plate, well, i, green_thresh, red_thresh))

    fig,ax = plt.subplots()

    ax.plot(slice_range/3, coeff, 'o-')
    ax.set_xlabel("Time (hours)")
    ax.set_ylabel("Manders' coefficient")
    ax.set_ylim([0,1])

    plot_folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/plots/colocalization/"
    plt.savefig(plot_folder + str(plate) + '_' + well + '.png', bbox_inches='tight')
    plt.show()

plate = 1737
well = "C2"

slice_range = np.arange(0, 220, 10)

plot_manders_vs_time(plate, well, slice_range)
