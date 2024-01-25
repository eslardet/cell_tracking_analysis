import os, sys, inspect
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
sys.path.insert(0, parentdir)

from analysis_functions_image import *
from analysis_functions_tracks import *
import tifffile as tiff
from tifffile import imread, imwrite
# from skimage.io import imread, imsave

def save_scaled_intensity(plate, well, phase):
    folder = "/Volumes/T7 Shield/Incucyte_data/processed_data_rescaled/" + str(plate) + "/"
    if not os.path.exists(folder):
        os.makedirs(folder)
    filename = folder + "VID" + str(plate) + '_' + str(phase) + "_" + well + '_1.tif'
    save_path = os.path.join(folder, filename)

    with tiff.TiffFile(filename) as tif:
        imagej_metadata = tif.imagej_metadata

    im = get_image(plate, well, phase)
    # im = im[:5]
    # im_scaled = np.zeros(im.shape, dtype=np.uint16)
    # for i in range(im.shape[0]):
    #     im_thresh = scaled_threshold(im[i], threshold=0.75)
    #     for j in range(im_thresh.shape[0]):
    #         for k in range(im_thresh.shape[1]):
    #             if im_thresh[j][k]:
    #                 im_scaled[i][j][k] = im[i][j][k]

    im_scaled = scaled_threshold(im, threshold=0.75).astype(np.uint16)*65535
    # imwrite(save_path, im_scaled, imagej=True, metadata=imagej_metadata)
    imwrite(save_path, im_scaled, imagej=True, metadata={'axes' : 'TYX'})


save_scaled_intensity(1737, "G2", "green")


# plate = 1737
# well = "G2"
# phase = "green"

# folder = "/Volumes/T7 Shield/Incucyte_data/processed_data/" + str(plate) + "/"
# filename = folder + "VID" + str(plate) + '_' + str(phase) + "_" + well + '_1.tif'
# # im = tiff.imread(filename)

# with tiff.TiffFile(filename) as tif:
#     im = tif.asarray()
#     imagej_metadata = tif.imagej_metadata
# print(imagej_metadata)

