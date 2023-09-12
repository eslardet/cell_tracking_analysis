# from PIL import Image, ImageEnhance
from skimage import io
from matplotlib import pyplot as plt
import numpy as np
from well_plate_dictionary import *
from image_analysis_functions import *
import os
import time
import freud

plate = 1728
well = "G4"
phase = "red"

im = get_image(plate, well, phase)[-1]
h,w = im.shape
box = freud.Box.from_box([h, w])
freud.data.Unit(box)