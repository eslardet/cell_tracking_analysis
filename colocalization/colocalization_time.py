import os, sys, inspect
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
sys.path.insert(0, parentdir)

import numpy as np
# from analysis_functions_xml import *
from analysis_functions_image import *


plate = 1723
well = "C2"

slice_range = np.arange(5, 216, 10)

t0 = time.time()
plot_manders_vs_time(plate, well, slice_range, show_fit=False)
print(time.time()-t0)

# print(get_mander_mult_time(plate, well, slice_range, mult_thresh=2, green_thresh=0.75, red_thresh=0.15))

# for well in all_wells:
#     t0 = time.time()
#     plot_manders_vs_time(plate, well, slice_range)
#     # print(get_manders_coeff(plate, well, -1))
#     print(time.time()-t0)
