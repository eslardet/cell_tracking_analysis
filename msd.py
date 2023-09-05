import os
import analysis_functions_xml as fun
import numpy as np
import time
from matplotlib import pyplot as plt

plate_dict = {'VID1711': "Astrocytes & T-cells", 'VID1712': "Microglia & T-cells", 'VID1713': "T-cells only"}
number_dict = {2:'low PVL, not stim. ', 3: 'low PVL, stim. ', 4: 'high PVL1, not stim.', 5: 'high PVL1, stim.', 6: 'high PVL2, not stim.', 7: 'high PVL2, stim.', 8: 'HAM1, not stim.', 9: 'HAM1, stim.', 10: 'HAM2, not stim. ', 11: 'HAM2, stim. '}
cell_dict = {'B': 'non-specific CD4', 'C': 'specific CD4', 'D': 'not-specific CD8', 'E': 'specific CD8'}
cell_dict_tcell = {'B': 'non-specific CD4', 'C': 'specific CD4 + non-specific CD8', 'D': 'specific CD8'}

plate = "1727"
well = "B2"
phase = "green"
pic = "1"

# file_path = "tracks/VID" + plate + "_day_red_" + well + "_" + pic + "_Tracks_bright.xml"

fig, ax = plt.subplots()

xml_folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/xml/" + plate
# plate = 'VID1711'
# phase = 'red'
# # well = 'B10'
# pic = '1-1'

filename = "VID" + plate + '_' + phase + '_' + well + '_' + pic
# filename = "VID1713_red_B10_1"
file_path = os.path.join(xml_folder, filename + '.xml')

# tracks = fun.read_xml(file_path)

msd = fun.get_msd(file_path, plot=True)
# print(msd)
# fun.plot_sdt_vs_t0(file_path, tau_range=[2,3,4], plot_name=filename, max_frames=25, max_plot=25, log=False)


# print(msd[1])
# ax.loglog(range(len(msd)), msd, label=well)
# ax.loglog(np.arange(1, 3, 1), 500*np.arange(1, 3, 1), label=r"$\propto t$")

# ax.legend()
# ax.set_xlabel("Time (hours)")
# ax.set_ylabel("MSD (pixels)")
# plt.show()
