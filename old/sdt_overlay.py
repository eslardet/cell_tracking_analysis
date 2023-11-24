import os
import analysis_functions_tracks as fun
import numpy as np
import time
from matplotlib import pyplot as plt



plate = "VID1712"
cell = "E"
number = 11
phase = "red"
pic = "1"

plot_name = plate + "_HAM2.png"

plate_dict = {'VID1711': "Astrocytes & T-cells", 'VID1712': "Microglia & T-cells", 'VID1713': "T-cells only"}
number_dict = {2:'low PVL, not stim. ', 3: 'low PVL, stim. ', 4: 'high PVL1, not stim.', 5: 'high PVL1, stim.', 6: 'high PVL2, not stim.', 7: 'high PVL2, stim.', 8: 'HAM1, not stim.', 9: 'HAM1, stim.', 10: 'HAM2, not stim. ', 11: 'HAM2, stim. '}
if plate == "1713":
    cell_dict = {'B': 'non-specific CD4', 'C': 'specific CD4 + non-specific CD8', 'D': 'specific CD8'}
else:
    cell_dict = {'B': 'non-specific CD4', 'C': 'specific CD4', 'D': 'not-specific CD8', 'E': 'specific CD8'}


max_frames = 25
max_plot = 25

# file_path = "tracks/VID" + plate + "_day_red_" + well + "_" + pic + "_Tracks_bright.xml"

fig, ax = plt.subplots()

xml_folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/xml_short"

filename = plate + '_' + phase + '_' + cell + str(number) + '_' + pic
file_path = os.path.join(xml_folder, filename + '.xml')

for cell in ['B', 'D']:
    for number in np.arange(10,12):
        well = cell + str(number)

        xml_folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/xml_short"
        phase = 'red'
        # well = 'B10'
        pic = '1'

        filename = plate + '_' + phase + '_' + well + '_' + pic
        file_path = os.path.join(xml_folder, filename + '.xml')

        # tracks = fun.read_xml(file_path)

        sdt = fun.get_sdt_vs_t0(file_path, tau=2, max_frames=25)
        ax.plot(range(max_frames), sdt, label=number_dict[number] + cell_dict[cell])

ax.set_xlabel(r"$t_0$ (hours)")
ax.set_ylabel("SDT (pixels)")
ax.set_xlim(right=max_plot)
ax.set_ylim(bottom=0)
ax.legend()
ax.set_title(plate_dict[plate])

plot_folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/plots/sdt_vs_t0"
plt.savefig(os.path.join(plot_folder, plot_name))
