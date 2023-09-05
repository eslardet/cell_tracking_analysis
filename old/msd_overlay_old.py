import os
import analysis_functions_xml as fun
import numpy as np
import time
from matplotlib import pyplot as plt
from well_plate_dictionary import *

save_plot=False
plate = '1712'
plot_name = plate + "_HAM2.png"
show_plot=True
log = False
max_t=15

# plate = "1712"
# well = "B2"
# pic = "1"

# file_path = "tracks/VID" + plate + "_day_red_" + well + "_" + pic + "_Tracks_bright.xml"

plate_dict = {'VID1711': "Astrocytes & T-cells", 'VID1712': "Microglia & T-cells", 'VID1713': "T-cells only"}
number_dict = {2:'low PVL, not stim. ', 3: 'low PVL, stim. ', 4: 'high PVL1, not stim.', 5: 'high PVL1, stim.', 6: 'high PVL2, not stim.', 7: 'high PVL2, stim.', 8: 'HAM1, not stim.', 9: 'HAM1, stim.', 10: 'HAM2, not stim. ', 11: 'HAM2, stim. '}
cell_dict = {'B': 'non-specific CD4', 'C': 'specific CD4', 'D': 'not-specific CD8', 'E': 'specific CD8'}
cell_dict_tcell = {'B': 'non-specific CD4', 'C': 'specific CD4 + non-specific CD8', 'D': 'specific CD8'}


fig, ax = plt.subplots()

for cell in ['B', 'C', 'D', 'E']:
    for number in np.arange(10,12):
        well = cell + str(number)

        xml_folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/xml/" + plate
        phase = 'red'
        # well = 'B10'
        pic = '1'

        filename = 'VID' + plate + '_' + phase + '_' + well + '_' + pic
        file_path = os.path.join(xml_folder, filename + '.xml')

        # tracks = fun.read_xml(file_path)

        msd = fun.get_msd(file_path, plot=False)[:max_t]
        # msd = fun.get_sdt(file_path, plot=False)
        if plate == "VID1713":
            if log == True:
                ax.loglog(range(len(msd)), msd, label=number_dict[number] + cell_dict_tcell[cell])
            else:
                ax.plot(range(len(msd)), msd, label=number_dict[number] + cell_dict_tcell[cell])
        else:
            if log == True:
                ax.loglog(range(len(msd)), msd, label=number_dict[number] + cell_dict[cell])
            else:
                ax.plot(range(len(msd)), msd, label=number_dict[number] + cell_dict[cell])

if log == True:
    ax.loglog(np.arange(1, 3, 1), 500*np.arange(1, 3, 1), '--', color="gray", label=r"$\propto t$")
    # ax.set_ylim(bottom=10**(2))
else:
    ax.set_ylim(bottom=0)
# ax.loglog(np.arange(1, 3, 1), 500*np.arange(1, 3, 1), 'r--', label=r"$\propto t$")
ax.set_title(plate_dict[plate])
ax.legend(fontsize=8)
ax.set_xlabel(r"$\tau$ (hours)")
ax.set_ylabel("MSD (pixels)")
# ax.set_ylabel("SDT (pixels)")
ax.set_xlim(right=max_t)


if save_plot == True:
    plot_folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/plots/msd"
    plt.savefig(os.path.join(plot_folder, plot_name))

if show_plot == True:
    plt.show()